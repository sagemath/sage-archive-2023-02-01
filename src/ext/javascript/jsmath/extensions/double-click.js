/*
 *  extensions/double-click.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file allows users to double click on typeset mathematics
 *  to view the TeX source for the given expression.  It will be loaded
 *  automatically when needed, or can be loaded by
 *
 *    jsMath.Extension.Require('double-click');
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2005-2006 by Davide P. Cervone
 *
 *  Licensed under the Apache License, Version 2.0 (the "License");
 *  you may not use this file except in compliance with the License.
 *  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 *  Unless required by applicable law or agreed to in writing, software
 *  distributed under the License is distributed on an "AS IS" BASIS,
 *  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *  See the License for the specific language governing permissions and
 *  limitations under the License.
 */

/********************************************************************/

jsMath.Add(jsMath.Click,{

  dragging: 0,

  /*
   *  Create the hidden DIV used for the tex source window
   */
  Init: function () {
    this.source = jsMath.Setup.DIV("float",{display:'none'});
    this.source.innerHTML =
        '<div class="drag"><div class="close"></div></div>'
      + '<div class="source"><span></span></div>';
    this.drag = this.source.firstChild;
    this.tex  = this.drag.nextSibling.firstChild;
    this.drag.firstChild.onclick = jsMath.Click.CloseSource;
    this.drag.onmousedown = jsMath.Click.StartDragging;
    this.drag.ondragstart = jsMath.Click.False;
    this.drag.onselectstart = jsMath.Click.False;
    this.source.onclick = jsMath.Click.CheckClose;
  },
  False: function () {return false},

  /*
   *  Handle a double-click on an equation
   */
  DblClick: function (data) {
    var event = data[0]; var TeX = data[1];
    var event = jsMath.Click.Event(event);

    var source = jsMath.Click.source
    var tex = jsMath.Click.tex;

    source.style.visibility = 'hidden';
    source.style.display = ''; source.style.width = '';
    source.style.left = ''; source.style.top = '';
    tex.innerHTML = '';

    TeX = TeX.replace(/^\s+|\s+$/g,'');
    TeX = TeX.replace(/&/g,'&amp;');
    TeX = TeX.replace(/</g,'&lt;');
    TeX = TeX.replace(/>/g,'&gt;');
    TeX = TeX.replace(/\n/g,'<br/>');
    tex.innerHTML = TeX;

    var h = source.offsetHeight; var w;
    if (jsMath.Browser.msieDivWidthBug) {
      tex.className = 'source';          // Work around MSIE bug where
      w = tex.offsetWidth + 5;           // DIV's don't collapse to
      tex.className = '';                // their natural widths
    } else {
      w = source.offsetWidth;
    }
    w = Math.max(50,Math.min(w,.8*event.W,event.W-40));
    var x = Math.floor(event.x-w/2); var y = Math.floor(event.y-h/2);
    x = event.X + Math.max(Math.min(x,event.W-w-20),20);
    y = event.Y + Math.max(Math.min(y,event.H-h-5),5);

    source.style.left = x+'px'; source.style.top = y+'px';
    source.style.width = w+'px';
    source.style.visibility = '';
    jsMath.Click.left = x + event.X; jsMath.Click.top = y + event.Y;
    jsMath.Click.w = w; jsMath.Click.h = source.offsetHeight;

    jsMath.Click.DeselectText(x,y);
    return false;
  },

  /*
   *  Get window width, height, and offsets plus
   *  position of pointer relative to the window
   */
  Event: function (event) {
    var W = jsMath.window.innerWidth  || jsMath.document.body.clientWidth;
    var H = jsMath.window.innerHeight || jsMath.document.body.clientHeight;
    var X = jsMath.window.pageXOffset; var Y = jsMath.window.pageYOffset;
    if (X == null) {
      X = jsMath.document.body.clientLeft;
      Y = jsMath.document.body.clientTop;
    }
    var x = event.pageX; var y = event.pageY;
    if (x == null) {
      x = event.clientX; y = event.clientY;
      if (jsMath.browser == 'MSIE' && jsMath.document.compatMode == 'CSS1Compat') {
        X = jsMath.document.documentElement.scrollLeft;
        Y = jsMath.document.documentElement.scrollTop;
        W = jsMath.document.documentElement.clientWidth;
        H = jsMath.document.documentElement.clientHeight;
      } else {
        X = jsMath.document.body.scrollLeft;
        Y = jsMath.document.body.scrollTop;
      }
    } else {x -= X; y -= Y}

    return {x: x, y: y, W: W, H: H, X: X, Y: Y};
  },

  /*
   *  Unselect whatever text is selected (since double-clicking
   *  usually selects something)
   */
  DeselectText: function (x,y) {
    if (jsMath.window.getSelection && jsMath.window.getSelection().removeAllRanges)
      {jsMath.window.getSelection().removeAllRanges()}
    else if (jsMath.document.getSelection && jsMath.document.getSelection().removeAllRanges)
      {jsMath.document.getSelection().removeAllRanges()}
    else if (jsMath.document.selection && jsMath.document.selection.empty)
      {jsMath.document.selection.empty()}
    else {
      /* Hack to deselect the text in Opera and Safari */
      if (jsMath.browser == 'MSIE') return;  // don't try it if MISE on Mac
      jsMath.hiddenTop.innerHTML =
        '<textarea style="visibility:hidden" rows="1" cols="1">a</textarea>';
      jsMath.hiddenTop.firstChild.style.position = 'absolute';
      jsMath.hiddenTop.firstChild.style.left = x+'px';
      jsMath.hiddenTop.firstChild.style.top  = y+'px';
      setTimeout(jsMath.Click.SelectHidden,1);
    }
  },
  SelectHidden: function () {
    jsMath.hiddenTop.firstChild.focus();
    jsMath.hiddenTop.firstChild.select();
    jsMath.hiddenTop.innerHTML = '';
  },

  /*
   *  Close the TeX source window
   */
  CloseSource: function () {
    jsMath.Click.tex.innerHTML = '';
    jsMath.Click.source.style.display = 'none';
    jsMath.Click.source.style.visibility = 'hidden';
    jsMath.Click.StopDragging();
    return false;
  },
  CheckClose: function (event) {
    if (!event) {event = jsMath.window.event}
    if (event.altKey) {jsMath.Click.CloseSource(); return false}
  },

  /*
   *  Set up for dragging the source panel
   */
  StartDragging: function (event) {
    if (!event) {event = jsMath.window.event}
    if (jsMath.Click.dragging) {jsMath.Click.StopDragging(event)}
    var event = jsMath.Click.Event(event);
    jsMath.Click.dragging = 1;
    jsMath.Click.x = event.x + 2*event.X - jsMath.Click.left;
    jsMath.Click.y = event.y + 2*event.Y - jsMath.Click.top;
    jsMath.Click.oldonmousemove = jsMath.document.body.onmousemove;
    jsMath.Click.oldonmouseup = jsMath.document.body.onmouseup;
    jsMath.document.body.onmousemove = jsMath.Click.DragSource;
    jsMath.document.body.onmouseup = jsMath.Click.StopDragging;
    return false;
  },

  /*
   *  Stop dragging the source window
   */
  StopDragging: function (event) {
    if (jsMath.Click.dragging) {
      jsMath.document.body.onmousemove = jsMath.Click.oldonmousemove;
      jsMath.document.body.onmouseup   = jsMath.Click.oldonmouseup;
      jsMath.Click.oldonmousemove = null;
      jsMath.Click.oldonmouseup   = null;
      jsMath.Click.dragging = 0;
    }
    return false;
  },

  /*
   *  Move the source window (but stay within the browser window)
   */
  DragSource: function (event) {
    if (!event) {event = jsMath.window.event}
    if (jsMath.Browser.buttonCheck && !event.button) {return jsMath.Click.StopDragging(event)}
    event = jsMath.Click.Event(event);
    var x = event.x + event.X - jsMath.Click.x;
    var y = event.y + event.Y - jsMath.Click.y;
    x = Math.max(event.X,Math.min(event.W+event.X-jsMath.Click.w,x));
    y = Math.max(event.Y,Math.min(event.H+event.Y-jsMath.Click.h,y));
    jsMath.Click.source.style.left = x + 'px';
    jsMath.Click.source.style.top  = y + 'px';
    jsMath.Click.left = x + event.X; jsMath.Click.top = y + event.Y;
    return false;
  }

});

jsMath.Click.Init();
