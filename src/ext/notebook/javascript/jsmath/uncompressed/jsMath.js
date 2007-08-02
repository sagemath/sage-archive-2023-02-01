/*****************************************************************************
 *
 *  jsMath: Mathematics on the Web
 *
 *  This jsMath package makes it possible to display mathematics in HTML pages
 *  that are viewable by a wide range of browsers on both the Mac and the IBM PC,
 *  including browsers that don't process MathML.  See
 *
 *            http://www.math.union.edu/locate/jsMath
 *
 *  for the latest version, and for documentation on how to use jsMath.
 *
 *  Copyright 2004-2007 by Davide P. Cervone
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
 *
 *****************************************************************************/

/*
 *  Prevent running everything again if this file is loaded twice
 */
if (!window.jsMath || !window.jsMath.loaded) {

var jsMath_old = window.jsMath;  // save user customizations

//
// debugging routine
//
/*
 * function ShowObject (obj,spaces) {
 *   var s = ''; if (!spaces) {spaces = ""}
 *   for (var i in obj) {
 *     if (obj[i] != null) {
 *       if (typeof(obj[i]) == "object") {
 *         s += spaces + i + ": {\n"
 *           + ShowObject(obj[i],spaces + '  ')
 *           + spaces + "}\n";
 *       } else if (typeof(obj[i]) != "function") {
 *         s += spaces + i + ': ' + obj[i] + "\n";
 *       }
 *     }
 *   }
 *   return s;
 * }
 */

/***************************************************************************/
//
//  Check for DOM support
//
if (!document.getElementById || !document.childNodes || !document.createElement) {
  alert('The mathematics on this page requires W3C DOM support in its JavaScript. '
      + 'Unfortunately, your browser doesn\'t seem to have this.');
} else {

/***************************************************************************/

window.jsMath = {

  version: "3.4c",  // change this if you edit the file, but don't edit this file

  document: document,  // the document loading jsMath
  window: window,      // the window of the of loading document

  platform: (navigator.platform.match(/Mac/) ? "mac" :
             navigator.platform.match(/Win/) ? "pc" : "unix"),

  // Font sizes for \tiny, \small, etc. (must match styles below)
  sizes: [50, 60, 70, 85, 100, 120, 144, 173, 207, 249],

  //
  //  The styles needed for the TeX fonts
  //
  styles: {
    '.math':              'font-family:serif; font-style:normal; font-weight:normal',

    '.typeset':           'font-family:serif; font-style:normal; font-weight:normal; line-height:normal',
    'div.typeset':        'text-align:center; margin:1em 0px;',
    'span.typeset':       'text-align:left',
    '.typeset span':      'text-align:left; border:0px; margin:0px; padding:0px',

    '.typeset .normal':   'font-family:serif; font-style:normal; font-weight:normal',

    '.typeset .size0':    'font-size:50%',  // tiny (\scriptscriptsize)
    '.typeset .size1':    'font-size:60%',  //       (50% of \large for consistency)
    '.typeset .size2':    'font-size:70%',  // scriptsize
    '.typeset .size3':    'font-size:85%',  // small (70% of \large for consistency)
    '.typeset .size4':    'font-size:100%', // normalsize
    '.typeset .size5':    'font-size:120%', // large
    '.typeset .size6':    'font-size:144%', // Large
    '.typeset .size7':    'font-size:173%', // LARGE
    '.typeset .size8':    'font-size:207%', // huge
    '.typeset .size9':    'font-size:249%', // Huge

    '.typeset .cmr10':    'font-family:jsMath-cmr10, serif',
    '.typeset .cmbx10':   'font-family:jsMath-cmbx10, jsMath-cmr10',
    '.typeset .cmti10':   'font-family:jsMath-cmti10, jsMath-cmr10',
    '.typeset .cmmi10':   'font-family:jsMath-cmmi10',
    '.typeset .cmsy10':   'font-family:jsMath-cmsy10',
    '.typeset .cmex10':   'font-family:jsMath-cmex10',

    '.typeset .textit':   'font-family:serif; font-style:italic',
    '.typeset .textbf':   'font-family:serif; font-weight:bold',

    '.typeset .link':     'text-decoration:none',
    '.typeset .error':    'font-size:10pt; font-style:italic; '
                             + 'background-color:#FFFFCC; padding:1px; '
                             + 'border:1px solid #CC0000',

    '.typeset .blank':    'display:inline-block; overflow:hidden; border:0px none; width:0px; height:0px;',
    '.typeset .spacer':   'display:inline-block',

    '#jsMath_hiddenSpan':      'visibility:hidden; position:absolute; top:0px;left:0px; '
                                  + 'line-height:normal; text-indent:0px',

    '#jsMath_message':         'position:fixed; bottom:1px; left:2px; background-color:#E6E6E6; '
                                 + 'border:solid 1px #959595; margin:0px; padding:1px 8px; '
                                 + 'z-index:102; color: black; font-size:small; width:auto;',
    '#jsMath_panel':           'position:fixed; bottom:1.5em; right:1.5em; padding:.8em 1.6em; '
                                 + 'background-color:#DDDDDD; border:outset 2px; '
                                 + 'z-index:103; width:auto; color:black; font-size:10pt; font-style:normal',
    '#jsMath_panel .disabled': 'color:#888888',
    '#jsMath_panel .infoLink': 'font-size:85%',
    '#jsMath_panel *':         'font-size:inherit; font-style:inherit; font-family:inherit; line-height:normal',
    '#jsMath_panel div':       'background-color:inherit; color:inherit;',
    '#jsMath_panel span':      'background-color:inherit; color:inherit;',
    '#jsMath_panel td':        'border:0px; padding:0px; margin:0px; background-color:inherit; color:inherit;',
    '#jsMath_panel tr':        'border:0px; padding:0px; margin:0px; background-color:inherit; color:inherit;',
    '#jsMath_panel table':     'border:0px; padding:0px; margin:0px; background-color:inherit; color:inherit;',

    '#jsMath_button':          'position:fixed; bottom:1px; right:2px; background-color:white; '
                                 + 'border:solid 1px #959595; margin:0px; padding:0px 3px 1px 3px; '
                                 + 'z-index:102; color:black; text-decoration:none; font-size:x-small; '
                                 + 'width:auto; cursor:hand;',
    '#jsMath_button *':        'padding:0px; border:0px; margin:0px; line-height:normal; '
                                 + 'font-size:inherit; font-style:inherit; font-family:inherit',

    '#jsMath_global':          'font-style:italic;',
    '#jsMath_float':           'position:absolute; top:0px; left:0px; max-width:80%; '
                                 + 'z-index:101; width:auto; height:auto;',
    '#jsMath_float .drag':     'background-color:#DDDDDD; border:outset 1px; height:12px; font-size:1px;',
    '#jsMath_float .close':    'background-color:#E6E6E6; border:inset 1px; width:8px; height:8px; margin:1px 2px;',
    '#jsMath_float .source':   'background-color:#E2E2E2; border:outset 1px; '
                                 + 'width:auto; height:auto; padding:8px 15px; '
                                 + 'font-family:courier, fixed; font-size:90%',

    '#jsMath_noFont .message': 'text-align: center; padding:.8em 1.6em; border:3px solid #DD0000; '
                                  + 'background-color:#FFF8F8; color: #AA0000; font-size:small; width:auto;',
    '#jsMath_noFont .link':    'padding:0px 5px 2px 5px; border:2px outset; background-color:#E8E8E8; '
                                  + 'color:black; font-size:80%; width:auto; cursor:hand;',

    '#jsMath_PrintWarning .message':
                                 'text-align:center; padding:.8em 1.6em; border:3px solid #DD0000; '
                                   + 'background-color: #FFF8F8; color: #AA0000; font-size:x-small; width:auto;',

    '@media print':      '#jsMath_button {display:none}\n' +
                         '#jsMath_Warning {display:none}',

    '@media screen':     '#jsMath_PrintWarning {display:none}'

  },


  /***************************************************************************/

  /*
   *  Get a jsMath DOM element
   */
  Element: function (name) {return jsMath.document.getElementById('jsMath_'+name)},

  /*
   *  Get the width and height (in pixels) of an HTML string
   */
  BBoxFor: function (s) {
    this.hidden.innerHTML =
      '<nobr><span class="typeset"><span class="scale">'+s+'</span></span></nobr>';
    var bbox = {w: this.hidden.offsetWidth, h: this.hidden.offsetHeight};
    this.hidden.innerHTML = '';
    return bbox;
  },

  /*
   *  Get the width and height (in ems) of an HTML string.
   *  Check the cache first to see if we've already measured it.
   */
  EmBoxFor: function (s) {
    var cache = jsMath.Global.cache.R;
    if (!cache[this.em]) {cache[this.em] = {}}
    if (!cache[this.em][s]) {
      var bbox = this.BBoxFor(s);
      cache[this.em][s] = {w: bbox.w/this.em, h: bbox.h/this.em};
    }
    return cache[this.em][s];
  },

  /*
   *  For browsers that don't handle sizes of italics properly (MSIE).
   *  Check the cache first to see if we've already measured it.
   */
  EmBoxForItalics: function (s) {
    var cache = jsMath.Global.cache.R;
    if (!cache[this.em]) {cache[this.em] = {}}
    if (!cache[this.em][s]) {
      var bbox = this.BBoxFor(s);
      if (s.match(/<i>|class=\"(icm|italic|igreek|iaccent)/i)) {
        bbox.w = this.BBoxFor(s+jsMath.Browser.italicString).w
                  - jsMath.Browser.italicCorrection;
      }
      cache[this.em][s] = {w: bbox.w/this.em, h: bbox.h/this.em};
    }
    return cache[this.em][s];
  },

  /*
   *  Initialize jsMath.  This determines the em size, and a variety
   *  of other parameters used throughout jsMath.
   */
  Init: function () {
    if (jsMath.Setup.inited != 1) {
      if (!jsMath.Setup.inited) {jsMath.Setup.Body()}
      if (jsMath.Setup.inited != 1) {
        if (jsMath.Setup.inited == -100) return;
        alert("It looks like jsMath failed to set up properly (error code "
               + jsMath.Setup.inited + ").  "
               + "I will try to keep going, but it could get ugly.");
        jsMath.Setup.inited = 1;
      }
    }
    this.em = this.BBoxFor('<span style="'+jsMath.Browser.block+';width:13em;height:1em"></span>').w/13;
    if (this.em == 0) {
      // handle older browsers
      this.em = this.BBoxFor('<img src="'+jsMath.blank+'" style="width:13em;height:1em"/>').w/13;
    }
    var cache = jsMath.Global.cache.B;
    if (!cache[this.em]) {
      cache[this.em] = {};
      cache[this.em].bb = this.BBoxFor('x'); var hh = cache[this.em].bb.h;
      cache[this.em].d = this.BBoxFor('x'+jsMath.HTML.Rule(1,hh/jsMath.em)).h - hh;
      if (jsMath.Browser.italicString)
        {cache[this.em].ic = jsMath.BBoxFor(jsMath.Browser.italicString).w}
    }
    jsMath.Browser.italicCorrection = cache[this.em].ic;
    var bb = cache[this.em].bb; var h = bb.h; var d = cache[this.em].d
    this.h = (h-d)/this.em; this.d = d/this.em;
    this.hd = this.h + this.d;
    this.xWidth = bb.w;  // used to tell if scale has changed

    this.Setup.TeXfonts();

    var x_height = this.EmBoxFor('<span class="cmr10">M</span>').w/2;
    this.TeX.M_height = x_height*(26/14);
    this.TeX.h = this.h; this.TeX.d = this.d; this.TeX.hd = this.hd;

    this.Img.Scale();
    if (!this.initialized) {
      this.Setup.Sizes();
      this.Img.UpdateFonts();
    }

    // factor for \big and its brethren
    this.p_height = (this.TeX.cmex10[0].h + this.TeX.cmex10[0].d) / .85;

    this.initialized = 1;
  },

  /*
   *  Get the xWidth size and if it has changed, reinitialize the sizes
   */
  ReInit: function () {
    var w = this.BBoxFor('x').w;
    if (w != this.xWidth) {this.Init()}
  },

  /*
   *  Mark jsMath as loaded and copy any user-provided overrides
   */
  Loaded: function () {
    if (jsMath_old) {
      var override = ['Process', 'ProcessBeforeShowing',
        'ConvertTeX','ConvertTeX2','ConvertLaTeX','ConvertCustom',
        'CustomSearch', 'Synchronize', 'Macro', 'document'];
      for (var i = 0; i < override.length; i++) {
        if (jsMath_old[override[i]]) {delete jsMath_old[override[i]]}
      }
    }
    if (jsMath_old) {this.Insert(jsMath,jsMath_old)}
    jsMath_old = null;
    jsMath.loaded = 1;
  },

  /*
   *  Manage JavaScript objects:
   *
   *      Add:        add/replace items in an object
   *      Insert:     add items to an object
   *      Package:    add items to an object prototype
   */
  Add: function (dst,src) {for (var id in src) {dst[id] = src[id]}},
  Insert: function (dst,src) {
    for (var id in src) {
      if (dst[id] && typeof(src[id]) == 'object'
                  && (typeof(dst[id]) == 'object'
                  ||  typeof(dst[id]) == 'function')) {
        this.Insert(dst[id],src[id]);
      } else {
        dst[id] = src[id];
      }
    }
  },
  Package: function (obj,def) {this.Insert(obj.prototype,def)}

};


/***************************************************************************/

  /*
   *  Implements items associated with the global cache.
   *
   *  This object will be replaced by a global version when
   *  (and if) jsMath-global.html is loaded.
   */
jsMath.Global = {
    isLocal: 1,  // a local copy if jsMath-global.html hasn't been loaded
    cache: {T: {}, D: {}, R: {}, B: {}},

    /*
     *  Clear the global (or local) cache
     */
    ClearCache: function () {jsMath.Global.cache = {T: {}, D: {}, R: {}, B: {}}},

    /*
     *  Initiate global mode
     */
    GoGlobal: function (cookie) {
      var url = String(jsMath.window.location);
      var c = (jsMath.isCHMmode ? '#' : '?');
      if (cookie) {url = url.replace(/\?.*/,'') + '?' + cookie}
      jsMath.Controls.Reload(jsMath.root + "jsMath-global.html" + c +escape(url));
    },

    /*
     *  Check if we need to go to global mode
     */
    Init: function () {
      if (jsMath.Controls.cookie.global == "always" && !jsMath.noGoGlobal) {
        if (navigator.accentColorName) return; // OmniWeb crashes on GoGlobal
        if (!jsMath.window) {jsMath.window = window}
        jsMath.Controls.loaded = 1;
        jsMath.Controls.defaults.hiddenGlobal = null;
        this.GoGlobal(jsMath.Controls.SetCookie(2));
      }
    },

    /*
     *  Try to register with a global.html window that contains us
     */
    Register: function () {
      var parent = jsMath.window.parent;
      if (!jsMath.isCHMode)
        {jsMath.isCHMmode = (jsMath.window.location.protocol == 'mk:')}
      try {
        if (!jsMath.isCHMmode) this.Domain();
        if (parent.jsMath && parent.jsMath.isGlobal)
          {parent.jsMath.Register(jsMath.window)}
      } catch (err) {jsMath.noGoGlobal = 1}
    },

    /*
     *  If we're not the parent window, try to set the domain to
     *  match the parent's domain (so we can use the Global data
     *  if the surrounding frame is a Global frame.
     */
    Domain: function () {
      // MSIE/Mac can't do domain changes, so don't bother trying
      if (navigator.appName == 'Microsoft Internet Explorer' &&
          jsMath.platform == 'mac' && navigator.userProfile != null) return;
      if (window == parent) return;
      var oldDomain = jsMath.document.domain;
      try {
        while (true) {
          try {if (parent.document.title != null) return} catch (err) {}
          if (!document.domain.match(/\..*\./)) break;
          jsMath.document.domain = jsMath.document.domain.replace(/^[^.]*\./,'');
        }
      } catch (err) {}
      jsMath.document.domain = oldDomain;
    }

};



/***************************************************************************/

/*
 *
 *  Implement loading of remote scripts using XMLHttpRequest, if
 *  possible, otherwise use a hidden IFRAME and fake it.  That
 *  method runs asynchronously, which causes lots of headaches.
 *  Solve these using Push command, which queues actions
 *  until files have loaded.
 */

jsMath.Script = {

  request: null, // the XMLHttpRequest object

  /*
   *  Create the XMLHttpRequest object, if we can.
   *  Otherwise, use the iframe-based fallback method.
   */
  Init: function () {
    if (!(jsMath.Controls.cookie.asynch && jsMath.Controls.cookie.progress)) {
      if (window.XMLHttpRequest) {try {this.request = new XMLHttpRequest} catch (err) {}}
      if (!this.request && window.ActiveXObject) {
        var xml = ["MSXML2.XMLHTTP.5.0","MSXML2.XMLHTTP.4.0","MSXML2.XMLHTTP.3.0",
                   "MSXML2.XMLHTTP","Microsoft.XMLHTTP"];
        for (var i = 0; i < xml.length && !this.request; i++) {
          try {this.request = new ActiveXObject(xml[i])} catch (err) {}
        }
      }
    }
    //
    //  Use the delayed-script fallback for MSIE/Mac and old versions
    //  of several browsers (Opera 7.5, OmniWeb 4.5).
    //
    if (!this.request || jsMath.Setup.domainChanged)
      {this.Load = this.delayedLoad; this.needsBody = 1}
  },

  /*
   *  Load a script and evaluate it in the window's context
   */
  Load: function (url,show) {
    if (show) {
      jsMath.Message.Set("Loading "+url);
      jsMath.Script.Delay(1);
      jsMath.Script.Push(this,'xmlRequest',url);
      jsMath.Script.Push(jsMath.Message,'Clear');
    } else {
      jsMath.Script.Push(this,'xmlRequest',url);
    }
  },

  /*
   *  Load a URL and run the contents of the file
   */
  xmlRequest: function (url) {
    this.blocking = 1;
//    this.debug('xmlRequest: '+url);
    try {
      this.request.open("GET",url,false);
      this.request.send(null);
    } catch (err) {
      this.blocking = 0;
      if (jsMath.Translate.restart && jsMath.Translate.asynchronous) {return ""}
      throw "jsMath can't load the file '"+url+"'\n"
          + "Message: "+err.message;
    }
    if (this.request.status && this.request.status >= 400) {
      // Do we need to deal with redirected links?
      this.blocking = 0;
      if (jsMath.Translate.restart && jsMath.Translate.asynchronous) {return ""}
      throw "jsMath can't load the file '"+url+"'\n"
          + "Error status: "+this.request.status;
    }
    if (!url.match(/\.js$/)) {return(this.request.responseText)}
    var tmpQueue = this.queue; this.queue = [];
//    this.debug('xml Eval ['+tmpQueue.length+']');
    jsMath.window.eval(this.request.responseText);
//    this.debug('xml Done ['+this.queue.length+' + '+tmpQueue.length+']');
    this.blocking = 0; this.queue = this.queue.concat(tmpQueue);
    this.Process();
    return "";
  },

  /********************************************************************
   *
   *  Implement asynchronous loading and execution of scripts
   *  (via hidden IFRAME) interleved with other JavaScript commands
   *  that must be synchronized with the file loading.  (Basically, this
   *  is for MSIE/Mac and Opera 7.5, which don't have XMLHttpRequest.)
   */

  cancelTimeout: 30*1000,   // delay for canceling load (30 sec)

  iframe: null,      // the hidden iframe
  blocking: 0,       // true when an asynchronous action is being performed
  cancelTimer: null, // timer to cancel load if it takes too long
  needsBody: 0,      // true if loading files requires BODY to be present

  queue: [],         // the stack of pending actions

  /*
   *  Provide mechanism for synchronizing with the asynchronous jsMath
   *  file-loading mechanism.  'code' can be a string or a function.
   */
  Synchronize: function (code,data) {
    if (typeof(code) != 'string') {jsMath.Script.Push(null,code,data)}
      else {jsMath.Script.Push(jsMath.window,'eval',code)}
  },

  /*
   *  Queue a function to be processed.
   *  If nothing is being loaded, do the pending commands.
   */
  Push: function (object,method,data) {
//    this.debug('Pushing: '+method+' at '+this.queue.length);
    this.queue[this.queue.length] = [object,method,data];
    if (!(this.blocking || (this.needsBody && !jsMath.document.body))) this.Process();
  },

  /*
   *  Do any pending functions (stopping if a file load is started)
   */
  Process: function () {
    while (this.queue.length && !this.blocking) {
      var call = this.queue[0]; var tmpQueue = this.queue.slice(1); this.queue = [];
      var object = call[0]; var method = call[1]; var data = call[2];
//      this.debug('Calling: '+method+' ['+tmpQueue.length+']');
      if (object) {object[method](data)} else if (method) {method(data)}
//      this.debug('Done:    '+method+' ['+this.queue.length+' + '+tmpQueue.length+']');
      this.queue = this.queue.concat(tmpQueue);
    }
  },

  /*
   *  Handle loading of scripts that run asynchronously
   */
  delayedLoad: function (url) {
//    this.debug('Loading: '+url);
    this.Push(this,'startLoad',url);
  },
  startLoad: function (url) {
    this.iframe = jsMath.document.createElement('iframe');
    this.iframe.style.visibility = 'hidden';
    this.iframe.style.position = 'absolute';
    this.iframe.style.width = '0px';
    this.iframe.style.height = '0px';
    if (jsMath.document.body.firstChild) {
      jsMath.document.body.insertBefore(this.iframe,jsMath.document.body.firstChild);
    } else {
      jsMath.document.body.appendChild(this.iframe);
    }
    this.blocking = 1; this.url = url;
    if (!url.match(/\.js$/)) {this.iframe.src = url}
                        else {this.iframe.src = jsMath.root+"jsMath-loader.html"}
    if (url.substr(0,jsMath.root.length) == jsMath.root)
      {url = url.substr(jsMath.root.length)}
    jsMath.Message.Set("Loading "+url);
    this.cancelTimer = setTimeout('jsMath.Script.cancelLoad()',this.cancelTimeout);
  },
  endLoad: function (action) {
    if (this.cancelTimer) {clearTimeout(this.cancelTimer); this.cancelTimer = null;}
    jsMath.Message.Clear();
    if (action != 'cancel') {this.blocking = 0; this.Process()}
  },

  Start: function () {
//    this.debug('Starting: ['+this.queue.length+'] '+this.url);
    this.tmpQueue = this.queue; this.queue = [];
  },
  End:   function () {
//    this.debug('Ending:   ['+this.queue.length+' + '+this.tmpQueue.length+'] '+this.url);
    this.queue = this.queue.concat(this.tmpQueue); delete this.tmpQueue;
  },

  /*
   *  If the loading takes too long, cancel it and end the load.
   */
  cancelLoad: function () {
    this.cancelTimer = null;
    jsMath.Message.Set("Can't load file");
    this.endLoad("cancel");
  },

  /*
   *  Perform a delay (to let the browser catch up)
   */
  Delay: function (time) {
    this.blocking = 1;
    setTimeout('jsMath.Script.endDelay()',time);
  },
  endDelay: function () {
//    this.debug('endDelay');
    this.blocking = 0;
    this.Process();
  },

  /*
   *  Load an image and wait for it
   *  (so MSIE won't load extra copies of it)
   */
  imageCount: 0,
  WaitForImage: function (file) {
    this.blocking = 1; this.imageCount++;
    if (this.img == null) {this.img = []}
    var img = new Image(); this.img[this.img.length] = img;
    img.onload = function () {if (--jsMath.Script.imageCount == 0) jsMath.Script.endDelay()}
    img.onerror = img.onload; img.onabort = img.onload;
    img.src = file;
  },

  /*
   *  The code uncompressor
   */
  Uncompress: function (data) {
    for (var k = 0; k <  data.length; k++) {
      var d = data[k]; var n = d.length;
      for (var i = 0; i < n; i++) {if (typeof(d[i]) == 'number') {d[i] = d[d[i]]}}
      data[k] = d.join('');
    }
    window.eval(data.join(''));
  }//,

  /*
   *  for debugging the event queue
   */
//  debug: function (message) {
//    if (jsMath.document.body && jsMath.window.debug) {jsMath.window.debug(message)}
//      else {alert(message)}
//  }

};

/***************************************************************************/

/*
 *  Message and screen blanking facility
 */

jsMath.Message = {

  blank: null,    // the div to blank out the screen
  message: null,  // the div for the messages
  text: null,     // the text node for messages
  clear: null,    // timer for clearing message

  /*
   *  Create the elements needed for the message box
   */
  Init: function () {
    if (!jsMath.document.body || !jsMath.Controls.cookie.progress) return;
    if (jsMath.Setup.stylesReady) {
      this.message = jsMath.Setup.DIV('message',{visibility:'hidden'});
    } else {
      this.message = jsMath.Setup.DIV('message',{
        visibility:'hidden', position:'absolute', bottom:'1px', left:'2px',
        backgroundColor:'#E6E6E6', border:'solid 1px #959595',
        margin:'0px', padding:'1px 8px', zIndex:102,
        color:'black', fontSize:'small', width:'auto'
      });
    }
    this.text = jsMath.document.createTextNode('');
    this.message.appendChild(this.text);
    this.message.onmousedown = jsMath.Translate.Cancel;
  },

  /*
   *  Set the contents of the message box, or use the window status line
   */
  Set: function (text,canCancel) {
    if (this.clear) {clearTimeout(this.clear); this.clear = null}
    if (jsMath.Controls.cookie.progress) {
      if (!this.text) {this.Init(); if (!this.text) return}
      if (jsMath.Browser.textNodeBug) {this.message.innerHTML = text}
        else {this.text.nodeValue = text}
      this.message.style.visibility = 'visible';
      if (canCancel) {
        this.message.style.cursor = 'pointer';
        if (!this.message.style.cursor) {this.message.style.cursor = 'hand'}
        this.message.title = ' Cancel Processing of Math ';
      } else {
        this.message.style.cursor = '';
        this.message.title = '';
      }
    } else {
      if (text.substr(0,8) != "Loading ") {jsMath.window.status = text}
    }
  },

  /*
   *  Clear the message box or status line
   */
  Clear: function () {
    if (this.clear) {clearTimeout(this.clear)}
    this.clear = setTimeout("jsMath.Message.doClear()",1000);
  },
  doClear: function () {
    if (this.clear) {
      this.clear = null;
      jsMath.window.status = '';
      if (this.text) {this.text.nodeValue = ''}
      if (this.message) {this.message.style.visibility = 'hidden'}
    }
  },


  /*
   *  Put up a DIV that covers the window so that the
   *  "flicker" of processing the mathematics will not be visible
   */
  Blank: function () {
    if (this.blank || !jsMath.document.body) return;
    this.blank = jsMath.Setup.DIV("blank",{
      position:(jsMath.Browser.msiePositionFixedBug? 'absolute': 'fixed'),
      top:'0px', left:'0px', bottom:'0px', right:'0px',
      zIndex:101, backgroundColor:'white'
    });
    if (jsMath.Browser.msieBlankBug) {
      this.blank.innerHTML = '&nbsp;';
      this.blank.style.width = "110%";
      this.blank.style.height = "110%";
    }
  },

  UnBlank: function () {
    if (this.blank) {jsMath.document.body.removeChild(this.blank)}
    this.blank = null;
  }
};


/***************************************************************************/

/*
 *  Miscellaneous setup and initialization
 */
jsMath.Setup = {

  loaded: [],  // array of files already loaded

  /*
   *  Insert a DIV at the top of the page with given ID,
   *  attributes, and style settings
   */
  DIV: function (id,styles) {
    var div = jsMath.document.createElement('div');
    div.id = 'jsMath_'+id;
    for (var i in styles) {div.style[i]= styles[i]}
    if (!jsMath.document.body.hasChildNodes) {jsMath.document.body.appendChild(div)}
      else {jsMath.document.body.insertBefore(div,jsMath.document.body.firstChild)}
    return div;
  },

  /*
   *  Source a jsMath JavaScript file (only load any given file once)
   */
  Script: function (file,show) {
    if (this.loaded[file]) {return} else {this.loaded[file] = 1}
    if (!file.match('^([a-zA-Z]+:/?)?/')) {file = jsMath.root + file}
    jsMath.Script.Load(file,show);
  },

  /*
   *  Use a hidden <DIV> for measuring the BBoxes of things
   */
  Hidden: function () {
    jsMath.hidden = this.DIV("Hidden",{
      visibility: 'hidden', position:"absolute",
      top:0, left:0, border:0, padding:0, margin:0
    });
    jsMath.hiddenTop = jsMath.hidden;
    return;
  },

  /*
   *  Find the root URL for the jsMath files (so we can load
   *  the other .js and .gif files)
   */
  Source: function () {
    if (jsMath.Autoload && jsMath.Autoload.root) {
      jsMath.root = jsMath.Autoload.root;
    } else {
      jsMath.root = '';
      var script = jsMath.document.getElementsByTagName('script');
      if (script) {
        for (var i = 0; i < script.length; i++) {
          var src = script[i].src;
          if (src && src.match('(^|/)jsMath.js$')) {
            jsMath.root = src.replace(/jsMath.js$/,'');
            i = script.length;
          }
        }
      }
    }
    if (jsMath.root.charAt(0) == '/') {
      jsMath.root = jsMath.document.location.protocol + '//'
                  + jsMath.document.location.host + jsMath.root;
    } else if (!jsMath.root.match(/^[a-z]+:/i)) {
      src = new String(jsMath.document.location);
      jsMath.root = src.replace(new RegExp('[^/]*$'),'') + jsMath.root;
      while (jsMath.root.match('/[^/]*/\\.\\./')) {
        jsMath.root = jsMath.root.replace(new RegExp('/[^/]*/\\.\\./'),'/');
      }
    }
    jsMath.Img.root = jsMath.root + "fonts/";
    jsMath.blank = jsMath.root + "blank.gif";
    this.Domain();
  },

  /*
   *  Find the most restricted common domain for the main
   *  page and jsMath.  Report an error if jsMath is outside
   *  the domain of the calling page.
   */
  Domain: function () {
    try {jsMath.document.domain} catch (err) {return}
    var jsDomain = ''; var pageDomain = jsMath.document.domain;
    if (jsMath.root.match('://([^/]*)/')) {jsDomain = RegExp.$1}
    jsDomain = jsDomain.replace(/:\d+$/,'');
    if (jsDomain == "" || jsDomain == pageDomain) return;
    //
    // MSIE on the Mac can't change jsMath.document.domain and 'try' won't
    //   catch the error (Grrr!), so exit for them.
    //
    if (navigator.appName == 'Microsoft Internet Explorer' &&
        jsMath.platform == 'mac' && navigator.onLine &&
        navigator.userProfile && jsMath.document.all) return;
    jsDomain = jsDomain.split(/\./); pageDomain = pageDomain.split(/\./);
    if (jsDomain.length < 2 || pageDomain.length < 2 ||
        jsDomain[jsDomain.length-1] != pageDomain[pageDomain.length-1] ||
        jsDomain[jsDomain.length-2] != pageDomain[pageDomain.length-2]) {
      this.DomainWarning();
      return;
    }
    var domain = jsDomain[jsDomain.length-2] + '.' + jsDomain[jsDomain.length-1];
    for (var i = 3; i <= jsDomain.length && i <= pageDomain.length; i++) {
      if (jsDomain[jsDomain.length-i] != pageDomain[pageDomain.length-i]) break;
      domain = jsDomain[jsDomain.length-i] + '.' + domain;
    }
    jsMath.document.domain = domain;
    this.domainChanged = 1;
  },

  DomainWarning: function () {
    alert("In order for jsMath to be able to load the additional "
        + "components that it may need, the jsMath.js file must be "
        + "loaded from a server in the same domain as the page that "
        + "contains it.  Because that is not the case for this page, "
        + "the mathematics displayed here may not appear correctly.");
  },

  /*
   *  Initialize a font's encoding array
   */
  EncodeFont: function (name) {
    var font = jsMath.TeX[name];
    if (font[0].c != null) return;
    for (var k = 0; k < 128; k++) {
      var data = font[k]; font[k] = data[3];
      if (font[k] == null) {font[k] = {}};
      font[k].w = data[0]; font[k].h = data[1];
      if (data[2] != null) {font[k].d = data[2]}
      font[k].c = jsMath.TeX.encoding[k];
    }
  },

  /*
   *  Initialize the encodings for all fonts
   */
  Fonts: function () {
    for (var i = 0; i < jsMath.TeX.fam.length; i++) {
      var name = jsMath.TeX.fam[i];
      if (name) {this.EncodeFont(name)}
    }
  },

  /*
   *  Look up the default height and depth for a TeX font
   *  and set the skewchar
   */
  TeXfont: function (name) {
    var font = jsMath.TeX[name]; if (font == null) return;
    var WH = jsMath.EmBoxFor('<span class="'+name+'">'+font[65].c+'</span>');
    font.hd = WH.h; font.dh = .05;
    font.d = jsMath.EmBoxFor('<span class="'+name+'">'+ font[65].c +
      jsMath.HTML.Rule(1,font.hd) + '</span>').h - font.hd;
    font.h = font.hd - font.d;
    if (name == 'cmmi10') {font.skewchar = 0177}
    else if (name == 'cmsy10') {font.skewchar = 060}
  },

  /*
   *  Init all the TeX fonts
   */
  TeXfonts: function () {
    for (var i = 0; i < jsMath.TeX.fam.length; i++)
      {if (jsMath.TeX.fam[i]) {this.TeXfont(jsMath.TeX.fam[i])}}
  },

  /*
   *  Compute font parameters for various sizes
   */
  Sizes: function () {
    jsMath.TeXparams = [];
    for (var j=0; j < jsMath.sizes.length; j++) {jsMath.TeXparams[j] = {}}
    for (var i in jsMath.TeX) {
      if (typeof(jsMath.TeX[i]) != 'object') {
        for (var j=0; j < jsMath.sizes.length; j++) {
          jsMath.TeXparams[j][i] = jsMath.sizes[j]*jsMath.TeX[i]/100;
        }
      }
    }
  },

  /*
   *  Send the style definitions to the browser (these may be adjusted
   *  by the browser-specific code)
   */
  Styles: function (styles) {
    if (!styles) {
      styles = jsMath.styles;
      styles['.typeset .scale'] = 'font-size:'+jsMath.Controls.cookie.scale+'%';
      this.stylesReady = 1;
    }
    jsMath.Script.Push(this,'AddStyleSheet',styles);
    if (jsMath.Browser.styleChangeDelay) {jsMath.Script.Push(jsMath.Script,'Delay',1)}
  },

  AddStyleSheet: function (styles) {
    var head = jsMath.document.getElementsByTagName('head')[0];
    var string = ''; for (var id in styles) {string += id + ' {'+styles[id]+"}\n"}
    if (jsMath.document.createStyleSheet) {// check for MSIE
      head.insertAdjacentHTML('beforeEnd',
          '<span style="display:none">x</span>'  // MSIE needs this for some reason
          + '<style type="text/css">'+string+'</style>');
    } else {
      var style = jsMath.document.createElement('style'); style.type = "text/css";
      style.appendChild(jsMath.document.createTextNode(string));
      head.appendChild(style);
    }
  },

  /*
   *  Do the initialization that requires the <body> to be in place.
   */
  Body: function () {
    if (this.inited) return;

    this.inited = -1;

    jsMath.Setup.Hidden(); this.inited = -2;
    jsMath.Browser.Init(); this.inited = -3;

    // blank screen if necessary
    if (jsMath.Controls.cookie.blank) {jsMath.Message.Blank()}; this.inited = -4;

    jsMath.Setup.Styles(); this.inited = -5;
    jsMath.Controls.Init(); this.inited = -6;

    // do user-specific initialization
    jsMath.Script.Push(jsMath.Setup,'User','pre-font'); this.inited = -7;

    // make sure browser-specific loads are done before this
    jsMath.Script.Push(jsMath.Font,'Check');
    if (jsMath.Font.register.length)
      {jsMath.Script.Push(jsMath.Font,'LoadRegistered')}

    this.inited = 1;
  },

  /*
   *  Web page author can override the entries to the UserEvent hash
   *  functions that will be run at various times during jsMath's setup
   *  process.
   */
  User: function (when) {
    if (jsMath.Setup.UserEvent[when]) {(jsMath.Setup.UserEvent[when])()}
  },

  UserEvent: {
    "pre-font": null,  // after browser is set up but before fonts are tested
    "onload": null     // after jsMath.js is loaded and finished running
  }

};

jsMath.Update = {

  /*
   *  Update specific parameters for a limited number of font entries
   */
  TeXfonts: function (change) {
    for (var font in change) {
      for (var code in change[font]) {
        for (var id in change[font][code]) {
          jsMath.TeX[font][code][id] = change[font][code][id];
        }
      }
    }
  },

  /*
   *  Update the character code for every character in a list
   *  of fonts
   */
  TeXfontCodes: function (change) {
    for (var font in change) {
      for (var i = 0; i < change[font].length; i++) {
        jsMath.TeX[font][i].c = change[font][i];
      }
    }
  }

};

/***************************************************************************/

/*
 *  Implement browser-specific checks
 */

jsMath.Browser = {

  allowAbsolute: 1,           // tells if browser can nest absolutely positioned
                              //   SPANs inside relative SPANs
  allowAbsoluteDelim: 0,      // OK to use absolute placement for building delims?
  separateSkips: 0,           // MSIE doesn't do negative left margins, and
                              //   Netscape doesn't combine skips well

  valignBug: 0,               // Konqueror doesn't nest vertical-align
  operaHiddenFix: '',         // for Opera to fix bug with math in tables
  msieCenterBugFix: '',       // for MSIE centering bug with image fonts
  msieInlineBlockFix: '',     // for MSIE alignment bug in non-quirks mode
  msieSpaceFix: '',           // for MSIE to avoid dropping empty spans
  imgScale: 1,                // MSI scales images for 120dpi screens, so compensate

  renameOK: 1,                // tells if brower will find a tag whose name
                              //   has been set via setAttributes
  styleChangeDelay: 0,        // true if style changes need a delay in order
                              //   for them to be available

  delay: 1,                   // delay for asynchronous math processing

  version: 0,                 // browser version number (when needed)

  /*
   *  Determine if the "top" of a <SPAN> is always at the same height
   *  or varies with the height of the rest of the line (MSIE).
   */
  TestSpanHeight: function () {
    jsMath.hidden.innerHTML = '<span><span style="'+this.block+';height:2em;width:1px"></span></span>';
    var span = jsMath.hidden.firstChild;
    var img  = span.firstChild;
    this.spanHeightVaries = (span.offsetHeight >= img.offsetHeight && span.offsetHeight > 0);
    this.spanHeightTooBig = (span.offsetHeight > img.offsetHeight);
    jsMath.hidden.innerHTML = '';
  },

  /*
   *  Determine if an inline-block with 0 width is OK or not
   *  and decide whether to use spans or images for spacing
   */
  TestInlineBlock: function () {
    this.block = "display:-moz-inline-box";
    this.hasInlineBlock = jsMath.BBoxFor('<span style="'+this.block+';width:10px;height:5px"></span>').w > 0;
    if (this.hasInlineBlock) {
      jsMath.styles['.typeset .blank']  = jsMath.styles['.typeset .blank'].replace(/display:inline-block/,this.block);
      jsMath.styles['.typeset .spacer'] = jsMath.styles['.typeset .spacer'].replace(/display:inline-block/,'');
    } else {
      this.block = "display:inline-block";
      this.hasInlineBlock = jsMath.BBoxFor('<span style="'+this.block+';width:10px;height:5px"></span>').w > 0;
      if (!this.hasInlineBlock) return;
    }
    this.block += ';overflow:hidden';
    var h = jsMath.BBoxFor('x').h;
    this.mozInlineBlockBug = jsMath.BBoxFor(
       '<span style="'+this.block+';height:'+h+'px;width:1px"></span>x'+
       '<span style="'+this.block+';height:'+h+'px;width:1px;vertical-align:-'+h+'px"></span>').h > 2*h;
    this.widthAddsBorder = jsMath.BBoxFor('<span style="'+this.block+
        ';overflow:hidden;height:1px;width:10px;border-left:10px solid"></span>').w > 10;
    this.msieBorderBug =
      jsMath.BBoxFor('<span style="'+this.block+';height:'+h+'px;width:1px"></span>x').h !=
      jsMath.BBoxFor('<span style="'+this.block+';height:'+h+'px;width:1px;border-left:1px solid"></span>x').h;
    this.blankWidthBug = this.msieBorderBug ||
      jsMath.BBoxFor('<span style="'+this.block+';height:2em;width:0px"></span>').h == 0;
  },

  /*
   *  Determine if the NAME attribute of a tag can be changed
   *  using the setAttribute function, and then be properly
   *  returned by getElementByName.
   */
  TestRenameOK: function () {
    jsMath.hidden.innerHTML = '<span></span>';
    var test = jsMath.hidden.firstChild;
    test.setAttribute('name','jsMath_test');
    this.renameOK = (jsMath.document.getElementsByName('jsMath_test').length > 0);
    jsMath.hidden.innerHTML = '';
  },

  /*
   *  See if style changes occur immediately, or if we need to delay
   *  in order to let them take effect.
   */
  TestStyleChange: function () {
    jsMath.hidden.innerHTML = '<span ID="jsMath_test">x</span>';
    var span = jsMath.hidden.firstChild;
    var w = span.offsetWidth;
    jsMath.Setup.AddStyleSheet({'#jsMath_test': 'font-size:200%'});
    this.styleChangeDelay = (span.offsetWidth == w);
    jsMath.hidden.innerHTML = '';
  },

  /*
   *  Perform a version check on a standard version string
   */
  VersionAtLeast: function (v) {
    var bv = new String(this.version).split('.');
    v = new String(v).split('.'); if (v[1] == null) {v[1] = '0'}
    return bv[0] > v[0] || (bv[0] == v[0] && bv[1] >= v[1]);
  },

  /*
   *  Test for browser characteristics, and adjust things
   *  to overcome specific browser bugs
   */
  Init: function () {
    jsMath.browser = 'unknown';
    this.TestInlineBlock();
    this.TestSpanHeight();
    this.TestRenameOK();
    this.TestStyleChange();

    this.MSIE();
    this.Mozilla();
    this.Opera();
    this.OmniWeb();
    this.Safari();
    this.Konqueror();

    //
    // Change some routines depending on the browser
    //
    if (this.allowAbsoluteDelim) {
      jsMath.Box.DelimExtend = jsMath.Box.DelimExtendAbsolute;
      jsMath.Box.Layout = jsMath.Box.LayoutAbsolute;
    } else {
      jsMath.Box.DelimExtend = jsMath.Box.DelimExtendRelative;
      jsMath.Box.Layout = jsMath.Box.LayoutRelative;
    }

    if (this.separateSkips) {
      jsMath.HTML.Place = jsMath.HTML.PlaceSeparateSkips;
      jsMath.Typeset.prototype.Place = jsMath.Typeset.prototype.PlaceSeparateSkips;
    }
  },

  //
  //  Handle bug-filled Internet Explorer
  //
  MSIE: function () {
    if (this.spanHeightVaries && !this.spanHeightTooBig) {
      jsMath.browser = 'MSIE';
      if (jsMath.platform == 'pc') {
        this.IE7 = (window.XMLHttpRequest != null);
        this.quirks = (jsMath.document.compatMode == "BackCompat");
        this.allowAbsoluteDelim = 1; this.separateSkips = 1;
        this.buttonCheck = 1; this.msieBlankBug = 1;
        this.msieDivWidthBug = 1; this.msiePositionFixedBug = 1;
        this.msieIntegralBug = 1; this.waitForImages = 1;
        this.msieAlphaBug = !this.IE7; this.alphaPrintBug = !this.IE7;
        this.msieCenterBugFix = 'position:relative; ';
        this.msieInlineBlockFix = ' display:inline-block;';
        if (!this.IE7) {this.msieSpaceFix = '<span style="display:inline-block"></span>'}
        jsMath.Macro('joinrel','\\mathrel{\\kern-5mu}'),
        jsMath.styles['.typeset .arial'] = "font-family: 'Arial unicode MS'";
        if (!this.IE7 || this.quirks) {
          // MSIE doesn't implement fixed positioning, so use absolute
          jsMath.styles['#jsMath_message'] =
              jsMath.styles['#jsMath_message'].replace(/position:fixed/,"position:absolute").replace(/width:auto/,"");
          jsMath.styles['#jsMath_panel'] =
              jsMath.styles['#jsMath_panel'].replace(/position:fixed/,"position:absolute").replace(/width:auto/,"");
          jsMath.styles['#jsMath_button'] = 'width:1px; '
            + jsMath.styles['#jsMath_button'].replace(/position:fixed/,"position:absolute").replace(/width:auto/,"");
          jsMath.window.attachEvent("onscroll",jsMath.Controls.MoveButton);
          if (this.IE7) jsMath.window.attachEvent("onresize",jsMath.Controls.MoveButton);
	  this.msieMoveButtonHack = this.IE7;
	}
        // Make MSIE put borders around the whole button
        jsMath.styles['#jsMath_noFont .link'] += " display: inline-block;";
        // MSIE needs this NOT to be inline-block
        jsMath.styles['.typeset .spacer'] =
              jsMath.styles['.typeset .spacer'].replace(/display:inline-block/,'');
        // MSIE can't insert DIV's into text nodes, so tex2math must use SPAN's to fake DIV's
        jsMath.styles['.tex2math_div'] = jsMath.styles['div.typeset'] + '; width: 100%; display: inline-block';
        // MSIE will rescale images if the DPIs differ
        if (screen.deviceXDPI && screen.logicalXDPI
             && screen.deviceXDPI != screen.logicalXDPI) {
          this.imgScale *= screen.logicalXDPI/screen.deviceXDPI;
          jsMath.Controls.cookie.alpha = 0;
        }
        // Handle bug with getting width of italic text
        this.italicString = '<i>x</i>';
        jsMath.EmBoxFor = jsMath.EmBoxForItalics;
      } else if (jsMath.platform == 'mac') {
        this.msieAbsoluteBug = 1; this.msieButtonBug = 1;
        this.msieDivWidthBug = 1; this.msieBlankBug = 1;
        this.quirks = 1;
        jsMath.Setup.Script('jsMath-msie-mac.js');
        jsMath.Parser.prototype.macros.angle = ['Replace','ord','<font face="Symbol">&#x8B;</font>','normal'];
        jsMath.styles['#jsMath_panel']  = 'width:42em; ' + jsMath.styles['#jsMath_panel'].replace(/width:auto/,"");
        jsMath.Controls.cookie.printwarn = 0; // MSIE/Mac doesn't handle '@media screen'
      }
      jsMath.Macro('not','\\mathrel{\\rlap{\\kern3mu/}}');
    }
  },

  //
  //  Handle Netscape/Mozilla (any flavor)
  //
  Mozilla: function () {
    if (jsMath.hidden.ATTRIBUTE_NODE) {
      jsMath.browser = 'Mozilla';
      if (jsMath.platform == 'pc') {this.alphaPrintBug = 1}
      this.allowAbsoluteDelim = 1;
      jsMath.styles['#jsMath_button'] = jsMath.styles['#jsMath_button'].replace(/cursor:hand/,'cursor:pointer');
      jsMath.styles['#jsMath_noFont .link'] = jsMath.styles['#jsMath_noFont .link'].replace(/cursor:hand/,'cursor:pointer');
      jsMath.Macro('not','\\mathrel{\\rlap{\\kern3mu/}}');
      if (navigator.vendor == 'Firefox') {
        this.version = navigator.vendorSub;
      } else if (navigator.userAgent.match(' Firefox/([0-9.]+)( |$)')) {
        this.version = RegExp.$1;
      }
    }
  },

  //
  //  Handle OmniWeb
  //
  OmniWeb: function () {
    if (navigator.accentColorName) {
      jsMath.browser = 'OmniWeb';
      this.allowAbsolute = this.hasInlineBlock;
      this.allowAbsoluteDelim = this.allowAbsolute;
      this.valignBug = !this.allowAbsolute;
      this.buttonCheck = 1; this.textNodeBug = 1;
      jsMath.noChangeGlobal = 1; // OmniWeb craches on GoGlobal
      if (!this.hasInlineBlock) {jsMath.Setup.Script('jsMath-old-browsers.js')}
    }
  },

  //
  //  Handle Opera
  //
  Opera: function () {
    if (this.spanHeightTooBig) {
      jsMath.browser = 'Opera';
      var isOld = navigator.userAgent.match("Opera 7");
      this.allowAbsolute = 0;
      this.delay = 10;
      this.operaHiddenFix = '[Processing]';
      if (isOld) {jsMath.Setup.Script('jsMath-old-browsers.js')}
    }
  },

  //
  //  Handle Safari
  //
  Safari: function () {
    if (navigator.appVersion.match(/Safari\//)) {
      jsMath.browser = 'Safari';
      var version = navigator.userAgent.match("Safari/([0-9]+)");
      version = (version)? version[1] : 400;
      for (var i = 0; i < jsMath.TeX.fam.length; i++) {
        if (jsMath.TeX.fam[i] && jsMath.TeX[jsMath.TeX.fam[i]])
          {jsMath.TeX[jsMath.TeX.fam[i]].dh = .1}
      }
      jsMath.TeX.axis_height += .05;
      jsMath.TeX.default_rule_thickness += .025;
      this.allowAbsoluteDelim = version >= 125;
      this.safariIFRAMEbug = version >= 312 && version < 412;
      this.safariButtonBug = version < 412;
      this.safariImgBug = 1; this.textNodeBug = 1;
      this.buttonCheck = 1;
      this.styleChangeDelay = 1;
    }
  },

  //
  //  Handle Konqueror
  //
  Konqueror: function () {
    if (navigator.product && navigator.product.match("Konqueror")) {
      jsMath.browser = 'Konqueror';
      this.allowAbsolute = 0;
      this.allowAbsoluteDelim = 0;
      if (navigator.userAgent.match(/Konqueror\/(\d+)\.(\d+)/)) {
        if (RegExp.$1 < 3 || (RegExp.$1 == 3 && RegExp.$2 < 3)) {
          this.separateSkips = 1;
          this.valignBug = 1;
          jsMath.Setup.Script('jsMath-old-browsers.js');
        }
      }
      //  Apparently, Konqueror wants the names without the hyphen
      jsMath.Add(jsMath.styles,{
        '.typeset .cmr10':    'font-family: jsMath-cmr10, jsMath cmr10, serif',
        '.typeset .cmbx10':   'font-family: jsMath-cmbx10, jsMath cmbx10, jsMath-cmr10, jsMath cmr10',
        '.typeset .cmti10':   'font-family: jsMath-cmti10, jsMath cmti10, jsMath-cmr10, jsMath cmr10',
        '.typeset .cmmi10':   'font-family: jsMath-cmmi10, jsMath cmmi10',
        '.typeset .cmsy10':   'font-family: jsMath-cmsy10, jsMath cmsy10',
        '.typeset .cmex10':   'font-family: jsMath-cmex10, jsMath cmex10'
      });
      jsMath.Font.testFont = "jsMath-cmex10, jsMath cmex10";
    }
  }

};

/***************************************************************************/

/*
 *  Implement font check and messages
 */
jsMath.Font = {

  testFont: "jsMath-cmex10",
  fallback: "symbol", // the default fallback method
  register: [],       // list of fonts registered before jsMath.Init()

  // the HTML for the missing font message
  message:
    '<b>No jsMath TeX fonts found</b> -- using image fonts instead.<br/>\n'
      + 'These may be slow and might not print well.<br/>\n'
      + 'Use the jsMath control panel to get additional information.',

  extra_message:
    'Extra TeX fonts not found: <b><span id="jsMath_ExtraFonts"></span></b><br/>'
      + 'Using image fonts instead.  This may be slow and might not print well.<br/>\n'
      + 'Use the jsMath control panel to get additional information.',

  print_message:
    'To print higher-resolution math symbols, click the<br/>\n'
       + '<b>Hi-Res Fonts for Printing</b> button on the jsMath control panel.<br/>\n',

  alpha_message:
    'If the math symbols print as black boxes, turn off <b>image alpha channels</b><br/>\n'
       + 'using the <B>Options</B> pane of the jsMath control panel.<br/>\n',

  /*
   *  Look to see if a font is found.
   *  Check the character in a given position, and see if it is
   *  wider than the usual one in that position.
   */
  Test1: function (name,n,factor,prefix) {
    if (n == null) {n = 0x7C}; if (factor == null) {factor = 2}; if (prefix == null) {prefix = ''}
    var wh1 = jsMath.BBoxFor('<span style="font-family: '+prefix+name+', serif">'+jsMath.TeX[name][n].c+'</span>');
    var wh2 = jsMath.BBoxFor('<span style="font-family: serif">'+jsMath.TeX[name][n].c+'</span>');
    //alert([wh1.w,wh2.w,wh1.h,factor*wh2.w]);
    return (wh1.w > factor*wh2.w && wh1.h != 0);
  },

  Test2: function (name,n,factor,prefix) {
    if (n == null) {n = 0x7C}; if (factor == null) {factor = 2}; if (prefix == null) {prefix = ''}
    var wh1 = jsMath.BBoxFor('<span style="font-family: '+prefix+name+', serif">'+jsMath.TeX[name][n].c+'</span>');
    var wh2 = jsMath.BBoxFor('<span style="font-family: serif">'+jsMath.TeX[name][n].c+'</span>');
    //alert([wh2.w,wh1.w,wh1.h,factor*wh1.w]);
    return (wh2.w > factor*wh1.w && wh1.h != 0);
  },

  /*
   *  Check for the new jsMath versions of the fonts (blacker with
   *  different encoding) and if not found, look for old-style fonts.
   *  If they are found, load the BaKoMa encoding information.
   */
  CheckTeX: function () {
    var wh = jsMath.BBoxFor('<span style="font-family: '+jsMath.Font.testFont+', serif">'+jsMath.TeX.cmex10[1].c+'</span>');
    jsMath.nofonts = ((wh.w*3 > wh.h || wh.h == 0) && !this.Test1('cmr10',null,null,'jsMath-'));
    if (jsMath.nofonts && (jsMath.platform != "mac" ||
        jsMath.browser != 'Mozilla' || !jsMath.Browser.VersionAtLeast(1.5))) {
      wh = jsMath.BBoxFor('<span style="font-family: cmex10, serif">'+jsMath.TeX.cmex10[1].c+'</span>');
      jsMath.nofonts = ((wh.w*3 > wh.h || wh.h == 0) && !this.Test1('cmr10'));
      if (!jsMath.nofonts) {jsMath.Setup.Script("jsMath-BaKoMa-fonts.js")}
    }
  },

  /*
   *  Check for the availability of TeX fonts.  We do this by looking at
   *  the width and height of a character in the cmex10 font.  The cmex10
   *  font has depth considerably greater than most characters' widths (the
   *  whole font has the depth of the character with greatest depth).  This
   *  is not the case for most fonts, so if we can access cmex10, the
   *  height of a character should be much bigger than the width.
   *  Otherwise, if we don't have cmex10, we'll get a character in another
   *  font with normal height and width.  In this case, we insert a message
   *  pointing the user to the jsMath site, and load one of the fallback
   *  definitions.
   *
   */
  Check: function () {
    var cookie = jsMath.Controls.cookie;
    this.CheckTeX();
    if (jsMath.nofonts) {
      if (cookie.autofont || cookie.font == 'tex') {
        cookie.font = this.fallback;
        if (cookie.warn) {
          jsMath.nofontMessage = 1;
          cookie.warn = 0; jsMath.Controls.SetCookie(0);
          if (jsMath.window.NoFontMessage) {jsMath.window.NoFontMessage()}
                               else {this.Message(this.message)}
        }
      }
    } else {
      if (cookie.autofont) {cookie.font = 'tex'}
      if (cookie.font == 'tex') return;
    }
    if (jsMath.noImgFonts) {cookie.font = 'unicode'}
    if (cookie.font == 'unicode') {
      jsMath.Setup.Script('jsMath-fallback-'+jsMath.platform+'.js');
      jsMath.Box.TeXnonfallback = jsMath.Box.TeX;
      jsMath.Box.TeX = jsMath.Box.TeXfallback;
      return;
    }
    if (!cookie.print && cookie.printwarn) {
      this.PrintMessage(
        (jsMath.Browser.alphaPrintBug && jsMath.Controls.cookie.alpha) ?
          this.print_message + this.alpha_message : this.print_message);
    }
    if (jsMath.Browser.waitForImages) {
      jsMath.Script.Push(jsMath.Script,"WaitForImage",jsMath.blank);
    }
    if (cookie.font == 'symbol') {
      jsMath.Setup.Script('jsMath-fallback-symbols.js');
      jsMath.Box.TeXnonfallback = jsMath.Box.TeX;
      jsMath.Box.TeX = jsMath.Box.TeXfallback;
      return;
    }
    jsMath.Img.SetFont({
      cmr10:  ['all'], cmmi10: ['all'], cmsy10: ['all'],
      cmex10: ['all'], cmbx10: ['all'], cmti10: ['all']
    });
    jsMath.Img.LoadFont('cm-fonts');
  },

  /*
   *  The message for when no TeX fonts.  You can eliminate this message
   *  by including
   *
   *      <script>jsMath = {Font: {Message: function () {}}}</script>
   *
   *  in your HTML file, before loading jsMath.js, if you want.  But this
   *  means the user may not know that he or she can get a better version
   *  of your page.
   */
  Message: function (message) {
    if (jsMath.Element("Warning")) return;
    var div = jsMath.Setup.DIV("Warning",{});
    div.innerHTML =
      '<center><table><tr><td>'
      + '<div id="jsMath_noFont"><div class="message">' + message
      + '<div style="text-align:left"><span style="float:left; margin: 8px 0px 0px 20px">'
      + '<span onclick="jsMath.Controls.Panel()" title=" Open the jsMath Control Panel " class="link">jsMath Control Panel</span>'
      + '</span><span style="margin: 8px 20px 0px 0px; float:right">'
      + '<span onclick="jsMath.Font.HideMessage()" title=" Remove this font warning message " class="link">Hide this Message</span>'
      + '</span></div><div style="height:6px"></div><br clear="all"/></div></div>'
      + '<div style="width:22em; height:1px"></div>'
      + '</td></tr></table></center><hr/>';
  },

  HideMessage: function () {
    var message = jsMath.Element("Warning");
    if (message) {message.style.display = "none"}
  },

  PrintMessage: function (message) {
    if (jsMath.Element("PrintWarning")) return;
    var div = jsMath.Setup.DIV("PrintWarning",{});
    div.innerHTML =
      '<center><table><tr><td>'
      + '<div class="message">' + message + '</div>'
      + '<div style="width:22em; height:1px"></div>'
      + '</td></tr></table></center><hr/>';
  },

  /*
   *  Register an extra font so jsMath knows about it
   */
  Register: function (data,force) {
    if (typeof(data) == 'string') {data = {name: data}}
    if (!jsMath.Setup.inited && !force) {
      this.register[this.register.length] = data;
      return;
    }
    var fontname = data.name; var name = fontname.replace(/10$/,'');
    var fontfam = jsMath.TeX.fam.length;
    if (data.prefix == null) {data.prefix = ""}
    if (!data.style) {data.style = "font-family: "+data.prefix+fontname+", serif"}
    if (!data.styles) {data.styles = {}}
    if (!data.macros) {data.macros = {}}
    /*
     *  Register font family
     */
    jsMath.TeX.fam[fontfam] = fontname;
    jsMath.TeX.famName[fontname] = fontfam;
    data.macros[name] = ['HandleFont',fontfam];
    jsMath.Add(jsMath.Parser.prototype.macros,data.macros);
    /*
     *  Set up styles
     */
    data.styles['.typeset .'+fontname] = data.style;
    jsMath.Setup.Styles(data.styles);
    if (jsMath.initialized) {jsMath.Script.Push(jsMath.Setup,'TeXfont',fontname)}
    /*
     *  Check for font and give message if missing
     */
    var cookie = jsMath.Controls.cookie;
    var hasTeXfont = !jsMath.nofonts &&
                      data.test(fontname,data.testChar,data.testFactor,data.prefix);
    if (hasTeXfont && cookie.font == 'tex') {
      if (data.tex) {data.tex(fontname,fontfam,data)}
      return;
    }
    if (!hasTeXfont && cookie.warn && cookie.font == 'tex' && !jsMath.nofonts) {
      if (!cookie.fonts.match("/"+fontname+"/")) {
        cookie.fonts += fontname + "/"; jsMath.Controls.SetCookie(0);
        if (!jsMath.Element("Warning")) this.Message(this.extra_message);
        var extra = jsMath.Element("ExtraFonts");
        if (extra) {
          if (extra.innerHTML != "") {extra.innerHTML += ','}
          extra.innerHTML += " " + data.prefix+fontname;
        }
      }
    }
    if (cookie.font == 'unicode' || jsMath.noImgFonts) {
      if (data.fallback) {data.fallback(fontname,fontfam,data)}
      return;
    }
    //  Image fonts
    var font = {}; font[fontname] = ['all'];
    jsMath.Img.SetFont(font);
    jsMath.Img.LoadFont(fontname);
    if (jsMath.initialized) {
      jsMath.Script.Push(jsMath.Img,'Scale');
      jsMath.Script.Push(jsMath.Img,'UpdateFonts');
    }
  },
  /*
   *  If fonts are registered before jsMath.Init() is called, jsMath.em
   *  will not be available, so they need to be delayed.
   */
  LoadRegistered: function () {
    var i = 0;
    while (i < this.register.length) {this.Register(this.register[i++],1)}
    this.register = [];
  },

  /*
   *  Load a font
   */
  Load: function (name) {jsMath.Setup.Script(this.URL(name))},
  URL: function (name) {return jsMath.Img.root+name+'/def.js'}

};

/***************************************************************************/

/*
 *  Implements the jsMath control panel.
 *  Much of the code is in jsMath-controls.html, which is
 *  loaded into a hidden IFRAME on demand
 */
jsMath.Controls = {

  //  Data stored in the jsMath cookie
  cookie: {
    scale: 100,
    font: 'tex', autofont: 1, scaleImg: 0, alpha: 1,
    warn: 1, fonts: '/', printwarn: 1, stayhires: 0,
    button: 1, progress: 1, asynch: 0, blank: 0,
    print: 0, keep: '0D', global: 'auto', hiddenGlobal: 1
  },

  cookiePath: '/',  // can also set cookieDomain
  noCookiePattern: /^(file|mk):$/,  // pattern for handling cookies locally


  /*
   *  Create the HTML needed for control panel
   */
  Init: function () {
    this.panel = jsMath.Setup.DIV("panel",{display:'none'});
    if (!jsMath.Browser.msieButtonBug) {this.Button()}
      else {setTimeout("jsMath.Controls.Button()",500)}
  },

  /*
   *  Load the control panel
   */
  Panel: function () {
    jsMath.Translate.Cancel();
    if (this.loaded) {this.Main()}
      else {jsMath.Script.delayedLoad(jsMath.root+"jsMath-controls.html")}
  },

  /*
   *  Create the control panel button
   */
  Button: function () {
    var button = jsMath.Setup.DIV("button",{});
    button.title = ' Open jsMath Control Panel ';
    button.innerHTML =
      '<span onclick="jsMath.Controls.Panel()">jsMath</span>';
    if (!jsMath.Global.isLocal && !jsMath.noShowGlobal) {
      button.innerHTML +=
        '<span id="jsMath_global" title=" Open jsMath Global Panel " '
          + 'onclick="jsMath.Global.Show(1)">Global&nbsp;</span>';
    }
    if (button.offsetWidth < 30) {button.style.width = "auto"}
    if (!this.cookie.button) {button.style.display = "none"}
  },

  /*
   *  Since MSIE doesn't handle position:float, we need to have the
   *  window repositioned every time the window scrolls.  We do that
   *  by hiding then showing the window, which apparently causes MSIE
   *  to recompute its location.  In MSIE7, that doesn't work anymore,
   *  so we have to move the window by hand.
   */
  MoveButton: function () {
    var controls = jsMath.Controls;
    if (!controls.button) {controls.button = jsMath.Element("button")}
    if (controls.button) controls.MoveElement(controls.button,3,2);
    var dx = 20; var dy = 20;
    if (controls.button) {dy = controls.button.offsetHeight + 6; dx = dy + 5}
    if (controls.panel)  controls.MoveElement(controls.panel,dx,dy);
  },
  MoveElement: function (obj,dx,dy) {
    if (jsMath.Browser.IE7) {
      var body = document.body;
      obj.style.right = "auto";
      obj.style.bottom = "auto";
      //
      // This position can't be overridden by CSS (grr)
      // (Perhaps we can look up the current position and which sides it's
      // attached to and use those.  What a pain.)
      //
      obj.style.left = body.clientWidth + body.scrollLeft - obj.offsetWidth - dx + "px";
      obj.style.top = body.clientHeight + body.scrollTop -  obj.offsetHeight - dy + "px";
    } else {
      obj.style.visibility = "hidden";
      obj.style.visibility = "visible";
    }
  },

  /*
   *  Get the cookie data from the browser
   *  (for file: references, use url '?' syntax)
   */
  GetCookie: function () {
    // save the current cookie settings as the defaults
    if (this.defaults == null) {this.defaults = {}}
    jsMath.Add(this.defaults,this.cookie); this.userSet = {};
    // get the browser's cookie data
    var cookies = jsMath.document.cookie;
    if (jsMath.window.location.protocol.match(this.noCookiePattern)) {
      cookies = this.localGetCookie();
      this.isLocalCookie = 1;
    }
    if (cookies.match(/jsMath=([^;]+)/)) {
      var data = unescape(RegExp.$1).split(/,/);
      for (var i = 0; i < data.length; i++) {
        var x = data[i].match(/(.*):(.*)/);
        if (x[2].match(/^\d+$/)) {x[2] = 1*x[2]} // convert from string
        this.cookie[x[1]] = x[2];
        this.userSet[x[1]] = 1;
      }
    }
  },
  localGetCookie: function () {
    return jsMath.window.location.search.substr(1);
  },

  /*
   *  Save the cookie data in the browser
   *  (for file: urls, append data like CGI reference)
   */
  SetCookie: function (warn) {
    var cookie = [];
    for (var id in this.cookie) {
      if (this.defaults[id] == null || this.cookie[id] != this.defaults[id])
        {cookie[cookie.length] = id + ':' + this.cookie[id]}
    }
    cookie = cookie.join(',');
    if (this.isLocalCookie) {
      if (warn == 2) {return 'jsMath='+escape(cookie)}
      this.localSetCookie(cookie,warn);
    } else {
      cookie = escape(cookie);
      if (cookie == '') {warn = 0}
      if (this.cookiePath) {cookie += '; path='+this.cookiePath}
      if (this.cookieDomain) {cookie += '; domain='+this.cookieDomain}
      if (this.cookie.keep != '0D') {
        var ms = {
          D: 1000*60*60*24,
          W: 1000*60*60*24*7,
          M: 1000*60*60*24*30,
          Y: 1000*60*60*24*365
        };
        var exp = new Date;
        exp.setTime(exp.getTime() +
            this.cookie.keep.substr(0,1) * ms[this.cookie.keep.substr(1,1)]);
        cookie += '; expires=' + exp.toGMTString();
      }
      if (cookie != '') {
        jsMath.document.cookie = 'jsMath='+cookie;
        var cookies = jsMath.document.cookie;
        if (warn && !cookies.match(/jsMath=/))
          {alert("Cookies must be enabled in order to save jsMath options")}
      }
    }
    return null;
  },
  localSetCookie: function (cookie,warn) {
    if (!warn) return;
    var href = String(jsMath.window.location).replace(/\?.*/,"");
    if (cookie != '') {href += '?jsMath=' + escape(cookie)}
    if (href != jsMath.window.location.href) {this.Reload(href)}
  },

  /*
   *  Reload the page (with the given URL)
   */
  Reload: function (url) {
    if (!this.loaded) return;
    this.loaded = 0; jsMath.Setup.inited = -100;
    jsMath.Global.ClearCache();
    if (url) {jsMath.window.location.replace(url)}
        else {jsMath.window.location.reload()}
  }

};

/***************************************************************************/

/*
 *  Implements the actions for clicking and double-clicking
 *  on math formulas
 */
jsMath.Click = {

  /*
   *  Handle clicking on math to get control panel
   */
  CheckClick: function (event) {
    if (!event) {event = jsMath.window.event}
    if (event.altKey) jsMath.Controls.Panel();
  },

  /*
   *  Handle double-click for seeing TeX code
   */

  CheckDblClick: function (event) {
    if (!event) {event = jsMath.window.event}
    if (!jsMath.Click.DblClick) {
      jsMath.Extension.Require('double-click',1);
      // Firefox clears the event, so copy it
      var tmpEvent = event; event = {};
      for (var id in tmpEvent) {event[id] = tmpEvent[id]}
    }
    jsMath.Script.Push(jsMath.Click,'DblClick',[event,this.alt]);
  }

};

/***************************************************************************/

/*
 *  The TeX font information
 */
jsMath.TeX = {

  //
  //  The TeX font parameters
  //
  thinmuskip:   3/18,
  medmuskip:    4/18,
  thickmuskip:  5/18,

  x_height:    .430554,
  quad:        1,
  num1:        .676508,
  num2:        .393732,
  num3:        .44373,
  denom1:      .685951,
  denom2:      .344841,
  sup1:        .412892,
  sup2:        .362892,
  sup3:        .288888,
  sub1:        .15,
  sub2:        .247217,
  sup_drop:    .386108,
  sub_drop:    .05,
  delim1:     2.39,
  delim2:     1.0,
  axis_height: .25,
  default_rule_thickness: .06,
  big_op_spacing1:  .111111,
  big_op_spacing2:  .166666,
  big_op_spacing3:  .2,
  big_op_spacing4:  .6,
  big_op_spacing5:  .1,

  integer:          6553.6,     // conversion of em's to TeX internal integer
  scriptspace:         .05,
  nulldelimiterspace:  .12,
  delimiterfactor:     901,
  delimitershortfall:   .5,
  scale:                 1,     //  scaling factor for font dimensions

  //  The TeX math atom types (see Appendix G of the TeXbook)
  atom: ['ord', 'op', 'bin', 'rel', 'open', 'close', 'punct', 'ord'],

  //  The TeX font families
  fam: ['cmr10','cmmi10','cmsy10','cmex10','cmti10','','cmbx10',''],
  famName: {cmr10:0, cmmi10:1, cmsy10:2, cmex10:3, cmti10:4, cmbx10:6},

  //  Encoding used by jsMath fonts
  encoding: [
    '&#xC0;', '&#xC1;', '&#xC2;', '&#xC3;', '&#xC4;', '&#xC5;', '&#xC6;', '&#xC7;',
    '&#xC8;', '&#xC9;', '&#xCA;', '&#xCB;', '&#xCC;', '&#xCD;', '&#xCE;', '&#xCF;',

    '&#xB0;', '&#xD1;', '&#xD2;', '&#xD3;', '&#xD4;', '&#xD5;', '&#xD6;', '&#xB7;',
    '&#xD8;', '&#xD9;', '&#xDA;', '&#xDB;', '&#xDC;', '&#xB5;', '&#xB6;', '&#xDF;',

    '&#xEF;', '!', '&#x22;', '#', '$', '%', '&#x26;', '&#x27;',
    '(', ')', '*', '+', ',', '-', '.', '/',

    '0', '1', '2', '3', '4', '5', '6', '7',
    '8', '9', ':', ';', '&#x3C;', '=', '&#x3E;', '?',

    '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
    'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',

    'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
    'X', 'Y', 'Z', '[', '&#x5C;', ']', '^', '_',

    '`', 'a', 'b', 'c', 'd', 'e', 'f', 'g',
    'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',

    'p', 'q', 'r', 's', 't', 'u', 'v', 'w',
    'x', 'y', 'z', '{', '|', '}', '&#x7E;', '&#xFF;'
  ],

  /*
   *  The following are the TeX font mappings and metrics.  The metric
   *  information comes directly from the TeX .tfm files.  Browser-specific
   *  adjustments are made to these tables in the Browser.Init() routine
   */
  cmr10: [
    [0.625,0.683], [0.833,0.683], [0.778,0.683], [0.694,0.683],
    [0.667,0.683], [0.75,0.683], [0.722,0.683], [0.778,0.683],
    [0.722,0.683], [0.778,0.683], [0.722,0.683],
    [0.583,0.694,0,{ic: 0.0778, krn: {'39': 0.0778, '63': 0.0778, '33': 0.0778, '41': 0.0778, '93': 0.0778}, lig: {'105': 14, '108': 15}}],
    [0.556,0.694], [0.556,0.694], [0.833,0.694], [0.833,0.694],

    [0.278,0.431], [0.306,0.431,0.194], [0.5,0.694], [0.5,0.694],
    [0.5,0.628], [0.5,0.694], [0.5,0.568], [0.75,0.694],
    [0.444,0,0.17], [0.5,0.694], [0.722,0.431], [0.778,0.431],
    [0.5,0.528,0.0972], [0.903,0.683], [1.01,0.683], [0.778,0.732,0.0486],

    [0.278,0.431,0,{krn: {'108': -0.278, '76': -0.319}}],
    [0.278,0.694,0,{lig: {'96': 60}}],
    [0.5,0.694], [0.833,0.694,0.194], [0.5,0.75,0.0556],
    [0.833,0.75,0.0556], [0.778,0.694],
    [0.278,0.694,0,{krn: {'63': 0.111, '33': 0.111}, lig: {'39': 34}}],
    [0.389,0.75,0.25], [0.389,0.75,0.25], [0.5,0.75],
    [0.778,0.583,0.0833], [0.278,0.106,0.194],
    [0.333,0.431,0,{lig: {'45': 123}}],
    [0.278,0.106], [0.5,0.75,0.25],

    [0.5,0.644], [0.5,0.644], [0.5,0.644], [0.5,0.644],
    [0.5,0.644], [0.5,0.644], [0.5,0.644], [0.5,0.644],
    [0.5,0.644], [0.5,0.644], [0.278,0.431], [0.278,0.431,0.194],
    [0.278,0.5,0.194], [0.778,0.367,-0.133], [0.472,0.5,0.194],
    [0.472,0.694,0,{lig: {'96': 62}}],

    [0.778,0.694],
    [0.75,0.683,0,{krn: {'116': -0.0278, '67': -0.0278, '79': -0.0278, '71': -0.0278, '85': -0.0278, '81': -0.0278, '84': -0.0833, '89': -0.0833, '86': -0.111, '87': -0.111}}],
    [0.708,0.683], [0.722,0.683],
    [0.764,0.683,0,{krn: {'88': -0.0278, '87': -0.0278, '65': -0.0278, '86': -0.0278, '89': -0.0278}}],
    [0.681,0.683],
    [0.653,0.683,0,{krn: {'111': -0.0833, '101': -0.0833, '117': -0.0833, '114': -0.0833, '97': -0.0833, '65': -0.111, '79': -0.0278, '67': -0.0278, '71': -0.0278, '81': -0.0278}}],
    [0.785,0.683], [0.75,0.683], [0.361,0.683,0,{krn: {'73': 0.0278}}],
    [0.514,0.683],
    [0.778,0.683,0,{krn: {'79': -0.0278, '67': -0.0278, '71': -0.0278, '81': -0.0278}}],
    [0.625,0.683,0,{krn: {'84': -0.0833, '89': -0.0833, '86': -0.111, '87': -0.111}}],
    [0.917,0.683], [0.75,0.683],
    [0.778,0.683,0,{krn: {'88': -0.0278, '87': -0.0278, '65': -0.0278, '86': -0.0278, '89': -0.0278}}],

    [0.681,0.683,0,{krn: {'65': -0.0833, '111': -0.0278, '101': -0.0278, '97': -0.0278, '46': -0.0833, '44': -0.0833}}],
    [0.778,0.683,0.194],
    [0.736,0.683,0,{krn: {'116': -0.0278, '67': -0.0278, '79': -0.0278, '71': -0.0278, '85': -0.0278, '81': -0.0278, '84': -0.0833, '89': -0.0833, '86': -0.111, '87': -0.111}}],
    [0.556,0.683],
    [0.722,0.683,0,{krn: {'121': -0.0278, '101': -0.0833, '111': -0.0833, '114': -0.0833, '97': -0.0833, '65': -0.0833, '117': -0.0833}}],
    [0.75,0.683],
    [0.75,0.683,0,{ic: 0.0139, krn: {'111': -0.0833, '101': -0.0833, '117': -0.0833, '114': -0.0833, '97': -0.0833, '65': -0.111, '79': -0.0278, '67': -0.0278, '71': -0.0278, '81': -0.0278}}],
    [1.03,0.683,0,{ic: 0.0139, krn: {'111': -0.0833, '101': -0.0833, '117': -0.0833, '114': -0.0833, '97': -0.0833, '65': -0.111, '79': -0.0278, '67': -0.0278, '71': -0.0278, '81': -0.0278}}],
    [0.75,0.683,0,{krn: {'79': -0.0278, '67': -0.0278, '71': -0.0278, '81': -0.0278}}],
    [0.75,0.683,0,{ic: 0.025, krn: {'101': -0.0833, '111': -0.0833, '114': -0.0833, '97': -0.0833, '65': -0.0833, '117': -0.0833}}],
    [0.611,0.683], [0.278,0.75,0.25], [0.5,0.694],
    [0.278,0.75,0.25], [0.5,0.694], [0.278,0.668],

    [0.278,0.694,0,{lig: {'96': 92}}],
    [0.5,0.431,0,{krn: {'118': -0.0278, '106': 0.0556, '121': -0.0278, '119': -0.0278}}],
    [0.556,0.694,0,{krn: {'101': 0.0278, '111': 0.0278, '120': -0.0278, '100': 0.0278, '99': 0.0278, '113': 0.0278, '118': -0.0278, '106': 0.0556, '121': -0.0278, '119': -0.0278}}],
    [0.444,0.431,0,{krn: {'104': -0.0278, '107': -0.0278}}],
    [0.556,0.694], [0.444,0.431],
    [0.306,0.694,0,{ic: 0.0778, krn: {'39': 0.0778, '63': 0.0778, '33': 0.0778, '41': 0.0778, '93': 0.0778}, lig: {'105': 12, '102': 11, '108': 13}}],
    [0.5,0.431,0.194,{ic: 0.0139, krn: {'106': 0.0278}}],
    [0.556,0.694,0,{krn: {'116': -0.0278, '117': -0.0278, '98': -0.0278, '121': -0.0278, '118': -0.0278, '119': -0.0278}}],
    [0.278,0.668], [0.306,0.668,0.194],
    [0.528,0.694,0,{krn: {'97': -0.0556, '101': -0.0278, '97': -0.0278, '111': -0.0278, '99': -0.0278}}],
    [0.278,0.694],
    [0.833,0.431,0,{krn: {'116': -0.0278, '117': -0.0278, '98': -0.0278, '121': -0.0278, '118': -0.0278, '119': -0.0278}}],
    [0.556,0.431,0,{krn: {'116': -0.0278, '117': -0.0278, '98': -0.0278, '121': -0.0278, '118': -0.0278, '119': -0.0278}}],
    [0.5,0.431,0,{krn: {'101': 0.0278, '111': 0.0278, '120': -0.0278, '100': 0.0278, '99': 0.0278, '113': 0.0278, '118': -0.0278, '106': 0.0556, '121': -0.0278, '119': -0.0278}}],

    [0.556,0.431,0.194,{krn: {'101': 0.0278, '111': 0.0278, '120': -0.0278, '100': 0.0278, '99': 0.0278, '113': 0.0278, '118': -0.0278, '106': 0.0556, '121': -0.0278, '119': -0.0278}}],
    [0.528,0.431,0.194], [0.392,0.431], [0.394,0.431],
    [0.389,0.615,0,{krn: {'121': -0.0278, '119': -0.0278}}],
    [0.556,0.431,0,{krn: {'119': -0.0278}}],
    [0.528,0.431,0,{ic: 0.0139, krn: {'97': -0.0556, '101': -0.0278, '97': -0.0278, '111': -0.0278, '99': -0.0278}}],
    [0.722,0.431,0,{ic: 0.0139, krn: {'101': -0.0278, '97': -0.0278, '111': -0.0278, '99': -0.0278}}],
    [0.528,0.431],
    [0.528,0.431,0.194,{ic: 0.0139, krn: {'111': -0.0278, '101': -0.0278, '97': -0.0278, '46': -0.0833, '44': -0.0833}}],
    [0.444,0.431], [0.5,0.431,0,{ic: 0.0278, lig: {'45': 124}}],
    [1,0.431,0,{ic: 0.0278}], [0.5,0.694], [0.5,0.668], [0.5,0.668]
  ],

  cmmi10: [
    [0.615,0.683,0,{ic: 0.139, krn: {'61': -0.0556, '59': -0.111, '58': -0.111, '127': 0.0833}}],
    [0.833,0.683,0,{krn: {'127': 0.167}}],
    [0.763,0.683,0,{ic: 0.0278, krn: {'127': 0.0833}}],
    [0.694,0.683,0,{krn: {'127': 0.167}}],
    [0.742,0.683,0,{ic: 0.0757, krn: {'127': 0.0833}}],
    [0.831,0.683,0,{ic: 0.0812, krn: {'61': -0.0556, '59': -0.0556, '58': -0.0556, '127': 0.0556}}],
    [0.78,0.683,0,{ic: 0.0576, krn: {'127': 0.0833}}],
    [0.583,0.683,0,{ic: 0.139, krn: {'61': -0.0556, '59': -0.111, '58': -0.111, '127': 0.0556}}],
    [0.667,0.683,0,{krn: {'127': 0.0833}}],
    [0.612,0.683,0,{ic: 0.11, krn: {'61': -0.0556, '59': -0.0556, '58': -0.0556, '127': 0.0556}}],
    [0.772,0.683,0,{ic: 0.0502, krn: {'127': 0.0833}}],
    [0.64,0.431,0,{ic: 0.0037, krn: {'127': 0.0278}}],
    [0.566,0.694,0.194,{ic: 0.0528, krn: {'127': 0.0833}}],
    [0.518,0.431,0.194,{ic: 0.0556}],
    [0.444,0.694,0,{ic: 0.0378, krn: {'59': -0.0556, '58': -0.0556, '127': 0.0556}}],
    [0.406,0.431,0,{krn: {'127': 0.0556}}],

    [0.438,0.694,0.194,{ic: 0.0738, krn: {'127': 0.0833}}],
    [0.497,0.431,0.194,{ic: 0.0359, krn: {'127': 0.0556}}],
    [0.469,0.694,0,{ic: 0.0278, krn: {'127': 0.0833}}],
    [0.354,0.431,0,{krn: {'127': 0.0556}}],
    [0.576,0.431], [0.583,0.694],
    [0.603,0.431,0.194,{krn: {'127': 0.0278}}],
    [0.494,0.431,0,{ic: 0.0637, krn: {'59': -0.0556, '58': -0.0556, '127': 0.0278}}],
    [0.438,0.694,0.194,{ic: 0.046, krn: {'127': 0.111}}],
    [0.57,0.431,0,{ic: 0.0359}],
    [0.517,0.431,0.194,{krn: {'127': 0.0833}}],
    [0.571,0.431,0,{ic: 0.0359, krn: {'59': -0.0556, '58': -0.0556}}],
    [0.437,0.431,0,{ic: 0.113, krn: {'59': -0.0556, '58': -0.0556, '127': 0.0278}}],
    [0.54,0.431,0,{ic: 0.0359, krn: {'127': 0.0278}}],
    [0.596,0.694,0.194,{krn: {'127': 0.0833}}],
    [0.626,0.431,0.194,{krn: {'127': 0.0556}}],

    [0.651,0.694,0.194,{ic: 0.0359, krn: {'127': 0.111}}],
    [0.622,0.431,0,{ic: 0.0359}],
    [0.466,0.431,0,{krn: {'127': 0.0833}}],
    [0.591,0.694,0,{krn: {'127': 0.0833}}],
    [0.828,0.431,0,{ic: 0.0278}],
    [0.517,0.431,0.194,{krn: {'127': 0.0833}}],
    [0.363,0.431,0.0972,{ic: 0.0799, krn: {'127': 0.0833}}],
    [0.654,0.431,0.194,{krn: {'127': 0.0833}}],
    [1,0.367,-0.133], [1,0.367,-0.133], [1,0.367,-0.133], [1,0.367,-0.133],
    [0.278,0.464,-0.0363], [0.278,0.464,-0.0363], [0.5,0.465,-0.0347], [0.5,0.465,-0.0347],

    [0.5,0.431], [0.5,0.431], [0.5,0.431], [0.5,0.431,0.194],
    [0.5,0.431,0.194], [0.5,0.431,0.194], [0.5,0.644], [0.5,0.431,0.194],
    [0.5,0.644], [0.5,0.431,0.194], [0.278,0.106], [0.278,0.106,0.194],
    [0.778,0.539,0.0391],
    [0.5,0.75,0.25,{krn: {'1': -0.0556, '65': -0.0556, '77': -0.0556, '78': -0.0556, '89': 0.0556, '90': -0.0556}}],
    [0.778,0.539,0.0391], [0.5,0.465,-0.0347],

    [0.531,0.694,0,{ic: 0.0556, krn: {'127': 0.0833}}],
    [0.75,0.683,0,{krn: {'127': 0.139}}],
    [0.759,0.683,0,{ic: 0.0502, krn: {'127': 0.0833}}],
    [0.715,0.683,0,{ic: 0.0715, krn: {'61': -0.0278, '59': -0.0556, '58': -0.0556, '127': 0.0833}}],
    [0.828,0.683,0,{ic: 0.0278, krn: {'127': 0.0556}}],
    [0.738,0.683,0,{ic: 0.0576, krn: {'127': 0.0833}}],
    [0.643,0.683,0,{ic: 0.139, krn: {'61': -0.0556, '59': -0.111, '58': -0.111, '127': 0.0833}}],
    [0.786,0.683,0,{krn: {'127': 0.0833}}],
    [0.831,0.683,0,{ic: 0.0812, krn: {'61': -0.0556, '59': -0.0556, '58': -0.0556, '127': 0.0556}}],
    [0.44,0.683,0,{ic: 0.0785, krn: {'127': 0.111}}],
    [0.555,0.683,0,{ic: 0.0962, krn: {'61': -0.0556, '59': -0.111, '58': -0.111, '127': 0.167}}],
    [0.849,0.683,0,{ic: 0.0715, krn: {'61': -0.0556, '59': -0.0556, '58': -0.0556, '127': 0.0556}}],
    [0.681,0.683,0,{krn: {'127': 0.0278}}],
    [0.97,0.683,0,{ic: 0.109, krn: {'61': -0.0556, '59': -0.0556, '58': -0.0556, '127': 0.0833}}],
    [0.803,0.683,0,{ic: 0.109, krn: {'61': -0.0833, '61': -0.0278, '59': -0.0556, '58': -0.0556, '127': 0.0833}}],
    [0.763,0.683,0,{ic: 0.0278, krn: {'127': 0.0833}}],

    [0.642,0.683,0,{ic: 0.139, krn: {'61': -0.0556, '59': -0.111, '58': -0.111, '127': 0.0833}}],
    [0.791,0.683,0.194,{krn: {'127': 0.0833}}],
    [0.759,0.683,0,{ic: 0.00773, krn: {'127': 0.0833}}],
    [0.613,0.683,0,{ic: 0.0576, krn: {'61': -0.0556, '59': -0.0556, '58': -0.0556, '127': 0.0833}}],
    [0.584,0.683,0,{ic: 0.139, krn: {'61': -0.0278, '59': -0.0556, '58': -0.0556, '127': 0.0833}}],
    [0.683,0.683,0,{ic: 0.109, krn: {'59': -0.111, '58': -0.111, '61': -0.0556, '127': 0.0278}}],
    [0.583,0.683,0,{ic: 0.222, krn: {'59': -0.167, '58': -0.167, '61': -0.111}}],
    [0.944,0.683,0,{ic: 0.139, krn: {'59': -0.167, '58': -0.167, '61': -0.111}}],
    [0.828,0.683,0,{ic: 0.0785, krn: {'61': -0.0833, '61': -0.0278, '59': -0.0556, '58': -0.0556, '127': 0.0833}}],
    [0.581,0.683,0,{ic: 0.222, krn: {'59': -0.167, '58': -0.167, '61': -0.111}}],
    [0.683,0.683,0,{ic: 0.0715, krn: {'61': -0.0556, '59': -0.0556, '58': -0.0556, '127': 0.0833}}],
    [0.389,0.75], [0.389,0.694,0.194], [0.389,0.694,0.194],
    [1,0.358,-0.142], [1,0.358,-0.142],

    [0.417,0.694,0,{krn: {'127': 0.111}}],
    [0.529,0.431], [0.429,0.694], [0.433,0.431,0,{krn: {'127': 0.0556}}],
    [0.52,0.694,0,{krn: {'89': 0.0556, '90': -0.0556, '106': -0.111, '102': -0.167, '127': 0.167}}],
    [0.466,0.431,0,{krn: {'127': 0.0556}}],
    [0.49,0.694,0.194,{ic: 0.108, krn: {'59': -0.0556, '58': -0.0556, '127': 0.167}}],
    [0.477,0.431,0.194,{ic: 0.0359, krn: {'127': 0.0278}}],
    [0.576,0.694,0,{krn: {'127': -0.0278}}], [0.345,0.66],
    [0.412,0.66,0.194,{ic: 0.0572, krn: {'59': -0.0556, '58': -0.0556}}],
    [0.521,0.694,0,{ic: 0.0315}], [0.298,0.694,0,{ic: 0.0197, krn: {'127': 0.0833}}],
    [0.878,0.431], [0.6,0.431], [0.485,0.431,0,{krn: {'127': 0.0556}}],

    [0.503,0.431,0.194,{krn: {'127': 0.0833}}],
    [0.446,0.431,0.194,{ic: 0.0359, krn: {'127': 0.0833}}],
    [0.451,0.431,0,{ic: 0.0278, krn: {'59': -0.0556, '58': -0.0556, '127': 0.0556}}],
    [0.469,0.431,0,{krn: {'127': 0.0556}}], [0.361,0.615,0,{krn: {'127': 0.0833}}],
    [0.572,0.431,0,{krn: {'127': 0.0278}}],
    [0.485,0.431,0,{ic: 0.0359, krn: {'127': 0.0278}}],
    [0.716,0.431,0,{ic: 0.0269, krn: {'127': 0.0833}}],
    [0.572,0.431,0,{krn: {'127': 0.0278}}],
    [0.49,0.431,0.194,{ic: 0.0359, krn: {'127': 0.0556}}],
    [0.465,0.431,0,{ic: 0.044, krn: {'127': 0.0556}}],
    [0.322,0.431,0,{krn: {'127': 0.0278}}],
    [0.384,0.431,0.194,{krn: {'127': 0.0833}}],
    [0.636,0.431,0.194,{krn: {'127': 0.111}}],
    [0.5,0.714,0,{ic: 0.154}], [0.278,0.694,0,{ic: 0.399}]
  ],

  cmsy10: [
    [0.778,0.583,0.0833], [0.278,0.444,-0.0556], [0.778,0.583,0.0833],
    [0.5,0.465,-0.0347], [0.778,0.583,0.0833], [0.5,0.444,-0.0556],
    [0.778,0.583,0.0833], [0.778,0.583,0.0833], [0.778,0.583,0.0833],
    [0.778,0.583,0.0833], [0.778,0.583,0.0833], [0.778,0.583,0.0833],
    [0.778,0.583,0.0833], [1,0.694,0.194], [0.5,0.444,-0.0556], [0.5,0.444,-0.0556],

    [0.778,0.464,-0.0363], [0.778,0.464,-0.0363], [0.778,0.636,0.136],
    [0.778,0.636,0.136], [0.778,0.636,0.136], [0.778,0.636,0.136],
    [0.778,0.636,0.136], [0.778,0.636,0.136], [0.778,0.367,-0.133],
    [0.778,0.483,-0.0169], [0.778,0.539,0.0391], [0.778,0.539,0.0391],
    [1,0.539,0.0391], [1,0.539,0.0391], [0.778,0.539,0.0391], [0.778,0.539,0.0391],

    [1,0.367,-0.133], [1,0.367,-0.133], [0.5,0.694,0.194], [0.5,0.694,0.194],
    [1,0.367,-0.133], [1,0.694,0.194], [1,0.694,0.194], [0.778,0.464,-0.0363],
    [1,0.367,-0.133], [1,0.367,-0.133], [0.611,0.694,0.194], [0.611,0.694,0.194],
    [1,0.367,-0.133], [1,0.694,0.194], [1,0.694,0.194], [0.778,0.431],

    [0.275,0.556], [1,0.431], [0.667,0.539,0.0391], [0.667,0.539,0.0391],
    [0.889,0.694,0.194], [0.889,0.694,0.194], [0,0.694,0.194], [0,0.367,-0.133],
    [0.556,0.694], [0.556,0.694], [0.667,0.431], [0.5,0.75,0.0556],
    [0.722,0.694], [0.722,0.694], [0.778,0.694], [0.778,0.694],

    [0.611,0.694], [0.798,0.683,0,{krn: {'48': 0.194}}],
    [0.657,0.683,0,{ic: 0.0304, krn: {'48': 0.139}}],
    [0.527,0.683,0,{ic: 0.0583, krn: {'48': 0.139}}],
    [0.771,0.683,0,{ic: 0.0278, krn: {'48': 0.0833}}],
    [0.528,0.683,0,{ic: 0.0894, krn: {'48': 0.111}}],
    [0.719,0.683,0,{ic: 0.0993, krn: {'48': 0.111}}],
    [0.595,0.683,0.0972,{ic: 0.0593, krn: {'48': 0.111}}],
    [0.845,0.683,0,{ic: 0.00965, krn: {'48': 0.111}}],
    [0.545,0.683,0,{ic: 0.0738, krn: {'48': 0.0278}}],
    [0.678,0.683,0.0972,{ic: 0.185, krn: {'48': 0.167}}],
    [0.762,0.683,0,{ic: 0.0144, krn: {'48': 0.0556}}],
    [0.69,0.683,0,{krn: {'48': 0.139}}], [1.2,0.683,0,{krn: {'48': 0.139}}],
    [0.82,0.683,0,{ic: 0.147, krn: {'48': 0.0833}}],
    [0.796,0.683,0,{ic: 0.0278, krn: {'48': 0.111}}],

    [0.696,0.683,0,{ic: 0.0822, krn: {'48': 0.0833}}],
    [0.817,0.683,0.0972,{krn: {'48': 0.111}}],
    [0.848,0.683,0,{krn: {'48': 0.0833}}],
    [0.606,0.683,0,{ic: 0.075, krn: {'48': 0.139}}],
    [0.545,0.683,0,{ic: 0.254, krn: {'48': 0.0278}}],
    [0.626,0.683,0,{ic: 0.0993, krn: {'48': 0.0833}}],
    [0.613,0.683,0,{ic: 0.0822, krn: {'48': 0.0278}}],
    [0.988,0.683,0,{ic: 0.0822, krn: {'48': 0.0833}}],
    [0.713,0.683,0,{ic: 0.146, krn: {'48': 0.139}}],
    [0.668,0.683,0.0972,{ic: 0.0822, krn: {'48': 0.0833}}],
    [0.725,0.683,0,{ic: 0.0794, krn: {'48': 0.139}}],
    [0.667,0.556], [0.667,0.556], [0.667,0.556], [0.667,0.556], [0.667,0.556],

    [0.611,0.694], [0.611,0.694], [0.444,0.75,0.25], [0.444,0.75,0.25],
    [0.444,0.75,0.25], [0.444,0.75,0.25], [0.5,0.75,0.25], [0.5,0.75,0.25],
    [0.389,0.75,0.25], [0.389,0.75,0.25], [0.278,0.75,0.25], [0.5,0.75,0.25],
    [0.5,0.75,0.25], [0.611,0.75,0.25], [0.5,0.75,0.25], [0.278,0.694,0.194],

    [0.833,0.04,0.96], [0.75,0.683], [0.833,0.683], [0.417,0.694,0.194,{ic: 0.111}],
    [0.667,0.556], [0.667,0.556], [0.778,0.636,0.136], [0.778,0.636,0.136],
    [0.444,0.694,0.194], [0.444,0.694,0.194], [0.444,0.694,0.194],
    [0.611,0.694,0.194], [0.778,0.694,0.13], [0.778,0.694,0.13],
    [0.778,0.694,0.13], [0.778,0.694,0.13]
  ],

  cmex10: [
    [0.458,0.04,1.16,{n: 16}], [0.458,0.04,1.16,{n: 17}],
    [0.417,0.04,1.16,{n: 104}], [0.417,0.04,1.16,{n: 105}],
    [0.472,0.04,1.16,{n: 106}], [0.472,0.04,1.16,{n: 107}],
    [0.472,0.04,1.16,{n: 108}], [0.472,0.04,1.16,{n: 109}],
    [0.583,0.04,1.16,{n: 110}], [0.583,0.04,1.16,{n: 111}],
    [0.472,0.04,1.16,{n: 68}], [0.472,0.04,1.16,{n: 69}],
    [0.333,0,0.6,{delim: {rep: 12}}], [0.556,0,0.6,{delim: {rep: 13}}],
    [0.578,0.04,1.16,{n: 46}], [0.578,0.04,1.16,{n: 47}],

    [0.597,0.04,1.76,{n: 18}], [0.597,0.04,1.76,{n: 19}],
    [0.736,0.04,2.36,{n: 32}], [0.736,0.04,2.36,{n: 33}],
    [0.528,0.04,2.36,{n: 34}], [0.528,0.04,2.36,{n: 35}],
    [0.583,0.04,2.36,{n: 36}], [0.583,0.04,2.36,{n: 37}],
    [0.583,0.04,2.36,{n: 38}], [0.583,0.04,2.36,{n: 39}],
    [0.75,0.04,2.36,{n: 40}], [0.75,0.04,2.36,{n: 41}],
    [0.75,0.04,2.36,{n: 42}], [0.75,0.04,2.36,{n: 43}],
    [1.04,0.04,2.36,{n: 44}], [1.04,0.04,2.36,{n: 45}],

    [0.792,0.04,2.96,{n: 48}], [0.792,0.04,2.96,{n: 49}],
    [0.583,0.04,2.96,{n: 50}], [0.583,0.04,2.96,{n: 51}],
    [0.639,0.04,2.96,{n: 52}], [0.639,0.04,2.96,{n: 53}],
    [0.639,0.04,2.96,{n: 54}], [0.639,0.04,2.96,{n: 55}],
    [0.806,0.04,2.96,{n: 56}], [0.806,0.04,2.96,{n: 57}],
    [0.806,0.04,2.96], [0.806,0.04,2.96],
    [1.28,0.04,2.96], [1.28,0.04,2.96],
    [0.811,0.04,1.76,{n: 30}], [0.811,0.04,1.76,{n: 31}],

    [0.875,0.04,1.76,{delim: {top: 48, bot: 64, rep: 66}}],
    [0.875,0.04,1.76,{delim: {top: 49, bot: 65, rep: 67}}],
    [0.667,0.04,1.76,{delim: {top: 50, bot: 52, rep: 54}}],
    [0.667,0.04,1.76,{delim: {top: 51, bot: 53, rep: 55}}],
    [0.667,0.04,1.76,{delim: {bot: 52, rep: 54}}],
    [0.667,0.04,1.76,{delim: {bot: 53, rep: 55}}],
    [0.667,0,0.6,{delim: {top: 50, rep: 54}}],
    [0.667,0,0.6,{delim: {top: 51, rep: 55}}],
    [0.889,0,0.9,{delim: {top: 56, mid: 60, bot: 58, rep: 62}}],
    [0.889,0,0.9,{delim: {top: 57, mid: 61, bot: 59, rep: 62}}],
    [0.889,0,0.9,{delim: {top: 56, bot: 58, rep: 62}}],
    [0.889,0,0.9,{delim: {top: 57, bot: 59, rep: 62}}],
    [0.889,0,1.8,{delim: {rep: 63}}],
    [0.889,0,1.8,{delim: {rep: 119}}],
    [0.889,0,0.3,{delim: {rep: 62}}],
    [0.667,0,0.6,{delim: {top: 120, bot: 121, rep: 63}}],

    [0.875,0.04,1.76,{delim: {top: 56, bot: 59, rep: 62}}],
    [0.875,0.04,1.76,{delim: {top: 57, bot: 58, rep: 62}}],
    [0.875,0,0.6,{delim: {rep: 66}}], [0.875,0,0.6,{delim: {rep: 67}}],
    [0.611,0.04,1.76,{n: 28}], [0.611,0.04,1.76,{n: 29}],
    [0.833,0,1,{n: 71}], [1.11,0.1,1.5], [0.472,0,1.11,{ic: 0.194, n: 73}],
    [0.556,0,2.22,{ic: 0.444}], [1.11,0,1,{n: 75}], [1.51,0.1,1.5],
    [1.11,0,1,{n: 77}], [1.51,0.1,1.5], [1.11,0,1,{n: 79}], [1.51,0.1,1.5],

    [1.06,0,1,{n: 88}], [0.944,0,1,{n: 89}], [0.472,0,1.11,{ic: 0.194, n: 90}],
    [0.833,0,1,{n: 91}], [0.833,0,1,{n: 92}], [0.833,0,1,{n: 93}],
    [0.833,0,1,{n: 94}], [0.833,0,1,{n: 95}], [1.44,0.1,1.5],
    [1.28,0.1,1.5], [0.556,0,2.22,{ic: 0.444}], [1.11,0.1,1.5],
    [1.11,0.1,1.5], [1.11,0.1,1.5], [1.11,0.1,1.5], [1.11,0.1,1.5],

    [0.944,0,1,{n: 97}], [1.28,0.1,1.5], [0.556,0.722,0,{n: 99}],
    [1,0.75,0,{n: 100}], [1.44,0.75], [0.556,0.722,0,{n: 102}],
    [1,0.75,0,{n: 103}], [1.44,0.75], [0.472,0.04,1.76,{n: 20}],
    [0.472,0.04,1.76,{n: 21}], [0.528,0.04,1.76,{n: 22}],
    [0.528,0.04,1.76,{n: 23}], [0.528,0.04,1.76,{n: 24}],
    [0.528,0.04,1.76,{n: 25}], [0.667,0.04,1.76,{n: 26}],
    [0.667,0.04,1.76,{n: 27}],

    [1,0.04,1.16,{n: 113}], [1,0.04,1.76,{n: 114}], [1,0.04,2.36,{n: 115}],
    [1,0.04,2.96,{n: 116}], [1.06,0,1.8,{delim: {top: 118, bot: 116, rep: 117}}],
    [1.06,0,0.6], [1.06,0.04,0.56],
    [0.778,0,0.6,{delim: {top: 126, bot: 127, rep: 119}}],
    [0.667,0,0.6,{delim: {top: 120, rep: 63}}],
    [0.667,0,0.6,{delim: {bot: 121, rep: 63}}],
    [0.45,0.12], [0.45,0.12], [0.45,0.12], [0.45,0.12],
    [0.778,0,0.6,{delim: {top: 126, rep: 119}}],
    [0.778,0,0.6,{delim: {bot: 127, rep: 119}}]
  ],

  cmti10: [
    [0.627,0.683,0,{ic: 0.133}], [0.818,0.683], [0.767,0.683,0,{ic: 0.094}],
    [0.692,0.683], [0.664,0.683,0,{ic: 0.153}], [0.743,0.683,0,{ic: 0.164}],
    [0.716,0.683,0,{ic: 0.12}], [0.767,0.683,0,{ic: 0.111}],
    [0.716,0.683,0,{ic: 0.0599}], [0.767,0.683,0,{ic: 0.111}],
    [0.716,0.683,0,{ic: 0.103}],
    [0.613,0.694,0.194,{ic: 0.212, krn: {'39': 0.104, '63': 0.104, '33': 0.104, '41': 0.104, '93': 0.104}, lig: {'105': 14, '108': 15}}],
    [0.562,0.694,0.194,{ic: 0.103}], [0.588,0.694,0.194,{ic: 0.103}],
    [0.882,0.694,0.194,{ic: 0.103}], [0.894,0.694,0.194,{ic: 0.103}],

    [0.307,0.431,0,{ic: 0.0767}], [0.332,0.431,0.194,{ic: 0.0374}],
    [0.511,0.694], [0.511,0.694,0,{ic: 0.0969}], [0.511,0.628,0,{ic: 0.083}],
    [0.511,0.694,0,{ic: 0.108}], [0.511,0.562,0,{ic: 0.103}], [0.831,0.694],
    [0.46,0,0.17], [0.537,0.694,0.194,{ic: 0.105}], [0.716,0.431,0,{ic: 0.0751}],
    [0.716,0.431,0,{ic: 0.0751}], [0.511,0.528,0.0972,{ic: 0.0919}],
    [0.883,0.683,0,{ic: 0.12}], [0.985,0.683,0,{ic: 0.12}],
    [0.767,0.732,0.0486,{ic: 0.094}],

    [0.256,0.431,0,{krn: {'108': -0.256, '76': -0.321}}],
    [0.307,0.694,0,{ic: 0.124, lig: {'96': 60}}],
    [0.514,0.694,0,{ic: 0.0696}], [0.818,0.694,0.194,{ic: 0.0662}],
    [0.769,0.694], [0.818,0.75,0.0556,{ic: 0.136}],
    [0.767,0.694,0,{ic: 0.0969}],
    [0.307,0.694,0,{ic: 0.124, krn: {'63': 0.102, '33': 0.102}, lig: {'39': 34}}],
    [0.409,0.75,0.25,{ic: 0.162}], [0.409,0.75,0.25,{ic: 0.0369}],
    [0.511,0.75,0,{ic: 0.149}], [0.767,0.562,0.0567,{ic: 0.0369}],
    [0.307,0.106,0.194], [0.358,0.431,0,{ic: 0.0283, lig: {'45': 123}}],
    [0.307,0.106], [0.511,0.75,0.25,{ic: 0.162}],

    [0.511,0.644,0,{ic: 0.136}], [0.511,0.644,0,{ic: 0.136}],
    [0.511,0.644,0,{ic: 0.136}], [0.511,0.644,0,{ic: 0.136}],
    [0.511,0.644,0.194,{ic: 0.136}], [0.511,0.644,0,{ic: 0.136}],
    [0.511,0.644,0,{ic: 0.136}], [0.511,0.644,0.194,{ic: 0.136}],
    [0.511,0.644,0,{ic: 0.136}], [0.511,0.644,0,{ic: 0.136}],
    [0.307,0.431,0,{ic: 0.0582}], [0.307,0.431,0.194,{ic: 0.0582}],
    [0.307,0.5,0.194,{ic: 0.0756}], [0.767,0.367,-0.133,{ic: 0.0662}],
    [0.511,0.5,0.194], [0.511,0.694,0,{ic: 0.122, lig: {'96': 62}}],

    [0.767,0.694,0,{ic: 0.096}],
    [0.743,0.683,0,{krn: {'110': -0.0256, '108': -0.0256, '114': -0.0256, '117': -0.0256, '109': -0.0256, '116': -0.0256, '105': -0.0256, '67': -0.0256, '79': -0.0256, '71': -0.0256, '104': -0.0256, '98': -0.0256, '85': -0.0256, '107': -0.0256, '118': -0.0256, '119': -0.0256, '81': -0.0256, '84': -0.0767, '89': -0.0767, '86': -0.102, '87': -0.102, '101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.704,0.683,0,{ic: 0.103}], [0.716,0.683,0,{ic: 0.145}],
    [0.755,0.683,0,{ic: 0.094, krn: {'88': -0.0256, '87': -0.0256, '65': -0.0256, '86': -0.0256, '89': -0.0256}}],
    [0.678,0.683,0,{ic: 0.12}],
    [0.653,0.683,0,{ic: 0.133, krn: {'111': -0.0767, '101': -0.0767, '117': -0.0767, '114': -0.0767, '97': -0.0767, '65': -0.102, '79': -0.0256, '67': -0.0256, '71': -0.0256, '81': -0.0256}}],
    [0.774,0.683,0,{ic: 0.0872}], [0.743,0.683,0,{ic: 0.164}],
    [0.386,0.683,0,{ic: 0.158}], [0.525,0.683,0,{ic: 0.14}],
    [0.769,0.683,0,{ic: 0.145, krn: {'79': -0.0256, '67': -0.0256, '71': -0.0256, '81': -0.0256}}],
    [0.627,0.683,0,{krn: {'84': -0.0767, '89': -0.0767, '86': -0.102, '87': -0.102, '101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.897,0.683,0,{ic: 0.164}], [0.743,0.683,0,{ic: 0.164}],
    [0.767,0.683,0,{ic: 0.094, krn: {'88': -0.0256, '87': -0.0256, '65': -0.0256, '86': -0.0256, '89': -0.0256}}],

    [0.678,0.683,0,{ic: 0.103, krn: {'65': -0.0767}}],
    [0.767,0.683,0.194,{ic: 0.094}],
    [0.729,0.683,0,{ic: 0.0387, krn: {'110': -0.0256, '108': -0.0256, '114': -0.0256, '117': -0.0256, '109': -0.0256, '116': -0.0256, '105': -0.0256, '67': -0.0256, '79': -0.0256, '71': -0.0256, '104': -0.0256, '98': -0.0256, '85': -0.0256, '107': -0.0256, '118': -0.0256, '119': -0.0256, '81': -0.0256, '84': -0.0767, '89': -0.0767, '86': -0.102, '87': -0.102, '101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.562,0.683,0,{ic: 0.12}],
    [0.716,0.683,0,{ic: 0.133, krn: {'121': -0.0767, '101': -0.0767, '111': -0.0767, '114': -0.0767, '97': -0.0767, '117': -0.0767, '65': -0.0767}}],
    [0.743,0.683,0,{ic: 0.164}],
    [0.743,0.683,0,{ic: 0.184, krn: {'111': -0.0767, '101': -0.0767, '117': -0.0767, '114': -0.0767, '97': -0.0767, '65': -0.102, '79': -0.0256, '67': -0.0256, '71': -0.0256, '81': -0.0256}}],
    [0.999,0.683,0,{ic: 0.184, krn: {'65': -0.0767}}],
    [0.743,0.683,0,{ic: 0.158, krn: {'79': -0.0256, '67': -0.0256, '71': -0.0256, '81': -0.0256}}],
    [0.743,0.683,0,{ic: 0.194, krn: {'101': -0.0767, '111': -0.0767, '114': -0.0767, '97': -0.0767, '117': -0.0767, '65': -0.0767}}],
    [0.613,0.683,0,{ic: 0.145}], [0.307,0.75,0.25,{ic: 0.188}],
    [0.514,0.694,0,{ic: 0.169}], [0.307,0.75,0.25,{ic: 0.105}],
    [0.511,0.694,0,{ic: 0.0665}], [0.307,0.668,0,{ic: 0.118}],

    [0.307,0.694,0,{ic: 0.124, lig: {'96': 92}}], [0.511,0.431,0,{ic: 0.0767}],
    [0.46,0.694,0,{ic: 0.0631, krn: {'101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.46,0.431,0,{ic: 0.0565, krn: {'101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.511,0.694,0,{ic: 0.103, krn: {'108': 0.0511}}],
    [0.46,0.431,0,{ic: 0.0751, krn: {'101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.307,0.694,0.194,{ic: 0.212, krn: {'39': 0.104, '63': 0.104, '33': 0.104, '41': 0.104, '93': 0.104}, lig: {'105': 12, '102': 11, '108': 13}}],
    [0.46,0.431,0.194,{ic: 0.0885}], [0.511,0.694,0,{ic: 0.0767}],
    [0.307,0.655,0,{ic: 0.102}], [0.307,0.655,0.194,{ic: 0.145}],
    [0.46,0.694,0,{ic: 0.108}], [0.256,0.694,0,{ic: 0.103, krn: {'108': 0.0511}}],
    [0.818,0.431,0,{ic: 0.0767}], [0.562,0.431,0,{ic: 0.0767, krn: {'39': -0.102}}],
    [0.511,0.431,0,{ic: 0.0631, krn: {'101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],

    [0.511,0.431,0.194,{ic: 0.0631, krn: {'101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.46,0.431,0.194,{ic: 0.0885}],
    [0.422,0.431,0,{ic: 0.108, krn: {'101': -0.0511, '97': -0.0511, '111': -0.0511, '100': -0.0511, '99': -0.0511, '103': -0.0511, '113': -0.0511}}],
    [0.409,0.431,0,{ic: 0.0821}], [0.332,0.615,0,{ic: 0.0949}],
    [0.537,0.431,0,{ic: 0.0767}], [0.46,0.431,0,{ic: 0.108}],
    [0.664,0.431,0,{ic: 0.108, krn: {'108': 0.0511}}],
    [0.464,0.431,0,{ic: 0.12}], [0.486,0.431,0.194,{ic: 0.0885}],
    [0.409,0.431,0,{ic: 0.123}], [0.511,0.431,0,{ic: 0.0921, lig: {'45': 124}}],
    [1.02,0.431,0,{ic: 0.0921}], [0.511,0.694,0,{ic: 0.122}],
    [0.511,0.668,0,{ic: 0.116}], [0.511,0.668,0,{ic: 0.105}]
  ],

  cmbx10: [
    [0.692,0.686], [0.958,0.686], [0.894,0.686], [0.806,0.686],
    [0.767,0.686], [0.9,0.686], [0.831,0.686], [0.894,0.686],
    [0.831,0.686], [0.894,0.686], [0.831,0.686],
    [0.671,0.694,0,{ic: 0.109, krn: {'39': 0.109, '63': 0.109, '33': 0.109, '41': 0.109, '93': 0.109}, lig: {'105': 14, '108': 15}}],
    [0.639,0.694], [0.639,0.694], [0.958,0.694], [0.958,0.694],

    [0.319,0.444], [0.351,0.444,0.194], [0.575,0.694], [0.575,0.694],
    [0.575,0.632], [0.575,0.694], [0.575,0.596], [0.869,0.694],
    [0.511,0,0.17], [0.597,0.694], [0.831,0.444], [0.894,0.444],
    [0.575,0.542,0.0972], [1.04,0.686], [1.17,0.686], [0.894,0.735,0.0486],

    [0.319,0.444,0,{krn: {'108': -0.319, '76': -0.378}}],
    [0.35,0.694,0,{lig: {'96': 60}}], [0.603,0.694], [0.958,0.694,0.194],
    [0.575,0.75,0.0556], [0.958,0.75,0.0556], [0.894,0.694],
    [0.319,0.694,0,{krn: {'63': 0.128, '33': 0.128}, lig: {'39': 34}}],
    [0.447,0.75,0.25], [0.447,0.75,0.25], [0.575,0.75], [0.894,0.633,0.133],
    [0.319,0.156,0.194], [0.383,0.444,0,{lig: {'45': 123}}],
    [0.319,0.156], [0.575,0.75,0.25],

    [0.575,0.644], [0.575,0.644], [0.575,0.644], [0.575,0.644],
    [0.575,0.644], [0.575,0.644], [0.575,0.644], [0.575,0.644],
    [0.575,0.644], [0.575,0.644], [0.319,0.444], [0.319,0.444,0.194],
    [0.35,0.5,0.194], [0.894,0.391,-0.109], [0.543,0.5,0.194],
    [0.543,0.694,0,{lig: {'96': 62}}],

    [0.894,0.694],
    [0.869,0.686,0,{krn: {'116': -0.0319, '67': -0.0319, '79': -0.0319, '71': -0.0319, '85': -0.0319, '81': -0.0319, '84': -0.0958, '89': -0.0958, '86': -0.128, '87': -0.128}}],
    [0.818,0.686], [0.831,0.686],
    [0.882,0.686,0,{krn: {'88': -0.0319, '87': -0.0319, '65': -0.0319, '86': -0.0319, '89': -0.0319}}],
    [0.756,0.686],
    [0.724,0.686,0,{krn: {'111': -0.0958, '101': -0.0958, '117': -0.0958, '114': -0.0958, '97': -0.0958, '65': -0.128, '79': -0.0319, '67': -0.0319, '71': -0.0319, '81': -0.0319}}],
    [0.904,0.686], [0.9,0.686], [0.436,0.686,0,{krn: {'73': 0.0319}}],
    [0.594,0.686],
    [0.901,0.686,0,{krn: {'79': -0.0319, '67': -0.0319, '71': -0.0319, '81': -0.0319}}],
    [0.692,0.686,0,{krn: {'84': -0.0958, '89': -0.0958, '86': -0.128, '87': -0.128}}],
    [1.09,0.686], [0.9,0.686],
    [0.864,0.686,0,{krn: {'88': -0.0319, '87': -0.0319, '65': -0.0319, '86': -0.0319, '89': -0.0319}}],

    [0.786,0.686,0,{krn: {'65': -0.0958, '111': -0.0319, '101': -0.0319, '97': -0.0319, '46': -0.0958, '44': -0.0958}}],
    [0.864,0.686,0.194],
    [0.862,0.686,0,{krn: {'116': -0.0319, '67': -0.0319, '79': -0.0319, '71': -0.0319, '85': -0.0319, '81': -0.0319, '84': -0.0958, '89': -0.0958, '86': -0.128, '87': -0.128}}],
    [0.639,0.686],
    [0.8,0.686,0,{krn: {'121': -0.0319, '101': -0.0958, '111': -0.0958, '114': -0.0958, '97': -0.0958, '65': -0.0958, '117': -0.0958}}],
    [0.885,0.686],
    [0.869,0.686,0,{ic: 0.016, krn: {'111': -0.0958, '101': -0.0958, '117': -0.0958, '114': -0.0958, '97': -0.0958, '65': -0.128, '79': -0.0319, '67': -0.0319, '71': -0.0319, '81': -0.0319}}],
    [1.19,0.686,0,{ic: 0.016, krn: {'111': -0.0958, '101': -0.0958, '117': -0.0958, '114': -0.0958, '97': -0.0958, '65': -0.128, '79': -0.0319, '67': -0.0319, '71': -0.0319, '81': -0.0319}}],
    [0.869,0.686,0,{krn: {'79': -0.0319, '67': -0.0319, '71': -0.0319, '81': -0.0319}}],
    [0.869,0.686,0,{ic: 0.0287, krn: {'101': -0.0958, '111': -0.0958, '114': -0.0958, '97': -0.0958, '65': -0.0958, '117': -0.0958}}],
    [0.703,0.686], [0.319,0.75,0.25], [0.603,0.694], [0.319,0.75,0.25],
    [0.575,0.694], [0.319,0.694],

    [0.319,0.694,0,{lig: {'96': 92}}],
    [0.559,0.444,0,{krn: {'118': -0.0319, '106': 0.0639, '121': -0.0319, '119': -0.0319}}],
    [0.639,0.694,0,{krn: {'101': 0.0319, '111': 0.0319, '120': -0.0319, '100': 0.0319, '99': 0.0319, '113': 0.0319, '118': -0.0319, '106': 0.0639, '121': -0.0319, '119': -0.0319}}],
    [0.511,0.444,0,{krn: {'104': -0.0319, '107': -0.0319}}],
    [0.639,0.694], [0.527,0.444],
    [0.351,0.694,0,{ic: 0.109, krn: {'39': 0.109, '63': 0.109, '33': 0.109, '41': 0.109, '93': 0.109}, lig: {'105': 12, '102': 11, '108': 13}}],
    [0.575,0.444,0.194,{ic: 0.016, krn: {'106': 0.0319}}],
    [0.639,0.694,0,{krn: {'116': -0.0319, '117': -0.0319, '98': -0.0319, '121': -0.0319, '118': -0.0319, '119': -0.0319}}],
    [0.319,0.694], [0.351,0.694,0.194],
    [0.607,0.694,0,{krn: {'97': -0.0639, '101': -0.0319, '97': -0.0319, '111': -0.0319, '99': -0.0319}}],
    [0.319,0.694],
    [0.958,0.444,0,{krn: {'116': -0.0319, '117': -0.0319, '98': -0.0319, '121': -0.0319, '118': -0.0319, '119': -0.0319}}],
    [0.639,0.444,0,{krn: {'116': -0.0319, '117': -0.0319, '98': -0.0319, '121': -0.0319, '118': -0.0319, '119': -0.0319}}],
    [0.575,0.444,0,{krn: {'101': 0.0319, '111': 0.0319, '120': -0.0319, '100': 0.0319, '99': 0.0319, '113': 0.0319, '118': -0.0319, '106': 0.0639, '121': -0.0319, '119': -0.0319}}],

    [0.639,0.444,0.194,{krn: {'101': 0.0319, '111': 0.0319, '120': -0.0319, '100': 0.0319, '99': 0.0319, '113': 0.0319, '118': -0.0319, '106': 0.0639, '121': -0.0319, '119': -0.0319}}],
    [0.607,0.444,0.194], [0.474,0.444], [0.454,0.444],
    [0.447,0.635,0,{krn: {'121': -0.0319, '119': -0.0319}}],
    [0.639,0.444,0,{krn: {'119': -0.0319}}],
    [0.607,0.444,0,{ic: 0.016, krn: {'97': -0.0639, '101': -0.0319, '97': -0.0319, '111': -0.0319, '99': -0.0319}}],
    [0.831,0.444,0,{ic: 0.016, krn: {'101': -0.0319, '97': -0.0319, '111': -0.0319, '99': -0.0319}}],
    [0.607,0.444],
    [0.607,0.444,0.194,{ic: 0.016, krn: {'111': -0.0319, '101': -0.0319, '97': -0.0319, '46': -0.0958, '44': -0.0958}}],
    [0.511,0.444], [0.575,0.444,0,{ic: 0.0319, lig: {'45': 124}}],
    [1.15,0.444,0,{ic: 0.0319}], [0.575,0.694], [0.575,0.694], [0.575,0.694]
  ]
};

/***************************************************************************/

/*
 *  Implement image-based fonts for fallback method
 */
jsMath.Img = {

  // font sizes available
  fonts: [50, 60, 70, 85, 100, 120, 144, 173, 207, 249, 298, 358, 430],

  // em widths for the various font size directories
  w: {'50': 6.9, '60': 8.3, '70': 9.7, '85': 11.8, '100': 13.9,
      '120': 16.7, '144': 20.0, '173': 24.0, '207': 28.8, '249': 34.6,
      '298': 41.4, '358': 49.8, '430': 59.8},

  best: 4,     // index of best font size in the fonts list
  update: {},  // fonts to update (see UpdateFonts below)
  factor: 1,   // factor by which to shrink images (for better printing)
  loaded: 0,   // image fonts are loaded

  // add characters to be drawn using images
  SetFont: function (change) {
    for (var font in change) {
      if (!this.update[font]) {this.update[font] = []}
      this.update[font] = this.update[font].concat(change[font]);
    }
  },

  /*
   *  Called by the exta-font definition files to add an image font
   *  into the mix
   */
  AddFont: function (size,def) {
    if (!jsMath.Img[size]) {jsMath.Img[size] = {}};
    jsMath.Add(jsMath.Img[size],def);
  },

  /*
   *  Update font(s) to use image data rather than native fonts
   *  It looks in the jsMath.Img.update array to find the names
   *  of the fonts to udpate, and the arrays of character codes
   *  to set (or 'all' to change every character);
   */
  UpdateFonts: function () {
    var change = this.update; if (!this.loaded) return;
    var best = this[jsMath.Img.fonts[this.best]];
    for (var font in change) {
      for (var i = 0; i < change[font].length; i++) {
        var c = change[font][i];
        if (c == 'all') {for (c in jsMath.TeX[font]) {jsMath.TeX[font][c].img = {}}}
          else {jsMath.TeX[font][c].img = {}}
      }
    }
    this.update = {};
  },

  /*
   *  Find the font size that best fits our current font
   *  (this is the directory name for the img files used
   *  in some fallback modes).
   */
  BestSize: function () {
    var w = jsMath.em * this.factor;
    var m = this.w[this.fonts[0]];
    for (var i = 1; i < this.fonts.length; i++) {
      if (w < (this.w[this.fonts[i]] + 2*m) / 3) {return i-1}
      m = this.w[this.fonts[i]];
    }
    return i-1;
  },

  /*
   *  Get the scaling factor for the image fonts
   */
  Scale: function () {
    if (!this.loaded) return;
    this.best = this.BestSize();
    this.em = jsMath.Img.w[this.fonts[this.best]];
    this.scale = (jsMath.em/this.em);
    if (Math.abs(this.scale - 1) < .12) {this.scale = 1}
  },

  /*
   *  Get URL to directory for given font and size, based on the
   *  user's alpha/plain setting
   */
  URL: function (name,size,C) {
    var type = (jsMath.Controls.cookie.alpha) ? '/alpha/': '/plain/';
    if (C == null) {C = "def.js"} else {C = 'char'+C+'.png'}
    if (size != "") {size += '/'}
    return this.root+name+type+size+C;
  },

  /*
   *  Laod the data for an image font
   */
  LoadFont: function (name) {
    if (!this.loaded) this.Init();
    jsMath.Setup.Script(this.URL(name,""));
  },

  /*
   *  Setup for print mode, and create the hex code table
   */
  Init: function () {
    if (jsMath.Controls.cookie.print || jsMath.Controls.cookie.stayhires) {
      jsMath.Controls.cookie.print = jsMath.Controls.cookie.stayhires;
      this.factor *= 3;
      if (!jsMath.Controls.isLocalCookie || !jsMath.Global.isLocal) {jsMath.Controls.SetCookie(0)}
      if (jsMath.Browser.alphaPrintBug) {jsMath.Controls.cookie.alpha = 0}
    }
    var codes = '0123456789ABCDEF';
    this.HexCode = [];
    for (var i = 0; i < 128; i++) {
      var h = Math.floor(i/16); var l = i - 16*h;
      this.HexCode[i] = codes.charAt(h)+codes.charAt(l);
    }
    this.loaded = 1;
  }

};

/***************************************************************************/

/*
 *  jsMath.HTML handles creation of most of the HTML needed for
 *  presenting mathematics in HTML pages.
 */

jsMath.HTML = {

  /*
   *  Produce a string version of a measurement in ems,
   *  showing only a limited number of digits, and
   *  using 0 when the value is near zero.
   */
  Em: function (m) {
    var n = 5; if (m < 0) {n++}
    if (Math.abs(m) < .000001) {m = 0}
    var s = String(m); s = s.replace(/(\.\d\d\d).+/,'$1');
    return s+'em'
  },

  /*
   *  Create a horizontal space of width w
   */
  Spacer: function (w) {
    if (w == 0) {return ''};
    return jsMath.Browser.msieSpaceFix+'<span class="spacer" style="margin-left:'+this.Em(w)+'"></span>';
  },

  /*
   *  Create a blank rectangle of the given size
   *  If the height is small, it is converted to pixels so that it
   *  will not disappear at small font sizes.
   */

  Blank: function (w,h,d,isRule) {
    var backspace = ''; var style = ''
    if (isRule) {
      style += 'border-left:'+this.Em(w)+' solid;';
      if (jsMath.Browser.widthAddsBorder) {w = 0};
    }
    if (w == 0) {
      if (jsMath.Browser.blankWidthBug) {
        style += 'width:1px;';
        backspace = '<span class="spacer" style="margin-right:-1px"></span>'
      }
    } else {style += 'width:'+this.Em(w)+';'}
    if (d == null) {d = 0}
    if (h) {
      var H = this.Em(h+d);
      if (isRule && h*jsMath.em < 1.5) {H = "1px"; h = 1/jsMath.em}
      style += 'height:'+H+';';
    }
    if (jsMath.Browser.mozInlineBlockBug) {d = -h}
    if (jsMath.Browser.msieBorderBug && !isRule) {d -= jsMath.d}
    if (d) {style += 'vertical-align:'+this.Em(-d)}
    return backspace+'<span class="blank" style="'+style+'"></span>';
  },

  /*
   *  Create a rule line for fractions, etc.
   */
  Rule: function (w,h) {
    if (h == null) {h = jsMath.TeX.default_rule_thickness}
    return this.Blank(w,h,0,1);
  },

  /*
   *  Add a <SPAN> tag to activate a specific CSS class
   */
  Class: function (tclass,html) {
    return '<span class="'+tclass+'">'+html+'</span>';
  },

  /*
   *  Use a <SPAN> to place some HTML at a specific position.
   *  (This can be replaced by the ones below to overcome
   *   some browser-specific bugs.)
   */
  Place: function (html,x,y) {
    if (Math.abs(x) < .0001) {x = 0}
    if (Math.abs(y) < .0001) {y = 0}
    if (x || y) {
      var span = '<span style="position: relative;';
      if (x) {span += ' margin-left:'+this.Em(x)+';'}
      if (y) {span += ' top:'+this.Em(-y)+';'}
      html = span + '">' + html + '</span>';
    }
    return html;
  },

  /*
   *  For MSIE on Windows, backspacing must be done in a separate
   *  <SPAN>, otherwise the contents will be clipped.  Netscape
   *  also doesn't combine vertical and horizontal spacing well.
   *  Here the x and y positioning are done in separate <SPAN> tags
   */
  PlaceSeparateSkips: function (html,x,y) {
    if (Math.abs(x) < .0001) {x = 0}
    if (Math.abs(y) < .0001) {y = 0}
    if (y) {html = '<span style="position: relative; top:'+this.Em(-y)+';'
                       + '">' + html + '</span>'}
    if (x) {html = this.Spacer(x) + html}
    return html;
  },

  /*
   *  Place a SPAN with absolute coordinates
   */
  PlaceAbsolute: function (html,x,y) {
    if (Math.abs(x) < .0001) {x = 0}
    if (Math.abs(y) < .0001) {y = 0}
    html = '<span style="position:absolute; left:'+this.Em(x)+'; '
              + 'top:'+this.Em(y)+';">' + html + '&nbsp;</span>';
              //  space normalizes line height in script styles
    return html;
  },

  Absolute: function(html,w,h,d,y,H) {
    if (y != "none") {
      if (Math.abs(y) < .0001) {y = 0}
      html = '<span style="position:absolute; '
               + 'top:'+jsMath.HTML.Em(y)+'; left:0em;">'
               + html + '&nbsp;' // space normalizes line height in script styles
             + '</span>';
    }
    if (d == "none") {d = 0}
    html += this.Blank(w,h-d,d);
    if (jsMath.Browser.msieAbsoluteBug) {           // for MSIE (Mac)
      html = '<span style="position:relative;">' + html + '</span>';
    }
    if (jsMath.Browser.spanHeightVaries) {
      var style = '';
      style = jsMath.Browser.msieInlineBlockFix
            + ' width:'+jsMath.HTML.Em(w)+';';
      if (jsMath.Browser.quirks) {
        style += ' height:'+jsMath.HTML.Em(H)+';'
      } else {
        style += ' height: 0px;'
               + ' vertical-align:'+jsMath.HTML.Em(H)+';'
      }
      html = '<span style="position:relative;' + style + '">'
           +   html
           + '</span>';
    } else {
      html = '<span style="position:relative">' + html + '</span>';
    }
    return html;
  }

};


/***************************************************************************/

/*
 *  jsMath.Box handles TeX's math boxes and jsMath's equivalent of hboxes.
 */

jsMath.Box = function (format,text,w,h,d) {
  if (d == null) {d = jsMath.d}
  this.type = 'typeset';
  this.w = w; this.h = h; this.d = d; this.bh = h; this.bd = d;
  this.x = 0; this.y = 0;
  this.html = text; this.format = format;
};


jsMath.Add(jsMath.Box,{

  defaultH: 0, // default height for characters with none specified

  /*
   *  An empty box
   */
  Null: function () {return new jsMath.Box('null','',0,0,0)},

  /*
   *  A box containing only text whose class and style haven't been added
   *  yet (so that we can combine ones with the same styles).  It gets
   *  the text dimensions, if needed.  (In general, this has been
   *  replaced by TeX() below, but is still used in fallback mode.)
   */
  Text: function (text,tclass,style,size,a,d) {
    var html = jsMath.Typeset.AddClass(tclass,text);
        html = jsMath.Typeset.AddStyle(style,size,html);
    var BB = jsMath.EmBoxFor(html); var TeX = jsMath.Typeset.TeX(style,size);
    var bd = ((tclass == 'cmsy10' || tclass == 'cmex10')? BB.h-TeX.h: TeX.d*BB.h/TeX.hd);
    var box = new jsMath.Box('text',text,BB.w,BB.h-bd,bd);
    box.style = style; box.size = size; box.tclass = tclass;
    if (d != null) {box.d = d*TeX.scale} else {box.d = 0}
    if (a == null || a == 1) {box.h = .9*TeX.M_height}
      else {box.h = 1.1*TeX.x_height + TeX.scale*a}
    return box;
  },

  /*
   *  Produce a box containing a given TeX character from a given font.
   *  The box is a text box (like the ones above), so that characters from
   *  the same font can be combined.
   */
  TeX: function (C,font,style,size) {
    var c = jsMath.TeX[font][C];
    if (c.d == null) {c.d = 0}; if (c.h == null) {c.h = 0}
    if (c.img != null && c.c != '') this.TeXIMG(font,C,jsMath.Typeset.StyleSize(style,size));
    var scale = jsMath.Typeset.TeX(style,size).scale;
    var h = c.h + jsMath.TeX[font].dh
    var box = new jsMath.Box('text',c.c,c.w*scale,h*scale,c.d*scale);
    box.style = style; box.size = size;
    if (c.tclass) {
      box.tclass = c.tclass;
      if (c.img) {box.bh = c.img.bh; box.bd = c.img.bd}
            else {box.bh = scale*jsMath.h; box.bd = scale*jsMath.d}
    } else {
      box.tclass = font;
      box.bh = scale*jsMath.TeX[font].h;
      box.bd = scale*jsMath.TeX[font].d;
      if (jsMath.Browser.msieFontBug && box.html.match(/&#/)) {
        // hack to avoid font changing back to the default
        // font when a unicode reference is not followed
        // by a letter or number
        box.html += '<span style="display:none">x</span>';
      }
    }
    return box;
  },

  /*
   *  In fallback modes, handle the fact that we don't have the
   *  sizes of the characters precomputed
   */
  TeXfallback: function (C,font,style,size) {
    var c = jsMath.TeX[font][C]; if (!c.tclass) {c.tclass = font}
    if (c.img != null) {return this.TeXnonfallback(C,font,style,size)}
    if (c.h != null && c.a == null) {c.a = c.h-1.1*jsMath.TeX.x_height}
    var a = c.a; var d = c.d; // avoid Firefox warnings
    var box = this.Text(c.c,c.tclass,style,size,a,d);
    var scale = jsMath.Typeset.TeX(style,size).scale;
    if (c.bh != null) {
      box.bh = c.bh*scale;
      box.bd = c.bd*scale;
    } else {
      var h = box.bd+box.bh;
      var html = jsMath.Typeset.AddClass(box.tclass,box.html);
          html = jsMath.Typeset.AddStyle(style,size,html);
      box.bd = jsMath.EmBoxFor(html+jsMath.HTML.Blank(1,h)).h - h;
      box.bh = h - box.bd;
      if (scale == 1) {c.bh = box.bh; c.bd = box.bd}
    }
    if (jsMath.msieFontBug && box.html.match(/&#/))
      {box.html += '<span style="display:none">x</span>'}
    return box;
  },

  /*
   *  Set the character's string to the appropriate image file
   */
  TeXIMG: function (font,C,size) {
    var c = jsMath.TeX[font][C];
    if (c.img.size != null && c.img.size == size &&
        c.img.best != null && c.img.best == jsMath.Img.best) return;
    var mustScale = (jsMath.Img.scale != 1);
    var id = jsMath.Img.best + size - 4;
    if (id < 0) {id = 0; mustScale = 1} else
    if (id >= jsMath.Img.fonts.length) {id = jsMath.Img.fonts.length-1; mustScale = 1}
    var imgFont = jsMath.Img[jsMath.Img.fonts[id]];
    var img = imgFont[font][C];
    var scale = 1/jsMath.Img.w[jsMath.Img.fonts[id]];
    if (id != jsMath.Img.best + size - 4) {
      if (c.w != null) {scale = c.w/img[0]} else {
        scale *= jsMath.Img.fonts[size]/jsMath.Img.fonts[4]
              *  jsMath.Img.fonts[jsMath.Img.best]/jsMath.Img.fonts[id];
      }
    }
    var w = img[0]*scale; var h = img[1]*scale; var d = -img[2]*scale; var v;
    var wadjust = (c.w == null || Math.abs(c.w-w) < .01)? "" : " margin-right:"+jsMath.HTML.Em(c.w-w)+';';
    var resize = ""; C = jsMath.Img.HexCode[C];
    if (!mustScale && !jsMath.Controls.cookie.scaleImg) {
      if (2*w < h || (jsMath.Browser.msieAlphaBug && jsMath.Controls.cookie.alpha))
         {resize = "height:"+(img[1]*jsMath.Browser.imgScale)+'px;'}
      resize += " width:"+(img[0]*jsMath.Browser.imgScale)+'px;'
      v = -img[2]+'px';
    } else {
      if (2*w < h || (jsMath.Browser.msieAlphaBug && jsMath.Controls.cookie.alpha))
         {resize = "height:"+jsMath.HTML.Em(h*jsMath.Browser.imgScale)+';'}
      resize += " width:"+jsMath.HTML.Em(w*jsMath.Browser.imgScale)+';'
      v = jsMath.HTML.Em(d);
    }
    var vadjust = (Math.abs(d) < .01 && !jsMath.Browser.valignBug)?
                         "": " vertical-align:"+v+';';
    var URL = jsMath.Img.URL(font,jsMath.Img.fonts[id],C);
    if (jsMath.Browser.msieAlphaBug && jsMath.Controls.cookie.alpha) {
      c.c = '<img src="'+jsMath.blank+'" '
               + 'style="'+jsMath.Browser.msieCenterBugFix
               + resize + vadjust + wadjust
               + ' filter:progid:DXImageTransform.Microsoft.AlphaImageLoader(src=' + "'"
               + URL + "', sizingMethod='scale'" + ');" />';
    } else {
      c.c = '<img src="'+URL+'" style="'+jsMath.Browser.msieCenterBugFix
                  + resize + vadjust + wadjust + '" />';
    }
    c.tclass = "normal";
    c.img.bh = h+d; c.img.bd = -d;
    c.img.size = size; c.img.best = jsMath.Img.best;
  },

  /*
   *  A box containing a spacer of a specific width
   */
  Space: function (w) {
    return new jsMath.Box('html',jsMath.HTML.Spacer(w),w,0,0);
  },

  /*
   *  A box containing a horizontal rule
   */
  Rule: function (w,h) {
    if (h == null) {h = jsMath.TeX.default_rule_thickness}
    var html = jsMath.HTML.Rule(w,h);
    return new jsMath.Box('html',html,w,h,0);
  },

  /*
   *  Get a character from a TeX font, and make sure that it has
   *  its metrics specified.
   */
  GetChar: function (code,font) {
    var c = jsMath.TeX[font][code];
    if (c.img != null) {this.TeXIMG(font,code,4)}
    if (c.tclass == null) {c.tclass = font}
    if (!c.computedW) {
      c.w = jsMath.EmBoxFor(jsMath.Typeset.AddClass(c.tclass,c.c)).w;
      if (c.h == null) {c.h = jsMath.Box.defaultH}; if (c.d == null) {c.d = 0}
      c.computedW = 1;
    }
    return c;
  },

  /*
   *  Locate the TeX delimiter character that matches a given height.
   *  Return the character, font, style and actual height used.
   */
  DelimBestFit: function (H,c,font,style) {
    if (c == 0 && font == 0) return null;
    var C; var h; font = jsMath.TeX.fam[font];
    var isSS = (style.charAt(1) == 'S');
    var isS  = (style.charAt(0) == 'S');
    while (c != null) {
      C = jsMath.TeX[font][c];
      if (C.h == null) {C.h = jsMath.Box.defaultH}; if (C.d == null) {C.d = 0}
      h = C.h+C.d;
      if (C.delim) {return [c,font,'',H]}
      if (isSS && .5*h >= H) {return [c,font,'SS',.5*h]}
      if (isS  && .7*h >= H) {return [c,font,'S',.7*h]}
      if (h >= H || C.n == null) {return [c,font,'T',h]}
      c = C.n;
    }
    return null;
  },

  /*
   *  Create the HTML needed for a stretchable delimiter of a given height,
   *  either centered or not.  This version uses relative placement (i.e.,
   *  backspaces, not line-breaks).  This works with more browsers, but
   *  if the font size changes, the backspacing may not be right, so the
   *  delimiters may become jagged.
   */
  DelimExtendRelative: function (H,c,font,a,nocenter) {
    var C = jsMath.TeX[font][c];
    var top = this.GetChar(C.delim.top? C.delim.top: C.delim.rep,font);
    var rep = this.GetChar(C.delim.rep,font);
    var bot = this.GetChar(C.delim.bot? C.delim.bot: C.delim.rep,font);
    var ext = jsMath.Typeset.AddClass(rep.tclass,rep.c);
    var w = rep.w; var h = rep.h+rep.d
    var y; var dx;
    if (C.delim.mid) {// braces
      var mid = this.GetChar(C.delim.mid,font);
      var n = Math.ceil((H-(top.h+top.d)-(mid.h+mid.d)-(bot.h+bot.d))/(2*(rep.h+rep.d)));
      H = 2*n*(rep.h+rep.d) + (top.h+top.d) + (mid.h+mid.d) + (bot.h+bot.d);
      if (nocenter) {y = 0} else {y = H/2+a}; var Y = y;
      var html = jsMath.HTML.Place(jsMath.Typeset.AddClass(top.tclass,top.c),0,y-top.h)
               + jsMath.HTML.Place(jsMath.Typeset.AddClass(bot.tclass,bot.c),-(top.w+bot.w)/2,y-(H-bot.d))
               + jsMath.HTML.Place(jsMath.Typeset.AddClass(mid.tclass,mid.c),-(bot.w+mid.w)/2,y-(H+mid.h-mid.d)/2);
      dx = (w-mid.w)/2; if (Math.abs(dx) < .0001) {dx = 0}
      if (dx) {html += jsMath.HTML.Spacer(dx)}
      y -= top.h+top.d + rep.h;
      for (var i = 0; i < n; i++) {html += jsMath.HTML.Place(ext,-w,y-i*h)}
      y -= H/2 - rep.h/2;
      for (var i = 0; i < n; i++) {html += jsMath.HTML.Place(ext,-w,y-i*h)}
    } else {// everything else
      var n = Math.ceil((H - (top.h+top.d) - (bot.h+bot.d))/(rep.h+rep.d));
      // make sure two-headed arrows have an extender
      if (top.h+top.d < .9*(rep.h+rep.d)) {n = Math.max(1,n)}
      H = n*(rep.h+rep.d) + (top.h+top.d) + (bot.h+bot.d);
      if (nocenter) {y = 0} else {y = H/2+a}; var Y = y;
      var html = jsMath.HTML.Place(jsMath.Typeset.AddClass(top.tclass,top.c),0,y-top.h)
      dx = (w-top.w)/2; if (Math.abs(dx) < .0001) {dx = 0}
      if (dx) {html += jsMath.HTML.Spacer(dx)}
      y -= top.h+top.d + rep.h;
      for (var i = 0; i < n; i++) {html += jsMath.HTML.Place(ext,-w,y-i*h)}
      html += jsMath.HTML.Place(jsMath.Typeset.AddClass(bot.tclass,bot.c),-(w+bot.w)/2,Y-(H-bot.d));
    }
    if (nocenter) {h = top.h} else {h = H/2+a}
    var box = new jsMath.Box('html',html,rep.w,h,H-h);
    box.bh = jsMath.TeX[font].h; box.bd = jsMath.TeX[font].d;
    return box;
  },

  /*
   *  Create the HTML needed for a stretchable delimiter of a given height,
   *  either centered or not.  This version uses absolute placement (i.e.,
   *  line-breaks, not backspacing).  This gives more reliable results,
   *  but doesn't work with all browsers.
   */
  DelimExtendAbsolute: function (H,c,font,a,nocenter) {
    var Font = jsMath.TeX[font];
    var C = Font[c]; var html;
    var top = this.GetChar(C.delim.top? C.delim.top: C.delim.rep,font);
    var rep = this.GetChar(C.delim.rep,font);
    var bot = this.GetChar(C.delim.bot? C.delim.bot: C.delim.rep,font);

    if (C.delim.mid) {// braces
      var mid = this.GetChar(C.delim.mid,font);
      var n = Math.ceil((H-(top.h+top.d)-(mid.h+mid.d-.05)-(bot.h+bot.d-.05))/(2*(rep.h+rep.d-.05)));
      H = 2*n*(rep.h+rep.d-.05) + (top.h+top.d) + (mid.h+mid.d-.05) + (bot.h+bot.d-.05);

      html = jsMath.HTML.PlaceAbsolute(jsMath.Typeset.AddClass(top.tclass,top.c),0,0);
      var h = rep.h+rep.d - .05; var y = top.d-.05 + rep.h;
      var ext = jsMath.Typeset.AddClass(font,rep.c)
      for (var i = 0; i < n; i++) {html += jsMath.HTML.PlaceAbsolute(ext,0,y+i*h)}
      html += jsMath.HTML.PlaceAbsolute(jsMath.Typeset.AddClass(mid.tclass,mid.c),0,y+n*h-rep.h+mid.h);
      y += n*h + mid.h+mid.d - .05;
      for (var i = 0; i < n; i++) {html += jsMath.HTML.PlaceAbsolute(ext,0,y+i*h)}
      html += jsMath.HTML.PlaceAbsolute(jsMath.Typeset.AddClass(bot.tclass,bot.c),0,y+n*h-rep.h+bot.h);
    } else {// all others
      var n = Math.ceil((H - (top.h+top.d) - (bot.h+bot.d-.05))/(rep.h+rep.d-.05));
      H = n*(rep.h+rep.d-.05) + (top.h+top.d) + (bot.h+bot.d-.05);

      html = jsMath.HTML.PlaceAbsolute(jsMath.Typeset.AddClass(top.tclass,top.c),0,0);
      var h = rep.h+rep.d-.05; var y = top.d-.05 + rep.h;
      var ext = jsMath.Typeset.AddClass(rep.tclass,rep.c);
      for (var i = 0; i < n; i++) {html += jsMath.HTML.PlaceAbsolute(ext,0,y+i*h)}
      html += jsMath.HTML.PlaceAbsolute(jsMath.Typeset.AddClass(bot.tclass,bot.c),0,y+n*h-rep.h+bot.h);
    }

    var w = top.w;
    if (nocenter) {h = top.h; y = 0} else {h = H/2 + a; y = h - top.h}
//    html = jsMath.HTML.Absolute(html,w,Font.h,"none",-y,top.h);
    html = jsMath.HTML.Absolute(html,w,Font.h,"none",-y,jsMath.h);
    var box = new jsMath.Box('html',html,rep.w,h,H-h);
    box.bh = jsMath.TeX[font].h; box.bd = jsMath.TeX[font].d;
    return box;
  },

  /*
   *  Get the HTML for a given delimiter of a given height.
   *  It will return either a single character, if one exists, or the
   *  more complex HTML needed for a stretchable delimiter.
   */
  Delimiter: function (H,delim,style,nocenter) {
    var size = 4;  //### pass this?
    var TeX = jsMath.Typeset.TeX(style,size);
    if (!delim) {return this.Space(TeX.nulldelimiterspace)}
    var CFSH = this.DelimBestFit(H,delim[2],delim[1],style);
    if (CFSH == null || CFSH[3] < H)
      {CFSH = this.DelimBestFit(H,delim[4],delim[3],style)}
    if (CFSH == null) {return this.Space(TeX.nulldelimiterspace)}
    if (CFSH[2] == '')
      {return this.DelimExtend(H,CFSH[0],CFSH[1],TeX.axis_height,nocenter)}
    var box = jsMath.Box.TeX(CFSH[0],CFSH[1],CFSH[2],size).Styled();
    if (nocenter) {box.y = -jsMath.TeX[CFSH[1]].dh*TeX.scale}
      else {box.y = -((box.h+box.d)/2 - box.d - TeX.axis_height)}
    if (Math.abs(box.y) < .0001) {box.y = 0}
    if (box.y) {box = jsMath.Box.SetList([box],CFSH[2],size)}
    return box;
  },

  /*
   *  Get a character by its TeX charcode, and make sure its width
   *  is specified.
   */
  GetCharCode: function (code) {
    var font = jsMath.TeX.fam[code[0]];
    var Font = jsMath.TeX[font];
    var c = Font[code[1]];
    if (c.img != null) {this.TeXIMG(font,code[1],4)}
    if (c.w == null) {c.w = jsMath.EmBoxFor(jsMath.Typeset.AddClass(c.tclass,c.c)).w}
    if (c.font == null) {c.font = font}
    return c;
  },

  /*
   * Add the class to the html, and use the font if there isn't one
   * specified already
   */

  AddClass: function (tclass,html,font) {
    if (tclass == null) {tclass = font}
    return jsMath.Typeset.AddClass(tclass,html);
  },

  /*
   *  Create the HTML for an alignment (e.g., array or matrix)
   *  Since the widths are not really accurate (they are based on pixel
   *  widths not the sub-pixel widths of the actual characters), there
   *  is some drift involved.  We lay out the table column by column
   *  to help reduce the problem.
   *
   *  ###  still need to allow users to specify row and column attributes,
   *       and do things like \span and \multispan  ###
   */
  LayoutRelative: function (size,table,align,cspacing,rspacing,vspace,useStrut,addWidth) {
    if (align == null) {align = []}
    if (cspacing == null) {cspacing = []}
    if (rspacing == null) {rspacing = []}
    if (useStrut == null) {useStrut = 1}
    if (addWidth == null) {addWidth = 1}

    // get row and column maximum dimensions
    var scale = jsMath.sizes[size]/100;
    var W = []; var H = []; var D = [];
    var unset = -1000; var bh = unset; var bd = unset;
    var i; var j; var row;
    for (i = 0; i < table.length; i++) {
      if (rspacing[i] == null) {rspacing[i] = 0}
      row = table[i];
      H[i] = useStrut*jsMath.h*scale; D[i] = useStrut*jsMath.d*scale;
      for (j = 0; j < row.length; j++) {
        row[j] = row[j].Remeasured();
        if (row[j].h > H[i]) {H[i] = row[j].h}
        if (row[j].d > D[i]) {D[i] = row[j].d}
        if (j >= W.length) {W[j] = row[j].w}
        else if (row[j].w > W[j]) {W[j] = row[j].w}
        if (row[j].bh > bh) {bh = row[j].bh}
        if (row[j].bd > bd) {bd = row[j].bd}
      }
    }
    if (rspacing[table.length] == null) {rspacing[table.length] = 0}
    if (bh == unset) {bh = 0}; if (bd == unset) {bd = 0}

    // lay out the columns
    var HD = useStrut*(jsMath.hd-.01)*scale;
    var dy = (vspace || 1) * scale/6;
    var html = ''; var pW = 0; var cW = 0;
    var w; var h; var y;
    var box; var mlist; var entry;
    for (j = 0; j < W.length; j++) {
      mlist = []; y = -H[0]-rspacing[0]; pW = 0;
      for (i = 0; i < table.length; i++) {
        entry = table[i][j];
        if (entry && entry.format != 'null') {
          if (align[j] == 'l') {w = 0} else
          if (align[j] == 'r') {w = W[j] - entry.w} else
            {w = (W[j] - entry.w)/2}
          entry.x = w - pW; pW = entry.w + w; entry.y = y;
          mlist[mlist.length] = entry;
        }
        if (i+1 < table.length) {y -= Math.max(HD,D[i]+H[i+1]) + dy + rspacing[i+1]}
      }
      if (cspacing[j] == null) cspacing[j] = scale;
      if (mlist.length > 0) {
        box = jsMath.Box.SetList(mlist,'T',size);
        html += jsMath.HTML.Place(box.html,cW,0);
        cW = W[j] - box.w + cspacing[j];
      } else {cW += cspacing[j]}
    }

    // get the full width and height
    w = -cspacing[W.length-1]; y = (H.length-1)*dy + rspacing[0];
    for (i = 0; i < W.length; i++) {w += W[i] + cspacing[i]}
    for (i = 0; i < H.length; i++) {y += Math.max(HD,H[i]+D[i]) + rspacing[i+1]}
    h = y/2 + jsMath.TeX.axis_height; var d = y-h;

    // adjust the final row width, and vcenter the table
    //   (add 1/6em at each side for the \,)
    html += jsMath.HTML.Spacer(cW-cspacing[W.length-1] + addWidth*scale/6);
    html = jsMath.HTML.Place(html,addWidth*scale/6,h);
    box = new jsMath.Box('html',html,w+addWidth*scale/3,h,d);
    box.bh = bh; box.bd = bd;
    return box;
  },

  /*
   *  Create the HTML for an alignment (e.g., array or matrix)
   *  Use absolute position for elements in the array.
   *
   *  ###  still need to allow users to specify row and column attributes,
   *       and do things like \span and \multispan  ###
   */
  LayoutAbsolute: function (size,table,align,cspacing,rspacing,vspace,useStrut,addWidth) {
    if (align == null) {align = []}
    if (cspacing == null) {cspacing = []}
    if (rspacing == null) {rspacing = []}
    if (useStrut == null) {useStrut = 1}
    if (addWidth == null) {addWidth = 1}

    // get row and column maximum dimensions
    var scale = jsMath.sizes[size]/100;
    var HD = useStrut*(jsMath.hd-.01)*scale;
    var dy = (vspace || 1) * scale/6;
    var W = []; var H = []; var D = [];
    var w = 0; var h; var x; var y;
    var i; var j; var row;
    for (i = 0; i < table.length; i++) {
      if (rspacing[i] == null) {rspacing[i] = 0}
      row = table[i];
      H[i] = useStrut*jsMath.h*scale; D[i] = useStrut*jsMath.d*scale;
      for (j = 0; j < row.length; j++) {
        row[j] = row[j].Remeasured();
        if (row[j].h > H[i]) {H[i] = row[j].h}
        if (row[j].d > D[i]) {D[i] = row[j].d}
        if (j >= W.length) {W[j] = row[j].w}
        else if (row[j].w > W[j]) {W[j] = row[j].w}
      }
    }
    if (rspacing[table.length] == null) {rspacing[table.length] = 0}

    // get the height and depth of the centered table
    y = (H.length-1)*dy + rspacing[0];
    for (i = 0; i < H.length; i++) {y += Math.max(HD,H[i]+D[i]) + rspacing[i+1]}
    h = y/2 + jsMath.TeX.axis_height; var d = y - h;

    // lay out the columns
    var html = ''; var entry; w = addWidth*scale/6;
    for (j = 0; j < W.length; j++) {
      y = H[0]-h + rspacing[0];
      for (i = 0; i < table.length; i++) {
        entry = table[i][j];
        if (entry && entry.format != 'null') {
          if (align[j] && align[j] == 'l') {x = 0} else
          if (align[j] && align[j] == 'r') {x = W[j] - entry.w} else
            {x = (W[j] - entry.w)/2}
          html += jsMath.HTML.PlaceAbsolute(entry.html,w+x,
                    y-Math.max(0,entry.bh-jsMath.h*scale));
        }
        if (i+1 < table.length) {y += Math.max(HD,D[i]+H[i+1]) + dy + rspacing[i+1]}
      }
      if (cspacing[j] == null) cspacing[j] = scale;
      w += W[j] + cspacing[j];
    }

    // get the full width
    w = -cspacing[W.length-1]+addWidth*scale/3;
    for (i = 0; i < W.length; i++) {w += W[i] + cspacing[i]}

    html = jsMath.HTML.Spacer(addWidth*scale/6)+html+jsMath.HTML.Spacer(addWidth*scale/6);
    if (jsMath.Browser.spanHeightVaries) {y = h-jsMath.h} else {y = 0}
    html = jsMath.HTML.Absolute(html,w,h+d,d,y,h);
    var box = new jsMath.Box('html',html,w+addWidth*scale/3,h,d);
    return box;
  },

  /*
   *  Look for math within \hbox and other non-math text
   */
  InternalMath: function (text,size) {
    text = text.replace(/@\(([^)]*)\)/g,'<$1>');
    if (!text.match(/\$|\\\(/)) {return this.Text(text,'normal','T',size).Styled()}

    var i = 0; var k = 0; var c; var match = '';
    var mlist = []; var parse; var html; var box;
    while (i < text.length) {
      c = text.charAt(i++);
      if (c == '$') {
        if (match == '$') {
          parse = jsMath.Parse(text.slice(k,i-1),null,size);
          if (parse.error) {
            mlist[mlist.length] = this.Text(parse.error,'error','T',size,1,.2);
          } else {
            parse.Atomize();
            mlist[mlist.length] = parse.mlist.Typeset('T',size).Styled();
          }
          match = ''; k = i;
        } else {
          mlist[mlist.length] = this.Text(text.slice(k,i-1),'normal','T',size,1,.2);
          match = '$'; k = i;
        }
      } else if (c == '\\') {
        c = text.charAt(i++);
        if (c == '(' && match == '') {
          mlist[mlist.length] = this.Text(text.slice(k,i-2),'normal','T',size,1,.2);
          match = ')'; k = i;
        } else if (c == ')' && match == ')') {
          parse = jsMath.Parse(text.slice(k,i-2),null,size);
          if (parse.error) {
            mlist[mlist.length] = this.Text(parse.error,'error','T',size,1,.2);
          } else {
            parse.Atomize();
            mlist[mlist.length] = parse.mlist.Typeset('T',size).Styled();
          }
          match = ''; k = i;
        }
      }
    }
    mlist[mlist.length] = this.Text(text.slice(k),'normal','T',size,1,.2);
    return this.SetList(mlist,'T',size);
  },

  /*
   *  Convert an abitrary box to a typeset box.  I.e., make an
   *  HTML version of the contents of the box, at its desired (x,y)
   *  position.
   */
  Set: function (box,style,size,addstyle) {
    if (box && box.type) {
      if (box.type == 'typeset') {return box}
      if (box.type == 'mlist') {
        box.mlist.Atomize(style,size);
        return box.mlist.Typeset(style,size);
      }
      if (box.type == 'text') {
        box = this.Text(box.text,box.tclass,style,size,box.ascend||null,box.descend||null);
        if (addstyle != 0) {box.Styled()}
        return box;
      }
      box = this.TeX(box.c,box.font,style,size);
      if (addstyle != 0) {box.Styled()}
      return box;
    }
    return jsMath.Box.Null();
  },

  /*
   *  Convert a list of boxes to a single typeset box.  I.e., finalize
   *  the HTML for the list of boxes, properly spaced and positioned.
   */
  SetList: function (boxes,style,size) {
    var mlist = []; var box;
    for (var i = 0; i < boxes.length; i++) {
      box = boxes[i];
      if (box.type == 'typeset') {box = jsMath.mItem.Typeset(box)}
      mlist[mlist.length] = box;
    }
    var typeset = new jsMath.Typeset(mlist);
    return typeset.Typeset(style,size);
  }

});


jsMath.Package(jsMath.Box,{

  /*
   *  Add the class and style to a text box (i.e., finalize the
   *  unpositioned HTML for the box).
   */
  Styled: function () {
    if (this.format == 'text') {
      this.html = jsMath.Typeset.AddClass(this.tclass,this.html);
      this.html = jsMath.Typeset.AddStyle(this.style,this.size,this.html);
      delete this.tclass; delete this.style;
      this.format = 'html';
    }
    return this;
  },

  /*
   *  Recompute the box width to make it more accurate.
   */
  Remeasured: function () {
    if (this.w > 0) {this.w = jsMath.EmBoxFor(this.html).w}
    return this;
  }

});


/***************************************************************************/

/*
 *  mItems are the building blocks of mLists (math lists) used to
 *  store the information about a mathematical expression.  These are
 *  basically the items listed in the TeXbook in Appendix G (plus some
 *  minor extensions).
 */
jsMath.mItem = function (type,def) {
  this.type = type;
  jsMath.Add(this,def);
}

jsMath.Add(jsMath.mItem,{

  /*
   *  A general atom (given a nucleus for the atom)
   */
  Atom: function (type,nucleus) {
    return new jsMath.mItem(type,{atom: 1, nuc: nucleus});
  },

  /*
   *  An atom whose nucleus is a piece of text, in a given
   *  class, with a given additional height and depth
   */
  TextAtom: function (type,text,tclass,a,d) {
    var atom = new jsMath.mItem(type,{
      atom: 1,
      nuc: {
        type: 'text',
        text: text,
        tclass: tclass
      }
    });
    if (a != null) {atom.nuc.ascend = a}
    if (d != null) {atom.nuc.descend = d}
    return atom;
  },

  /*
   *  An atom whose nucleus is a TeX character in a specific font
   */
  TeXAtom: function (type,c,font) {
    return new jsMath.mItem(type,{
      atom: 1,
      nuc: {
        type: 'TeX',
        c: c,
        font: font
      }
    });
  },

  /*
   *  A generalized fraction atom, with given delimiters, rule
   *  thickness, and a numerator and denominator.
   */
  Fraction: function (name,num,den,thickness,left,right) {
    return new jsMath.mItem('fraction',{
      from: name, num: num, den: den,
      thickness: thickness, left: left, right: right
    });
  },

  /*
   *  An atom that inserts some glue
   */
  Space: function (w) {return new jsMath.mItem('space',{w: w})},

  /*
   *  An atom that contains a typeset box (like an hbox or vbox)
   */
  Typeset: function (box) {return new jsMath.mItem('ord',{atom:1, nuc: box})},

  /*
   *  An atom that contains some finished HTML (acts like a typeset box)
   */
  HTML: function (html) {return new jsMath.mItem('html',{html: html})}

});

/***************************************************************************/

/*
 *  mLists are lists of mItems, and encode the contents of
 *  mathematical expressions and sub-expressions.  They act as
 *  the expression "stack" as the mathematics is parsed, and
 *  contain some state information, like the position of the
 *  most recent open paren and \over command, and the current font.
 */
jsMath.mList = function (list,font,size,style) {
  if (list) {this.mlist = list} else {this.mlist = []}
  if (style == null) {style = 'T'}; if (size == null) {size = 4}
  this.data = {openI: null, overI: null, overF: null,
               font: font, size: size, style: style};
  this.init = {size: size, style: style};
}

jsMath.Package(jsMath.mList,{

  /*
   *  Add an mItem to the list
   */
  Add: function (box) {return (this.mlist[this.mlist.length] = box)},

  /*
   *  Get the i-th mItem from the list
   */
  Get: function (i) {return this.mlist[i]},

  /*
   *  Get the length of the list
   */
  Length: function() {return this.mlist.length},

  /*
   *  Get the tail mItem of the list
   */
  Last: function () {
    if (this.mlist.length == 0) {return null}
    return this.mlist[this.mlist.length-1]
  },

  /*
   *  Get a sublist of an mList
   */
  Range: function (i,j) {
    if (j == null) {j = this.mlist.length}
    return new jsMath.mList(this.mlist.slice(i,j+1));
  },

  /*
   *  Remove a range of mItems from the list.
   */
  Delete: function (i,j) {
    if (j == null) {j = i}
    if (this.mlist.splice) {this.mlist.splice(i,j-i+1)} else {
      var mlist = [];
      for (var k = 0; k < this.mlist.length; k++)
        {if (k < i || k > j) {mlist[mlist.length] = this.mlist[k]}}
      this.mlist = mlist;
    }
  },

  /*
   *  Add an open brace and maintain the stack information
   *  about the previous open brace so we can recover it
   *  when this one os closed.
   */
  Open: function (left) {
    var box = this.Add(new jsMath.mItem('boundary',{data: this.data}));
    var olddata = this.data;
    this.data = {}; for (var i in olddata) {this.data[i] = olddata[i]}
    delete this.data.overI; delete this.data.overF;
    this.data.openI = this.mlist.length-1;
    if (left != null) {box.left = left}
    return box;
  },

  /*
   *  Attempt to close a brace.  Recover the stack information
   *  about previous open braces and \over commands.  If there was an
   *  \over (or \above, etc) in this set of braces, create a fraction
   *  atom from the two halves, otherwise create an inner or ord
   *  from the contents of the braces.
   *  Remove the braced material from the list and add the newly
   *  created atom (the fraction, inner or ord).
   */
  Close: function (right) {
    if (right != null) {right = new jsMath.mItem('boundary',{right: right})}
    var atom; var open = this.data.openI;
    var over = this.data.overI; var from = this.data.overF;
    this.data  = this.mlist[open].data;
    if (over) {
      atom = jsMath.mItem.Fraction(from.name,
        {type: 'mlist', mlist: this.Range(open+1,over-1)},
        {type: 'mlist', mlist: this.Range(over)},
        from.thickness,from.left,from.right);
      if (right) {
        var mlist = new jsMath.mList([this.mlist[open],atom,right]);
        atom = jsMath.mItem.Atom('inner',{type: 'mlist', mlist: mlist});
      }
    } else {
      var openI = open+1; if (right) {this.Add(right); openI--}
      atom = jsMath.mItem.Atom((right)?'inner':'ord',
                  {type: 'mlist', mlist: this.Range(openI)});
    }
    this.Delete(open,this.Length());
    return this.Add(atom);
  },

  /*
   *  Create a generalized fraction from an mlist that
   *  contains an \over (or \above, etc).
   */
  Over: function () {
    var over = this.data.overI; var from = this.data.overF;
    var atom = jsMath.mItem.Fraction(from.name,
      {type: 'mlist', mlist: this.Range(open+1,over-1)},
      {type: 'mlist', mlist: this.Range(over)},
      from.thickness,from.left,from.right);
    this.mlist = [atom];
  },

  /*
   *  Take a raw mList (that has been produced by parsing some TeX
   *  expression), and perform the modifications outlined in
   *  Appendix G of the TeXbook.
   */
  Atomize: function (style,size) {
    var mitem; var prev = '';
    this.style = style; this.size = size;
    for (var i = 0; i < this.mlist.length; i++) {
      mitem = this.mlist[i]; mitem.delta = 0;
      if (mitem.type == 'choice')
        {this.mlist = this.Atomize.choice(this.style,mitem,i,this.mlist); i--}
      else if (this.Atomize[mitem.type]) {
        var f = this.Atomize[mitem.type]; // Opera needs separate name
        f(this.style,this.size,mitem,prev,this,i);
      }
      prev = mitem;
    }
    if (mitem && mitem.type == 'bin') {mitem.type = 'ord'}
    if (this.mlist.length >= 2 && mitem.type == 'boundary' &&
        this.mlist[0].type == 'boundary') {this.AddDelimiters(style,size)}
  },

  /*
   *  For a list that has boundary delimiters as its first and last
   *  entries, we replace the boundary atoms by open and close
   *  atoms whose nuclii are the specified delimiters properly sized
   *  for the contents of the list.  (Rule 19)
   */
  AddDelimiters: function(style,size) {
    var unset = -10000; var h = unset; var d = unset;
    for (var i = 0; i < this.mlist.length; i++) {
      var mitem = this.mlist[i];
      if (mitem.atom || mitem.type == 'box') {
        h = Math.max(h,mitem.nuc.h+mitem.nuc.y);
        d = Math.max(d,mitem.nuc.d-mitem.nuc.y);
      }
    }
    var TeX = jsMath.TeX; var a = jsMath.Typeset.TeX(style,size).axis_height;
    var delta = Math.max(h-a,d+a);
    var H =  Math.max(Math.floor(TeX.integer*delta/500)*TeX.delimiterfactor,
                      TeX.integer*(2*delta-TeX.delimitershortfall))/TeX.integer;
    var left = this.mlist[0]; var right = this.mlist[this.mlist.length-1];
    left.nuc = jsMath.Box.Delimiter(H,left.left,style);
    right.nuc = jsMath.Box.Delimiter(H,right.right,style);
    left.type = 'open'; left.atom = 1; delete left.left;
    right.type = 'close'; right.atom = 1; delete right.right;
  },

  /*
   *  Typeset a math list to produce final HTML for the list.
   */
  Typeset: function (style,size) {
    var typeset = new jsMath.Typeset(this.mlist);
    return typeset.Typeset(style,size);
  }

});


/*
 *  These routines implement the main rules given in Appendix G of the
 *  TeXbook
 */

jsMath.Add(jsMath.mList.prototype.Atomize,{

  /*
   *  Handle \displaystyle, \textstyle, etc.
   */
  style: function (style,size,mitem,prev,mlist) {
    mlist.style = mitem.style;
  },

  /*
   *  Handle \tiny, \small, etc.
   */
  size: function (style,size,mitem,prev,mlist) {
    mlist.size = mitem.size;
  },

  /*
   *  Create empty boxes of the proper sizes for the various
   *  phantom-type commands
   */
  phantom: function (style,size,mitem) {
    var box = mitem.nuc = jsMath.Box.Set(mitem.phantom,style,size);
    if (mitem.h) {box.Remeasured(); box.html = jsMath.HTML.Spacer(box.w)}
      else {box.html = '', box.w = 0}
    if (!mitem.v) {box.h = box.d = 0}
    box.bd = box.bh = 0;
    delete mitem.phantom;
    mitem.type = 'box';
  },

  /*
   *  Create a box of zero height and depth containing the
   *  contents of the atom
   */
  smash: function (style,size,mitem) {
    var box = mitem.nuc = jsMath.Box.Set(mitem.smash,style,size).Remeasured();
    box.h = box.d = 0;
    delete mitem.smash;
    mitem.type = 'box';
  },

  /*
   *  Move a box up or down vertically
   */
  raise: function (style,size,mitem) {
    mitem.nuc = jsMath.Box.Set(mitem.nuc,style,size);
    var y = mitem.raise;
    mitem.nuc.html = jsMath.HTML.Place(mitem.nuc.html,0,y);
    mitem.nuc.h += y; mitem.nuc.d -= y;
    mitem.type = 'ord'; mitem.atom = 1;
  },

  /*
   *  Hide the size of a box so that it laps to the left or right, or
   *  up or down.
   */
  lap: function (style,size,mitem) {
    var box = jsMath.Box.Set(mitem.nuc,style,size).Remeasured();
    var mlist = [box];
    if (mitem.lap == 'llap') {box.x = -box.w} else
    if (mitem.lap == 'rlap') {mlist[1] = jsMath.mItem.Space(-box.w)} else
    if (mitem.lap == 'ulap') {box.y = box.d; box.h = box.d = 0} else
    if (mitem.lap == 'dlap') {box.y = -box.h; box.h = box.d = 0}
    mitem.nuc = jsMath.Box.SetList(mlist,style,size);
    if (mitem.lap == 'ulap' || mitem.lap == 'dlap') {mitem.nuc.h = mitem.nuc.d = 0}
    mitem.type = 'box'; delete mitem.atom;
  },

  /*
   *  Handle a Bin atom. (Rule 5)
   */
  bin: function (style,size,mitem,prev) {
    if (prev && prev.type) {
      var type  = prev.type;
      if (type == 'bin' || type == 'op' || type == 'rel' ||
          type == 'open' || type == 'punct' || type == '' ||
          (type == 'boundary' && prev.left != '')) {mitem.type = 'ord'}
    } else {mitem.type = 'ord'}
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle a Rel atom.  (Rule 6)
   */
  rel: function (style,size,mitem,prev) {
    if (prev.type && prev.type == 'bin') {prev.type = 'ord'}
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle a Close atom.  (Rule 6)
   */
  close: function (style,size,mitem,prev) {
    if (prev.type && prev.type == 'bin') {prev.type = 'ord'}
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle a Punct atom.  (Rule 6)
   */
  punct: function (style,size,mitem,prev) {
    if (prev.type && prev.type == 'bin') {prev.type = 'ord'}
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle an Open atom.  (Rule 7)
   */
  open: function (style,size,mitem) {
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle an Inner atom.  (Rule 7)
   */
  inner: function (style,size,mitem) {
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle a Vcent atom.  (Rule 8)
   */
  vcenter: function (style,size,mitem) {
    var box = jsMath.Box.Set(mitem.nuc,style,size);
    var TeX = jsMath.Typeset.TeX(style,size);
    box.y = TeX.axis_height - (box.h-box.d)/2;
    mitem.nuc = box; mitem.type = 'ord';
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle an Over atom.  (Rule 9)
   */
  overline: function (style,size,mitem) {
    var TeX = jsMath.Typeset.TeX(style,size);
    var box = jsMath.Box.Set(mitem.nuc,jsMath.Typeset.PrimeStyle(style),size).Remeasured();
    var t = TeX.default_rule_thickness;
    var rule = jsMath.Box.Rule(box.w,t);
    rule.x = -rule.w; rule.y = box.h + 3*t;
    mitem.nuc = jsMath.Box.SetList([box,rule],style,size);
    mitem.nuc.h += t;
    mitem.type = 'ord';
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle an Under atom.  (Rule 10)
   */
  underline: function (style,size,mitem) {
    var TeX = jsMath.Typeset.TeX(style,size);
    var box = jsMath.Box.Set(mitem.nuc,jsMath.Typeset.PrimeStyle(style),size).Remeasured();
    var t = TeX.default_rule_thickness;
    var rule = jsMath.Box.Rule(box.w,t);
    rule.x = -rule.w; rule.y = -box.d - 3*t - t;
    mitem.nuc = jsMath.Box.SetList([box,rule],style,size);
    mitem.nuc.d += t;
    mitem.type = 'ord';
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle a Rad atom.  (Rule 11 plus stuff for \root..\of)
   */
  radical: function (style,size,mitem) {
    var TeX = jsMath.Typeset.TeX(style,size);
    var Cp = jsMath.Typeset.PrimeStyle(style);
    var box = jsMath.Box.Set(mitem.nuc,Cp,size).Remeasured();
    var t = TeX.default_rule_thickness;
    var p = t; if (style == 'D' || style == "D'") {p = TeX.x_height}
    var r = t + p/4;
    var surd = jsMath.Box.Delimiter(box.h+box.d+r+t,[0,2,0x70,3,0x70],style,1);
//    if (surd.h > 0) {t = surd.h} // thickness of rule is height of surd character
    if (surd.d > box.h+box.d+r) {r = (r+surd.d-box.h-box.d)/2}
    surd.y = box.h+r;
    var rule = jsMath.Box.Rule(box.w,t);
    rule.y = surd.y-t/2; rule.h += 3*t/2; box.x = -box.w;
    var Cr = jsMath.Typeset.UpStyle(jsMath.Typeset.UpStyle(style));
    var root = jsMath.Box.Set(mitem.root || null,Cr,size).Remeasured();
    if (mitem.root) {
      root.y = .55*(box.h+box.d+3*t+r)-box.d;
      surd.x = Math.max(root.w-(11/18)*surd.w,0);
      rule.x = (7/18)*surd.w;
      root.x = -(root.w+rule.x);
    }
    mitem.nuc = jsMath.Box.SetList([surd,root,rule,box],style,size);
    mitem.type = 'ord';
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle an Acc atom.  (Rule 12)
   */
  accent: function (style,size,mitem) {
    var TeX = jsMath.Typeset.TeX(style,size);
    var Cp = jsMath.Typeset.PrimeStyle(style);
    var box = jsMath.Box.Set(mitem.nuc,Cp,size);
    var u = box.w; var s; var Font;
    if (mitem.nuc.type == 'TeX') {
      Font = jsMath.TeX[mitem.nuc.font];
      if (Font[mitem.nuc.c].krn && Font.skewchar)
        {s = Font[mitem.nuc.c].krn[Font.skewchar]}
    }
    if (s == null) {s = 0}

    var c = mitem.accent[2];
    var font = jsMath.TeX.fam[mitem.accent[1]]; Font = jsMath.TeX[font];
    while (Font[c].n && Font[Font[c].n].w <= u) {c = Font[c].n}

    var delta = Math.min(box.h,TeX.x_height);
    if (mitem.nuc.type == 'TeX') {
      var nitem = jsMath.mItem.Atom('ord',mitem.nuc);
      nitem.sup = mitem.sup; nitem.sub = mitem.sub; nitem.delta = 0;
      jsMath.mList.prototype.Atomize.SupSub(style,size,nitem);
      delta += (nitem.nuc.h - box.h);
      box = mitem.nuc = nitem.nuc;
      delete mitem.sup; delete mitem.sub;
    }
    var acc = jsMath.Box.TeX(c,font,style,size);
    acc.y = box.h - delta; acc.x = -box.w + s + (u-acc.w)/2;
    if (Font[c].ic) {acc.x -= Font[c].ic * TeX.scale}

    mitem.nuc = jsMath.Box.SetList([box,acc],style,size);
    if (mitem.nuc.w != box.w) {
      var space = jsMath.mItem.Space(box.w-mitem.nuc.w);
      mitem.nuc = jsMath.Box.SetList([mitem.nuc,space],style,size);
    }
    mitem.type = 'ord';
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle an Op atom.  (Rules 13 and 13a)
   */
  op: function (style,size,mitem) {
    var TeX = jsMath.Typeset.TeX(style,size); var box;
    mitem.delta = 0; var isD = (style.charAt(0) == 'D');
    if (mitem.limits == null && isD) {mitem.limits = 1}

    if (mitem.nuc.type == 'TeX') {
      var C = jsMath.TeX[mitem.nuc.font][mitem.nuc.c];
      if (isD && C.n) {mitem.nuc.c = C.n; C = jsMath.TeX[mitem.nuc.font][C.n]}
      box = mitem.nuc = jsMath.Box.Set(mitem.nuc,style,size);
      if (C.ic) {
        mitem.delta = C.ic * TeX.scale;
        if (mitem.limits || !mitem.sub || jsMath.Browser.msieIntegralBug) {
          box = mitem.nuc = jsMath.Box.SetList([box,jsMath.mItem.Space(mitem.delta)],style,size);
        }
      }
      box.y = -((box.h+box.d)/2 - box.d - TeX.axis_height);
      if (Math.abs(box.y) < .0001) {box.y = 0}
    }

    if (!box) {box = mitem.nuc = jsMath.Box.Set(mitem.nuc,style,size).Remeasured()}
    if (mitem.limits) {
      var W = box.w; var x = box.w;
      var mlist = [box]; var dh = 0; var dd = 0;
      if (mitem.sup) {
        var sup = jsMath.Box.Set(mitem.sup,jsMath.Typeset.UpStyle(style),size).Remeasured();
        sup.x = ((box.w-sup.w)/2 + mitem.delta/2) - x; dh = TeX.big_op_spacing5;
        W = Math.max(W,sup.w); x += sup.x + sup.w;
        sup.y = box.h+sup.d + box.y +
                    Math.max(TeX.big_op_spacing1,TeX.big_op_spacing3-sup.d);
        mlist[mlist.length] = sup; delete mitem.sup;
      }
      if (mitem.sub) {
        var sub = jsMath.Box.Set(mitem.sub,jsMath.Typeset.DownStyle(style),size).Remeasured();
        sub.x = ((box.w-sub.w)/2 - mitem.delta/2) - x; dd = TeX.big_op_spacing5;
        W = Math.max(W,sub.w); x += sub.x + sub.w;
        sub.y = -box.d-sub.h + box.y -
                   Math.max(TeX.big_op_spacing2,TeX.big_op_spacing4-sub.h);
        mlist[mlist.length] = sub; delete mitem.sub;
      }
      if (W > box.w) {box.x = (W-box.w)/2; x += box.x}
      if (x < W) {mlist[mlist.length] = jsMath.mItem.Space(W-x)}
      mitem.nuc = jsMath.Box.SetList(mlist,style,size);
      mitem.nuc.h += dh; mitem.nuc.d += dd;
    } else {
      if (jsMath.Browser.msieIntegralBug && mitem.sub && C && C.ic)
        {mitem.nuc = jsMath.Box.SetList([box,jsMath.Box.Space(-C.ic*TeX.scale)],style,size)}
      else if (box.y) {mitem.nuc = jsMath.Box.SetList([box],style,size)}
      jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
    }
  },

  /*
   *  Handle an Ord atom.  (Rule 14)
   */
  ord: function (style,size,mitem,prev,mList,i) {
    if (mitem.nuc.type == 'TeX' && !mitem.sup && !mitem.sub) {
      var nitem = mList.mlist[i+1];
      if (nitem && nitem.atom && nitem.type &&
          (nitem.type == 'ord' || nitem.type == 'op' || nitem.type == 'bin' ||
           nitem.type == 'rel' || nitem.type == 'open' ||
           nitem.type == 'close' || nitem.type == 'punct')) {
        if (nitem.nuc.type == 'TeX' && nitem.nuc.font == mitem.nuc.font) {
          mitem.textsymbol = 1;
          var krn = jsMath.TeX[mitem.nuc.font][mitem.nuc.c].krn;
          krn *= jsMath.Typeset.TeX(style,size).scale;
          if (krn && krn[nitem.nuc.c]) {
            for (var k = mList.mlist.length-1; k > i; k--)
              {mList.mlist[k+1] = mList.mlist[k]}
            mList.mlist[i+1] = jsMath.mItem.Space(krn[nitem.nuc.c]);
          }
        }
      }
    }
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Handle a generalized fraction.  (Rules 15 to 15e)
   */
  fraction: function (style,size,mitem) {
    var TeX = jsMath.Typeset.TeX(style,size); var t = 0;
    if (mitem.thickness != null) {t = mitem.thickness}
    else if (mitem.from.match(/over/)) {t = TeX.default_rule_thickness}
    var isD = (style.charAt(0) == 'D');
    var Cn = (style == 'D')? 'T': (style == "D'")? "T'": jsMath.Typeset.UpStyle(style);
    var Cd = (isD)? "T'": jsMath.Typeset.DownStyle(style);
    var num = jsMath.Box.Set(mitem.num,Cn,size).Remeasured();
    var den = jsMath.Box.Set(mitem.den,Cd,size).Remeasured();

    var u; var v; var w;
    var H = (isD)? TeX.delim1 : TeX.delim2;
    var mlist = [jsMath.Box.Delimiter(H,mitem.left,style)]
    var right = jsMath.Box.Delimiter(H,mitem.right,style);

    if (num.w < den.w) {
      num.x = (den.w-num.w)/2;
      den.x = -(num.w + num.x);
      w = den.w; mlist[1] = num; mlist[2] = den;
    } else {
      den.x = (num.w-den.w)/2;
      num.x = -(den.w + den.x);
      w = num.w; mlist[1] = den; mlist[2] = num;
    }
    if (isD) {u = TeX.num1; v = TeX.denom1} else {
      u = (t != 0)? TeX.num2: TeX.num3;
      v = TeX.denom2;
    }
    if (t == 0) {// atop
      var p = (isD)? 7*TeX.default_rule_thickness: 3*TeX.default_rule_thickness;
      var r = (u - num.d) - (den.h - v);
      if (r < p) {u += (p-r)/2; v += (p-r)/2}
    } else {// over
      var p = (isD)? 3*t: t; var a = TeX.axis_height;
      var r = (u-num.d)-(a+t/2); if (r < p) {u += p-r}
          r = (a-t/2)-(den.h-v); if (r < p) {v += p-r}
      var rule = jsMath.Box.Rule(w,t); rule.x = -w; rule.y = a - t/2;
      mlist[mlist.length] = rule;
    }
    num.y = u; den.y = -v;

    mlist[mlist.length] = right;
    mitem.nuc = jsMath.Box.SetList(mlist,style,size);
    mitem.type = 'ord'; mitem.atom = 1;
    delete mitem.num; delete mitem.den;
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  },

  /*
   *  Add subscripts and superscripts.  (Rules 17-18f)
   */
  SupSub: function (style,size,mitem) {
    var TeX = jsMath.Typeset.TeX(style,size);
    var nuc = mitem.nuc;
    var box = mitem.nuc = jsMath.Box.Set(mitem.nuc,style,size,0);
    if (box.format == 'null')
      {box = mitem.nuc = jsMath.Box.Text('','normal',style,size)}

    if (nuc.type == 'TeX') {
      if (!mitem.textsymbol) {
        var C = jsMath.TeX[nuc.font][nuc.c];
        if (C.ic) {
          mitem.delta = C.ic * TeX.scale;
          if (!mitem.sub) {
            box = mitem.nuc = jsMath.Box.SetList([box,jsMath.Box.Space(mitem.delta)],style,size);
            mitem.delta = 0;
          }
        }
      } else {mitem.delta = 0}
    }

    if (!mitem.sup && !mitem.sub) return;
    mitem.nuc.Styled();

    var Cd = jsMath.Typeset.DownStyle(style);
    var Cu = jsMath.Typeset.UpStyle(style);
    var q = jsMath.Typeset.TeX(Cu,size).sup_drop;
    var r = jsMath.Typeset.TeX(Cd,size).sub_drop;
    var u = 0; var v = 0; var p;
    if (nuc.type && nuc.type != 'text' && nuc.type != 'TeX' && nuc.type != 'null')
      {u = box.h - q; v = box.d + r}

    if (mitem.sub) {
      var sub = jsMath.Box.Set(mitem.sub,Cd,size);
      sub = jsMath.Box.SetList([sub,jsMath.mItem.Space(TeX.scriptspace)],style,size);
    }

    if (!mitem.sup) {
      sub.y = -Math.max(v,TeX.sub1,sub.h-(4/5)*jsMath.Typeset.TeX(Cd,size).x_height);
      mitem.nuc = jsMath.Box.SetList([box,sub],style,size).Styled(); delete mitem.sub;
      return;
    }

    var sup = jsMath.Box.Set(mitem.sup,Cu,size);
    sup = jsMath.Box.SetList([sup,jsMath.mItem.Space(TeX.scriptspace)],style,size);
    if (style == 'D') {p = TeX.sup1}
    else if (style.charAt(style.length-1) == "'") {p = TeX.sup3}
    else {p = TeX.sup2}
    u = Math.max(u,p,sup.d+jsMath.Typeset.TeX(Cu,size).x_height/4);

    if (!mitem.sub) {
      sup.y = u;
      mitem.nuc = jsMath.Box.SetList([box,sup],style,size); delete mitem.sup;
      return;
    }

    v = Math.max(v,jsMath.Typeset.TeX(Cd,size).sub2);
    var t = TeX.default_rule_thickness;
    if ((u-sup.d) - (sub.h -v) < 4*t) {
      v = 4*t + sub.h - (u-sup.d);
      p = (4/5)*TeX.x_height - (u-sup.d);
      if (p > 0) {u += p; v -= p}
    }
    sup.Remeasured(); sub.Remeasured();
    sup.y = u; sub.y = -v; sup.x = mitem.delta;
    if (sup.w+sup.x > sub.w)
      {sup.x -= sub.w; mitem.nuc = jsMath.Box.SetList([box,sub,sup],style,size)} else
      {sub.x -= (sup.w+sup.x); mitem.nuc = jsMath.Box.SetList([box,sup,sub],style,size)}

    delete mitem.sup; delete mitem.sub;
  }

});


/***************************************************************************/

/*
 *  The Typeset object handles most of the TeX-specific processing
 */

jsMath.Typeset = function (mlist) {
  this.type = 'typeset';
  this.mlist = mlist;
}

jsMath.Add(jsMath.Typeset,{

  /*
   *  The "C-uparrow" style table (TeXbook, p. 441)
   */
  upStyle: {
    D: "S", T: "S",  "D'": "S'", "T'": "S'",
    S: "SS",  SS: "SS",  "S'": "SS'", "SS'": "SS'"
  },

  /*
   *  The "C-downarrow" style table (TeXbook, p. 441)
   */
  downStyle: {
    D: "S'", T: "S'",  "D'": "S'", "T'": "S'",
    S: "SS'",  SS: "SS'",  "S'": "SS'", "SS'": "SS'"
  },

  /*
   *  Get the various styles given the current style
   *  (see TeXbook, p. 441)
   */
  UpStyle: function (style) {return this.upStyle[style]},
  DownStyle: function (style) {return this.downStyle[style]},
  PrimeStyle: function (style) {
    if (style.charAt(style.length-1) == "'") {return style}
    return style + "'"
  },

  /*
   *  A value scaled to the appropriate size for scripts
   */
  StyleValue: function (style,v) {
    if (style == "S" || style == "S'")   {return .7*v}
    if (style == "SS" || style == "SS'") {return .5*v}
    return v;
  },

  /*
   *  Return the size associated with a given style and size
   */
  StyleSize: function (style,size) {
    if      (style == "S" || style == "S'")   {size = Math.max(0,size-2)}
    else if (style == "SS" || style == "SS'") {size = Math.max(0,size-4)}
    return size;
  },

  /*
   *  Return the font parameter table for the given style
   */
  TeX: function (style,size) {
    if      (style == "S" || style == "S'")   {size = Math.max(0,size-2)}
    else if (style == "SS" || style == "SS'") {size = Math.max(0,size-4)}
    return jsMath.TeXparams[size];
  },


  /*
   *  Add the CSS class for the given TeX style
   */
  AddStyle: function (style,size,html) {
    if      (style == "S" || style == "S'")   {size = Math.max(0,size-2)}
    else if (style == "SS" || style == "SS'") {size = Math.max(0,size-4)}
    if (size != 4) {html = '<span class="size'+size+'">' + html + '</span>'}
    return html;
  },

  /*
   *  Add the font class, if needed
   */
  AddClass: function (tclass,html) {
    if (tclass != '' && tclass != 'normal') {html = jsMath.HTML.Class(tclass,html)}
    return html;
  }

});


jsMath.Package(jsMath.Typeset,{

  /*
   *  The spacing tables for inter-atom spacing
   *  (See rule 20, and Chapter 18, p 170)
   */
  DTsep: {
    ord: {op: 1, bin: 2, rel: 3, inner: 1},
    op:  {ord: 1, op: 1, rel: 3, inner: 1},
    bin: {ord: 2, op: 2, open: 2, inner: 2},
    rel: {ord: 3, op: 3, open: 3, inner: 3},
    open: {},
    close: {op: 1, bin:2, rel: 3, inner: 1},
    punct: {ord: 1, op: 1, rel: 1, open: 1, close: 1, punct: 1, inner: 1},
    inner: {ord: 1, op: 1, bin: 2, rel: 3, open: 1, punct: 1, inner: 1}
  },

  SSsep: {
    ord: {op: 1},
    op:  {ord: 1, op: 1},
    bin: {},
    rel: {},
    open: {},
    close: {op: 1},
    punct: {},
    inner: {op: 1}
  },

  /*
   *  The sizes used in the tables above
   */
  sepW: ['','thinmuskip','medmuskip','thickmuskip'],


  /*
   *  Find the amount of separation to use between two adjacent
   *  atoms in the given style
   */
  GetSeparation: function (l,r,style) {
    if (l && l.atom && r.atom) {
      var table = this.DTsep; if (style.charAt(0) == "S") {table = this.SSsep}
      var row = table[l.type];
      if (row && row[r.type] != null) {return jsMath.TeX[this.sepW[row[r.type]]]}
    }
    return 0;
  },

  /*
   *  Typeset an mlist (i.e., turn it into HTML).
   *  Here, text items of the same class and style are combined
   *  to reduce the number of <SPAN> tags used (though it is still
   *  huge).  Spaces are combined, when possible.
   *  ###  More needs to be done with that.  ###
   *  The width of the final box is recomputed at the end, since
   *  the final width is not necessarily the sum of the widths of
   *  the individual parts (widths are in pixels, but the browsers
   *  puts pieces together using sub-pixel accuracy).
   */
  Typeset: function (style,size) {
    this.style = style; this.size = size; var unset = -10000
    this.w = 0; this.h = unset; this.d = unset;
    this.bh = this.h; this.bd = this.d;
    this.tbuf = ''; this.tx = 0; this.tclass = '';
    this.cbuf = ''; this.hbuf = ''; this.hx = 0;
    var mitem = null; var prev; this.x = 0; this.dx = 0;

    for (var i = 0; i < this.mlist.length; i++) {
      prev = mitem; mitem = this.mlist[i];
      switch (mitem.type) {

        case 'size':
          this.FlushClassed();
          this.size = mitem.size;
          mitem = prev; // hide this from TeX
          break;

        case 'style':
          this.FlushClassed();
          if (this.style.charAt(this.style.length-1) == "'")
            {this.style = mitem.style + "'"} else {this.style = mitem.style}
          mitem = prev; // hide this from TeX
          break;

        case 'space':
          if (typeof(mitem.w) == 'object') {
            if (this.style.charAt(1) == 'S') {mitem.w = .5*mitem.w[0]/18}
            else if (this.style.charAt(0) == 'S') {mitem.w = .7*mitem.w[0]/18}
            else {mitem.w = mitem.w[0]/18}
          }
          this.dx += mitem.w-0; // mitem.w is sometimes a string?
          mitem = prev; // hide this from TeX
          break;

        case 'html':
          this.FlushClassed();
          if (this.hbuf == '') {this.hx = this.x}
          this.hbuf += mitem.html;
          mitem = prev; // hide this from TeX
          break;

        default:   // atom
          if (!mitem.atom && mitem.type != 'box') break;
          mitem.nuc.x += this.dx + this.GetSeparation(prev,mitem,this.style);
          if (mitem.nuc.y || mitem.nuc.x) mitem.nuc.Styled();
          this.dx = 0; this.x = this.x + this.w;
          this.w += mitem.nuc.w + mitem.nuc.x;
          if (mitem.nuc.format == 'text') {
            if (this.tclass != mitem.nuc.tclass && this.tclass != '') this.FlushText();
            if (this.tbuf == '' && this.cbuf == '') {this.tx = this.x}
            this.tbuf += mitem.nuc.html; this.tclass = mitem.nuc.tclass;
          } else  {
            this.FlushClassed();
            if (mitem.nuc.x || mitem.nuc.y) this.Place(mitem.nuc);
            if (this.hbuf == '') {this.hx = this.x}
            this.hbuf += mitem.nuc.html;
          }
          this.h = Math.max(this.h,mitem.nuc.h+mitem.nuc.y); this.bh = Math.max(this.bh,mitem.nuc.bh);
          this.d = Math.max(this.d,mitem.nuc.d-mitem.nuc.y); this.bd = Math.max(this.bd,mitem.nuc.bd);
          break;
      }
    }

    this.FlushClassed(); // make sure scaling is included
    if (this.dx) {this.hbuf += jsMath.HTML.Spacer(this.dx); this.w += this.dx}
    if (this.hbuf == '') {return jsMath.Box.Null()}
    if (this.h == unset) {this.h = 0}
    if (this.d == unset) {this.d = 0}
    var box = new jsMath.Box('html',this.hbuf,this.w,this.h,this.d);
    box.bh = this.bh; box.bd = this.bd;
    return box;
  },

  /*
   *  Add the font to the buffered text and move it to the
   *  classed-text buffer.
   */
  FlushText: function () {
    if (this.tbuf == '') return;
    this.cbuf += jsMath.Typeset.AddClass(this.tclass,this.tbuf);
    this.tbuf = ''; this.tclass = '';
  },

  /*
   *  Add the script or scriptscript style to the text and
   *  move it to the HTML buffer
   */
  FlushClassed: function () {
    this.FlushText();
    if (this.cbuf == '') return;
    if (this.hbuf == '') {this.hx = this.tx}
    this.hbuf += jsMath.Typeset.AddStyle(this.style,this.size,this.cbuf);
    this.cbuf = '';
  },

  /*
   *  Add a <SPAN> to position an item's HTML, and
   *  adjust the item's height and depth.
   *  (This may be replaced buy one of the following browser-specific
   *   versions by Browser.Init().)
   */
  Place: function (item) {
    var html = '<span style="position: relative;';
    if (item.x) {html += ' margin-left:'+jsMath.HTML.Em(item.x)+';'}
    if (item.y) {html += ' top:'+jsMath.HTML.Em(-item.y)+';'}
    item.html = html + '">' + item.html + '</span>';
    item.h += item.y; item.d -= item.y;
    item.x = 0; item.y = 0;
  },

  /*
   *  For MSIE on Windows, backspacing must be done in a separate
   *  <SPAN>, otherwise the contents will be clipped.  Netscape
   *  also doesn't combine vertical and horizontal spacing well.
   *  Here, the horizontal and vertical spacing are done separately.
   */
  PlaceSeparateSkips: function (item) {
    if (item.y) {
      item.html = '<span style="position: relative; '
                     + 'top:'+jsMath.HTML.Em(-item.y)+';'
                     + '">' + item.html + '</span>'
    }
    if (item.x) {item.html = jsMath.HTML.Spacer(item.x) + item.html}
    item.h += item.y; item.d -= item.y;
    item.x = 0; item.y = 0;
  }

});



/***************************************************************************/

/*
 *  The Parse object handles the parsing of the TeX input string, and creates
 *  the mList to be typeset by the Typeset object above.
 */

jsMath.Parse = function (s,font,size,style) {
  var parse = new jsMath.Parser(s,font,size,style);
  parse.Parse();
  return parse;
}

jsMath.Parser = function (s,font,size,style) {
  this.string = s; this.i = 0;
  this.mlist = new jsMath.mList(null,font,size,style);
}

jsMath.Package(jsMath.Parser,{

  // special characters
  cmd:   '\\',
  open:  '{',
  close: '}',

  // patterns for letters and numbers
  letter:  /[a-z]/i,
  number:  /[0-9]/,
  //  pattern for macros to ^ and _ that should be read with arguments
  scriptargs: /^((math|text)..|mathcal|[hm]box)$/,

  //  the \mathchar definitions (see Appendix B of the TeXbook).
  mathchar: {
    '!': [5,0,0x21],
    '(': [4,0,0x28],
    ')': [5,0,0x29],
    '*': [2,2,0x03], // \ast
    '+': [2,0,0x2B],
    ',': [6,1,0x3B],
    '-': [2,2,0x00],
    '.': [0,1,0x3A],
    '/': [0,1,0x3D],
    ':': [3,0,0x3A],
    ';': [6,0,0x3B],
    '<': [3,1,0x3C],
    '=': [3,0,0x3D],
    '>': [3,1,0x3E],
    '?': [5,0,0x3F],
    '[': [4,0,0x5B],
    ']': [5,0,0x5D],
//  '{': [4,2,0x66],
//  '}': [5,2,0x67],
    '|': [0,2,0x6A]
  },

  //  handle special \catcode characters
  special: {
    '~':   'Tilde',
    '^':   'HandleSuperscript',
    '_':   'HandleSubscript',
    ' ':   'Space',
    '\01': 'Space',
    "\t":  'Space',
    "\r":  'Space',
    "\n":  'Space',
    "'":   'Prime',
    '%':   'HandleComment',
    '&':   'HandleEntry',
    '#':   'Hash'
  },

  // the \mathchardef table (see Appendix B of the TeXbook).
  mathchardef: {
  // brace parts
    braceld:      [0,3,0x7A],
    bracerd:      [0,3,0x7B],
    bracelu:      [0,3,0x7C],
    braceru:      [0,3,0x7D],

  // Greek letters
    alpha:        [0,1,0x0B],
    beta:         [0,1,0x0C],
    gamma:        [0,1,0x0D],
    delta:        [0,1,0x0E],
    epsilon:      [0,1,0x0F],
    zeta:         [0,1,0x10],
    eta:          [0,1,0x11],
    theta:        [0,1,0x12],
    iota:         [0,1,0x13],
    kappa:        [0,1,0x14],
    lambda:       [0,1,0x15],
    mu:           [0,1,0x16],
    nu:           [0,1,0x17],
    xi:           [0,1,0x18],
    pi:           [0,1,0x19],
    rho:          [0,1,0x1A],
    sigma:        [0,1,0x1B],
    tau:          [0,1,0x1C],
    upsilon:      [0,1,0x1D],
    phi:          [0,1,0x1E],
    chi:          [0,1,0x1F],
    psi:          [0,1,0x20],
    omega:        [0,1,0x21],
    varepsilon:   [0,1,0x22],
    vartheta:     [0,1,0x23],
    varpi:        [0,1,0x24],
    varrho:       [0,1,0x25],
    varsigma:     [0,1,0x26],
    varphi:       [0,1,0x27],

    Gamma:        [7,0,0x00],
    Delta:        [7,0,0x01],
    Theta:        [7,0,0x02],
    Lambda:       [7,0,0x03],
    Xi:           [7,0,0x04],
    Pi:           [7,0,0x05],
    Sigma:        [7,0,0x06],
    Upsilon:      [7,0,0x07],
    Phi:          [7,0,0x08],
    Psi:          [7,0,0x09],
    Omega:        [7,0,0x0A],

  // Ord symbols
    aleph:        [0,2,0x40],
    imath:        [0,1,0x7B],
    jmath:        [0,1,0x7C],
    ell:          [0,1,0x60],
    wp:           [0,1,0x7D],
    Re:           [0,2,0x3C],
    Im:           [0,2,0x3D],
    partial:      [0,1,0x40],
    infty:        [0,2,0x31],
    prime:        [0,2,0x30],
    emptyset:     [0,2,0x3B],
    nabla:        [0,2,0x72],
    surd:         [1,2,0x70],
    top:          [0,2,0x3E],
    bot:          [0,2,0x3F],
    triangle:     [0,2,0x34],
    forall:       [0,2,0x38],
    exists:       [0,2,0x39],
    neg:          [0,2,0x3A],
    lnot:         [0,2,0x3A],
    flat:         [0,1,0x5B],
    natural:      [0,1,0x5C],
    sharp:        [0,1,0x5D],
    clubsuit:     [0,2,0x7C],
    diamondsuit:  [0,2,0x7D],
    heartsuit:    [0,2,0x7E],
    spadesuit:    [0,2,0x7F],

  // big ops
    coprod:      [1,3,0x60],
    bigvee:      [1,3,0x57],
    bigwedge:    [1,3,0x56],
    biguplus:    [1,3,0x55],
    bigcap:      [1,3,0x54],
    bigcup:      [1,3,0x53],
    intop:       [1,3,0x52],
    prod:        [1,3,0x51],
    sum:         [1,3,0x50],
    bigotimes:   [1,3,0x4E],
    bigoplus:    [1,3,0x4C],
    bigodot:     [1,3,0x4A],
    ointop:      [1,3,0x48],
    bigsqcup:    [1,3,0x46],
    smallint:    [1,2,0x73],

  // binary operations
    triangleleft:      [2,1,0x2F],
    triangleright:     [2,1,0x2E],
    bigtriangleup:     [2,2,0x34],
    bigtriangledown:   [2,2,0x35],
    wedge:       [2,2,0x5E],
    land:        [2,2,0x5E],
    vee:         [2,2,0x5F],
    lor:         [2,2,0x5F],
    cap:         [2,2,0x5C],
    cup:         [2,2,0x5B],
    ddagger:     [2,2,0x7A],
    dagger:      [2,2,0x79],
    sqcap:       [2,2,0x75],
    sqcup:       [2,2,0x74],
    uplus:       [2,2,0x5D],
    amalg:       [2,2,0x71],
    diamond:     [2,2,0x05],
    bullet:      [2,2,0x0F],
    wr:          [2,2,0x6F],
    div:         [2,2,0x04],
    odot:        [2,2,0x0C],
    oslash:      [2,2,0x0B],
    otimes:      [2,2,0x0A],
    ominus:      [2,2,0x09],
    oplus:       [2,2,0x08],
    mp:          [2,2,0x07],
    pm:          [2,2,0x06],
    circ:        [2,2,0x0E],
    bigcirc:     [2,2,0x0D],
    setminus:    [2,2,0x6E], // for set difference A\setminus B
    cdot:        [2,2,0x01],
    ast:         [2,2,0x03],
    times:       [2,2,0x02],
    star:        [2,1,0x3F],

  // Relations
    propto:      [3,2,0x2F],
    sqsubseteq:  [3,2,0x76],
    sqsupseteq:  [3,2,0x77],
    parallel:    [3,2,0x6B],
    mid:         [3,2,0x6A],
    dashv:       [3,2,0x61],
    vdash:       [3,2,0x60],
    leq:         [3,2,0x14],
    le:          [3,2,0x14],
    geq:         [3,2,0x15],
    ge:          [3,2,0x15],
    lt:          [3,1,0x3C],  // extra since < and > are hard
    gt:          [3,1,0x3E],  //   to get in HTML
    succ:        [3,2,0x1F],
    prec:        [3,2,0x1E],
    approx:      [3,2,0x19],
    succeq:      [3,2,0x17],
    preceq:      [3,2,0x16],
    supset:      [3,2,0x1B],
    subset:      [3,2,0x1A],
    supseteq:    [3,2,0x13],
    subseteq:    [3,2,0x12],
    'in':        [3,2,0x32],
    ni:          [3,2,0x33],
    owns:        [3,2,0x33],
    gg:          [3,2,0x1D],
    ll:          [3,2,0x1C],
    not:         [3,2,0x36],
    sim:         [3,2,0x18],
    simeq:       [3,2,0x27],
    perp:        [3,2,0x3F],
    equiv:       [3,2,0x11],
    asymp:       [3,2,0x10],
    smile:       [3,1,0x5E],
    frown:       [3,1,0x5F],

  // Arrows
    Leftrightarrow:   [3,2,0x2C],
    Leftarrow:        [3,2,0x28],
    Rightarrow:       [3,2,0x29],
    leftrightarrow:   [3,2,0x24],
    leftarrow:        [3,2,0x20],
    gets:             [3,2,0x20],
    rightarrow:       [3,2,0x21],
    to:               [3,2,0x21],
    mapstochar:       [3,2,0x37],
    leftharpoonup:    [3,1,0x28],
    leftharpoondown:  [3,1,0x29],
    rightharpoonup:   [3,1,0x2A],
    rightharpoondown: [3,1,0x2B],
    nearrow:          [3,2,0x25],
    searrow:          [3,2,0x26],
    nwarrow:          [3,2,0x2D],
    swarrow:          [3,2,0x2E],

    minuschar:  [3,2,0x00], // for longmapsto
    hbarchar:   [0,0,0x16], // for \hbar
    lhook:      [3,1,0x2C],
    rhook:      [3,1,0x2D],

    ldotp:      [6,1,0x3A], // ldot as a punctuation mark
    cdotp:      [6,2,0x01], // cdot as a punctuation mark
    colon:      [6,0,0x3A], // colon as a punctuation mark

    '#':        [7,0,0x23],
    '$':        [7,0,0x24],
    '%':        [7,0,0x25],
    '&':        [7,0,0x26]
  },

  // The delimiter table (see Appendix B of the TeXbook)
  delimiter: {
    '(':                [0,0,0x28,3,0x00],
    ')':                [0,0,0x29,3,0x01],
    '[':                [0,0,0x5B,3,0x02],
    ']':                [0,0,0x5D,3,0x03],
    '<':                [0,2,0x68,3,0x0A],
    '>':                [0,2,0x69,3,0x0B],
    '\\lt':             [0,2,0x68,3,0x0A],  // extra since < and > are
    '\\gt':             [0,2,0x69,3,0x0B],  //  hard to get in HTML
    '/':                [0,0,0x2F,3,0x0E],
    '|':                [0,2,0x6A,3,0x0C],
    '.':                [0,0,0x00,0,0x00],
    '\\':               [0,2,0x6E,3,0x0F],
    '\\lmoustache':     [4,3,0x7A,3,0x40],  // top from (, bottom from )
    '\\rmoustache':     [5,3,0x7B,3,0x41],  // top from ), bottom from (
    '\\lgroup':         [4,6,0x28,3,0x3A],  // extensible ( with sharper tips
    '\\rgroup':         [5,6,0x29,3,0x3B],  // extensible ) with sharper tips
    '\\arrowvert':      [0,2,0x6A,3,0x3C],  // arrow without arrowheads
    '\\Arrowvert':      [0,2,0x6B,3,0x3D],  // double arrow without arrowheads
//  '\\bracevert':      [0,7,0x7C,3,0x3E],  // the vertical bar that extends braces
    '\\bracevert':      [0,2,0x6A,3,0x3E],  // we don't load tt, so use | instead
    '\\Vert':           [0,2,0x6B,3,0x0D],
    '\\|':              [0,2,0x6B,3,0x0D],
    '\\vert':           [0,2,0x6A,3,0x0C],
    '\\uparrow':        [3,2,0x22,3,0x78],
    '\\downarrow':      [3,2,0x23,3,0x79],
    '\\updownarrow':    [3,2,0x6C,3,0x3F],
    '\\Uparrow':        [3,2,0x2A,3,0x7E],
    '\\Downarrow':      [3,2,0x2B,3,0x7F],
    '\\Updownarrow':    [3,2,0x6D,3,0x77],
    '\\backslash':      [0,2,0x6E,3,0x0F],  // for double coset G\backslash H
    '\\rangle':         [5,2,0x69,3,0x0B],
    '\\langle':         [4,2,0x68,3,0x0A],
    '\\rbrace':         [5,2,0x67,3,0x09],
    '\\lbrace':         [4,2,0x66,3,0x08],
    '\\}':              [5,2,0x67,3,0x09],
    '\\{':              [4,2,0x66,3,0x08],
    '\\rceil':          [5,2,0x65,3,0x07],
    '\\lceil':          [4,2,0x64,3,0x06],
    '\\rfloor':         [5,2,0x63,3,0x05],
    '\\lfloor':         [4,2,0x62,3,0x04],
    '\\lbrack':         [0,0,0x5B,3,0x02],
    '\\rbrack':         [0,0,0x5D,3,0x03]
  },

  /*
   *  The basic macros for plain TeX.
   *
   *  When the control sequence on the left is called, the JavaScript
   *  funtion on the right is called, with the name of the control sequence
   *  as its first parameter (this way, the same function can be called by
   *  several different control sequences to do similar actions, and the
   *  function can still tell which TeX command was issued).  If the right
   *  is an array, the first entry is the routine to call, and the
   *  remaining entries in the array are parameters to pass to the function
   *  as the second parameter (they are in an array reference).
   *
   *  Note:  TeX macros as defined by the user are discussed below.
   */
  macros: {
    displaystyle:      ['HandleStyle','D'],
    textstyle:         ['HandleStyle','T'],
    scriptstyle:       ['HandleStyle','S'],
    scriptscriptstyle: ['HandleStyle','SS'],

    rm:                ['HandleFont',0],
    mit:               ['HandleFont',1],
    oldstyle:          ['HandleFont',1],
    cal:               ['HandleFont',2],
    it:                ['HandleFont',4],
    bf:                ['HandleFont',6],

    font:              ['Extension','font'],

    left:              'HandleLeft',
    right:             'HandleRight',

    arcsin:       ['NamedOp',0],
    arccos:       ['NamedOp',0],
    arctan:       ['NamedOp',0],
    arg:          ['NamedOp',0],
    cos:          ['NamedOp',0],
    cosh:         ['NamedOp',0],
    cot:          ['NamedOp',0],
    coth:         ['NamedOp',0],
    csc:          ['NamedOp',0],
    deg:          ['NamedOp',0],
    det:           'NamedOp',
    dim:          ['NamedOp',0],
    exp:          ['NamedOp',0],
    gcd:           'NamedOp',
    hom:          ['NamedOp',0],
    inf:           'NamedOp',
    ker:          ['NamedOp',0],
    lg:           ['NamedOp',0],
    lim:           'NamedOp',
    liminf:       ['NamedOp',null,'lim<span style="margin-left: '+1/6+'em"></span>inf'],
    limsup:       ['NamedOp',null,'lim<span style="margin-left: '+1/6+'em"></span>sup'],
    ln:           ['NamedOp',0],
    log:          ['NamedOp',0],
    max:           'NamedOp',
    min:           'NamedOp',
    Pr:            'NamedOp',
    sec:          ['NamedOp',0],
    sin:          ['NamedOp',0],
    sinh:         ['NamedOp',0],
    sup:           'NamedOp',
    tan:          ['NamedOp',0],
    tanh:         ['NamedOp',0],

    vcenter:        ['HandleAtom','vcenter'],
    overline:       ['HandleAtom','overline'],
    underline:      ['HandleAtom','underline'],
    over:            'HandleOver',
    overwithdelims:  'HandleOver',
    atop:            'HandleOver',
    atopwithdelims:  'HandleOver',
    above:           'HandleOver',
    abovewithdelims: 'HandleOver',
    brace:           ['HandleOver','\\{','\\}'],
    brack:           ['HandleOver','[',']'],
    choose:          ['HandleOver','(',')'],

    overbrace:       ['Extension','leaders'],
    underbrace:      ['Extension','leaders'],
    overrightarrow:  ['Extension','leaders'],
    underrightarrow: ['Extension','leaders'],
    overleftarrow:   ['Extension','leaders'],
    underleftarrow:  ['Extension','leaders'],
    overleftrightarrow:  ['Extension','leaders'],
    underleftrightarrow: ['Extension','leaders'],
    overset:         ['Extension','underset-overset'],
    underset:        ['Extension','underset-overset'],

    llap:            'HandleLap',
    rlap:            'HandleLap',
    ulap:            'HandleLap',
    dlap:            'HandleLap',
    raise:           'RaiseLower',
    lower:           'RaiseLower',
    moveleft:        'MoveLeftRight',
    moveright:       'MoveLeftRight',

    frac:            'Frac',
    root:            'Root',
    sqrt:            'Sqrt',

    //  TeX substitution macros
    hbar:               ['Macro','\\hbarchar\\kern-.5em h'],
    ne:                 ['Macro','\\not='],
    neq:                ['Macro','\\not='],
    notin:              ['Macro','\\mathrel{\\rlap{\\kern2mu/}}\\in'],
    cong:               ['Macro','\\mathrel{\\lower2mu{\\mathrel{{\\rlap{=}\\raise6mu\\sim}}}}'],
    bmod:               ['Macro','\\mathbin{\\rm mod}'],
    pmod:               ['Macro','\\kern 18mu ({\\rm mod}\\,\\,#1)',1],
    'int':              ['Macro','\\intop\\nolimits'],
    oint:               ['Macro','\\ointop\\nolimits'],
    doteq:              ['Macro','\\buildrel\\textstyle.\\over='],
    ldots:              ['Macro','\\mathinner{\\ldotp\\ldotp\\ldotp}'],
    cdots:              ['Macro','\\mathinner{\\cdotp\\cdotp\\cdotp}'],
    vdots:              ['Macro','\\mathinner{\\rlap{\\raise8pt{.\\rule 0pt 6pt 0pt}}\\rlap{\\raise4pt{.}}.}'],
    ddots:              ['Macro','\\mathinner{\\kern1mu\\raise7pt{\\rule 0pt 7pt 0pt .}\\kern2mu\\raise4pt{.}\\kern2mu\\raise1pt{.}\\kern1mu}'],
    joinrel:            ['Macro','\\mathrel{\\kern-4mu}'],
    relbar:             ['Macro','\\mathrel{\\smash-}'], // \smash, because - has the same height as +
    Relbar:             ['Macro','\\mathrel='],
    bowtie:             ['Macro','\\mathrel\\triangleright\\joinrel\\mathrel\\triangleleft'],
    models:             ['Macro','\\mathrel|\\joinrel='],
    mapsto:             ['Macro','\\mathrel{\\mapstochar\\rightarrow}'],
    rightleftharpoons:  ['Macro','\\vcenter{\\mathrel{\\rlap{\\raise3mu{\\rightharpoonup}}}\\leftharpoondown}'],
    hookrightarrow:     ['Macro','\\lhook\\joinrel\\rightarrow'],
    hookleftarrow:      ['Macro','\\leftarrow\\joinrel\\rhook'],
    Longrightarrow:     ['Macro','\\Relbar\\joinrel\\Rightarrow'],
    longrightarrow:     ['Macro','\\relbar\\joinrel\\rightarrow'],
    longleftarrow:      ['Macro','\\leftarrow\\joinrel\\relbar'],
    Longleftarrow:      ['Macro','\\Leftarrow\\joinrel\\Relbar'],
    longmapsto:         ['Macro','\\mathrel{\\mapstochar\\minuschar\\joinrel\\rightarrow}'],
    longleftrightarrow: ['Macro','\\leftarrow\\joinrel\\rightarrow'],
    Longleftrightarrow: ['Macro','\\Leftarrow\\joinrel\\Rightarrow'],
    iff:                ['Macro','\\;\\Longleftrightarrow\\;'],
    mathcal:            ['Macro','{\\cal #1}',1],
    mathrm:             ['Macro','{\\rm #1}',1],
    mathbf:             ['Macro','{\\bf #1}',1],
    mathbb:             ['Macro','{\\bf #1}',1],
    mathit:             ['Macro','{\\it #1}',1],
    textrm:             ['Macro','\\mathord{\\hbox{#1}}',1],
    textit:             ['Macro','\\mathord{\\class{textit}{\\hbox{#1}}}',1],
    textbf:             ['Macro','\\mathord{\\class{textbf}{\\hbox{#1}}}',1],
    pmb:                ['Macro','\\rlap{#1}\\kern1px{#1}',1],

    TeX:                ['Macro','T\\kern-.1667em\\lower.5ex{E}\\kern-.125em X'],

    limits:       ['Limits',1],
    nolimits:     ['Limits',0],

    ',':          ['Spacer',1/6],
    ':':          ['Spacer',1/6],  // for LaTeX
    '>':          ['Spacer',2/9],
    ';':          ['Spacer',5/18],
    '!':          ['Spacer',-1/6],
    enspace:      ['Spacer',1/2],
    quad:         ['Spacer',1],
    qquad:        ['Spacer',2],
    thinspace:    ['Spacer',1/6],
    negthinspace: ['Spacer',-1/6],

    hskip:         'Hskip',
    kern:          'Hskip',
    rule:          ['Rule','colored'],
    space:         ['Rule','blank'],

    big:        ['MakeBig','ord',0.85],
    Big:        ['MakeBig','ord',1.15],
    bigg:       ['MakeBig','ord',1.45],
    Bigg:       ['MakeBig','ord',1.75],
    bigl:       ['MakeBig','open',0.85],
    Bigl:       ['MakeBig','open',1.15],
    biggl:      ['MakeBig','open',1.45],
    Biggl:      ['MakeBig','open',1.75],
    bigr:       ['MakeBig','close',0.85],
    Bigr:       ['MakeBig','close',1.15],
    biggr:      ['MakeBig','close',1.45],
    Biggr:      ['MakeBig','close',1.75],
    bigm:       ['MakeBig','rel',0.85],
    Bigm:       ['MakeBig','rel',1.15],
    biggm:      ['MakeBig','rel',1.45],
    Biggm:      ['MakeBig','rel',1.75],

    mathord:    ['HandleAtom','ord'],
    mathop:     ['HandleAtom','op'],
    mathopen:   ['HandleAtom','open'],
    mathclose:  ['HandleAtom','close'],
    mathbin:    ['HandleAtom','bin'],
    mathrel:    ['HandleAtom','rel'],
    mathpunct:  ['HandleAtom','punct'],
    mathinner:  ['HandleAtom','inner'],

    mathchoice: ['Extension','mathchoice'],
    buildrel:   'BuildRel',

    hbox:       'HBox',
    text:       'HBox',
    mbox:       'HBox',
    fbox:       ['Extension','fbox'],

    strut:      'Strut',
    mathstrut:  ['Macro','\\vphantom{(}'],
    phantom:    ['Phantom',1,1],
    vphantom:   ['Phantom',1,0],
    hphantom:   ['Phantom',0,1],
    smash:      'Smash',

    acute:      ['MathAccent', [7,0,0x13]],
    grave:      ['MathAccent', [7,0,0x12]],
    ddot:       ['MathAccent', [7,0,0x7F]],
    tilde:      ['MathAccent', [7,0,0x7E]],
    bar:        ['MathAccent', [7,0,0x16]],
    breve:      ['MathAccent', [7,0,0x15]],
    check:      ['MathAccent', [7,0,0x14]],
    hat:        ['MathAccent', [7,0,0x5E]],
    vec:        ['MathAccent', [0,1,0x7E]],
    dot:        ['MathAccent', [7,0,0x5F]],
    widetilde:  ['MathAccent', [0,3,0x65]],
    widehat:    ['MathAccent', [0,3,0x62]],

    '_':        ['Replace','ord','_','normal',-.4,.1],
    ' ':        ['Replace','ord','&nbsp;','normal'],
    angle:      ['Replace','ord','&#x2220;','normal'],

    matrix:     'Matrix',
    array:      'Matrix',  // ### still need to do alignment options ###
    pmatrix:    ['Matrix','(',')','c'],
    cases:      ['Matrix','\\{','.',['l','l'],null,2],
    eqalign:    ['Matrix',null,null,['r','l'],[5/18],3,'D'],
    displaylines: ['Matrix',null,null,['c'],null,3,'D'],
    cr:         'HandleRow',
    '\\':       'HandleRow',
    newline:    'HandleRow',
    noalign:    'HandleNoAlign',
    eqalignno:  ['Matrix',null,null,['r','l','r'],[5/8,3],3,'D'],
    leqalignno: ['Matrix',null,null,['r','l','r'],[5/8,3],3,'D'],

    //  LaTeX
    begin:      'Begin',
    end:        'End',
    tiny:       ['HandleSize',0],
    Tiny:       ['HandleSize',1],  // non-standard
    scriptsize: ['HandleSize',2],
    small:      ['HandleSize',3],
    normalsize: ['HandleSize',4],
    large:      ['HandleSize',5],
    Large:      ['HandleSize',6],
    LARGE:      ['HandleSize',7],
    huge:       ['HandleSize',8],
    Huge:       ['HandleSize',9],
    dots:       ['Macro','\\ldots'],

    newcommand:     ['Extension','newcommand'],
    newenvironment: ['Extension','newcommand'],
    def:            ['Extension','newcommand'],

    //  Extensions to TeX
    color:      ['Extension','HTML'],
    href:       ['Extension','HTML'],
    'class':    ['Extension','HTML'],
    style:      ['Extension','HTML'],
    cssId:      ['Extension','HTML'],
    unicode:    ['Extension','HTML'],
    bbox:       ['Extension','bbox'],

    require:    'Require',

    //  debugging and test routines
    'char':     'Char'
  },

  /*
   *  LaTeX environments
   */
  environments: {
    array:        'Array',
    matrix:       ['Array',null,null,'c'],
    pmatrix:      ['Array','(',')','c'],
    bmatrix:      ['Array','[',']','c'],
    Bmatrix:      ['Array','\\{','\\}','c'],
    vmatrix:      ['Array','\\vert','\\vert','c'],
    Vmatrix:      ['Array','\\Vert','\\Vert','c'],
    cases:        ['Array','\\{','.','ll',null,2],
    eqnarray:     ['Array',null,null,'rcl',[5/18,5/18],3,'D'],

    align:        ['Extension','AMSmath'],
    'align*':     ['Extension','AMSmath'],
    multline:     ['Extension','AMSmath'],
    'multline*':  ['Extension','AMSmath'],
    split:        ['Extension','AMSmath'],
    gather:       ['Extension','AMSmath'],
    'gather*':    ['Extension','AMSmath']
  },


  /***************************************************************************/

  /*
   *  Add special characters to list above.  (This makes it possible
   *  to define them in a variable that the user can change.)
   */
  AddSpecial: function (obj) {
    for (var id in obj) {
      jsMath.Parser.prototype.special[jsMath.Parser.prototype[id]] = obj[id];
    }
  },

  /*
   *  Throw an error
   */
  Error: function (s) {
    this.i = this.string.length;
    if (s.error) {this.error = s.error} else {
      if (!this.error) {this.error = s}
    }
  },

  /***************************************************************************/

  /*
   *  Check if the next character is a space
   */
  nextIsSpace: function () {
    return this.string.charAt(this.i) == ' ';
  },

  /*
   *  Trim spaces from a string
   */
  trimSpaces: function (text) {
    if (typeof(text) != 'string') {return text}
    return text.replace(/^\s+|\s+/g,'');
  },

  /*
   *  Parse a substring to get its mList, and return it.
   *  Check that no errors occured
   */
  Process: function (arg) {
    var data = this.mlist.data;
    arg = jsMath.Parse(arg,data.font,data.size,data.style);
      if (arg.error) {this.Error(arg); return null}
    if (arg.mlist.Length() == 0) {return null}
    if (arg.mlist.Length() == 1) {
      var atom = arg.mlist.Last();
      if (atom.atom && atom.type == 'ord' && atom.nuc &&
         !atom.sub && !atom.sup && (atom.nuc.type == 'text' || atom.nuc.type == 'TeX'))
             {return atom.nuc}
    }
    return {type: 'mlist', mlist: arg.mlist};
  },

  /*
   *  Get and return a control-sequence name from the TeX string
   */
  GetCommand: function () {
    var letter = /^([a-z]+|.) ?/i;
    var cmd = letter.exec(this.string.slice(this.i));
    if (cmd) {this.i += cmd[1].length; return cmd[1]}
    this.i++; return " ";
  },

  /*
   *  Get and return a TeX argument (either a single character or control sequence,
   *  or the contents of the next set of braces).
   */
  GetArgument: function (name,noneOK) {
    while (this.nextIsSpace()) {this.i++}
    if (this.i >= this.string.length) {if (!noneOK) this.Error("Missing argument for "+name); return null}
    if (this.string.charAt(this.i) == this.close) {if (!noneOK) this.Error("Extra close brace"); return null}
    if (this.string.charAt(this.i) == this.cmd) {this.i++; return this.cmd+this.GetCommand()}
    if (this.string.charAt(this.i) != this.open) {return this.string.charAt(this.i++)}
    var j = ++this.i; var pcount = 1; var c = '';
    while (this.i < this.string.length) {
      c = this.string.charAt(this.i++);
      if (c == this.cmd) {this.i++}
      else if (c == this.open) {pcount++}
      else if (c == this.close) {
        if (pcount == 0) {this.Error("Extra close brace"); return null}
        if (--pcount == 0) {return this.string.slice(j,this.i-1)}
      }
    }
    this.Error("Missing close brace");
    return null;
  },

  /*
   *  Get an argument and process it into an mList
   */
  ProcessArg: function (name) {
    var arg = this.GetArgument(name); if (this.error) {return null}
    return this.Process(arg);
  },

  /*
   *  Get and process an argument for a super- or subscript.
   *  (read extra args for \frac, \sqrt, \mathrm, etc.)
   *  This handles these macros as special cases, so is really
   *  rather a hack.  A more general method for indicating
   *  how to handle macros in scripts needs to be developed.
   */
  ProcessScriptArg: function (name) {
    var arg = this.GetArgument(name); if (this.error) {return null}
    if (arg.charAt(0) == this.cmd) {
      var csname = arg.substr(1);
      if (csname == "frac") {
        arg += '{'+this.GetArgument(csname)+'}'; if (this.error) {return null}
        arg += '{'+this.GetArgument(csname)+'}'; if (this.error) {return null}
      } else if (csname == "sqrt") {
        arg += '['+this.GetBrackets(csname)+']'; if (this.error) {return null}
        arg += '{'+this.GetArgument(csname)+'}'; if (this.error) {return null}
      } else if (csname.match(this.scriptargs)) {
        arg += '{'+this.GetArgument(csname)+'}'; if (this.error) {return null}
      }
    }
    return this.Process(arg);
  },

  /*
   *  Get the name of a delimiter (check it in the delimiter list).
   */
  GetDelimiter: function (name) {
    while (this.nextIsSpace()) {this.i++}
    var c = this.string.charAt(this.i);
    if (this.i < this.string.length) {
      this.i++;
      if (c == this.cmd) {c = '\\'+this.GetCommand(name); if (this.error) return null}
      if (this.delimiter[c] != null) {return this.delimiter[c]}
    }
    this.Error("Missing or unrecognized delimiter for "+name);
    return null;
  },

  /*
   *  Get a dimension (including its units).
   *  Convert the dimen to em's, except for mu's, which must be
   *  converted when typeset.
   */
  GetDimen: function (name,nomu) {
    var rest; var advance = 0;
    if (this.nextIsSpace()) {this.i++}
    if (this.string.charAt(this.i) == '{') {
      rest = this.GetArgument(name);
    } else {
      rest = this.string.slice(this.i);
      advance = 1;
    }
    return this.ParseDimen(rest,name,advance,nomu);
  },

  ParseDimen: function (dimen,name,advance,nomu) {
    var match = dimen.match(/^\s*([-+]?(\.\d+|\d+(\.\d*)?))(pt|em|ex|mu|px)/);
    if (!match) {this.Error("Missing dimension or its units for "+name); return null}
    if (advance) {
      this.i += match[0].length;
      if (this.nextIsSpace()) {this.i++}
    }
    var d = match[1]-0;
    if (match[4] == 'px') {d /= jsMath.em}
    else if (match[4] == 'pt') {d /= 10}
    else if (match[4] == 'ex') {d *= jsMath.TeX.x_height}
    else if (match[4] == 'mu') {if (nomu) {d = d/18} else {d = [d,'mu']}}
    return d;
  },

  /*
   *  Get the next non-space character
   */
  GetNext: function () {
    while (this.nextIsSpace()) {this.i++}
    return this.string.charAt(this.i);
  },

  /*
   *  Get an optional LaTeX argument in brackets
   */
  GetBrackets: function (name) {
    var c = this.GetNext(); if (c != '[') return '';
    var start = ++this.i; var pcount = 0;
    while (this.i < this.string.length) {
      var c = this.string.charAt(this.i++);
      if (c == '{') {pcount++}
      else if (c == '}') {
        if (pcount == 0)
          {this.Error("Extra close brace while looking for ']'"); return null}
        pcount --;
      } else if (c == this.cmd) {
        this.i++;
      } else if (c == ']') {
        if (pcount == 0) {return this.string.slice(start,this.i-1)}
      }
    }
    this.Error("Couldn't find closing ']' for argument to "+this.cmd+name);
    return null;
  },

  /*
   *  Get everything up to the given control sequence name (token)
   */
  GetUpto: function (name,token) {
    while (this.nextIsSpace()) {this.i++}
    var start = this.i; var pcount = 0;
    while (this.i < this.string.length) {
      var c = this.string.charAt(this.i++);
      if (c == '{') {pcount++}
      else if (c == '}') {
        if (pcount == 0)
          {this.Error("Extra close brace while looking for "+this.cmd+token); return null}
        pcount --;
      } else if (c == this.cmd) {
        // really need separate counter for begin/end
        // and it should really be a stack (new pcount for each begin)
        if (this.string.slice(this.i,this.i+5) == "begin") {pcount++; this.i+=4}
        else if (this.string.slice(this.i,this.i+3) == "end") {
          if (pcount > 0) {pcount--; this.i += 2}
        }
        if (pcount == 0)  {
          if (this.string.slice(this.i,this.i+token.length) == token) {
            c = this.string.charAt(this.i+token.length);
            if (c.match(/[^a-z]/i) || !token.match(/[a-z]/i)) {
              var arg = this.string.slice(start,this.i-1);
              this.i += token.length;
              return arg;
            }
          }
        }
        this.i++;
      }
    }
    this.Error("Couldn't find "+this.cmd+token+" for "+name);
    return null;
  },

  /*
   *  Get a parameter delimited by a control sequence, and
   *  process it to get its mlist
   */
  ProcessUpto: function (name,token) {
    var arg = this.GetUpto(name,token); if (this.error) return null;
    return this.Process(arg);
  },

  /*
   *  Get everything up to \end{env}
   */
  GetEnd: function (env) {
    var body = ''; var name = '';
    while (name != env) {
      body += this.GetUpto('begin{'+env+'}','end'); if (this.error) return null;
      name = this.GetArgument(this.cmd+'end'); if (this.error) return null;
    }
    return body;
  },


  /***************************************************************************/


  /*
   *  Ignore spaces
   */
  Space: function () {},

  /*
   *  Collect together any primes and convert them to a superscript
   */
  Prime: function (c) {
    var base = this.mlist.Last();
    if (base == null || (!base.atom && base.type != 'box' && base.type != 'frac'))
       {base = this.mlist.Add(jsMath.mItem.Atom('ord',{type:null}))}
    if (base.sup) {this.Error("Prime causes double exponent: use braces to clarify"); return}
    var sup = '';
    while (c == "'") {sup += this.cmd+'prime'; c = this.GetNext(); if (c == "'") {this.i++}}
    base.sup = this.Process(sup);
    base.sup.isPrime = 1;
  },

  /*
   *  Raise or lower its parameter by a given amount
   *  @@@ Note that this is different from TeX, which requires an \hbox @@@
   *  ### make this work with mu's ###
   */
  RaiseLower: function (name) {
    var h = this.GetDimen(this.cmd+name,1); if (this.error) return;
    var box = this.ProcessScriptArg(this.cmd+name); if (this.error) return;
    if (name == 'lower') {h = -h}
    this.mlist.Add(new jsMath.mItem('raise',{nuc: box, raise: h}));
  },

  /*
   *  Shift an expression to the right or left
   *  @@@ Note that this is different from TeX, which requires a \vbox @@@
   *  ### make this work with mu's ###
   */
  MoveLeftRight: function (name) {
    var x = this.GetDimen(this.cmd+name,1); if (this.error) return;
    var box = this.ProcessArg(this.cmd+name); if (this.error) return;
    if (name == 'moveleft') {x = -x}
    this.mlist.Add(jsMath.mItem.Space(x));
    this.mlist.Add(jsMath.mItem.Atom('ord',box));
    this.mlist.Add(jsMath.mItem.Space(-x));
  },

  /*
   *  Load an extension if it has not already been loaded
   */
  Require: function (name) {
    var file = this.GetArgument(this.cmd+name); if (this.error) return;
    file = jsMath.Extension.URL(file);
    if (jsMath.Setup.loaded[file]) return;
    this.Extension(null,[file]);
  },

  /*
   *  Load an extension file and restart processing the math
   */
  Extension: function (name,data) {
    jsMath.Translate.restart = 1;
    if (name != null) {delete jsMath.Parser.prototype[data[1]||'macros'][name]}
    jsMath.Extension.Require(data[0],jsMath.Translate.asynchronous);
    throw "restart";
  },

  /*
   *  Implements \frac{num}{den}
   */
  Frac: function (name) {
    var num = this.ProcessArg(this.cmd+name); if (this.error) return;
    var den = this.ProcessArg(this.cmd+name); if (this.error) return;
    this.mlist.Add(jsMath.mItem.Fraction('over',num,den));
  },

  /*
   *  Implements \sqrt[n]{...}
   */
  Sqrt: function (name) {
    var n = this.GetBrackets(this.cmd+name); if (this.error) return;
    var arg = this.ProcessArg(this.cmd+name); if (this.error) return;
    var box = jsMath.mItem.Atom('radical',arg);
    if (n != '') {box.root = this.Process(n); if (this.error) return}
    this.mlist.Add(box);
  },

  /*
   *  Implements \root...\of{...}
   */
  Root: function (name) {
    var n = this.ProcessUpto(this.cmd+name,'of'); if (this.error) return;
    var arg = this.ProcessArg(this.cmd+name); if (this.error) return;
    var box = jsMath.mItem.Atom('radical',arg);
    box.root = n; this.mlist.Add(box);
  },


  /*
   *  Implements \buildrel...\over{...}
   */
  BuildRel: function (name) {
    var top = this.ProcessUpto(this.cmd+name,'over'); if (this.error) return;
    var bot = this.ProcessArg(this.cmd+name); if (this.error) return;
    var op = jsMath.mItem.Atom('op',bot);
    op.limits = 1; op.sup = top;
    this.mlist.Add(op);
  },

  /*
   *  Create a delimiter of the type and size specified in the parameters
   */
  MakeBig: function (name,data) {
    var type = data[0]; var h = data[1] * jsMath.p_height;
    var delim = this.GetDelimiter(this.cmd+name); if (this.error) return;
    this.mlist.Add(jsMath.mItem.Atom(type,jsMath.Box.Delimiter(h,delim,'T')));
  },

  /*
   *  Insert the specified character in the given font.
   *  (Try to load the font if it is not already available.)
   */
  Char: function (name) {
    var font = this.GetArgument(this.cmd+name); if (this.error) return;
    var n = this.GetArgument(this.cmd+name); if (this.error) return;
    if (!jsMath.TeX[font]) {
      jsMath.TeX[font] = [];
      this.Extension(null,[jsMath.Font.URL(font)]);
    } else {
      this.mlist.Add(jsMath.mItem.Typeset(jsMath.Box.TeX(n-0,font,'T',this.mlist.data.size)));
    }
  },

  /*
   *  Create an array or matrix.
   */
  Matrix: function (name,delim) {
    var data = this.mlist.data;
    var arg = this.GetArgument(this.cmd+name); if (this.error) return;
    var parse = new jsMath.Parser(arg+this.cmd+'\\',null,data.size,delim[5] || 'T');
    parse.matrix = name; parse.row = []; parse.table = []; parse.rspacing = [];
    parse.Parse(); if (parse.error) {this.Error(parse); return}
    parse.HandleRow(name,1);  // be sure the last row is recorded
    var box = jsMath.Box.Layout(data.size,parse.table,delim[2]||null,delim[3]||null,parse.rspacing,delim[4]||null);
    // Add parentheses, if needed
    if (delim[0] && delim[1]) {
      var left  = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[0]],'T');
      var right = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[1]],'T');
      box = jsMath.Box.SetList([left,box,right],data.style,data.size);
    }
    this.mlist.Add(jsMath.mItem.Atom((delim[0]? 'inner': 'ord'),box));
  },

  /*
   *  When we see an '&', try to add a matrix entry to the row data.
   *  (Use all the data in the current mList, and then clear it)
   */
  HandleEntry: function (name) {
    if (!this.matrix)
      {this.Error(name+" can only appear in a matrix or array"); return}
    if (this.mlist.data.openI != null) {
      var open = this.mlist.Get(this.mlist.data.openI);
      if (open.left) {this.Error("Missing "+this.cmd+"right")}
        else {this.Error("Missing close brace")}
    }
    if (this.mlist.data.overI != null) {this.mlist.Over()}
    var data = this.mlist.data;
    this.mlist.Atomize(data.style,data.size);
    var box = this.mlist.Typeset(data.style,data.size);
    box.entry = data.entry; delete data.entry; if (!box.entry) {box.entry = {}};
    this.row[this.row.length] = box;
    this.mlist = new jsMath.mList(null,null,data.size,data.style);
  },

  /*
   *  When we see a \cr or \\, try to add a row to the table
   */
  HandleRow: function (name,last) {
    var dimen;
    if (!this.matrix) {this.Error(this.cmd+name+" can only appear in a matrix or array"); return}
    if (name == "\\") {
      dimen = this.GetBrackets(this.cmd+name); if (this.error) return;
      if (dimen) {dimen = this.ParseDimen(dimen,this.cmd+name,0,1)}
    }
    this.HandleEntry(name);
    if (!last || this.row.length > 1 || this.row[0].format != 'null')
      {this.table[this.table.length] = this.row}
    if (dimen) {this.rspacing[this.table.length] = dimen}
    this.row = [];
  },

  /*
   *  Look for \vskip or \vspace in \noalign parameters
   */
  HandleNoAlign: function (name) {
    var arg = this.GetArgument(this.cmd+name); if (this.error) return;
    var skip = arg.replace(/^.*(vskip|vspace)([^a-z])/i,'$2');
    if (skip.length == arg.length) return;
    var d = this.ParseDimen(skip,this.cmd+RegExp.$1,0,1); if (this.error) return;
    this.rspacing[this.table.length] = (this.rspacing[this.table.length] || 0) + d;
  },

  /*
   *  LaTeX array environment
   */
  Array: function (name,delim) {
    var columns = delim[2]; var cspacing = delim[3];
    if (!columns) {
      columns = this.GetArgument(this.cmd+'begin{'+name+'}');
      if (this.error) return;
    }
    columns = columns.replace(/[^clr]/g,'');
    columns = columns.split('');
    var data = this.mlist.data; var style = delim[5] || 'T';
    var arg = this.GetEnd(name); if (this.error) return;
    var parse = new jsMath.Parser(arg+this.cmd+'\\',null,data.size,style);
    parse.matrix = name; parse.row = []; parse.table = []; parse.rspacing = [];
    parse.Parse(); if (parse.error) {this.Error(parse); return}
    parse.HandleRow(name,1);  // be sure the last row is recorded
    var box = jsMath.Box.Layout(data.size,parse.table,columns,cspacing,parse.rspacing,delim[4],delim[6],delim[7]);
    // Add parentheses, if needed
    if (delim[0] && delim[1]) {
      var left  = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[0]],'T');
      var right = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[1]],'T');
      box = jsMath.Box.SetList([left,box,right],data.style,data.size);
    }
    this.mlist.Add(jsMath.mItem.Atom((delim[0]? 'inner': 'ord'),box));
  },

  /*
   *  LaTeX \begin{env}
   */
  Begin: function (name) {
    var env = this.GetArgument(this.cmd+name); if (this.error) return;
    if (env.match(/[^a-z*]/i)) {this.Error('Invalid environment name "'+env+'"'); return}
    if (!this.environments[env]) {this.Error('Unknown environment "'+env+'"'); return}
    var cmd = this.environments[env];
    if (typeof(cmd) == "string") {cmd = [cmd]}
    this[cmd[0]](env,cmd.slice(1));
  },

  /*
   *  LaTeX \end{env}
   */
  End: function (name) {
    var env = this.GetArgument(this.cmd+name); if (this.error) return;
    this.Error(this.cmd+name+'{'+env+'} without matching '+this.cmd+'begin');
  },

  /*
   *  Add a fixed amount of horizontal space
   */
  Spacer: function (name,w) {
    this.mlist.Add(jsMath.mItem.Space(w-0));
  },

  /*
   *  Add horizontal space given by the argument
   */
  Hskip: function (name) {
    var w = this.GetDimen(this.cmd+name); if (this.error) return;
    this.mlist.Add(jsMath.mItem.Space(w));
  },

  /*
   *  Typeset the argument as plain text rather than math.
   */
  HBox: function (name) {
    var text = this.GetArgument(this.cmd+name); if (this.error) return;
    var box = jsMath.Box.InternalMath(text,this.mlist.data.size);
    this.mlist.Add(jsMath.mItem.Typeset(box));
  },

  /*
   *  Insert a rule of a particular width, height and depth
   *  This replaces \hrule and \vrule
   *  @@@ not a standard TeX command, and all three parameters must be given @@@
   */
  Rule: function (name,style) {
    var w = this.GetDimen(this.cmd+name,1); if (this.error) return;
    var h = this.GetDimen(this.cmd+name,1); if (this.error) return;
    var d = this.GetDimen(this.cmd+name,1); if (this.error) return;
    h += d; var html;
    if (h != 0) {h = Math.max(1.05/jsMath.em,h)}
    if (h == 0 || w == 0 || style == "blank")
      {html = jsMath.HTML.Blank(w,h)} else {html = jsMath.HTML.Rule(w,h)}
    if (d) {
      html = '<span style="vertical-align:'+jsMath.HTML.Em(-d)+'">'
           +  html + '</span>';
    }
    this.mlist.Add(jsMath.mItem.Typeset(new jsMath.Box('html',html,w,h-d,d)));
  },

  /*
   *  Inserts an empty box of a specific height and depth
   */
  Strut: function () {
    var size = this.mlist.data.size;
    var box = jsMath.Box.Text('','normal','T',size).Styled();
    box.bh = box.bd = 0; box.h = .8; box.d = .3; box.w = 0;
    this.mlist.Add(jsMath.mItem.Typeset(box));
  },

  /*
   *  Handles \phantom, \vphantom and \hphantom
   */
  Phantom: function (name,data) {
    var arg = this.ProcessArg(this.cmd+name); if (this.error) return;
    this.mlist.Add(new jsMath.mItem('phantom',{phantom: arg, v: data[0], h: data[1]}));
  },

  /*
   *  Implements \smash
   */
  Smash: function (name,data) {
    var arg = this.ProcessArg(this.cmd+name); if (this.error) return;
    this.mlist.Add(new jsMath.mItem('smash',{smash: arg}));
  },

  /*
   *  Puts an accent on the following argument
   */
  MathAccent: function (name,accent) {
    var c = this.ProcessArg(this.cmd+name); if (this.error) return;
    var atom = jsMath.mItem.Atom('accent',c); atom.accent = accent[0];
    this.mlist.Add(atom);
  },

  /*
   *  Handles functions and operators like sin, cos, sum, etc.
   */
  NamedOp: function (name,data) {
    var a = (name.match(/[^acegm-su-z]/)) ? 1: 0;
    var d = (name.match(/[gjpqy]/)) ? .2: 0;
    if (data[1]) {name = data[1]}
    var box = jsMath.mItem.TextAtom('op',name,'cmr10',a,d);
    if (data[0] != null) {box.limits = data[0]}
    this.mlist.Add(box);
  },

  /*
   *  Implements \limits
   */
  Limits: function (name,data) {
    var atom = this.mlist.Last();
    if (!atom || atom.type != 'op')
      {this.Error(this.cmd+name+" is allowed only on operators"); return}
    atom.limits = data[0];
  },

  /*
   *  Implements macros like those created by \def.  The named control
   *  sequence is replaced by the string given as the first data value.
   *  If there is a second data value, this specifies how many arguments
   *  the macro uses, and in this case, those arguments are substituted
   *  for #1, #2, etc. within the replacement string.
   *
   *  See the jsMath.Macro() command below for more details.
   *  The "newcommand" extension implements \newcommand and \def
   *  and are loaded automatically if needed.
   */
  Macro: function (name,data) {
    var text = data[0];
    if (data[1]) {
      var args = [];
      for (var i = 0; i < data[1]; i++)
        {args[args.length] = this.GetArgument(this.cmd+name); if (this.error) return}
      text = this.SubstituteArgs(args,text);
    }
    this.string = this.AddArgs(text,this.string.slice(this.i));
    this.i = 0;
  },

  /*
   *  Replace macro paramters with their values
   */
  SubstituteArgs: function (args,string) {
    var text = ''; var newstring = ''; var c; var i = 0;
    while (i < string.length) {
      c = string.charAt(i++);
      if (c == this.cmd) {text += c + string.charAt(i++)}
      else if (c == '#') {
        c = string.charAt(i++);
        if (c == "#") {text += c} else {
          if (!c.match(/[1-9]/) || c > args.length)
            {this.Error("Illegal macro parameter reference"); return null}
          newstring = this.AddArgs(this.AddArgs(newstring,text),args[c-1]);
          text = '';
        }
      } else {text += c}
    }
    return this.AddArgs(newstring,text);
  },

  /*
   *  Make sure that macros are followed by a space if their names
   *  could accidentally be continued into the following text.
   */
  AddArgs: function (s1,s2) {
    if (s2.match(/^[a-z]/i) && s1.match(/(^|[^\\])(\\\\)*\\[a-z]+$/i)) {s1 += ' '}
    return s1+s2;
  },

  /*
   *  Replace the control sequence with the given text
   */
  Replace: function (name,data) {
    this.mlist.Add(jsMath.mItem.TextAtom(data[0],data[1],data[2],data[3]));
  },

  /*
   *  Error for # (must use \#)
   */
  Hash: function (name) {
    this.Error("You can't use 'macro parameter character #' in math mode");
  },

  /*
   *  Insert space for ~
   */
  Tilde: function (name) {
    this.mlist.Add(jsMath.mItem.TextAtom('ord','&nbsp;','normal'));
  },

  /*
   *  Implements \llap, \rlap, etc.
   */
  HandleLap: function (name) {
    var box = this.ProcessArg(); if (this.error) return;
    box = this.mlist.Add(new jsMath.mItem('lap',{nuc: box, lap: name}));
  },

  /*
   *  Adds the argument as a specific type of atom (for commands like
   *  \overline, etc.)
   */
  HandleAtom: function (name,data) {
    var arg = this.ProcessArg(this.cmd+name); if (this.error) return;
    this.mlist.Add(jsMath.mItem.Atom(data[0],arg));
  },


  /*
   *  Process the character associated with a specific \mathcharcode
   */
  HandleMathCode: function (name,code) {
    this.HandleTeXchar(code[0],code[1],code[2]);
  },

  /*
   *  Add a specific character from a TeX font (use the current
   *  font if the type is 7 (variable) or the font is not specified)
   *  Load the font if it is not already loaded.
   */
  HandleTeXchar: function (type,font,code) {
    if (type == 7 && this.mlist.data.font != null) {font = this.mlist.data.font}
    font = jsMath.TeX.fam[font];
    if (!jsMath.TeX[font]) {
      jsMath.TeX[font] = [];
      this.Extension(null,[jsMath.Font.URL(font)]);
    } else {
      this.mlist.Add(jsMath.mItem.TeXAtom(jsMath.TeX.atom[type],code,font));
    }
  },

  /*
   *  Add a TeX variable character or number
   */
  HandleVariable: function (c) {this.HandleTeXchar(7,1,c.charCodeAt(0))},
  HandleNumber: function (c) {this.HandleTeXchar(7,0,c.charCodeAt(0))},

  /*
   *  For unmapped characters, just add them in as normal
   *  (non-TeX) characters
   */
  HandleOther: function (c) {
    this.mlist.Add(jsMath.mItem.TextAtom('ord',c,'normal'));
  },

  /*
   *  Ignore comments in TeX data
   *  ### Some browsers remove the newlines, so this might cause
   *      extra stuff to be ignored; look into this ###
   */
  HandleComment: function () {
    var c;
    while (this.i < this.string.length) {
      c = this.string.charAt(this.i++);
      if (c == "\r" || c == "\n") return;
    }
  },

  /*
   *  Add a style change (e.g., \displaystyle, etc)
   */
  HandleStyle: function (name,style) {
    this.mlist.data.style = style[0];
    this.mlist.Add(new jsMath.mItem('style',{style: style[0]}));
  },

  /*
   *  Implements \small, \large, etc.
   */
  HandleSize: function (name,size) {
    this.mlist.data.size = size[0];
    this.mlist.Add(new jsMath.mItem('size',{size: size[0]}));
  },

  /*
   *  Set the current font (e.g., \rm, etc)
   */
  HandleFont: function (name,font) {
    this.mlist.data.font = font[0];
  },

  /*
   *  Look for and process a control sequence
   */
  HandleCS: function () {
    var cmd = this.GetCommand(); if (this.error) return;
    if (this.macros[cmd]) {
      var macro = this.macros[cmd];
      if (typeof(macro) == "string") {macro = [macro]}
      this[macro[0]](cmd,macro.slice(1)); return;
    }
    if (this.mathchardef[cmd]) {
      this.HandleMathCode(cmd,this.mathchardef[cmd]);
      return;
    }
    if (this.delimiter[this.cmd+cmd]) {
      this.HandleMathCode(cmd,this.delimiter[this.cmd+cmd].slice(0,3))
      return;
    }
    this.Error("Unknown control sequence '"+this.cmd+cmd+"'");
  },

  /*
   *  Process open and close braces
   */
  HandleOpen: function () {this.mlist.Open()},
  HandleClose: function () {
    if (this.mlist.data.openI == null) {this.Error("Extra close brace"); return}
    var open = this.mlist.Get(this.mlist.data.openI);
    if (!open || open.left == null) {this.mlist.Close()}
      else {this.Error("Extra close brace or missing "+this.cmd+"right"); return}
  },

  /*
   *  Implements \left
   */
  HandleLeft: function (name) {
    var left = this.GetDelimiter(this.cmd+name); if (this.error) return;
    this.mlist.Open(left);
  },

  /*
   *  Implements \right
   */
  HandleRight: function (name) {
    var right = this.GetDelimiter(this.cmd+name); if (this.error) return;
    var open = this.mlist.Get(this.mlist.data.openI);
    if (open && open.left != null) {this.mlist.Close(right)}
      else {this.Error("Extra open brace or missing "+this.cmd+"left");}
  },

  /*
   *  Implements generalized fractions (\over, \above, etc.)
   */
  HandleOver: function (name,data) {
    if (this.mlist.data.overI != null)
      {this.Error('Ambiguous use of '+this.cmd+name); return}
    this.mlist.data.overI = this.mlist.Length();
    this.mlist.data.overF = {name: name};
    if (data.length > 0) {
      this.mlist.data.overF.left  = this.delimiter[data[0]];
      this.mlist.data.overF.right = this.delimiter[data[1]];
    } else if (name.match(/withdelims$/)) {
      this.mlist.data.overF.left  = this.GetDelimiter(this.cmd+name); if (this.error) return;
      this.mlist.data.overF.right = this.GetDelimiter(this.cmd+name); if (this.error) return;
    } else {
      this.mlist.data.overF.left  = null;
      this.mlist.data.overF.right = null;
    }
    if (name.match(/^above/)) {
      this.mlist.data.overF.thickness = this.GetDimen(this.cmd+name,1);
      if (this.error) return;
    } else {
      this.mlist.data.overF.thickness = null;
    }
  },

  /*
   *  Add a superscript to the preceeding atom
   */
  HandleSuperscript: function () {
    var base = this.mlist.Last();
    if (this.mlist.data.overI == this.mlist.Length()) {base = null}
    if (base == null || (!base.atom && base.type != 'box' && base.type != 'frac'))
       {base = this.mlist.Add(jsMath.mItem.Atom('ord',{type:null}))}
    if (base.sup) {
      if (base.sup.isPrime) {base = this.mlist.Add(jsMath.mItem.Atom('ord',{type:null}))}
        else {this.Error("Double exponent: use braces to clarify"); return}
    }
    base.sup = this.ProcessScriptArg('superscript'); if (this.error) return;
  },

  /*
   *  Add a subscript to the preceeding atom
   */
  HandleSubscript: function () {
    var base = this.mlist.Last();
    if (this.mlist.data.overI == this.mlist.Length()) {base = null}
    if (base == null || (!base.atom && base.type != 'box' && base.type != 'frac'))
       {base = this.mlist.Add(jsMath.mItem.Atom('ord',{type:null}))}
    if (base.sub) {this.Error("Double subscripts: use braces to clarify"); return}
    base.sub = this.ProcessScriptArg('subscript'); if (this.error) return;
  },

  /*
   *  Parse a TeX math string, handling macros, etc.
   */
  Parse: function () {
    var c;
    while (this.i < this.string.length) {
      c = this.string.charAt(this.i++);
      if (this.mathchar[c]) {this.HandleMathCode(c,this.mathchar[c])}
      else if (this.special[c]) {this[this.special[c]](c)}
      else if (this.letter.test(c)) {this.HandleVariable(c)}
      else if (this.number.test(c)) {this.HandleNumber(c)}
      else {this.HandleOther(c)}
    }
    if (this.mlist.data.openI != null) {
      var open = this.mlist.Get(this.mlist.data.openI);
      if (open.left) {this.Error("Missing "+this.cmd+"right")}
        else {this.Error("Missing close brace")}
    }
    if (this.mlist.data.overI != null) {this.mlist.Over()}
  },

  /*
   *  Perform the processing of Appendix G
   */
  Atomize: function () {
    var data = this.mlist.init;
    if (!this.error) this.mlist.Atomize(data.style,data.size)
  },

  /*
   *  Produce the final HTML.
   *
   *  We have to wrap the HTML it appropriate <SPAN> tags to hide its
   *  actual dimensions when these don't match the TeX dimensions of the
   *  results.  We also include an image to force the results to take up
   *  the right amount of space.  The results may need to be vertically
   *  adjusted to make the baseline appear in the correct place.
   */
  Typeset: function () {
    var data = this.mlist.init;
    var box = this.typeset = this.mlist.Typeset(data.style,data.size);
    if (this.error) {return '<span class="error">'+this.error+'</span>'}
    if (box.format == 'null') {return ''};

    box.Styled().Remeasured(); var isSmall = 0; var isBig = 0;
    if (box.bh > box.h && box.bh > jsMath.h+.001) {isSmall = 1}
    if (box.bd > box.d && box.bd > jsMath.d+.001) {isSmall = 1}
    if (box.h > jsMath.h || box.d > jsMath.d) {isBig = 1}

    var html = box.html;
    if (isSmall) {// hide the extra size
      if (jsMath.Browser.allowAbsolute) {
        var y = 0;
        if (box.bh > jsMath.h+.001) {y = jsMath.h - box.bh}
        html = jsMath.HTML.Absolute(html,box.w,jsMath.h,0,y,jsMath.h);
      } else if (jsMath.Browser.valignBug) {
        // remove line height
        html = '<span style="line-height:'+jsMath.HTML.Em(jsMath.d)+';">'
             +    html + '</span>';
      } else {
        // remove line height and try to hide the depth
        var dy = jsMath.HTML.Em(Math.max(0,box.bd-jsMath.hd)/3);
        html = '<span style="line-height:'+jsMath.HTML.Em(jsMath.d)+';'
               + ' position:relative; top:'+dy+'; vertical-align:'+dy
               + '">' + html + '</span>';
      }
      isBig = 1;
    }
    if (isBig) {
      // add height and depth to the line
      //   (force a little extra to separate lines if needed)
      html += jsMath.HTML.Blank(0,box.h+.05,box.d+.05);
    }
    return '<nobr><span class="scale">'+html+'</span></nobr>';
  }

});

/*
 *  Make these characters special (and call the given routines)
 */
jsMath.Parser.prototype.AddSpecial({
  cmd:   'HandleCS',
  open:  'HandleOpen',
  close: 'HandleClose'
});


/*
 *  The web-page author can call jsMath.Macro to create additional
 *  TeX macros for use within his or her mathematics.  See the
 *  author's documentation for more details.
 */

jsMath.Add(jsMath,{
  Macro: function (name) {
    var macro = jsMath.Parser.prototype.macros;
    macro[name] = ['Macro'];
    for (var i = 1; i < arguments.length; i++)
      {macro[name][macro[name].length] = arguments[i]}
  }
});

/*
 *  Use these commands to create macros that load
 *  JavaScript files and reprocess the mathematics when
 *  the file is loaded.  This lets you to have macros or
 *  LaTeX environments that autoload their own definitions
 *  only when they are needed, saving initial download time
 *  on pages where they are not used.  See the author's
 *  documentation for more details.
 *
 */

jsMath.Extension = {

  safeRequire: 1,   // disables access to files outside of jsMath/extensions

  Macro: function (name,file) {
    var macro = jsMath.Parser.prototype.macros;
    if (file == null) {file = name}
    macro[name] = ['Extension',file];
  },

  LaTeX: function (env,file) {
    var latex = jsMath.Parser.prototype.environments;
    latex[env] = ['Extension',file,'environments'];
  },

  Font: function (name,font) {
    if (font == null) {font = name + "10"}
    var macro = jsMath.Parser.prototype.macros;
    macro[name] = ['Extension',jsMath.Font.URL(font)];
  },

  MathChar: function (font,defs) {
    var fam = jsMath.TeX.famName[font];
    if (fam == null) {
      fam = jsMath.TeX.fam.length;
      jsMath.TeX.fam[fam] = font;
      jsMath.TeX.famName[font] = fam;
    }
    var mathchardef = jsMath.Parser.prototype.mathchardef;
    for (var c in defs) {mathchardef[c] = [defs[c][0],fam,defs[c][1]]}
  },

  Require: function (file,show) {
    if (this.safeRequire && (file.match(/\.\.\/|[^-a-z0-9.\/:_+=%~]/i) ||
         (file.match(/:/) && file.substr(0,jsMath.root.length) != jsMath.root))) {
      jsMath.Setup.loaded[file] = 1;
      return;
    }
    jsMath.Setup.Script(this.URL(file),show);
  },

  URL: function (file) {
    file = file.replace(/^\s+|\s+$/g,'');
    if (!file.match(/^([a-z]+:|\/|fonts|extensions\/)/i)) {file = 'extensions/'+file}
    if (!file.match(/\.js$/)) {file += '.js'}
    return file;
  }
}


/***************************************************************************/

/*
 *  These routines look through the web page for math elements to process.
 *  There are two main entry points you can call:
 *
 *      <script> jsMath.Process() </script>
 *  or
 *      <script> jsMath.ProcessBeforeShowing() </script>
 *
 *  The first will process the page asynchronously (so the user can start
 *  reading the top of the file while jsMath is still processing the bottom)
 *  while the second does not update until all the mathematics is typeset.
 */

jsMath.Add(jsMath,{
  /*
   *  Call this at the bottom of your HTML page to have the
   *  mathematics typeset asynchronously.  This lets the user
   *  start reading the mathematics while the rest of the page
   *  is being processed.
   */
  Process: function (obj) {
    jsMath.Setup.Body();
    jsMath.Script.Push(jsMath.Translate,'Asynchronous',obj);
  },

  /*
   *  Call this at the bottom of your HTML page to have the
   *  mathematics typeset before the page is displayed.
   *  This can take a long time, so the user could cancel the
   *  page before it is complete; use it with caution, and only
   *  when there is a relatively small amount of math on the page.
   */
  ProcessBeforeShowing: function (obj) {
    jsMath.Setup.Body();
    var method = (jsMath.Controls.cookie.asynch ? "Asynchronous": "Synchronous");
    jsMath.Script.Push(jsMath.Translate,method,obj);
  }

});

jsMath.Translate = {

  element: [],  // the list of math elements on the page
  cancel: 0,    // set to 1 to cancel asynchronous processing

  /*
   *  Parse a TeX string in Text or Display mode and return
   *  the HTML for it (taking it from the cache, if available)
   */
  Parse: function (style,s,noCache) {
    var cache = jsMath.Global.cache[style];
    if (!cache[jsMath.em]) {cache[jsMath.em] = {}}
    var HTML = cache[jsMath.em][s];
    if (!HTML || noCache) {
      var parse = jsMath.Parse(s,null,null,style);
      parse.Atomize(); HTML = parse.Typeset();
      if (!noCache) {cache[jsMath.em][s] = HTML}
    }
    return HTML;
  },

  TextMode:    function (s,noCache) {this.Parse('T',s,noCache)},
  DisplayMode: function (s,noCache) {this.Parse('D',s,noCache)},

  /*
   *  Return the text of a given DOM element
   */
  GetElementText: function (element) {
    var text = this.recursiveElementText(element);
    element.alt = text;
    if (text.search('&') >= 0) {
      text = text.replace(/&lt;/g,'<');
      text = text.replace(/&gt;/g,'>');
      text = text.replace(/&quot;/g,'"');
      text = text.replace(/&amp;/g,'&');
    }
    return text;
  },
  recursiveElementText: function (element) {
    if (element.nodeValue != null) {return element.nodeValue}
    if (element.childNodes.length == 0) {return " "}
    var text = '';
    for (var i = 0; i < element.childNodes.length; i++)
      {text += this.recursiveElementText(element.childNodes[i])}
    return text;
  },

  /*
   *  Move hidden to the location of the math element to be
   *  processed and reinitialize sizes for that location.
   */
  ResetHidden: function (element) {
    element.innerHTML =
      '<span class="jsMath_hiddenSpan" style="position:absolute"></span>'
        + jsMath.Browser.operaHiddenFix; // needed by Opera in tables
    element.className = '';
    jsMath.hidden = element.firstChild;
    if (!jsMath.BBoxFor("x").w) {jsMath.hidden = jsMath.hiddenTop}
    jsMath.ReInit();
  },


  /*
   *  Typeset the contents of an element in \textstyle
   */
  ConvertText: function (element,noCache) {
    var text = this.GetElementText(element);
    this.ResetHidden(element);
    if (text.match(/^\s*\\nocache([^a-zA-Z])/))
      {noCache = true; text = text.replace(/\s*\\nocache/,'')}
    text = this.Parse('T',text,noCache);
    element.className = 'typeset';
    element.innerHTML = text;
  },

  /*
   *  Typeset the contents of an element in \displaystyle
   */
  ConvertDisplay: function (element,noCache) {
    var text = this.GetElementText(element);
    this.ResetHidden(element);
    if (text.match(/^\s*\\nocache([^a-zA-Z])/))
      {noCache = true; text = text.replace(/\s*\\nocache/,'')}
    text = this.Parse('D',text,noCache);
    element.className = 'typeset';
    element.innerHTML = text;
  },

  /*
   *  Process a math element
   */
  ProcessElement: function (element) {
    this.restart = 0;
    var noCache = (element.className.toLowerCase().match(/(^| )nocache( |$)/) != null);
    try {
      if (element.tagName.toLowerCase() == 'div') {
        this.ConvertDisplay(element,noCache);
        element.onclick = jsMath.Click.CheckClick;
        element.ondblclick = jsMath.Click.CheckDblClick;
      } else if (element.tagName.toLowerCase() == 'span') {
        this.ConvertText(element,noCache);
        element.onclick = jsMath.Click.CheckClick;
        element.ondblclick = jsMath.Click.CheckDblClick;
      }
    } catch (err) {
      if (element.alt) {
        var tex = element.alt;
        tex = tex.replace(/&/g,'&amp;');
        tex = tex.replace(/</g,'&lt;');
        tex = tex.replace(/>/g,'&gt;');
        element.innerHTML = tex;
        element.className = 'math';
        if (noCache) {element.className += ' nocache'}
      }
      jsMath.hidden = jsMath.hiddenTop;
    }
  },

  /*
   *  Asynchronously process all the math elements starting with
   *  the k-th one
   */
  ProcessElements: function (k) {
    jsMath.Script.blocking = 1;
    if (k >= this.element.length || this.cancel) {
      this.ProcessComplete();
      if (this.cancel) {
        jsMath.Message.Set("Process Math: Canceled");
        jsMath.Message.Clear()
      }
      jsMath.Script.blocking = 0;
      jsMath.Script.Process();
    } else {
      this.ProcessElement(this.element[k]);
      if (this.restart) {
        jsMath.Script.Push(this,'ProcessElements',k);
        jsMath.Script.blocking = 0;
        setTimeout('jsMath.Script.Process()',jsMath.Browser.delay);
      } else {
        k++; var p = Math.floor(100 * k / this.element.length);
        jsMath.Message.Set('Processing Math: '+p+'%');
        setTimeout('jsMath.Translate.ProcessElements('+k+')',jsMath.Browser.delay);
      }
    }
  },

  /*
   *  Start the asynchronous processing of mathematics
   */
  Asynchronous: function (obj) {
    if (!jsMath.initialized) {jsMath.Init()}
    this.element = this.GetMathElements(obj);
    this.cancel = 0; this.asynchronous = 1;
    jsMath.Script.blocking = 1;
    jsMath.Message.Set('Processing Math: 0%',1);
    setTimeout('jsMath.Translate.ProcessElements(0)',jsMath.Browser.delay);
  },

  /*
   *  Do synchronous processing of mathematics
   */
  Synchronous: function (obj,i) {
    if (i == null) {
      if (!jsMath.initialized) {jsMath.Init()}
      this.element = this.GetMathElements(obj);
      i = 0;
    }
    this.asynchronous = 0;
    while (i < this.element.length) {
      this.ProcessElement(this.element[i]);
      if (this.restart) {
        jsMath.Synchronize('jsMath.Translate.Synchronous(null,'+i+')');
        jsMath.Script.Process();
        return;
      }
      i++;
    }
    this.ProcessComplete(1);
  },

  /*
   *  Look up all the math elements on the page and
   *  put them in a list sorted from top to bottom of the page
   */
  GetMathElements: function (obj) {
    var element = [];
    if (!obj) {obj = jsMath.document}
    if (typeof(obj) == 'string') {obj = jsMath.document.getElementById(obj)}
    if (!obj.getElementsByTagName) return null;
    var math = obj.getElementsByTagName('div');
    for (var k = 0; k < math.length; k++) {
      if (math[k].className && math[k].className.match(/(^| )math( |$)/)) {
        if (jsMath.Browser.renameOK && obj.getElementsByName)
               {math[k].setAttribute('name','_jsMath_')}
          else {element[element.length] = math[k]}
      }
    }
    math = obj.getElementsByTagName('span');
    for (var k = 0; k < math.length; k++) {
      if (math[k].className && math[k].className.match(/(^| )math( |$)/)) {
        if (jsMath.Browser.renameOK && obj.getElementsByName)
               {math[k].setAttribute('name','_jsMath_')}
          else {element[element.length] = math[k]}
      }
    }
    // this gets the SPAN and DIV elements interleaved in order
    if (jsMath.Browser.renameOK && obj.getElementsByName) {
      element = obj.getElementsByName('_jsMath_');
    } else if (jsMath.hidden.sourceIndex) {
      element.sort(function (a,b) {return a.sourceIndex - b.sourceIndex});
    }
    return element;
  },

  /*
   *  Remove the window message about processing math
   *  and clean up any marked <SPAN> or <DIV> tags
   */
  ProcessComplete: function (noMessage) {
    if (jsMath.Browser.renameOK) {
      var element = jsMath.document.getElementsByName('_jsMath_');
      for (var i = element.length-1; i >= 0; i--) {
        element[i].removeAttribute('name');
      }
    }
    jsMath.hidden = jsMath.hiddenTop;
    this.element = [];
    if (!noMessage) {
      jsMath.Message.Set('Processing Math: Done');
      jsMath.Message.Clear();
    }
    jsMath.Message.UnBlank();
    if (jsMath.Browser.safariImgBug &&
        (jsMath.Controls.cookie.font == 'symbol' ||
         jsMath.Controls.cookie.font == 'image')) {
      //
      //  For Safari, the images don't always finish
      //  updating, so nudge the window to cause a
      //  redraw.  (Hack!)
      //
      if (this.timeout) {clearTimeout(this.timeout)}
      this.timeout = setTimeout("jsMath.window.resizeBy(-1,0); "
                              + "jsMath.window.resizeBy(1,0); "
                              + "jsMath.Translate.timeout = null",2000);
    }
  },

  /*
   *  Cancel procesing elements
   */
  Cancel: function () {
    jsMath.Translate.cancel = 1;
    if (jsMath.Script.cancelTimer) {jsMath.Script.cancelLoad()}
  }

};

jsMath.Add(jsMath,{
  //
  //  Synchronize these with the loading of the tex2math plugin.
  //
  ConvertTeX: function (element) {jsMath.Script.Push(jsMath.tex2math,'ConvertTeX',element)},
  ConvertTeX2: function (element) {jsMath.Script.Push(jsMath.tex2math,'ConvertTeX2',element)},
  ConvertLaTeX: function (element) {jsMath.Script.Push(jsMath.tex2math,'ConvertLaTeX',element)},
  ConvertCustom: function (element) {jsMath.Script.Push(jsMath.tex2math,'ConvertCustom',element)},
  CustomSearch: function (om,cm,od,cd) {jsMath.Script.Push(null,function () {jsMath.tex2math.CustomSearch(om,cm,od,cd)})},
  tex2math: {
    ConvertTeX: function () {},
    ConvertTeX2: function () {},
    ConvertLaTeX: function () {},
    ConvertCustom: function () {},
    CustomSearch: function () {}
  }
});
jsMath.Synchronize = jsMath.Script.Synchronize;

/***************************************************************************/


/*
 *  Initialize things
 */
try {
  if (window.parent != window && window.jsMathAutoload) {
    window.parent.jsMath = jsMath;
    jsMath.document = window.parent.document;
    jsMath.window = window.parent;
  }
} catch (err) {}

jsMath.Global.Register();
jsMath.Loaded();
jsMath.Controls.GetCookie();
jsMath.Setup.Source();
jsMath.Global.Init();
jsMath.Script.Init();
jsMath.Setup.Fonts();
if (jsMath.document.body) {jsMath.Setup.Body()}
jsMath.Setup.User("onload");

}}
