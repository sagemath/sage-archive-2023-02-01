/*
 *  extensions/bbox.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements the \bbox macro, which creates an HTML box that
 *  can be styled (for background colors, and so on).  You can include
 *  an optional dimension that tells how much extra padding to include
 *  around the bounding box for the mathematics.  E.g.,
 *
 *    \bbox[2pt]{x+y}        %  an invisible box around x+y with 2pt of extra space
 *    \bbox[green]{x+y}      %  a green box around x+y
 *    \bbox[green,2pt]{x+y}  %  a green box with 2pt of extra space
 *    \bbox[yellow,2pt,border:1px solid red]{x+y}
 *                           %  a yellow box with a red border and 2pt space
 *
 *  This extension is loaded automatically when needed, or you can call
 *  it directly via
 *
 *    jsMath.Extension.Require('bbox');
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2006 by Davide P. Cervone
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

jsMath.Add(jsMath.HTML,{

  /*
   *  Create a colored bbounding box
   */
  BBox: function (w,h,d,c,s) {
    if (w <= 0) {return ''}
    if (d == null) {d = 0}
    var style = (jsMath.Browser.msieInlineBlockFix ? '' : 'overflow:visible;');
    style += 'width:'+this.Em(w)+'; height:'+this.Em(h+d)+';';
    if (jsMath.Browser.mozInlineBlockBug) {d = -h}
    if (jsMath.Browser.msieInlineBlockFix) {d -= jsMath.d}
    if (d) {style += ' vertical-align:'+this.Em(-d)+';'}
    if (c) {style += ' background-color:'+c+';'}
    var html = '<span class="blank" style="'+style+s+'"></span>';
    return html;
  }

});

jsMath.Add(jsMath.mList.prototype.Atomize,{
  /*
   *  Creates the box HTML
   */
  bbox: function (style,size,mitem,prev,mlist) {
    var box; var w; var h; var d;
    var nuc = mitem.nuc = jsMath.Box.Set(mitem.nuc,style,size).Remeasured();
    if (box == null) {w = nuc.w; h = nuc.h; d = nuc.d} // values before super/subs-cript
    var nuc = mitem.nuc; nuc.Styled(); var pad = mitem.pad;
    if (pad) {w += 2*pad; h += pad; d += pad; nuc.w += pad}
    if (jsMath.Browser.msieCenterBugFix)
      {nuc.html = '<span style="position:relative">'+nuc.html+'</span>'}
    nuc.html =
      jsMath.HTML.BBox(w,h,d,mitem.color,mitem.style) +
      jsMath.HTML.Spacer(pad-w) +
      nuc.html;
    nuc.Remeasured();
    if (pad && nuc.w < w) {
      nuc.html += jsMath.HTML.Spacer(w-nuc.w);
      nuc.w = w;
    }
    nuc.h  = Math.max(nuc.h,h);  nuc.d  = Math.max(nuc.d,d);
    nuc.bh = Math.max(nuc.bh,nuc.h); nuc.bd = Math.max(nuc.bd,nuc.d);
    mitem.type = 'ord';
    jsMath.mList.prototype.Atomize.SupSub(style,size,mitem);
  }
});

jsMath.Package(jsMath.Parser,{

  macros: {bbox: 'BBox'},

  /*
   *  Implement \bbox[...]{...}
   */
  BBox: function (name) {
    var extra = this.GetBrackets(this.cmd+name); if (this.error) return;
    var arg = this.GetArgument(this.cmd+name); if (this.error) return;
    var nuc = this.Process(arg); if (this.error) return;
    var color; var pad = 0; var style = '';
    if (extra != '') {
      var parts = extra.split(/,/);
      for (var i in parts) {
        if (parts[i].match(/^\s*([-+]?(\.\d+|\d+(\.\d*)?))(pt|em|ex|mu|px)\s*$/))
          {pad = this.ParseDimen(parts[i],'',0,1)}
        else if (parts[i].match(/:/)) {style = parts[i]}
        else {color = parts[i]}
      }
    }
    var atom = {nuc: nuc, atom: 1, pad: pad, color: color, style: style};
    this.mlist.Add(new jsMath.mItem('bbox',atom));
  }
});
