/*
 *  extensions/leaders.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements the \overbrace, \underbrace, \overrightarrow
 *  and \overleftarrow macros. It will be loaded automatically when needed,
 *  or can be loaded by
 *
 *    jsMath.Extension.Require('leaders');
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

jsMath.Add(jsMath.Box,{

  /*
   *  Create a horizontally stretchable "delimiter" (like over- and
   *  underbraces).
   */
//###  Add size?
  Leaders: function (W,leader) {
    var h; var d; var w; var html; var font;
    if (leader.lmid) {// braces
      font = jsMath.TeX.fam[leader.left[0]];
      var left = this.GetCharCode(leader.left);
      var right = this.GetCharCode(leader.right);
      var lmid = this.GetCharCode(leader.lmid);
      var rmid = this.GetCharCode(leader.rmid);
      w = (W - left.w - right.w - lmid.w - rmid.w)/2 - .1; h = .4; d = .3;
      if (w < 0) {w = 0}
      html = this.AddClass(left.tclass,left.c,left.font)
           + jsMath.HTML.Rule(w,left.h)
           + this.AddClass(lmid.tclass,lmid.c+rmid.c,lmid.font)
           + jsMath.HTML.Rule(w,right.h)
           + this.AddClass(right.tclass,right.c,right.font);
    } else { //arrows
      font = jsMath.TeX.fam[leader.rep[0]];
      var left = this.GetCharCode(leader.left? leader.left: leader.rep);
      var rep = this.GetCharCode(leader.rep);
      var right = this.GetCharCode(leader.right? leader.right: leader.rep);
      var n = Math.ceil((W - left.w - right.w + .4)/(rep.w - .3));
      w = (W - left.w - right.w + .4 - n*(rep.w - .3));
      if (leader.left) {h = left.h; d = left.d} else {h = right.h; d = right.d}
      if (d == null) {d = 0}; if (h == null) {h = 0}
      var html = this.AddClass(left.tclass,left.c,left.font); var m = Math.floor(n/2);
      var ext = jsMath.HTML.Place(rep.c,-.3,0);
      var ehtml = ''; for (var i = 0; i < m; i++) {ehtml += ext};
      html += this.AddClass(rep.tclass,ehtml,rep.font) + jsMath.HTML.Spacer(w);
      ehtml = ''; for (var i = m; i < n; i++) {ehtml += ext};
      html += this.AddClass(rep.tclass,ehtml,rep.font);
      if (jsMath.Browser.msieFontBug) {html += '<span style="display: none">x</span>'}
      html += jsMath.HTML.Place(this.AddClass(right.tclass,right.c,right.font),-.4,0);
    }
    w = jsMath.EmBoxFor(html).w;
    if (w != W) {
      w = jsMath.HTML.Spacer((W-w)/2);
      html = w + html + w;
    }
    var box = new jsMath.Box('html',html,W,h,d);
    box.bh = jsMath.TeX[font].h; box.bd = jsMath.TeX[font].d;
    return box;
  }

});

jsMath.Package(jsMath.Parser,{

  macros: {
    overbrace:       ['HandleLeaders','downbrace',1],
    underbrace:      ['HandleLeaders','upbrace',1,1],
    overrightarrow:  ['HandleLeaders','rightarrow'],
    overleftarrow:   ['HandleLeaders','leftarrow']
  },

  /*
   *  The horizontally stretchable delimiters
   */
  leaders: {
    downbrace:  {left: [3,0x7A], lmid: [3,0x7D], rmid: [3,0x7C], right: [3,0x7B]},
    upbrace:    {left: [3,0x7C], lmid: [3,0x7B], rmid: [3,0x7A], right: [3,0x7D]},
    leftarrow:  {left: [2,0x20], rep:   [2,0x00]},
    rightarrow: {rep:  [2,0x00], right: [2,0x21]}
  },

  /*
   *  Implements \overbrace, \underbrace, etc.
   */
  HandleLeaders: function (name,data) {
    var box = this.ProcessArg(this.cmd+name); if (this.error) return;
    box = jsMath.Box.Set(box,'D',this.mlist.data.size).Remeasured();
    var leader = jsMath.Box.Leaders(box.w,this.leaders[data[0]]);
    if (data[2]) {leader.y = -leader.h - box.d}
            else {leader.y = box.h + Math.max(0,leader.d)}
    box.x = -(leader.w + box.w)/2;
    var space = jsMath.Box.Space((leader.w-box.w)/2);
    box = jsMath.mItem.Atom(data[1]? 'op': 'inner',
      jsMath.Box.SetList([leader,box,space],'T',this.mlist.data.size));
    box.limits = (data[1]? 1: 0);
    this.mlist.Add(box);
  }

});
