/*
 *  extensions/mathchoice.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements the 4-way math choice.  It will be loaded
 *  automatically when needed, or can be loaded by
 *
 *    jsMath.Extension.Require('mathchoice');
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

jsMath.Add(jsMath.mList.prototype.Atomize,{

  /*
   *  Handle a 4-way choice atom.  (Rule 4)
   */
  choice: function (style,mitem,i,mlist) {
    if (style.charAt(style.length-1) == "'") {style = style.slice(0,style.length-1)}
    var nlist = []; var M = mitem[style];
    if (!M) {M = {type: 'mlist', mlist: {mlist: []}}}
    if (M.type == 'mlist') {
      M = M.mlist.mlist;
      for (var k = 0; k < i; k++) {nlist[k] = mlist[k]}
      for (k = 0; k < M.length; k++) {nlist[i+k] = M[k]}
      for (k = i+1; k < mlist.length; k++) {nlist[nlist.length] = mlist[k]}
      return nlist;
    } else {
      mlist[i] = jsMath.mItem.Atom('ord',M);
      return mlist;
    }
  }

});

jsMath.Package(jsMath.Parser,{

  macros: {mathchoice: 'MathChoice'},

  /*
   *  Implements \mathchoice{}{}{}{}
   */
  MathChoice: function (name) {
    var D  = this.ProcessArg(this.cmd+name); if (this.error) return;
    var T  = this.ProcessArg(this.cmd+name); if (this.error) return;
    var S  = this.ProcessArg(this.cmd+name); if (this.error) return;
    var SS = this.ProcessArg(this.cmd+name); if (this.error) return;
    var box = new jsMath.mItem('choice',{D: D, T: T, S: S, SS: SS});
    this.mlist.Add(new jsMath.mItem('choice',{D: D, T: T, S: S, SS: SS}));
  }

});

