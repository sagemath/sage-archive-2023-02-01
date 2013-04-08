/*
 *  extensions/moreArrows.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements additional arrow macros with under- and
 *  overset labels.  It can be loaded by
 *
 *    jsMath.Extension.Require('moreArrows');
 *
 *  or using \require{moreArrows} within a math formula.
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

jsMath.Extension.Require('leaders');

jsMath.Package(jsMath.Parser,{

  macros: {
    xrightarrow:       ['HandleArrows','rightarrow'],
    xleftarrow:        ['HandleArrows','leftarrow'],
    xuprightharpoon:   ['HandleArrows','uprightharpoon'],
    xupleftharpoon:    ['HandleArrows','upleftharpoon'],
    xdownrightharpoon: ['HandleArrows','downrightharpoon'],
    xdownleftharpoon:  ['HandleArrows','downleftharpoon']
  },

  leaders: {
    upleftharpoon:    {left: [1,0x28], rep:   [2,0x00]},
    uprightharpoon:   {rep:  [2,0x00], right: [1,0x2A]},
    downleftharpoon:  {left: [1,0x29], rep:   [2,0x00]},
    downrightharpoon: {rep:  [2,0x00], right: [1,0x2B]}
  },

  HandleArrows: function (name,data) {
    var bot = this.GetBrackets(this.cmd+name); if (this.error) return;
    var top = this.ProcessArg(this.cmd+name); if (this.error) return;
    var box = jsMath.Box.Set(top,'S',this.mlist.data.size).Remeasured();
    var w = box.w;
    if (bot) {
      bot = this.Process(bot); if (this.error) return;
      var box = jsMath.Box.Set(bot,'S',this.mlist.data.size).Remeasured();
      w = Math.max(w,box.w);
    }
    var leader = jsMath.Box.Leaders(w+.75,this.leaders[data[0]]);
    box = jsMath.mItem.Atom('op',jsMath.Box.SetList([leader],'T',this.mlist.data.size));
    box.limits = 1; box.sup = top; box.sub = bot;
    this.mlist.Add(box);
  }

});
