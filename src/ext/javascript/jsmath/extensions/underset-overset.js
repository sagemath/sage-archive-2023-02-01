/*
 *  extensions/underset-overset.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements \underset and \overset macros.  It will be loaded
 *  automatically when needed, or can be loaded by
 *
 *    jsMath.Extension.Require('underset-overset');
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 200-20065 by Davide P. Cervone
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

jsMath.Package(jsMath.Parser,{

  macros: {
    overset:            'Overset',
    underset:           'Underset'
  },

  Overset: function (name) {
    var top = this.ProcessArg(this.cmd+name); if (this.error) return;
    var bot = this.ProcessArg(this.cmd+name); if (this.error) return;
    var op = jsMath.mItem.Atom('op',bot);
    op.limits = 1; op.sup = top;
    this.mlist.Add(op);
  },

  Underset: function (name) {
    var bot = this.ProcessArg(this.cmd+name); if (this.error) return;
    var top = this.ProcessArg(this.cmd+name); if (this.error) return;
    var op = jsMath.mItem.Atom('op',top);
    op.limits = 1; op.sub = bot;
    this.mlist.Add(op);
  }

});
