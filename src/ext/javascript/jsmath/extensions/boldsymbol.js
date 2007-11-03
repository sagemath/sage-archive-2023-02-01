/*
 *  extensions/boldsymbol.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements the \boldsymbol macro.  You can activate it
 *  by calling
 *
 *    jsMath.Extension.Macro('boldsymbol');
 *
 *  which will cause the extension to be loaded only when it is
 *  needed, or you can force it to be loaded via
 *
 *    jsMath.Extension.Require('boldsymbol');
 *
 *  once jsMath.js has been loaded.
 *
 *  Note that you will need to install the cmmib10 and cmbsy10 fonts
 *  that are available from the jsMath extra font page at
 *
 *      http://www.math.union.edu/locate/jsMath/download/extra-fonts/
 *
 *  to make this work in image mode.  Note that there is no unicode
 *  fallback for these fonts at the moment.
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

jsMath.Package(jsMath.Parser,{

  macros: {boldsymbol: 'BoldSymbol'},

  /*
   *  Implement \boldsymbol{...}
   */
  BoldSymbol: function (name) {
    var fam = jsMath.TeX.fam
    var oldfam = [fam[0],fam[1],fam[2]];
    fam[0] = "cmbx10"; fam[1] = "cmmib10"; fam[2] = "cmbsy10";
    var box = this.ProcessArg(this.cmd+name);
    fam[0] = oldfam[0]; fam[1] = oldfam[1]; fam[2] = oldfam[2];
    if (this.error) return;
    this.mlist.Add(jsMath.mItem.Atom('ord',box));
  }

});
