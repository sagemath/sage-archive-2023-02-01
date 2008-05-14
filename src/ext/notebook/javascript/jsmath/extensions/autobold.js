/*
 *  extensions/autobold.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file causes jsMath to use \boldsymbol{...} around mathematics
 *  that appears within <B>...</B> tags or has font-weight:bold applied
 *  via CSS rule.  You can activate it by calling
 *
 *    jsMath.Extension.Require('autobold');
 *
 *  once jsMath.js has been loaded, or by adding "extensions/autobold.js"
 *  to the loadFiles array in jsMath/easy/load.js.
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
 *  Copyright 2008 by Davide P. Cervone
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

jsMath.Extension.Require("boldsymbol");

jsMath.Translate.OldParse = jsMath.Translate.Parse;
jsMath.Translate.Parse = function (style,text,noCache) {
  if (jsMath.BBoxFor('</SPAN></SPAN>MMMMMMMMMM<SPAN><SPAN>').w >
      jsMath.BBoxFor('</SPAN></SPAN><SPAN STYLE="font-weight:normal">MMMMMMMMMM</SPAN><SPAN><SPAN>').w) {
    text = '\\boldsymbol{' + text + '}';
  }
  return jsMath.Translate.OldParse(style,text,noCache);
}
