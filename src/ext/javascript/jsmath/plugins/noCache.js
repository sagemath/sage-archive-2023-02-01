/*
 *  noCache.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file disables the equation cache that jsMath uses to
 *  store the typeset versions of TeX code so that common expressions
 *  won't need to be re-typeset.
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

jsMath.Add(jsMath,{

  /*
   *  Get the width and height (in ems) of an HTML string
   */
  EmBoxFor: function (s) {
    var bbox = this.BBoxFor(s);
    return {w: bbox.w/this.em, h: bbox.h/this.em};
  },

  /*
   *  For browsers that don't handle sizes of italics properly (MSIE)
   */
  EmBoxForItalics: function (s) {
    var bbox = this.BBoxFor(s);
    if (s.match(/<i>|class=\"(icm|italic|igreek|iaccent)/i)) {
      bbox.w = this.BBoxFor(s+jsMath.Browser.italicString).w
                - jsMath.Browser.italicCorrection;
    }
    return {w: bbox.w/this.em, h: bbox.h/this.em};
  }

});

jsMath.Add(jsMath.Translate,{

  /*
   *  Typeset a string in \textstyle and return the HTML for it
   */
  TextMode: function (s) {
    var parse = jsMath.Parse(s,null,null,'T');
    parse.Atomize();
    return parse.Typeset();
  },

  /*
   *  Typeset a string in \displaystyle and return the HTML for it
   */
  DisplayMode: function (s) {
    var parse = jsMath.Parse(s,null,null,'D');
    parse.Atomize();
    return parse.Typeset();
  }

});
