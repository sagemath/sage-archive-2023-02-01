/*
 *  smallFonts.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file changes the sizes of fonts used in subscripts so that they
 *  are larger.  This can be helpful if jsMath is used on a page with a
 *  small font size, where the subscripts may tend to disappear.
 *  It should be loaded BEFORE jsMath.js.
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

if (!window.jsMath) {window.jsMath = {}}
if (!jsMath.styles) {jsMath.styles = {}}
if (!jsMath.Img) {jsMath.Img = {}}
if (!jsMath.Box) {jsMath.Box = {}}
if (!jsMath.Typeset) {jsMath.Typeset = {}}

/*
 *  Replace the smaller font sizes
 */
jsMath.sizes = [70, 77, 85, 92, 100, 120, 144, 173, 207, 249];
jsMath.Img.fonts = [70, 70, 85, 85, 100, 120, 144, 173, 207, 249, 298, 358, 430];

jsMath.styles['.typeset .size0'] = 'font-size: 70%';
jsMath.styles['.typeset .size1'] = 'font-size: 77%';
jsMath.styles['.typeset .size2'] = 'font-size: 85%';
jsMath.styles['.typeset .size3'] = 'font-size: 92%';

/*
 *  Fix multiplication factors in these routines
 */
jsMath.Typeset.StyleValue = function (style,v) {
  if (style == "S" || style == "S'")   {return .85*v}
  if (style == "SS" || style == "SS'") {return .70*v}
  return v;
};

jsMath.Box.DelimBestFit = function (H,c,font,style) {
  if (c == 0 && font == 0) return null;
  var C; var h; font = jsMath.TeX.fam[font];
  var isSS = (style.charAt(1) == 'S');
  var isS  = (style.charAt(0) == 'S');
  while (c != null) {
    C = jsMath.TeX[font][c];
    if (C.h == null) {C.h = jsMath.Box.defaultH}; if (C.d == null) {C.d = 0}
    h = C.h+C.d;
    if (C.delim) {return [c,font,'',H]}
    if (isSS && .70*h >= H) {return [c,font,'SS',.7*h]}
    if (isS  && .85*h >= H) {return [c,font,'S',.85*h]}
    if (h >= H || C.n == null) {return [c,font,'T',h]}
    c = C.n;
  }
  return null;
};
