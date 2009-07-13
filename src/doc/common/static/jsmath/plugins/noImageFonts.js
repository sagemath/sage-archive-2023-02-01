/*
 *  noImageFonts.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file indicates that no image fonts are available.
 *  It should be loaded BEFORE jsMath.js is loaded.
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2004-2006 by Davide P. Cervone
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
jsMath.noImgFonts = 1;

if (!jsMath.Font) {jsMath.Font = {}}
if (!jsMath.Font.extra_message) {
  jsMath.Font.extra_message =
    'Extra TeX fonts not found: <b><span id="jsMath_ExtraFonts"></span></b><br/>'
      + 'Using unicode fonts instead.  This may be slow and might not print well.<br/>\n'
      + 'Use the jsMath control panel to get additional information.';
}