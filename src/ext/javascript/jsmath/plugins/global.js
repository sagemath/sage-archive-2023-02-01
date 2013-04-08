/*
 *  plugins/global.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file call up the global frame, if it is not already in place.
 *
 *  This should be called BEFORE loading jsMath.js
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

if (!parent.jsMath || !parent.jsMath.isGlobal) {
  var cookie = []; var cookies = document.cookies;
  if (window.location.protocol == 'file:') {cookies = unescape(window.location.search.substr(1))}
  else if (window.location.protocol == 'mk:') {cookies = unescape(window.location.hash.substr(1))}
  if (cookies.match(/jsMath=([^;]+)/)) {
    var data = RegExp.$1.split(/,/);
    for (var i = 0; i < data.length; i++) {
      var x = data[i].match(/(.*):(.*)/);
      cookie[x[1]] = x[2];
    }
  }
  if (cookie.global != "never" && !navigator.accentColorName) {
    var script = document.getElementsByTagName('script');
    if (script) {
      for (var i = 0; i < script.length; i++) {
        src = script[i].src;
        if (src && src.match('(^|/)plugins/global.js$')) {
          src = src.replace(/plugins\/global.js$/,'jsMath-global.html');
          var sep = (window.location.protocol == 'mk:') ? '#' : '?';
          window.location.replace(src + sep + escape(window.location));
          break;
        }
      }
    }
  }
}