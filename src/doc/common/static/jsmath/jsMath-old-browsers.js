/*
 *  jsMath-old-browsers.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes changes needed by older versions of some browsers
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

jsMath.Add(jsMath.HTML,{
  /*
   *  Use the blank GIF image for spacing and rules
   */
  Blank: function (w,h,d,isRule) {
    var style = '';
    if (isRule) {
      if (h*jsMath.em < 1.5) {h = '1px'} else {h = jsMath.HTML.Em(h)}
      style = 'border-top:'+h+' solid;'; h = 0;
    }
    if (d == null) {d = 0}
    style += 'width:'+this.Em(w)+'; height:'+this.Em(h+d)+';';
    if (d) {style += 'vertical-align:'+this.Em(-d)}
    return '<img src="'+jsMath.blank+'" style="'+style+'" />';
  }
});

if (jsMath.browser == 'Konqueror') {

  jsMath.Package(jsMath.Box,{Remeasured: function() {return this}});

  jsMath.Add(jsMath.HTML,{
    Spacer: function (w) {
      if (w == 0) {return ''};
      return '<span style="margin-left:'+this.Em(w-jsMath.Browser.spaceWidth)+'">'
             + '&nbsp;</span>';
    }
  });

  jsMath.Browser.spaceWidth = this.EmBoxFor('&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;').w/5;

}

jsMath.styles['.typeset .spacer'] = '';
