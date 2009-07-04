/*
 *  extensions/fbox.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements the \fbox macro.  It will be loaded
 *  automatically when needed, or can be loaded by
 *
 *    jsMath.Extension.Require('fbox');
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

jsMath.Add(jsMath.HTML,{

  /*
   *  Create a colored frame
   */
  Frame: function (x,y,w,h,c,pos) {

    h -= 2/jsMath.em; // use 2 pixels to compensate for border size
    w -= 2/jsMath.em;
    y -= 1/jsMath.em;
    if (!c) {c = ''} else {c = ' '+c};
    if (pos) {pos = 'absolute;'} else
             {pos = 'relative; margin-right: '+this.Em(-(w+2/jsMath.em))+'; '}
    return '<img src="'+jsMath.blank+'" style="position:' + pos
             + 'vertical-align: '+this.Em(y)+'; left: '+this.Em(x)+'; '
             + 'width:' +this.Em(w*jsMath.Browser.imgScale)+'; '
             + 'height:'+this.Em(h*jsMath.Browser.imgScale)+'; '
             + 'border: 1px solid'+c+';" />';
  }

});

jsMath.Package(jsMath.Parser,{

  macros: {fbox: 'FBox'},

  /*
   *  Implement \fbox{...}
   */
  FBox: function (name) {
    var text = this.GetArgument(this.cmd+name); if (this.error) return;
    var arg = jsMath.Box.InternalMath(text,this.mlist.data.size);
    var f = 0.25 * jsMath.sizes[this.mlist.data.size]/100;
    var box = jsMath.Box.Set(arg,this.mlist.data.style,this.mlist.data.size,1).Remeasured();
    var frame = jsMath.HTML.Frame(-f,-box.d-f,box.w+2*f,box.h+box.d+2*f);
    box.html = frame + box.html + jsMath.HTML.Spacer(f);
    box.h += f; box.d += f; box.w +=2*f; box.x += f;
    box.bh = Math.max(box.bh,box.h); box.bd = Math.max(box.bd,box.d);
    this.mlist.Add(jsMath.mItem.Atom('ord',box));
  }

});
