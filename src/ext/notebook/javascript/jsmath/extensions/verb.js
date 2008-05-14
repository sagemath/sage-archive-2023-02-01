/*
 *  extensions/verb.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements the \verb macro.  You can activate it
 *  by calling
 *
 *    jsMath.Extension.Macro('verb');
 *
 *  which will cause the extension to be loaded only when it is
 *  needed, or you can force it to be loaded via
 *
 *    jsMath.Extension.Require('verb');
 *
 *  once jsMath.js has been loaded, or by adding "extensions/verb.js"
 *  to the loadFiles array in the easy/load.js file.
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

jsMath.Package(jsMath.Parser,{

  macros: {verb: 'Verb'},

  /*
   *  Implement \verb|...|
   */
  Verb: function (name) {
    var c = this.GetNext(); var start = ++this.i;
    if (c == "" ) {this.Error(this.cmd+name+" requires an argument"); return}
    while (this.i < this.string.length && this.string.charAt(this.i) != c) {this.i++}
    if (this.i == this.string.length)
      {this.Error("Can't find closing delimiter for "+this.cmd+name); return}
    var text = this.string.slice(start,this.i); this.i++;
    text = text.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/>/g,'&gt;');
    text = '<span style="font-family:monospace">'+text+'</span>';
    var box = jsMath.Box.Text(text,'normal','T',this.mlist.data.size).Styled();
    box.h = box.bh+box.bd -jsMath.d; box.d = jsMath.d;
    this.mlist.Add(jsMath.mItem.Typeset(box));
  }
});
