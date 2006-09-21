/*
 *  extensions/HTML.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file implements a number of HTML-specific extensions to TeX,
 *  including \color, \style, \class, \unicode, etc.  It will be loaded
 *  automatically when needed, or can be loaded by
 *
 *    jsMath.Extension.Require('HTML');
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

jsMath.Package(jsMath.Parser,{

  macros: {
    color:      'Color',
    href:       'Href',
    'class':    'Class',
    style:      'Style',
    cssId:      'CSSId',
    unicode:    'Unicode'
  },

  /*
   *  Show the argument in a particular color
   */
  Color: function (name) {
    var color = this.GetArgument(this.cmd+name); if (this.error) return;
    // check that it looks like a color?
    this.AddHTML(name,['<span style="color: '+color+'">','</span>']);
  },

  /*
   *  Make the argument be a link
   */
  Href: function (name) {
    var href = this.GetArgument(this.cmd+name); if (this.error) return;
    this.AddHTML(name,['<a class="link" href="'+href+'">','</a>']);
  },

  /*
   *  Apply a CSS class to the argument
   */
  Class: function (name) {
    var clss = this.GetArgument(this.cmd+name); if (this.error) return;
    this.AddHTML(name,['<span class="'+clss+'">','</span>']);
  },

  /*
   *  Apply a CSS style to the argument
   */
  Style: function (name) {
    var style = this.GetArgument(this.cmd+name); if (this.error) return;
    this.AddHTML(name,['<span style="'+style+'">','</span>']);
  },

  /*
   *  Add a CSS element ID to the argument
   */
  CSSId: function (name) {
    var id = this.GetArgument(this.cmd+name); if (this.error) return;
    this.AddHTML(name,['<span id="'+id+'">','</span>']);
  },

  /*
   *  Insert some raw HTML around the argument (this will not affect
   *  the spacing or other TeX features)
   */
  AddHTML: function (name,params) {
    var data = this.mlist.data;
    var arg = this.GetArgument(this.cmd+name); if (this.error) return;
    arg = jsMath.Parse(arg,data.font,data.size,data.style);
      if (arg.error) {this.Error(arg); return}
    this.mlist.Add(jsMath.mItem.HTML(params[0]));
    for (var i = 0; i < arg.mlist.Length(); i++) {this.mlist.Add(arg.mlist.Get(i))}
    this.mlist.Add(jsMath.mItem.HTML(params[1]));
  },

  /*
   *  Insert a unicode reference as an Ord atom.  Its argument should
   *  be the unicode code point, e.g. \unicode{8211}, or \unicode{x203F}.
   *  You can also specify the height (offset from the x height) and depth
   *  in ems, together with a CSS class for the character, e.g.,
   *  \unicode{8211,class,.2,-.3}
   */
  Unicode: function (name) {
    var arg = this.GetArgument(this.cmd+name); if (this.error) return;
    arg = arg.split(','); arg[0] = '&#'+arg[0]+';';
    if (!arg[1]) {arg[1] = 'normal'}
    this.mlist.Add(jsMath.mItem.TextAtom('ord',arg[0],arg[1],arg[2],arg[3]));
  }

});
