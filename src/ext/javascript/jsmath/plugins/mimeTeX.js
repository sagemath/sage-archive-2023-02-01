/*
 *  mimeTeX.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes jsMath more compatible with the mimeTeX program.
 *  It does not make everything work, but it goes a long way.
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


/*
 *  Treat ~ as space
 */
jsMath.Parser.prototype.nextIsSpace = function () {
  return this.string.charAt(this.i) == ' ' ||
         this.string.charAt(this.i) == '~';
}
jsMath.Parser.prototype.special['~'] = 'Space';

/*
 *  Implement \[ ... \], \( ... \), etc.
 */
jsMath.Macro('[','\\left[');  jsMath.Macro(']','\\right]');
jsMath.Macro('(','\\left(');  jsMath.Macro(')','\\right)');
jsMath.Macro('<','\\left<');  jsMath.Macro('>','\\right>');
// can't do \. in a reasonable way
jsMath.Parser.prototype.macros['|'] = ['HandleLR','|','|'];
jsMath.Parser.prototype.macros['='] = ['HandleLR','\\|','\\|'];

/*
 *  Make non-standard \left{ and \right} work
 */
jsMath.Parser.prototype.delimiter['}'] = [5,2,0x67,3,0x09];
jsMath.Parser.prototype.delimiter['{'] = [4,2,0x66,3,0x08];


/*
 *  Immitate mimeTeX \big... and \Big... ops
 */

// make the normal ones in text mode
jsMath.Macro('int','\\intop\\nolimits');
jsMath.Macro('oint','\\ointop\\nolimits');
jsMath.Macro('sum','\\sumop\\nolimits');
jsMath.Macro('prod','\\prodop\\nolimits');
jsMath.Macro('coprod','\\coprodop\\nolimits');

jsMath.Macro('bigint','\\bigintop\\nolimits');
jsMath.Macro('bigoint','\\bigointop\\nolimits');
jsMath.Macro('bigsum','\\bigsumop\\nolimits');
jsMath.Macro('bigprod','\\bigprodop\\nolimits');
jsMath.Macro('bigcoprod','\\bigcoprodop\\nolimits');

jsMath.Macro('Bigint','\\bigintop\\limits');
jsMath.Macro('Bigoint','\\bigointop\\limits');
jsMath.Macro('Bigsum','\\bigsumop\\limits');
jsMath.Macro('Bigprod','\\bigprodop\\limits');
jsMath.Macro('Bigcoprod','\\bigcoprod\\limits');

/*
 *  The characters needed for the macros above
 */
jsMath.Parser.prototype.mathchardef['coprodop'] = [1,3,0x60];
jsMath.Parser.prototype.mathchardef['prodop']   = [1,3,0x51];
jsMath.Parser.prototype.mathchardef['sumop']    = [1,3,0x50];

jsMath.Parser.prototype.mathchardef['bigintop']    = [1,3,0x5A];
jsMath.Parser.prototype.mathchardef['bigointop']   = [1,3,0x49];
jsMath.Parser.prototype.mathchardef['bigcoprodop'] = [1,3,0x61];
jsMath.Parser.prototype.mathchardef['bigprodop']   = [1,3,0x59];
jsMath.Parser.prototype.mathchardef['bigsumop']    = [1,3,0x58];

/*
 * Unlink the small versions so they don't enlarge in display mode
 */
jsMath.TeX['cmex10'][0x48].n = null;
jsMath.TeX['cmex10'][0x50].n = null;
jsMath.TeX['cmex10'][0x51].n = null;
jsMath.TeX['cmex10'][0x52].n = null;
jsMath.TeX['cmex10'][0x60].n = null;


/*
 *  Some other missing items
 */
jsMath.Macro('/','{}'); // insert an empty box \/
jsMath.Macro('raisebox','\\raise #1px ',1); // convert to \raise
jsMath.Macro('hfill','\\quad ',1); // punt
jsMath.Macro('fbox','\\oldfbox{$#1$}',1); // do fbox in math mode

/*
 *  These get new JavaScript routines
 */
jsMath.Parser.prototype.macros['unitlength'] = 'unitlength';
jsMath.Parser.prototype.macros['hspace']     = 'hspace';
jsMath.Parser.prototype.macros['fs']         = 'fs';
jsMath.Parser.prototype.macros['oldfbox']    = 'FBox';

/*
 *  Add some JavaScript functions to the parser
 */
jsMath.Package(jsMath.Parser,{

  /*
   *  Implement \left x ... \right x
   */
  HandleLR: function (name,data) {
    var arg = this.GetUpto(name,name); if (this.error) return;
    this.string = '\\left'+data[0]+arg+'\\right'+data[1];
    this.i = 0;
  },

  /*
   *  Hold the unit length in mlist.data
   */
  unitlength: function (name) {
    var n = this.GetArgument(this.cmd+name); if (this.error) return;
    if (!n.match(/^-?(\d+(\.\d*)?|\.\d+)$/)) {
      this.Error("Argument for "+this.cmd+name+" must be a number");
      return;
    }
    this.mlist.data['unitlength'] = n;
  },

  /*
   *  Get the length (converted to ems) and multiply by the unit length
   */
  hspace: function (name) {
    var w = this.GetArgument(this.cmd+name); if (this.error) return;
    if (!w.match(/^-?(\d+(\.\d*)?|\.\d+)$/)) {
      this.Error("Argument for "+this.cmd+name+" must be a number");
      return;
    }
    w /= jsMath.em
    if (this.mlist.data['unitlength']) {w *= this.mlist.data['unitlength']}
    this.mlist.Add(jsMath.mItem.Space(w));
  },

  /*
   *  Implement \fs{...} for font-size changing
   */
  fs: function (name) {
    var n = this.GetArgument(this.cmd+name); if (this.error) return;
    if (!n.match(/^[-+]?\d+$/)) {
      this.Error("Argument for "+this.cmd+name+" must be an integer");
      return;
    }
    if (n.match(/[-+]/)) {n = n - 0; n += this.mlist.data.size}
    this.mlist.data.size = n = Math.max(0,Math.min(9,n));
    this.mlist.Add(new jsMath.mItem('size',{size: n}));
  },

  /*
   *  Repalce the Array function by one that accepts an optional
   *  parameter for the column types, and that handle's mimeTeX's
   *  "preamble" format.
   */
  Array: function (name,delim) {
    var columns = delim[2]; var cspacing = delim[3];
    if (!columns && this.GetNext() == '{') {
      columns = this.GetArgument(this.cmd+'begin{'+name+'}');
      if (this.error) return;
    } else {
      columns = '';
    }
    columns = columns.replace(/[^clr]/g,'');
    columns = columns.split('');
    var data = this.mlist.data; var style = delim[5] || 'T';
    var arg = this.GetEnd(name); if (this.error) return;
    if (arg.match(/\$/)) {arg = arg.replace(/^([^$]*)\$/,''); columns = RegExp.$1}
    var parse = new jsMath.Parser(arg+this.cmd+'\\',null,data.size,style);
    parse.matrix = name; parse.row = []; parse.table = []; parse.rspacing = [];
    parse.Parse(); if (parse.error) {this.Error(parse); return}
    parse.HandleRow(name,1);  // be sure the last row is recorded
    var box = jsMath.Box.Layout(data.size,parse.table,columns,cspacing,parse.rspacing,delim[4]||null);
    // Add parentheses, if needed
    if (delim[0] && delim[1]) {
      var left  = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[0]],'T');
      var right = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[1]],'T');
      box = jsMath.Box.SetList([left,box,right],data.style,data.size);
    }
    this.mlist.Add(jsMath.mItem.Atom((delim[0]? 'inner': 'ord'),box));
  },

  /*
   *  Similarly for Matrix (used by \matrix and \array)
   */
  Matrix: function (name,delim) {
    var data = this.mlist.data;
    var arg = this.GetArgument(this.cmd+name); if (this.error) return;
    if (arg.match(/\$/)) {arg = arg.replace(/^([^$]*)\$/,''); delim[2] = RegExp.$1}
    var parse = new jsMath.Parser(arg+this.cmd+'\\',null,data.size,delim[5] || 'T');
    parse.matrix = name; parse.row = []; parse.table = []; parse.rspacing = [];
    parse.Parse(); if (parse.error) {this.Error(parse); return}
    parse.HandleRow(name,1);  // be sure the last row is recorded
    var box = jsMath.Box.Layout(data.size,parse.table,delim[2]||null,delim[3]||null,parse.rspacing,delim[4]||null);
    // Add parentheses, if needed
    if (delim[0] && delim[1]) {
      var left  = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[0]],'T');
      var right = jsMath.Box.Delimiter(box.h+box.d-jsMath.hd/4,this.delimiter[delim[1]],'T');
      box = jsMath.Box.SetList([left,box,right],data.style,data.size);
    }
    this.mlist.Add(jsMath.mItem.Atom((delim[0]? 'inner': 'ord'),box));
  }
});
