/*
 *  extensions/AMSmath.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file defines most of the macros and environments from
 *  the amsmath LaTeX package.  You can activate it by calling
 *
 *    jsMath.Extension.Require('AMSmath');
 *
 *  once jsMath.js has been loaded, or by adding "extensions/AMSmath.js"
 *  to the loadFiles array in jsMath/easy/load.js.
 *
 *  You may wish to load AMSsymbols.js as well, but note that it
 *  requires the extra msam10 and msb10 fonts that you will have
 *  to install on your server first.
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2007 by Davide P. Cervone
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

jsMath.Extension.Require("moreArrows");

jsMath.Package(jsMath.Parser,{
  macros: {
    intI:       ['Macro','\\mathchoice{\\!}{}{}{}\\!\\!\\int'],
    iint:       ['Macro','\\!\\!\\!\\mathop{\\,\\,\\,\\int\\intI}'],
    iiint:      ['Macro','\\!\\!\\!\\mathop{\\,\\,\\,\\int\\intI\\intI}'],
    iiiint:     ['Macro','\\!\\!\\!\\mathop{\\,\\,\\,\\int\\intI\\intI\\intI}'],
    idotsint:   ['Macro','\\!\\!\\mathop{\\,\\,\\int\\cdots\\int}'],

    dddot:      ['Macro','\\mathop{#1}\\limits^{\\textstyle ...}',1],
    ddddot:     ['Macro','\\mathop{#1}\\limits^{\\textstyle ....}',1],

    sideset:    ['Macro','\\mathop{\\rlap{\\phantom{#3}}}#1\\!{#3}#2',3],
    stackrel:   ['Macro','\\mathrel{\\mathop{#2}\\limits^{#1}}',2],

    boxed:      ['Macro','\\fbox{$\\displaystyle{#1}$}',1],

    tag:        'HandleTag',
    notag:      ['Macro',''],

    substack:   ['Macro','\\begin{subarray}{c}#1\\end{subarray}',1],

    varliminf:  ['Macro','\\mathop{\\underline{\\raise1.5pt{\\rule{0pt}{.6em}{0pt}\\smash{\\lower1.5pt{\\rm lim}}}}}'],
    varlimsup:  ['Macro','\\mathop{\\overline{\\rule{0pt}{.6em}{0pt}\\smash{\\rm lim}}}'],
    varinjlim:  ['Macro','\\mathop{\\underrightarrow{\\rm lim}}'],
    varprojlim: ['Macro','\\mathop{\\underleftarrow{\\rm lim}}'],

    DeclareMathOperator: 'HandleDeclareOp',

    genfrac:    'Genfrac',
    frac:       ['Genfrac',"","","",""],
    tfrac:      ['Genfrac',"","","","1"],
    dfrac:      ['Genfrac',"","","","0"],
    binom:      ['Genfrac',"(",")","0pt",""],
    tbinom:     ['Genfrac',"(",")","0pt","1"],
    dbinom:     ['Genfrac',"(",")","0pt","0"],

    cfrac:      'CFrac',

    shoveleft:  ['HandleShove','left'],
    shoveright: ['HandleShove','right']
  },

  environments: {
    align:         ['Array',null,null,'rlrlrlrlrlrl',[5/18,2,5/18,2,5/18,2,5/18,2,5/18,2,5/18],1,'D'],
    'align*':      ['Array',null,null,'rlrlrlrlrlrl',[5/18,2,5/18,2,5/18,2,5/18,2,5/18,2,5/18],1,'D'],
    multline:     'Multline',
    'multline*':  'Multline',
    split:         ['Array',null,null,'rl',[5/18],1,'D'],
    gather:        ['Array',null,null,'c',null,1,'D'],
    'gather*':     ['Array',null,null,'c',null,1,'D'],
    subarray:      ['Array',null,null,null,[0,0,0,0],1,'S',0,.25],
    smallmatrix:   ['Array',null,null,'cccccccccc',[1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3,1/3],1,'S',0]
  },

  delimiter: {
    '\\lvert':     [4,2,0x6A,3,0x0C],
    '\\rvert':     [5,2,0x6A,3,0x0C],
    '\\lVert':     [4,2,0x6B,3,0x0D],
    '\\rVert':     [5,2,0x6B,3,0x0D]
  },

  /*
   *  Ignore the tag for now
   */
  HandleTag: function (name) {
    var arg = this.trimSpaces(this.GetArgument(this.cmd+name)); if (this.error) return;
    if (arg == "*") this.GetArgument(this.cmd+name);
  },

  /*
   *  Handle \DeclareMathOperator
   */
  HandleDeclareOp: function (name) {
    var limits = "";
    var cs = this.trimSpaces(this.GetArgument(this.cmd+name)); if (this.error) return;
    if (cs == "*") {
      limits = "\\limits";
      cs = this.trimSpaces(this.GetArgument(this.cmd+name)); if (this.error) return;
    }
    if (cs.charAt(0) == "\\") {cs = cs.substr(1)}
    var op = this.GetArgument(this.cmd+name); if (this.error) return;
    op = op.replace(/\*/g,'\\char{cmr10}{0x2A}').replace(/-/g,'\\char{cmr10}{0x2D}');
    jsMath.Parser.prototype.macros[cs] = ['Macro','\\mathop{\\rm '+op+'}'+limits];
  },

  /*
   *  Record presence of \shoveleft and \shoveright
   */
  HandleShove: function (name,data) {
    if (this.mlist.data.entry == null) {this.mlist.data.entry = {}}
    this.mlist.data.entry.shove = data[0];
  },

  /*
   *  Handle \cfrac
   */
  CFrac: function (name) {
    var lr = this.GetBrackets(this.cmd+name); if (this.error) return;
    var num = this.GetArgument(this.cmd+name); if (this.error) return;
    var den = this.GetArgument(this.cmd+name); if (this.error) return;

    num = this.Process('\\strut\\textstyle{'+num+'}'); if (this.error) return;
    den = this.Process('\\strut\\textstyle{'+den+'}'); if (this.error) return;
    var data = this.mlist.data;
    var TeX = jsMath.Typeset.TeX(data.style,data.size);

    if (lr != "") {
      if (lr != 'l' && lr != 'r') {this.Error("Illegal alignment specified in "+this.cmd+name); return}
      num = jsMath.Box.Set(num,data.style,data.size);
      den = jsMath.Box.Set(den,data.style,data.size);
      if (num.w > den.w) {
        if (lr == 'l') {den.html += jsMath.HTML.Spacer(num.w-den.w)}
                  else {den.html = jsMath.HTML.Spacer(num.w-den.w) + den.html}
        den.w = num.w;
      } else if (num.w < den.w) {
        if (lr == 'l') {num.html += jsMath.HTML.Spacer(den.w-num.w)}
                  else {num.html = jsMath.HTML.Spacer(den.w-num.w) + num.html}
        num.w = den.w;
      }
    }

    this.mlist.Add(jsMath.mItem.Fraction(name,num,den,TeX.default_rule_thickness));
  },

  /*
   *  Implement AMS generalizes fraction
   */
  Genfrac: function (name,data) {
    var left = data[0]; var right = data[1];
    var thickness = data[2]; var style = data[3];

    if (left != null) {left = this.delimiter[left]} else
      {left = this.GetDelimiterArg(this.cmd+name); if (this.error) return}
    if (right != null) {right = this.delimiter[right]} else
      {right = this.GetDelimiterArg(this.cmd+name); if (this.error) return}
    if (thickness == null) {thickness = this.GetArgument(this.cmd+name); if (this.error) return}
    if (style == null) {style = this.GetArgument(this.cmd+name); if (this.error) return}

    var num = this.ProcessArg(this.cmd+name); if (this.error) return;
    var den = this.ProcessArg(this.cmd+name); if (this.error) return;

    if (left == "") {left = null}; if (right == "") {right = null}
    if (thickness == "") {
      var TeX =jsMath.Typeset.TeX(this.mlist.data.style,this.mlist.data.size);
      thickness = TeX.default_rule_thickness;
    } else {
      thickness = this.ParseDimen(thickness,this.cmd+name,0,0);
    }

    var frac = jsMath.mItem.Fraction(name,num,den,thickness,left,right);

    if (style != "") {
      style = (["D","T","S","SS"])[style];
      if (style == null) {this.Error("Bad math style for "+this.cmd+name); return}
      var mlist = new jsMath.mList([new jsMath.mItem('style',{style:style}),frac]);
      this.mlist.Add(jsMath.mItem.Atom('inner',{type:'mlist',mlist: mlist}));
    } else {
      this.mlist.Add(frac);
    }
  },

  /*
   *  Implements the multline environment
   */
  Multline: function (name,delim) {
    var data = this.mlist.data;
    var width = this.GetBrackets(this.cmd+'begin{'+name+'}'); if (this.error) return;
    var arg = this.GetEnd(name); if (this.error) return;

    var parse = new jsMath.Parser(arg+this.cmd+'\\',null,data.size,'D');
    parse.matrix = name; parse.row = []; parse.table = []; parse.rspacing = [];
    parse.Parse(); if (parse.error) {this.Error(parse); return}
    parse.HandleRow(name,1);  // be sure the last row is recorded

    //
    // check rows for extra columns and maximum width
    //
    var i; var row; var W = 0;
    for (i = 0; i < parse.table.length; i++) {
      row = parse.table[i];
      if (row.length > 1) {
        this.Error("Rows can contain only one equation in '"+name+"' environment");
        return;
      }
      if (row[0].w > W) {W = row[0].w}
    }

    //
    //  Determine width of display
    //
    if (width == "") {width = W+2} else {
      width = this.ParseDimen(width,name,0,0);
      if (width < W) {width = W}
    }

    //
    //  Shove the top and bottom lines
    //
    if (parse.table.length > 1) {
      parse.table[0][0].entry.shove = 'left';
      row = parse.table[parse.table.length-1];
      if (!row[0].entry.shove) {row[0].entry.shove = 'right'}
    }
    //
    //  Adjust widths of shoved lines
    //
    for (i = 0; i < parse.table.length; i++) {
      row = parse.table[i][0];
      if (row.entry.shove && row.w < width) {
        switch (row.entry.shove) {
          case 'left':
            row.html += jsMath.HTML.Spacer(width-row.w);
            break;

          case 'right':
            row.html = jsMath.HTML.Spacer(width-row.w)+row.html;
            break;
        }
        row.w = width;
      }
    }

    //
    //  Do the layout
    //
    var box = jsMath.Box.Layout(data.size,parse.table);
    this.mlist.Add(jsMath.mItem.Atom('ord',box));
  },

  /*
   *  Get a delimiter or empty argument
   */
  GetDelimiterArg: function (name) {
    var c = this.trimSpaces(this.GetArgument(name)); if (this.error) return null;
    if (c == "") return null;
    if (this.delimiter[c]) return this.delimiter[c];
    this.Error("Missing or unrecognized delimiter for "+name);
    return null;
  }
});