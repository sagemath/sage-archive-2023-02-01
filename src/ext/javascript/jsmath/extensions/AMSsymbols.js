/*
 *  extensions/AMSsymbol.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file defines the macros needed to access the AMS symbol fonts
 *  available in msam10 and msbm10.  You can activate it by calling
 *
 *    jsMath.Extension.Require('AMSsymbols');
 *
 *  once jsMath.js has been loaded.
 *
 *  Note that you will need to install the msam10 and msbm10 fonts
 *  that are available from the jsMath extra font page at
 *
 *      http://www.math.union.edu/locate/jsMath/download/extra-fonts/
 *
 *  in order to make this work in image mode.  Note that there is no
 *  unicode fallback mode for these fonts at this time.
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

/********************************************************************/


delete jsMath.Parser.prototype.macros['hbar'];
delete jsMath.Parser.prototype.macros['angle'];
delete jsMath.Parser.prototype.macros['rightleftharpoons'];

jsMath.Extension.MathChar("msam10",{
  // Miscellaneous symbols
  vartriangle:        [3,0x4D],
  triangledown:       [0,0x4F],
  square:             [0,0x03],
  lozenge:            [0,0x06],
  circledS:           [0,0x73],
  angle:              [0,0x5C],
  measuredangle:      [0,0x5D],
  backprime:          [0,0x38],
  blacktriangle:      [0,0x4E],
  blacktriangledown:  [0,0x48],
  blacksquare:        [0,0x04],
  blacklozenge:       [0,0x07],
  bigstar:            [0,0x46],
  sphericalangle:     [0,0x5E],
  complement:         [0,0x7B],

  // Binary operators
  dotplus:            [2,0x75],
  Cap:                [2,0x65],
  doublecap:          [2,0x65],
  Cup:                [2,0x64],
  doublecup:          [2,0x64],
  barwedge:           [2,0x5A],
  veebar:             [2,0x59],
  doublebarwedge:     [2,0x5B],
  boxminus:           [2,0x0C],
  boxtimes:           [2,0x02],
  boxdot:             [2,0x00],
  boxplus:            [2,0x01],
  leftthreetimes:     [2,0x68],
  rightthreetimes:    [2,0x69],
  curlywedge:         [2,0x66],
  curlyvee:           [2,0x67],
  circleddash:        [2,0x7F],
  circledast:         [2,0x7E],
  circledcirc:        [2,0x7D],
  centerdot:          [2,0x05],
  intercal:           [2,0x7C],

  // Binary relations
  leqq:               [3,0x35],
  leqslant:           [3,0x36],
  eqslantless:        [3,0x30],
  lesssim:            [3,0x2E],
  lessapprox:         [3,0x2F],
  lll:                [3,0x6E],
  llless:             [3,0x6E],
  lessgtr:            [3,0x37],
  lesseqgtr:          [3,0x51],
  lesseqqgtr:         [3,0x53],
  doteqdot:           [3,0x2B],
  Doteq:              [3,0x2B],
  risingdotseq:       [3,0x3A],
  fallingdotseq:      [3,0x3B],
  backsim:            [3,0x76],
  backsimeq:          [3,0x77],
  subseteqq:          [3,0x6A],
  Subset:             [3,0x62],
  sqsubset:           [3,0x40],
  preccurlyeq:        [3,0x34],
  curlyeqprec:        [3,0x32],
  precsim:            [3,0x2D],
  vartriangleleft:    [3,0x43],
  trianglelefteq:     [3,0x45],
  vDash:              [3,0x0F],
  Vvdash:             [3,0x0E],
  smallsmile:         [3,0x60],
  smallfrown:         [3,0x61],
  bumpeq:             [3,0x6C],
  Bumpeq:             [3,0x6D],
  varpropto:          [3,0x5F],
  blacktriangleleft:  [3,0x4A],
  therefore:          [3,0x29],
  geqq:               [3,0x3D],
  geqslant:           [3,0x3E],
  eqslantgtr:         [3,0x31],
  gtrsim:             [3,0x26],
  gtrapprox:          [3,0x27],
  ggg:                [3,0x6F],
  gggtr:              [3,0x6F],
  gtrless:            [3,0x3F],
  gtreqless:          [3,0x52],
  gtreqqless:         [3,0x54],
  eqcirc:             [3,0x50],
  circeq:             [3,0x24],
  triangleq:          [3,0x2C],
  supseteqq:          [3,0x6B],
  Supset:             [3,0x63],
  sqsupset:           [3,0x41],
  succcurlyeq:        [3,0x3C],
  curlyeqsucc:        [3,0x33],
  succsim:            [3,0x25],
  vartriangleright:   [3,0x42],
  trianglerighteq:    [3,0x44],
  Vdash:              [3,0x0D],
  between:            [3,0x47],
  pitchfork:          [3,0x74],
  blacktriangleright: [3,0x49],
  because:            [3,0x2A],

  // Arrows
  leftleftarrows:     [3,0x12],
  leftrightarrows:    [3,0x1C],
  Lleftarrow:         [3,0x57],
  twoheadleftarrow:   [3,0x11],
  leftarrowtail:      [3,0x1B],
  looparrowleft:      [3,0x22],
  leftrightharpoons:  [3,0x0B],
  circlearrowleft:    [3,0x09],
  Lsh:                [3,0x1E],
  upuparrows:         [3,0x14],
  upharpoonleft:      [3,0x18],
  downharpoonleft:    [3,0x19],
  multimap:           [3,0x28],
  leftrightsquigarrow:[3,0x21],
  rightrightarrows:   [3,0x13],
  rightleftarrows:    [3,0x1D],
  Rrightarrow:        [3,0x56],
  twoheadrightarrow:  [3,0x10],
  rightarrowtail:     [3,0x1A],
  looparrowright:     [3,0x23],
  rightleftharpoons:  [3,0x0A],
  circlearrowright:   [3,0x08],
  Rsh:                [3,0x1F],
  downdownarrows:     [3,0x15],
  upharpoonright:     [3,0x16],
  downharpoonright:   [3,0x17],
  rightsquigarrow:    [3,0x20]
});

jsMath.Extension.MathChar("msbm10",{
  // Lowercase Greek letters
  digamma:            [0,0x7A],
  varkappa:           [0,0x7B],

  // Hebrew letters
  beth:               [0,0x69],
  daleth:             [0,0x6B],
  gimel:              [0,0x6A],

  // Miscellaneous symbols
  hbar:               [0,0x7E],
  hslash:             [0,0x7D],
  nexists:            [0,0x40],
  mho:                [0,0x66],
  Finv:               [0,0x60],
  Game:               [0,0x61],
  Bbbk:               [0,0x7C],
  varnothing:         [0,0x3F],
  eth:                [0,0x67],
  diagup:             [0,0x1E],
  diagdown:           [0,0x1F],

  // Binary operators
  smallsetminus:      [2,0x72],
  divideontimes:      [2,0x3E],
  ltimes:             [2,0x6E],
  rtimes:             [2,0x6F],

  // Binary relations
  approxeq:           [3,0x75],
  lessdot:            [3,0x6C],
  precapprox:         [3,0x77],
  gtrdot:             [3,0x6D],
  thicksim:           [3,0x73],
  thickapprox:        [3,0x74],
  succapprox:         [3,0x76],
  shortmid:           [3,0x70],
  shortparallel:      [3,0x71],
  backepsilon:        [3,0x7F],

  // Negated relations
  nless:              [3,0x04],
  nleq:               [3,0x02],
  nleqslant:          [3,0x0A],
  nleqq:              [3,0x14],
  lneq:               [3,0x0C],
  lneqq:              [3,0x08],
  lvertneqq:          [3,0x00],
  lnsim:              [3,0x12],
  lnapprox:           [3,0x1A],
  nprec:              [3,0x06],
  npreceq:            [3,0x0E],
  precneqq:           [3,0x16],
  precnsim:           [3,0x10],
  precnapprox:        [3,0x18],
  nsim:               [3,0x1C],
  nshortmid:          [3,0x2E],
  nmid:               [3,0x2D],
  nvdash:             [3,0x30],
  nVdash:             [3,0x31],
  ntriangleleft:      [3,0x36],
  ntrianglelefteq:    [3,0x35],
  nsubseteq:          [3,0x2A],
  nsubseteqq:         [3,0x22],
  subsetneq:          [3,0x28],
  varsubsetneq:       [3,0x20],
  subsetneqq:         [3,0x24],
  varsubsetneqq:      [3,0x26],
  ngtr:               [3,0x05],
  ngeq:               [3,0x03],
  ngeqslant:          [3,0x0B],
  ngeqq:              [3,0x15],
  gneq:               [3,0x0D],
  gneqq:              [3,0x09],
  gvertneqq:          [3,0x01],
  gnsim:              [3,0x13],
  gnapprox:           [3,0x1B],
  nsucc:              [3,0x07],
  nsucceq:            [3,0x0F],
  succneqq:           [3,0x17],
  succnsim:           [3,0x11],
  succnapprox:        [3,0x19],
  ncong:              [3,0x1D],
  nshortparallel:     [3,0x2F],
  nparallel:          [3,0x2C],
  nvDash:             [3,0x32],
  nVDash:             [3,0x33],
  ntriangleright:     [3,0x37],
  ntrianglerighteq:   [3,0x34],
  nsupseteq:          [3,0x2B],
  nsupseteqq:         [3,0x23],
  supsetneq:          [3,0x29],
  varsupsetneq:       [3,0x21],
  supsetneqq:         [3,0x25],
  varsupsetneqq:      [3,0x27],

  // Arrows
  curvearrowleft:     [3,0x78],
  curvearrowright:    [3,0x79],

  // Negated arrows
  nleftarrow:         [3,0x38],
  nLeftarrow:         [3,0x3A],
  nleftrightarrow:    [3,0x3D],
  nrightarrow:        [3,0x39],
  nRightarrow:        [3,0x3B],
  nLeftrightarrow:    [3,0x3C]
});

jsMath.Macro('Bbb','{\\msbm #1}',1);
jsMath.Extension.Font('msbm');
jsMath.Extension.Font('msam');
