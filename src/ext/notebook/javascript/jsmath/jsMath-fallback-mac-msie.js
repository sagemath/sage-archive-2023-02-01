/*
 *  jsMath-fallback-mac-msie.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes changes needed by Internet Explorer on the Mac
 *  for when the TeX fonts are not available.
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



/********************************************************************
 *
 *  Fix the default non-TeX-font characters to work with MSIE
 *
 */

jsMath.Update.TeXfonts({
  cmr10: {
    '0':  {c: 'G', tclass: 'greek'},
    '1':  {c: 'D', tclass: 'greek'},
    '2':  {c: 'Q', tclass: 'greek'},
    '3':  {c: 'L', tclass: 'greek'},
    '4':  {c: 'X', tclass: 'greek'},
    '5':  {c: 'P', tclass: 'greek'},
    '6':  {c: 'S', tclass: 'greek'},
    '7':  {c: '&#161;', tclass: 'greek'},
    '8':  {c: 'F', tclass: 'greek'},
    '9':  {c: 'Y', tclass: 'greek'},
    '10': {c: 'W', tclass: 'greek'},
    '22': {c: '<span style="position:relative; top:.1em">&#96;</span>', tclass: 'symbol3'}
  },

  cmti10: {
    '0':  {c: 'G', tclass: 'igreek'},
    '1':  {c: 'D', tclass: 'igreek'},
    '2':  {c: 'Q', tclass: 'igreek'},
    '3':  {c: 'L', tclass: 'igreek'},
    '4':  {c: 'X', tclass: 'igreek'},
    '5':  {c: 'P', tclass: 'igreek'},
    '6':  {c: 'S', tclass: 'igreek'},
    '7':  {c: '&#161;', tclass: 'igreek'},
    '8':  {c: 'F', tclass: 'igreek'},
    '9':  {c: 'Y', tclass: 'igreek'},
    '10': {c: 'W', tclass: 'igreek'},
    '22': {c: '<span style="position:relative; top:.1em">&#96;</span>', tclass: 'symbol3'}
  },

  cmbx10: {
    '0':  {c: 'G', tclass: 'bgreek'},
    '1':  {c: 'D', tclass: 'bgreek'},
    '2':  {c: 'Q', tclass: 'bgreek'},
    '3':  {c: 'L', tclass: 'bgreek'},
    '4':  {c: 'X', tclass: 'bgreek'},
    '5':  {c: 'P', tclass: 'bgreek'},
    '6':  {c: 'S', tclass: 'bgreek'},
    '7':  {c: '&#161;', tclass: 'bgreek'},
    '8':  {c: 'F', tclass: 'bgreek'},
    '9':  {c: 'Y', tclass: 'bgreek'},
    '10': {c: 'W', tclass: 'bgreek'},
    '22': {c: '<span style="position:relative; top:.1em">&#96;</span>', tclass: 'symbol3'}
  },
  cmmi10: {
    '0':  {c: 'G', tclass: 'igreek'},
    '1':  {c: 'D', tclass: 'igreek'},
    '2':  {c: 'Q', tclass: 'igreek'},
    '3':  {c: 'L', tclass: 'igreek'},
    '4':  {c: 'X', tclass: 'igreek'},
    '5':  {c: 'P', tclass: 'igreek'},
    '6':  {c: 'S', tclass: 'igreek'},
    '7':  {c: '&#161;', tclass: 'igreek'},
    '8':  {c: 'F', tclass: 'igreek'},
    '9':  {c: 'Y', tclass: 'igreek'},
    '10': {c: 'W', tclass: 'igreek'},
    '11': {c: 'a', tclass: 'greek'},
    '12': {c: 'b', tclass: 'greek'},
    '13': {c: 'g', tclass: 'greek'},
    '14': {c: 'd', tclass: 'greek'},
    '15': {c: 'e', tclass: 'greek'},
    '16': {c: 'z', tclass: 'greek'},
    '17': {c: 'h', tclass: 'greek'},
    '18': {c: 'q', tclass: 'greek'},
    '19': {c: 'i', tclass: 'greek'},
    '20': {c: 'k', tclass: 'greek'},
    '21': {c: 'l', tclass: 'greek'},
    '22': {c: 'm', tclass: 'greek'},
    '23': {c: 'n', tclass: 'greek'},
    '24': {c: 'x', tclass: 'greek'},
    '25': {c: 'p', tclass: 'greek'},
    '26': {c: 'r', tclass: 'greek'},
    '27': {c: 's', tclass: 'greek'},
    '28': {c: 't', tclass: 'greek'},
    '29': {c: 'u', tclass: 'greek'},
    '30': {c: 'f', tclass: 'greek'},
    '31': {c: 'c', tclass: 'greek'},
    '32': {c: 'y', tclass: 'greek'},
    '33': {c: 'w', tclass: 'greek'},
//  '41':  // leftharpoondown
//  '43':  // rightharpoondown
//  '44':  // hook left
//  '45':  // hook right
//  '92':  // natural
    '94': {c: '<span style="position:relative; top:.3em">&#xFE36;</span>'},
    '95': {c: '<span style="position:relative; top:-.2em">&#xFE35;</span>'}
//  '127': // half-circle down accent?
  },

  cmsy10: {
    '0':  {c: '&ndash;', tclass: 'normal'},
    '11': {c: '<span style="font-size: 70%">&#x25EF;</span><span style="position:relative; margin-left:-.5em; top:.1em; margin-right:.3em">/</span>', tclass: 'normal'},
    '16': {c: '<span style="position:relative;top:-.1em; font-size: 67%">&#xFE35;</span><span style="position:relative;top:.1em;font-size:67%;margin-left:-1em">&#xFE36;</span>', tclass: 'normal'},
    '48': {c: '<span style="font-size: 133%; margin-left:-.1em; margin-right: -.6em; position: relative; top:.4em">&#x2032;</span>'},
    '93': {c: '&#x222A;<span style="font-size: 50%; margin-left:-1.3em; position: relative; top:-.3em; margin-right:.6em">+</span>'},
    '96': {c: '<span style="font-size:67%; position:relative; top:-.3em;">|</span><span style="position:relative; top:-.15em; margin-left:-.1em">&ndash;</span>', tclass: 'normal'},
    '104': {c: '<span style="position:relative; top:.2em; margin-left:-.6em">&#x3008;</span>'},
    '105': {c: '<span style="position:relative; top:.2em; margin-right:-.6em">&#x3009;</span>'},
    '109': {c: '&#x21D1;<span style="position:relative; top:.1em; margin-left:-1em">&#x21D3;</span>'},
    '110': {c: '\\', d:0, tclass: 'normal'}
//  '111': // wr
//, '113': // amalg
//  '116': // sqcup
//  '117': // sqcap
//  '118': // sqsubseteq
//  '119': // sqsupseteq
  },

  cmex10: {
    '10': {c: '<span style="position:relative; top:.1em; margin-left:-.6em">&#x3008;</span>'},
    '11': {c: '<span style="position:relative; top:.1em; margin-right:-.6em">&#x3009;</span>'},
    '14': {c: '/'}, '15': {c: '\\'},
    '28': {c: '<span style="position:relative; top:.05em; margin-left:-.6em">&#x3008;</span>'},
    '29': {c: '<span style="position:relative; top:.05em; margin-right:-.6em">&#x3009;</span>'},
    '30': {c: '/'}, '31': {c: '\\'},
    '42': {c: '<span style="margin-left:-.6em">&#x3008;</span>'},
    '43': {c: '<span style="margin-right:-.6em">&#x3009;</span>'},
    '44': {c: '/'}, '45': {c: '\\'},
    '46': {c: '/'}, '47': {c: '\\'},
    '68': {c: '<span style="margin-left:-.6em">&#x3008;</span>'},
    '69': {c: '<span style="margin-right:-.6em">&#x3009;</span>'},
//  '70':  // sqcup
//  '71':  // big sqcup
    '72': {ic: 0},  '73': {ic: 0},
    '82': {tclass: 'bigop1cx', ic: .15}, '90': {tclass: 'bigop2cx', ic:.6},
    '85': {c: '&#x222A;<span style="font-size: 50%; margin-left:-1.25em; position: relative; top:-.3em; margin-right:.6em">+</span>'},
    '93': {c: '&#x222A;<span style="font-size: 50%; margin-left:-1.25em; position: relative; top:-.3em; margin-right:.6em">+</span>'},
//  '96': // coprod
//  '97': // big coprod
    '98': {c: '&#xFE3F;', h: 0.722, w: .58, tclass: 'wide1'},
    '99': {c: '&#xFE3F;', h: 0.722, w: .58, tclass: 'wide2'},
    '100': {c: '&#xFE3F;', h: 0.722, w: .58, tclass: 'wide3'},
    '101': {c: '~', h: 0.722, w: .42, tclass: 'wide1a'},
    '102': {c: '~', h: 0.8, w: .73, tclass: 'wide2a'},
    '103': {c: '~', h: 0.8, w: 1.1, tclass: 'wide3a'}
  }

});

jsMath.Setup.Styles({
  '.typeset .arrow1':   "font-family: Osaka; position: relative; top: .125em; margin: -1px",
  '.typeset .arrow2':   "font-family: Osaka; position: relative; top: .1em; margin:-1px",
  '.typeset .bigop1':   "font-family: Symbol; font-size: 110%; position:relative; top: .8em; margin:-.05em",
  '.typeset .bigop1b':  "font-family: Symbol; font-size: 140%; position: relative; top: .8em; margin:-.1em",
  '.typeset .bigop1c':  "font-family: Osaka; font-size: 125%; position:relative; top: .85em; margin:-.3em",
  '.typeset .bigop1cx': "font-family: 'Apple Chancery'; font-size: 125%; position:relative; top: .7em; margin:-.1em",
  '.typeset .bigop2':   "font-family: Symbol; font-size: 175%; position:relative; top: .8em; margin:-.07em",
  '.typeset .bigop2a':  "font-family: Baskerville; font-size: 175%; position: relative; top: .65em",
  '.typeset .bigop2b':  "font-family: Symbol; font-size: 175%; position: relative; top: .8em; margin:-.07em",
  '.typeset .bigop2c':  "font-family: Osaka; font-size: 230%; position:relative; top: .85em; margin:-.35em",
  '.typeset .bigop2cx': "font-family: 'Apple Chancery'; font-size: 250%; position:relative; top: .6em; margin-left:-.1em; margin-right:-.2em",
  '.typeset .delim1b':  "font-family: Times; font-size: 150%; position:relative; top:.8em",
  '.typeset .delim2b':  "font-family: Times; font-size: 210%; position:relative; top:.75em;",
  '.typeset .delim3b':  "font-family: Times; font-size: 300%; position:relative; top:.7em;",
  '.typeset .delim4b':  "font-family: Times; font-size: 400%; position:relative; top:.65em;",
  '.typeset .symbol3':  "font-family: Symbol",
  '.typeset .wide1':    "font-size: 50%; position: relative; top:-1.1em",
  '.typeset .wide2':    "font-size: 80%; position: relative; top:-.7em",
  '.typeset .wide3':    "font-size: 125%; position: relative; top:-.5em",
  '.typeset .wide1a':   "font-size: 75%; position: relative; top:-.5em",
  '.typeset .wide2a':   "font-size: 133%; position: relative; top: -.15em",
  '.typeset .wide3a':   "font-size: 200%; position: relative; top: -.05em",
  '.typeset .greek':    "font-family: Symbol",
  '.typeset .igreek':   "font-family: Symbol; font-style:italic",
  '.typeset .bgreek':   "font-family: Symbol; font-weight:bold"
});
