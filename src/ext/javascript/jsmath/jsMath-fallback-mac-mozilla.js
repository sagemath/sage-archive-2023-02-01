/*
 *  jsMath-fallback-mac-mozilla.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes changes needed by Mozilla-based browsers on the Mac
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
 *  Fix the default non-TeX-font characters to work with Mozilla
 *
 */

jsMath.Update.TeXfonts({
  cmmi10: {
//  '41':  // leftharpoondown
//  '43':  // rightharpoondown
    '44': {c: '<span style="position:relative; top:.15em; margin-right:-.1em; margin-left:-.2em">&#x02D3;</span>'},
    '45': {c: '<span style="position:relative; top:.15em; margin-right:-.1em; margin-left:-.2em">&#x02D2;</span>'},
    '47': {c: '<span style="font-size:60%">&#x25C1;</span>'},
//  '92':  // natural
    '126': {c: '<span style="position:relative; left: .3em; top: -.7em; font-size: 50%">&#x2192;</span>'}
  },

  cmsy10: {
    '0':  {c: '&ndash;', tclass: 'normal'},
    '11': {c: '<span style="font-size: 70%">&#x25EF;</span><span style="position:relative; margin-left:-.5em; top:.1em; margin-right:.3em">/</span>', tclass: 'normal'},
    '42': {c: '&#x2963;'}, '43': {c: '&#x2965'},
    '48': {c: '<span style="font-size: 133%; margin-right: -.75em; position: relative; top:.4em">&#x2032;</span>', tclass: 'normal'},
    '93': {c: '&#x222A;<span style="font-size: 50%; margin-left:-1.3em; position: relative; top:-.3em; margin-right:.6em">+</span>'},
    '104': {c: '<span style="position:relative; top:.15em; margin-left:-.6em">&#x3008;</span>'},
    '105': {c: '<span style="position:relative; top:.15em; margin-right:-.6em">&#x3009;</span>'},
    '109': {c: '&#x2963;<span style="position:relative; top:.1em; margin-left:-1em">&#x2965;</span>'}
//, '116':  // sqcup
//  '117':  // sqcap
//  '118':  // sqsubseteq
//  '119':  // sqsupseteq
  },

  cmex10: {
    '10': {c: '<span style="position:relative; top:.1em; margin-left:-.6em">&#x3008;</span>'},
    '11': {c: '<span style="position:relative; top:.1em; margin-right:-.6em">&#x3009;</span>'},
    '14': {c: '/'}, '15': {c: '\\'},
    '28': {c: '<span style="position:relative; top:.1em; margin-left:-.6em">&#x3008;</span>'},
    '29': {c: '<span style="position:relative; top:.1em; margin-right:-.6em">&#x3009;</span>'},
    '30': {c: '/'}, '31': {c: '\\'},
    '42': {c: '<span style="position:relative; top:.1em; margin-left:-.6em">&#x3008;</span>'},
    '43': {c: '<span style="position:relative; top:.1em; margin-right:-.6em">&#x3009;</span>'},
    '44': {c: '/'}, '45': {c: '\\'},
    '46': {c: '/'}, '47': {c: '\\'},
    '68': {c: '<span style="position:relative; top:.1em; margin-left:-.6em">&#x3008;</span>'},
    '69': {c: '<span style="position:relative; top:.1em; margin-right:-.6em">&#x3009;</span>'},
//  '70':  // sqcup
//  '71':  // big sqcup
    '72': {ic: .194},  '73': {ic: .444},
    '82': {tclass: 'bigop1cx', ic: .15}, '90': {tclass: 'bigop2cx', ic:.6},
    '85': {c: '&#x222A;<span style="font-size: 50%; margin-left:-1.3em; position: relative; top:-.3em; margin-right:.6em">+</span>'},
    '93': {c: '&#x222A;<span style="font-size: 50%; margin-left:-1.3em; position: relative; top:-.3em; margin-right:.6em">+</span>'}
  }

});

jsMath.Setup.Styles({
  '.typeset .symbol':   "font-family: Osaka",
  '.typeset .arrow1':   "font-family: Osaka; position: relative; top: .125em; margin: -1px",
  '.typeset .arrow2':   "font-family: AppleGothic; font-size: 100%; position:relative; top: .11em; margin:-1px",
  '.typeset .bigop1':   "font-family: AppleGothic; font-size: 110%; position:relative; top: .9em; margin:-.05em",
  '.typeset .bigop1b':  "font-family: Osaka; font-size: 140%; position: relative; top: .8em; margin:-.1em",
  '.typeset .bigop1c':  "font-family: AppleGothic; font-size: 125%; position:relative; top: .85em; margin:-.3em",
  '.typeset .bigop1cx': "font-family: 'Apple Chancery'; font-size: 125%; position:relative; top: .7em; margin:-.1em",
  '.typeset .bigop2':   "font-family: AppleGothic; font-size: 175%; position:relative; top: .85em; margin:-.1em",
  '.typeset .bigop2b':  "font-family: Osaka; font-size: 200%; position: relative; top: .75em; margin:-.15em",
  '.typeset .bigop2c':  "font-family: AppleGothic; font-size: 300%; position:relative; top: .75em; margin:-.35em",
  '.typeset .bigop2cx': "font-family: 'Apple Chancery'; font-size: 250%; position:relative; top: .7em; margin-left:-.1em; margin-right:-.2em",
  '.typeset .delim1b':  "font-family: Times; font-size: 150%; position:relative; top:.8em; margin:.01em",
  '.typeset .delim2b':  "font-family: Times; font-size: 210%; position:relative; top:.8em; margin:.01em",
  '.typeset .delim3b':  "font-family: Times; font-size: 300%; position:relative; top:.75em; margin:.01em",
  '.typeset .delim4b':  "font-family: Times; font-size: 400%; position:relative; top:.725em; margin:.01em"
});


/*
 *  replace \not and \joinrel with better dimensions
 */

jsMath.Macro('not','\\mathrel{\\rlap{\\kern 3mu/}}');
jsMath.Macro('joinrel','\\mathrel{\\kern-3mu}');
