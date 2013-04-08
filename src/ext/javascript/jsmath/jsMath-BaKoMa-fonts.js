/*
 *  jsMath-BaKoMa-fonts.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes changes needed to use the BaKoMa fonts and
 *  their encoding.
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
 *  The BaKoMa fonts have a different encoding, so change the characters
 *  to correspond to the their encoding.
 */

if (jsMath.browser == "Mozilla" && navigator.platform != "MacPPC") {

  /*
   *  Mozilla/PC
   */
  jsMath.Update.TeXfontCodes({
    cmr10: [
      '&#x0393;', '&#x0394;', '&#x0398;', '&#x039B;',
      '&#x039E;', '&#x03A0;', '&#x03A3;', '&#x03A5;',
      '&#x03A6;', '&#x03A8;', '&#x03A9;', '&#xFB00;',
      '&#xFB01;', '&#xFB02;', '&#xFB03;', '&#xFB04;',

      '&#x0131;', '&#xEEF0;', '&#x0300;', '&#x0301;',
      '&#x030C;', '&#x0306;', '&#x0305;', '&#x030A;',
      '&#x0327;', '&#x00DF;', '&#x00E6;', '&#x0153;',
      '&#x00F8;', '&#x00C6;', '&#x0152;', '&#x00D8;',

      '&#x337;', '&#x21;', '&#x201D;', '&#x23;',
      '&#x24;', '&#x25;', '&#x26;', '&#x27;',
      '&#x28;', '&#x29;', '&#x2A;', '&#x2B;',
      '&#x2C;', '&#x2D;', '&#x2E;', '&#x2F;',

      '0', '1', '2', '3', '4', '5', '6', '7',
      '8', '9', ':', ';', '&#xA1;', '&#x3D;', '&#xBf;', '&#x3F;',

      '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
      'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',

      'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
      'X', 'Y', 'Z', '&#x5B;', '&#x201C;', '&#x5D;', '&#x302;', '&#x307;',

      '&#x2018;', 'a', 'b', 'c', 'd', 'e', 'f', 'g',
      'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',

      'p', 'q', 'r', 's', 't', 'u', 'v', 'w',
      'x', 'y', 'z', '&#x2013;', '&#x2014;', '&#x30B;', '&#x303;', '&#x308;'
    ],

    cmmi10: [
      '&#x393;', '&#x394;', '&#x398;', '&#x39B;',
      '&#x39E;', '&#x3A0;', '&#x3A3;', '&#x3A5;',
      '&#x3A6;', '&#x3A8;', '&#x3A9;', '&#x3B1;',
      '&#x3B2;', '&#x3B3;', '&#x3B4;', '&#x3B5;',

      '&#x3B6;', '&#x3B7;', '&#x3B8;', '&#x3B9;',
      '&#x3BA;', '&#x3BB;', '&#x3BC;', '&#x3BD;',
      '&#x3BE;', '&#x3C0;', '&#x3C1;', '&#x3C3;',
      '&#x3C4;', '&#x3C5;', '&#x3C6;', '&#x3C7;',

      '&#x3C8;', '&#x3C9;', '&#x25B;', '&#x3D1;',
      '&#x3D6;', '&#x3F1;', '&#x3C2;', '&#x3D5;',
      '&#x21BC;', '&#x21BD;', '&#x21C0;', '&#x21C1;',
      '&#xEFBA;', '&#xEFBB;', '&#x25B9;', '&#x25C3;',

      '0', '1', '2', '3', '4', '5', '6', '7',
      '8', '9', '&#x2E;', '&#x2C;', '&#x3C;', '&#x2F;', '&#x3E;', '&#x22C6;',

      '&#x2202;', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
      'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',

      'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
      'X', 'Y', 'Z', '&#x266D;', '&#x266E;', '&#x266F;', '&#x2323;', '&#x2322;',

      '&#x2113;', 'a', 'b', 'c', 'd', 'e', 'f', 'g',
      'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',

      'p', 'q', 'r', 's', 't', 'u', 'v', 'w',
      'x', 'y', 'z', '&#x131;', '&#xEEF0;', '&#x2118;', '&#x20D7;', '&#x0311;'
    ],

    cmsy10: [
      '&#x2212;', '&#xB7;',   '&#xD7;',   '&#x22C6;',
      '&#xF7;',   '&#x22C4;', '&#xB1;',   '&#x2213;',
      '&#x2295;', '&#x2296;', '&#x2297;', '&#x2298;',
      '&#x2299;', '&#x25CB;', '&#x2218;', '&#x2219;',

      '&#x2243;', '&#x224D;', '&#x2286;', '&#x2287;',
      '&#x2264;', '&#x2265;', '&#x227C;', '&#x227D;',
      '&#x223C;', '&#x2245;', '&#x2282;', '&#x2283;',
      '&#x226A;', '&#x226B;', '&#x227A;', '&#x227B;',

      '&#x2190;', '&#x2192;', '&#x2191;', '&#x2193;',
      '&#x2194;', '&#x2197;', '&#x2198;', '&#x2242;',
      '&#x21D0;', '&#x21D2;', '&#x21D1;', '&#x21D3;',
      '&#x21D4;', '&#x2196;', '&#x2199;', '&#x221D;',

      '&#x2032;', '&#x221E;', '&#x2208;', '&#x220B;',
      '&#x25B3;', '&#x25BD;', '&#x338;',  '&#xEFB9;',
      '&#x2200;', '&#x2203;', '&#xAC;',   '&#x2205;',
      '&#x211C;', '&#x2111;', '&#x22A4;', '&#x22A5;',

      '&#x2135;', '&#xEF35;', '&#x212C;', '&#xEF36;',
      '&#xEF37;', '&#x2130;', '&#x2131;', '&#xEF38;',
      '&#x210B;', '&#x2110;', '&#xEF39;', '&#xEF3A;',
      '&#x2112;', '&#x2133;', '&#xEF3B;', '&#xEF3C;',

      '&#xEF3D;', '&#xEF3E;', '&#x211B;', '&#xEF3F;',
      '&#xEF40;', '&#xEF41;', '&#xEF42;', '&#xEF43;',
      '&#xEF44;', '&#xEF45;', '&#xEF46;', '&#x222A;',
      '&#x2229;', '&#x228E;', '&#x2227;', '&#x2228;',

      '&#x22A2;', '&#x22A3;', '&#x230A;', '&#x230B;',
      '&#x2308;', '&#x2309;', '&#x7B;',   '&#x7D;',
      '&#x2329;', '&#x232A;', '&#x2223;', '&#x2225;',
      '&#x2195;', '&#x21D5;', '&#x2216;', '&#x2240;',

      '&#x221A;', '&#x2210;', '&#x2207;', '&#x222B;',
      '&#x2294;', '&#x2293;', '&#x2291;', '&#x2292;',
      '&#xA7;',   '&#x2020;', '&#x2021;', '&#xB6;',
      '&#x2663;', '&#x2662;', '&#x2661;', '&#x2660;'
    ],

    cmex10: [
      '&#xF006;', '&#xF007;', '&#xF008;', '&#xF009;',
      '&#xF00A;', '&#xF00B;', '&#xF00C;', '&#xF00D;',
      '&#xF00E;', '&#xF00F;', '&#xF02C;', '&#xF02D;',
      '&#xF012;', '&#xF013;', '&#xF014;', '&#xF015;',

      '&#xF016;', '&#xF017;', '&#xF018;', '&#xF019;',
      '&#xF01A;', '&#xF01B;', '&#xF01C;', '&#xF01D;',
      '&#xF01E;', '&#xF01F;', '&#xF020;', '&#xF021;',
      '&#xF02E;', '&#xF02F;', '&#xF024;', '&#xF025;',

      '&#xF026;', '&#xEFBC;', '&#xEFBD;', '&#xEFBE;',
      '&#xEFBF;', '&#xEFC0;', '&#xEFC1;', '&#xEFC2;',
      '&#xEFC3;', '&#xEFC4;', '&#xF028;', '&#xF029;',
      '&#xEFC7;', '&#xEFC8;', '&#xEFC9;', '&#xEFCA;',

      '&#xF8EB;', '&#xF8F6;', '&#xF8EE;', '&#xF8F9;',
      '&#xF8F0;', '&#xF8FB;', '&#xF8EF;', '&#xF8FA;',
      '&#xF8F1;', '&#xF8FC;', '&#xF8F3;', '&#xF8FE;',
      '&#xF8F2;', '&#xF8FD;', '&#xF8F4;', '&#xF8E6;',

      '&#xF8ED;', '&#xF8F8;', '&#xF8EC;', '&#xF8F7;',
      '&#xF02A;', '&#xF02B;', '&#xEFCD;', '&#xEFCE;',
      '&#xEFCF;', '&#xEFD0;', '&#xEFD1;', '&#xEFD2;',
      '&#xEFD3;', '&#xEFD4;', '&#xEFD5;', '&#xEFD6;',

      '&#xEFD7;', '&#xEFD8;', '&#xEFD9;', '&#xEFDA;',
      '&#xEFDB;', '&#xEFDC;', '&#xEFDD;', '&#xEFDE;',
      '&#xEFDF;', '&#xEFE0;', '&#xEFE1;', '&#xEFE2;',
      '&#xEFE3;', '&#xEFE4;', '&#xEFE5;', '&#xEFE6;',

      '&#xEFE7;', '&#xEFE8;', '&#xEFE9;', '&#xEFEA;',
      '&#xEFEB;', '&#xEFEC;', '&#xEFED;', '&#xEFEE;',
      '&#xEFEF;', '&#xEFF0;', '&#xEFF1;', '&#xEFF2;',
      '&#xEFF3;', '&#xEFF4;', '&#xEFF5;', '&#xEFF6;',

      '&#xEFF7;', '&#xEFF8;', '&#xEFF9;', '&#xEFFA;',
      '&#xEFFB;', '&#xEFFC;', '&#xEFFD;', '&#xEFFE;',
      '&#xEFFF;', '&#xF000;', '&#xF001;', '&#xF002;',
      '&#xF003;', '&#xF004;', '&#xF005;', '&#xF027;'
    ]
  });

  /*
   *  Adjust a few other characters as well
   */

  jsMath.Update.TeXfonts({
    cmr10:  {'20': {c: '&#x02C7;', tclass: 'normal', w: .3}},
    cmmi10: {
      '20': {c: '<i>&kappa</i>', tclass: 'normal'},
      '58': {c: '.', tclass: 'normal'},
      '59': {c: ',', tclass: 'normal'},
      '61': {c: '&#x2F;', tclass: 'cmr10'}
    },
    cmsy10: {
      '3':  {c: '*', tclass: 'normal'},
      '16': {c: '&#x224D;'},
      '17': {c: '&equiv;', tclass: 'normal'},
      '25': {c: '&#x2248;', tclass: 'normal'},
      '39': {c: '&#x2243;'},
      '20': {c: '&le;', tclass: 'normal'}
    },
    cmex10: {'20': {c: '<span style="font-size: 80%">&#xEFBD;</span>'}},
    cmti10: {'10': {c: '<i>&Omega;</i>', tclass: 'normal'}},
    cmbx10: {'10': {c: '<b>&Omega;</b>', tclass: 'normal'}}
  });

} else {

  jsMath.Font.BaKoMa = [
    '&#xA1;', '&#xA2;', '&#xA3;', '&#xA4;', '&#xA5;', '&#xA6;', '&#xA7;', '&#xA8;',
    '&#xA9;', '&#xAA;', '&#xAD;', '&#xAE;', '&#xAF;', '&#xB0;', '&#xB1;', '&#xB2;',

    '&#xB3;', '&#xB4;', '&#xB5;', '&#xB6;', '&#x2219;', '&#xB8;', '&#xB9;', '&#xBA;',
    '&#xBB;', '&#xBC;', '&#xBD;', '&#xBE;', '&#xBF;', '&#xC0;', '&#xC1;', '&#xC2;',

    '&#xC3;', '!', '"', '#', '$', '%', '&#x26;', '\'',
    '(', ')', '*', '+', ',', '-', '.', '/',

    '0', '1', '2', '3', '4', '5', '6', '7',
    '8', '9', ':', ';', '&#x3C;', '=', '&#x3E;', '?',

    '@', 'A', 'B', 'C', 'D', 'E', 'F', 'G',
    'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O',

    'P', 'Q', 'R', 'S', 'T', 'U', 'V', 'W',
    'X', 'Y', 'Z', '[', '\\', ']', '^', '_',

    '&#x60;', 'a', 'b', 'c', 'd', 'e', 'f', 'g',
    'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',

    'p', 'q', 'r', 's', 't', 'u', 'v', 'w',
    'x', 'y', 'z', '&#x7B;', '&#x7C;', '&#x7D;', '&#x7E;', '&#xC4;'
  ];

  jsMath.Update.TeXfontCodes({
    cmr10:  jsMath.Font.BaKoMa,
    cmmi10: jsMath.Font.BaKoMa,
    cmsy10: jsMath.Font.BaKoMa,
    cmex10: jsMath.Font.BaKoMa,
    cmti10: jsMath.Font.BaKoMa,
    cmbx10: jsMath.Font.BaKoMa
  });

  /*
   *  MSIE corrections
   */
  switch (jsMath.browser) {

    case "MSIE":
      if (navigator.platform == "Win32") {
        /*
         *  MSIE/PC
         */
        jsMath.Browser.msieFontBug = 1;
        jsMath.Update.TeXfonts({
          cmr10:  {'10': {c: '&Omega;', tclass: 'normal'}},
          cmmi10: {
            '10':  {c: '<i>&Omega;</i>', tclass: 'normal'},
            '126': {c: '&#x7E;<span style="margin-left:.1em"></span>'}
          },
          cmsy10: {
            '10': {c: '&#x2297;', tclass: 'arial'},
            '55': {c: '<span style="margin-right:-.54em">7</span>'}
          },
          cmex10: {'10': {c: '<span style="font-size: 67%">D</span>'}},
          cmti10: {'10': {c: '<i>&Omega;</i>', tclass: 'normal'}},
          cmbx10: {'10': {c: '<b>&Omega;</b>', tclass: 'normal'}}
        });
      } else {
        /*
         *  MSIE/Mac
         */
        jsMath.Update.TeXfonts({

          cmr10: {
            '3':  {c: '<font face="Symbol">L</font>', tclass: 'normal'},
            '5':  {c: '<font face="Symbol">P</font>', tclass: 'normal'},
            '10': {c: '<font face="Symbol">W</font>', tclass: 'normal'},
            '15': {c: 'ffl', tclass: 'normal'},
            '16': {c: '&#x0131;', tclass: 'normal'},
            '20': {c: '&#x02C7;', tclass: 'normal'},
            '22': {c: '&#xAF;', tclass: 'normal', w: .3},
            '25': {c: '&#xDF;', tclass: 'normal'},
            '26': {c: '&#xE6;', tclass: 'normal'},
            '27': {c: '&#x153;', tclass: 'normal'}
          },

          cmmi10: {
            '3':  {c: '<font face="Symbol">L</font>', tclass: 'italic'},
            '5':  {c: '<font face="Symbol">P</font>', tclass: 'italic'},
            '10': {c: '<font face="Symbol">W</font>', tclass: 'italic'},
            '15': {c: '<font face="Symbol">e</font>', tclass: 'italic'},
            '16': {c: '<font face="Symbol">z</font>', tclass: 'italic'},
            '20': {c: '<font face="Symbol">k</font>', tclass: 'italic'},
            '22': {c: '<font face="Symbol">m</font>', tclass: 'italic'},
            '25': {c: '<font face="Symbol">p</font>', tclass: 'italic'},
            '26': {c: '<font face="Symbol">r</font>', tclass: 'italic'},
            '27': {c: '<font face="Symbol">s</font>', tclass: 'italic'}
          },

          cmsy10: {
            '3':  {c: '<span style="vertical-align:-.3em">*</span>', tclass: 'normal'},
            '5':  {c: '&#x389;', tclass: 'normal'},
            '10': {c: '&otimes;', tclass: 'normal'},
            '15': {c: '&#x2022;', tclass: 'normal'},
            '16': {c: '&#x224D;', tclass: 'normal'},
            '20': {c: '&le;', tclass: 'normal'},
            '22': {c: '&le;', tclass: 'normal'},
            '25': {c: '&#x2248;', tclass: 'normal'},
            '26': {c: '<font face="Symbol">&#xCC;</font>', tclass: 'normal'},
            '27': {c: '<font face="Symbol">&#xC9;</font>', tclass: 'normal'}
          },

          cmex10: {
            '3':  {c: '<span style="font-size: 67%">&#x69;</span>'},
            '5':  {c: '<span style="font-size: 67%">&#x6B;</span>'},
            '10': {c: '<span style="font-size: 67%">&#x44;</span>'},
            '15': {c: '<span style="font-size: 55%">&#xC2;</span>'},
            '16': {c: '<span style="font-size: 83%">&#xB5;</span>'},
            '20': {c: '<span style="font-size: 83%">"</span>'},
            '22': {c: '<span style="font-size: 83%">$</span>'},
            '25': {c: '<span style="font-size: 83%">\'</span>'},
            '26': {c: '<span style="font-size: 83%">(</span>'},
            '27': {c: '<span style="font-size: 83%">)</span>'}
          },

          cmti10: {
            '3':  {c: '<font face="Symbol">L</font>', tclass: 'italic'},
            '5':  {c: '<font face="Symbol">P</font>', tclass: 'italic'},
            '10': {c: '<font face="Symbol">W</font>', tclass: 'italic'},
            '16': {c: '&#x0131;', tclass: 'italic'},
            '20': {c: '&#xAD;', tclass: 'italic'},
            '22': {c: '&#xAF;', tclass: 'italic', w: .3},
            '25': {c: '&#xDF;', tclass: 'italic'},
            '26': {c: '&#xE6;', tclass: 'italic'},
            '27': {c: '&#x153;', tclass: 'italic'}
          },

          cmbx10: {
            '3':  {c: '<font face="Symbol">L</font>', tclass: 'bold'},
            '5':  {c: '<font face="Symbol">P</font>', tclass: 'bold'},
            '10': {c: '<font face="Symbol">W</font>', tclass: 'bold'},
            '16': {c: '&#x0131;', tclass: 'bold'},
            '20': {c: '&#xAD;', tclass: 'bold'},
            '22': {c: '&#xAF;', tclass: 'bold', w: .3},
            '25': {c: '&#xDF;', tclass: 'bold'},
            '26': {c: '&#xE6;', tclass: 'bold'},
            '27': {c: '&#x153;', tclass: 'bold'}
          }
        });
      }
      break;

    case "Mozilla":
      if (navigator.platform == "MacPPC") {
        /*
         *  Mozilla/Mac
         */
        jsMath.Update.TeXfonts({
          cmr10:  {'10': {c: '&Omega;', tclass: 'normal'}},
          cmmi10: {'10': {c: '<i>&Omega;</i>', tclass: 'normal'}},
          cmsy10: {'10': {c: '&otimes;', tclass: 'normal'}},
          cmex10: {'10': {c: '<span style="font-size: 67%">D</span>'}},
          cmti10: {'10': {c: '<i>&Omega;</i>', tclass: 'normal'}},
          cmbx10: {'10': {c: '<b>&Omega;</b>', tclass: 'normal'}}
        });
      }
      break;

    case "Opera":
      jsMath.Update.TeXfonts({
        cmr10:  {
          '10': {c: '&Omega;', tclass: 'normal'},
          '20': {c: '&#x2C7;', tclass: 'normal'}
        },
        cmmi10: {
          '10': {c: '<i>&Omega;</i>', tclass: 'normal'},
          '20': {c: '&kappa;', tclass: 'normal'}
        },
        cmsy10: {
          '10': {c: '&otimes;', tclass: 'normal'},
          '20': {c: '&#x2264;', tclass: 'normal'}
        },
        cmex10: {
          '10': {c: '<span style="font-size: 67%">D</span>'},
          '20': {c: '<span style="font-size: 82%">"</span>'}
        },
        cmti10: {
          '10': {c: '<i>&Omega;</i>', tclass: 'normal'},
          '20': {c: '<i>&#x2C7;</i>', tclass: 'normal'}
        },
        cmbx10: {
          '10': {c: '<b>&Omega;</b>', tclass: 'normal'},
          '20': {c: '<b>&#x2C7;</b>', tclass: 'normal'}
        }
      });
      break;

    case "Konqueror":
      jsMath.Update.TeXfonts({
        cmr10:  {'20': {c: '&#x2C7;', tclass: 'normal'}},
        cmmi10: {'20': {c: '&kappa;', tclass: 'normal'}},
        cmsy10: {'20': {c: '&#x2264;', tclass: 'normal'}},
        cmex10: {'20': {c: '<span style="font-size: 84%">"</span>'}},
        cmti10: {'20': {c: '<i>&#x2C7;</i>', tclass: 'normal'}},
        cmbx10: {'20': {c: '<b>&#x2C7;</b>', tclass: 'normal'}}
      });
      break;
  }

}

jsMath.Setup.Styles({
  '.typeset .cmr10':  'font-family: cmr10, serif',
  '.typeset .cmbx10': 'font-family: cmbx10, cmr10',
  '.typeset .cmti10': 'font-family: cmti10, cmr10',
  '.typeset .cmmi10': 'font-family: cmmi10',
  '.typeset .cmsy10': 'font-family: cmsy10',
  '.typeset .cmex10': 'font-family: cmex10',
  '.typeset .arial':  "font-family: 'Arial unicode MS'"
});
