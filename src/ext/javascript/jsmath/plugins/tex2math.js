/*
 *  tex2math.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file is a plugin that searches text within a web page
 *  for \(...\), \[...\], $...$ and $$...$$ and converts them to
 *  the appropriate <SPAN CLASS="math">...</SPAN> or
 *  <DIV CLASS="math">...</DIV> tags.
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

if (!jsMath.tex2math) {jsMath.tex2math = {}}  // make sure jsMath.tex2math is defined
if (!jsMath.tex2math.loaded) {                // only load it once

if (!jsMath.Controls) {jsMath.Controls = {}}
if (!jsMath.Controls.cookie) {jsMath.Controls.cookie = {}}

jsMath.Add(jsMath.tex2math,{

  loaded: 1,
  window: window,

  /*
   *  Call the main conversion routine with appropriate flags
   */

  ConvertTeX: function (element) {
    this.Convert(element,{
      processSingleDollars: 1, processDoubleDollars: 1,
      processSlashParens: 1, processSlashBrackets: 1,
      custom: 0, fixEscapedDollars: 1
    });
  },

  ConvertTeX2: function (element) {
    this.Convert(element,{
      processSingleDollars: 0, processDoubleDollars: 1,
      processSlashParens: 1, processSlashBrackets: 1,
      custom: 0, fixEscapedDollars: 0
    });
  },

  ConvertLaTeX: function (element) {
    this.Convert(element,{
      processSingleDollars: 0, processDoubleDollars: 0,
      processSlashParens: 1, processSlashBrackets: 1,
      custom: 0, fixEscapedDollars: 0
    });
  },

  ConvertCustom: function (element) {
    this.Convert(element,{custom: 1, fixEscapedDollars: 0});
  },

  /*******************************************************************/

  /*
   *  Define a custom search by indicating the
   *  strings to use for starting and ending
   *  in-line and display mathematics
   */
  CustomSearch: function (iOpen,iClose,dOpen,dClose) {
    this.inLineOpen = iOpen; this.inLineClose = iClose;
    this.displayOpen = dOpen; this.displayClose = dClose;
    this.createPattern('customPattern',new RegExp(
      '('+this.patternQuote(dOpen)+'|'
         +this.patternQuote(iOpen)+'|'
         +this.patternQuote(dClose)+'|'
         +this.patternQuote(iClose)+'|\\\\.)','g'
    ));
  },

  patternQuote: function (s) {
    s = s.replace(/([\^$(){}+*?\-|\[\]\:\\])/g,'\\$1');
    return s;
  },

  /*
   *  MSIE on the Mac doesn't handle lastIndex correctly, so
   *  override it and implement it correctly.
   */
  createPattern: function (name,pattern) {
    jsMath.tex2math[name] = pattern;
    if (this.fixPatterns) {
      pattern.oldExec = pattern.exec;
      pattern.exec = this.msiePatternExec;
    }
  },
  msiePatternExec: function (string) {
    if (this.lastIndex == null) (this.lastIndex = 0);
    var match = this.oldExec(string.substr(this.lastIndex));
    if (match) {this.lastIndex += match.lastIndex}
          else {this.lastIndex = null}
    return match;
  },

  /*******************************************************************/

  /*
   *  Set up for the correct type of search, and recursively
   *  convert the mathematics.  Disable tex2math if the cookie
   *  isn't set, or of there is an element with ID of 'tex2math_off'.
   */
  Convert: function (element,flags) {
    this.Init();
    if (!element) {element = jsMath.document.body}
    if (typeof(element) == 'string') {element = jsMath.document.getElementById(element)}
    if (jsMath.Controls.cookie.tex2math &&
        (!jsMath.tex2math.allowDisableTag || !jsMath.document.getElementById('tex2math_off'))) {
      this.custom = 0; for (var i in flags) {this[i] = flags[i]}
      if (this.custom) {
        this.pattern = this.customPattern;
        this.ProcessMatch = this.customProcessMatch;
      } else {
        this.pattern = this.stdPattern;
        this.ProcessMatch = this.stdProcessMatch;
      }
      if (this.processDoubleDollars || this.processSingleDollars ||
          this.processSlashParens   || this.processSlashBrackets ||
          this.custom) this.ScanElement(element);
    }
  },

  /*
   *  Recursively look through a document for text nodes that could
   *  contain mathematics.
   */
  ScanElement: function (element,ignore) {
    if (!element) {element = jsMath.document.body}
    if (typeof(element) == 'string') {element = jsMath.document.getElementById(element)}
    while (element) {
      if (element.nodeName == '#text') {
        if (!ignore) {element = this.ScanText(element)}
      } else {
        if (element.className == null) {element.className = ''}
        if (element.firstChild && element.className != 'math') {
          var off = ignore || element.className == 'tex2math_ignore' ||
             (element.tagName && element.tagName.match(/^(script|noscript|style|textarea|pre)$/i));
          off = off && element.className != 'tex2math_process';
          this.ScanElement(element.firstChild,off);
        }
      }
      if (element) {element = element.nextSibling}
    }
  },

  /*
   *  Looks through a text element for math delimiters and
   *  process them.  If <BR> tags are found in the middle, they
   *  are ignored (this is for BBS systems that have editors
   *  that insert these automatically).
   */
  ScanText: function (element) {
    if (element.nodeValue.replace(/\s+/,'') == '') {return element}
    var match; var prev; this.search = {};
    while (element) {
      this.pattern.lastIndex = 0;
      while (element && element.nodeName == '#text' &&
            (match = this.pattern.exec(element.nodeValue))) {
        this.pattern.match = match;
        element = this.ProcessMatch(match[0],match.index,element);
      }
      if (this.search.matched) {element = this.EncloseMath(element)}
      if (!element) {return null}
      prev = element; element = element.nextSibling;
      while (element && element.nodeName.toLowerCase() == 'br')
        {prev = element; element = element.nextSibling}
      if (!element || element.nodeName != '#text') {return prev}
    }
    return element;
  },

  /*
   *  If a matching end tag has been found, process the mathematics.
   *  Otherwise, update the search data for the given delimiter,
   *  or ignore it, as the item dictates.
   */
  stdProcessMatch: function (match,index,element) {
    if (match == this.search.end) {
      this.search.close = element;
      this.search.clength = match.length;
      this.search.cpos = this.pattern.lastIndex;
      element = this.EncloseMath(element);
    } else {
      switch (match) {
        case '\\(':
          if (this.search.end == null ||
             (this.search.end != '$' && this.search.end != '$$') &&
              this.processSlashParens) {
            this.ScanMark('span',element,'\\)');
          }
          break;

        case '\\[':
          if (this.search.end == null ||
             (this.search.end != '$' && this.search.end != '$$') &&
              this.processSlashBrackets) {
            this.ScanMark('div',element,'\\]');
          }
          break;

        case '$$':
          if (this.processDoubleDollars) {
            var type = (this.doubleDollarsAreInLine? 'span': 'div');
            this.ScanMark(type,element,'$$');
          }
          break;

        case '$':
          if (this.search.end == null && this.processSingleDollars) {
            this.ScanMark('span',element,'$');
          }
          break;

        case '\\$':
          if (this.search.end == null && this.fixEscapedDollars) {
            element.nodeValue = element.nodeValue.substr(0,index)
                              + element.nodeValue.substr(index+1);
          }
          break;
      }
    }
    return element;
  },

  /*
   *  If a matching end tag has been found, process the mathematics.
   *  Otherwise, update the search data for the given delimiter,
   *  or ignore it, as the item dictates.
   */
  customProcessMatch: function (match,index,element) {
    if (match == this.search.end) {
      this.search.close = element;
      this.search.clength = match.length;
      this.search.cpos = this.pattern.lastIndex;
      this.search.matched = 1;
    } else if (match == this.inLineOpen) {
      if (this.search.matched) {element = this.EncloseMath(element)}
      this.ScanMark('span',element,this.inLineClose);
    } else if (match == this.displayOpen) {
      if (this.search.matched) {element = this.EncloseMath(element)}
      this.ScanMark('div',element,this.displayClose);
    }
    return element;
  },

  /*
   *  Return a structure that records the starting location
   *  for the math element, and the end delimiter we want to find.
   */
  ScanMark: function (type,element,end) {
    var len = this.pattern.match[1].length;
    this.search = {
      type: type, end: end, open: element, olength: len,
      pos: this.pattern.lastIndex - len
    };
  },

  /*******************************************************************/

  /*
   *  Surround the mathematics by an appropriate
   *  SPAN or DIV element marked as CLASS="math".
   */
  EncloseMath: function (element) {
    if (this.callback) {if (!this.callback()) {return null}}
    var search = this.search;
    var close = search.close;
    if (search.cpos == close.length) {close = close.nextSibling}
       else {close = close.splitText(search.cpos)}
    if (!close) {close = jsMath.document.createTextNode("")}
    if (element == search.close) {element = close}
    var math = search.open.splitText(search.pos);
    while (math.nextSibling && math.nextSibling != close) {
      if (math.nextSibling.nodeValue) {math.nodeValue += math.nextSibling.nodeValue}
        else {math.nodeValue += ' '}
      math.parentNode.removeChild(math.nextSibling);
    }
    var TeX = math.nodeValue.substr(search.olength,
      math.nodeValue.length-search.olength-search.clength);
    math.parentNode.removeChild(math);
    math = this.createMathTag(search.type,TeX);
    //
    //  This is where older, buggy browsers can fail under unpredicatble
    //  circumstances, so we trap errors and at least get to continue
    //  with the rest of the math.  (## should add error message ##)
    //
    try {
      if (close && close.parentNode) {
        close.parentNode.insertBefore(math,close);
      } else if (search.open.nextSibling) {
        search.open.parentNode.insertBefore(math,search.open.nextSibling);
      } else {
        search.open.parentNode.appendChild(math);
      }
    } catch (err) {}
    this.search = {}; this.pattern.lastIndex = 0;
    return math;
  },

  /*
   *  Create an element for the mathematics
   */
  createMathTag: function (type,text) {
    var tag = jsMath.document.createElement(type); tag.className = "math";
    var math = jsMath.document.createTextNode(text);
    tag.appendChild(math);
    return tag;
  },

  //
  //  MSIE won't let you insert a DIV within tags that are supposed to
  //  contain in-line data (like <P> or <SPAN>), so we have to fake it
  //  using SPAN tags that force the formatting to work like DIV.  We
  //  use a separate SPAN that is the full width of the containing
  //  item, and that has the margins and centering from the div.typeset
  //  style.
  //
  MSIEcreateMathTag: function (type,text) {
    var tag = jsMath.document.createElement("span");
    tag.className = "math";
    text = text.replace(/</g,'&lt;').replace(/>/g,'&gt;');
    if (type == 'div') {
      tag.className = "";
      tag.style.width = "100%"; tag.style.margin = jsMath.tex2math.margin;
      tag.style.display = "inline-block";
      text = '<span class="math">\\displaystyle{'+text+'}</span>';
      if (jsMath.tex2math.center) {
        tag.style.textAlign = "center";
        text = '<span style="text-align:left">'+text+'</span>'
      }
    }
    tag.innerHTML = text;
    return tag;
  },

  /*******************************************************************/

  Init: function () {

    if (this.inited || !jsMath.browser) return
    /*
     *  MSIE can't handle the DIV's properly, so we need to do it by
     *  hand.  Look up the style for typeset math to see if the user
     *  has changed it, and get whether it is centered or indented
     *  so we can mirror that using a SPAN
     */
    if (jsMath.browser == 'MSIE' && navigator.platform == 'Win32') {
      this.createMathTag = this.MSIEcreateMathTag;
      this.margin = ""; this.center = 0;
      for (var i = 0; i < jsMath.document.styleSheets.length; i++) {
        var rules = jsMath.document.styleSheets[i].cssRules;
        if (!rules) {rules = jsMath.document.styleSheets[i].rules}
        for (var j = 0; j < rules.length; j++) {
          if (rules[j].selectorText.toLowerCase() == 'div.typeset') {
            if (rules[j].style.margin != "") {this.margin = rules[j].style.margin}
            this.center = (rules[j].style.textAlign == 'center');
          }
        }
      }
    }
    this.inited = 1;
  },

  /*
   *  Test to see if we need to override the pattern exec() call
   *  (for MSIE on the Mac).
   */
  TestPatterns: function () {
    var pattern = /a/g;
    var match = pattern.exec("xax");
    this.fixPatterns = (pattern.lastIndex != 2 && match.lastIndex == 2);
  }

});

/*
 *  Initialize
 */
if (jsMath.Controls.cookie.tex2math == null) {jsMath.Controls.cookie.tex2math = 1}
if (jsMath.tex2math.allowDisableTag == null) {jsMath.tex2math.allowDisableTag = 1}
jsMath.tex2math.TestPatterns();
jsMath.tex2math.createPattern('stdPattern',/(\\[\(\)\[\]$]|\$\$|\$)/g);

}
