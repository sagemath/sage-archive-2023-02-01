/*
 *  autoload.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file is a plugin that checks if a page contains any math
 *  that must be processed by jsMath, and only loads jsMath.js
 *  when there is.
 *
 *  You can control the items to look for via the variables
 *
 *      jsMath.Autoload.findMathElements
 *      jsMath.Autoload.findTeXstrings
 *      jsMath.Autoload.findLaTeXstrings
 *      jsMath.Autoload.findCustomStrings
 *      jsMath.Autoload.findCustomSettings
 *
 *  which control whether to look for SPAN and DIV elements of class
 *  "math", TeX strings that will be converted by jsMath.ConvertTeX(), or
 *  LaTeX strings that will be converted by jsMath.ConvertLaTeX().  By
 *  default, the first is true and the second and third are false.  The
 *  findCustomStrings can be used to specify your own delimiters for
 *  in-line and display mathematics, e.g
 *
 *      jsMath.Autoload.findCustomStrings = [
 *         '[math],'[/math]',          // start and end in-line math
 *         '[display]','[/display]'    // start and end display math
 *      ];
 *
 *  Finally, findCustomSettings can be set to an object reference whose
 *  name:value pairs control the individual search settings for tex2math.
 *  (See the plugins/tex2math.js file for more details).
 *
 *  If any math strings are found, jsMath.js will be loaded automatically,
 *  but not loaded otherwise.  If any of the last four are set and TeX math
 *  strings are found, then plugins/tex2ath.js will be loaded
 *  automatically.  jsMath.Autoload.needsJsMath will be set to true or
 *  false depending on whether jsMath needed to be loaded.
 *
 *  The value of jsMath.Autoload.element controls the element to be
 *  searched by the autoload plug-in.  If unset, the complete document will
 *  be searched.  If set to a string, the element with that name will be
 *  searched.  If set to a DOM object, that object and its children will
 *  be searched.
 *
 *  Finally, there are two additional parameters that control files to
 *  be loaded after jsMath.js, should it be needed.  These are
 *
 *      jsMath.Autoload.loadFonts
 *      jsMath.Autoload.loadFiles
 *
 *  If jsMath.js is loaded, the fonts contained in the loadFonts array
 *  will be loaded, and the JavaScript files listed in the loadFiles array
 *  will be run.  Relative URL's are loaded based from the URL containing
 *  jsMath.js.
 *
 *  The autoload plugin can be loaded in the document HEAD or in the BODY.
 *  If it is loaded in the HEAD, you will need to call jsMath.Autoload.Check()
 *  at the end of the BODY (say in the window.onload handler) in order to
 *  get it to check the page for math that needs to be tagged, otherwise load
 *  the file at the bottom of the BODY and it will run the check automatically.
 *
 *  You can call jsMath.Autoload.Run() after the check has been performed
 *  in order to call the appropriate tex2math routines for the given Autoload
 *  settings.  You can call jsMath.Autoload.Run() even when jsMath isn't loaded.
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

/*************************************************************************/

/*
 *  Make sure jsMath.Autoload is available
 */
if (!window.jsMath) {window.jsMath = {}}
if (jsMath.Autoload == null) {jsMath.Autoload = {}}
jsMath.Add = function (dst,src) {for (var id in src) {dst[id] = src[id]}},
jsMath.document = document; // tex2math needs this

jsMath.Add(jsMath.Autoload,{

  Script: {

    request: null,  // XMLHttpRequest object (if we can get it)
    iframe: null,   // the hidden iframe (if not)

    /*
     *  Get XMLHttpRequest object, if possible, and look up the URL root
     */
    Init: function () {
      if (window.XMLHttpRequest) {try {this.request = new XMLHttpRequest} catch (err) {}}
      if (!this.request && window.ActiveXObject) {
        var xml = ["MSXML2.XMLHTTP.5.0","MSXML2.XMLHTTP.4.0","MSXML2.XMLHTTP.3.0",
                   "MSXML2.XMLHTTP","Microsoft.XMLHTTP"];
        for (var i = 0; i < xml.length && !this.request; i++) {
          try {this.request = new ActiveXObject(xml[i])} catch (err) {}
        }
      }
      this.Root();
    },

    /*
     *  Load an external JavaScript file
     */
    Load: function (url) {
      if (this.request) {
        setTimeout(function () {jsMath.Autoload.Script.xmlLoad(url)},1);
      } else {
        this.startLoad(url);
      }
    },

    /*
     *  Load an external JavaScript file via XMLHttpRequest
     */
    xmlLoad: function (url) {
      try {
        this.request.open("GET",jsMath.Autoload.root+url,false);
        this.request.send(null);
      } catch (err) {
        throw "autoload: can't load the file '"+url+"'\n"
            + "Message: "+err.message;
      }
      if (this.request.status && this.request.status >= 400) {
        throw "autoload: can't load the file '"+url+"'\n"
            + "Error status: "+this.request.status;
      }
      window.eval(this.request.responseText);
      this.endLoad();
    },

    /*
     *  Load an external JavaScript file via jsMath-autoload.html
     */
    startLoad: function (url) {
      this.iframe = document.createElement('iframe');
      this.iframe.style.visibility = 'hidden';
      this.iframe.style.position = 'absolute';
      this.iframe.style.width  = '0px';
      this.iframe.style.height = '0px';
      if (document.body.firstChild) {
        document.body.insertBefore(this.iframe,document.body.firstChild);
      } else {
        document.body.appendChild(this.iframe);
      }
      this.url = url; setTimeout('jsMath.Autoload.Script.setURL()',100);
    },
    endLoad: function () {setTimeout('jsMath.Autoload.Script.AfterLoad()',1)},

    /*
     *  Use location.replace() to avoid browsers placing the file in
     *  the history (and messing up the BACK button action).  Not
     *  completely effective in Firefox 1.0.x.  Safari won't handle
     *  replace() if the document is local (not sure if that's really
     *  the issue, but that's the only time I see it).
     */
    setURL: function () {
      var url = jsMath.Autoload.root+"jsMath-autoload.html";
      var doc = this.iframe.contentDocument;
      if (!doc && this.iframe.contentWindow) {doc = this.iframe.contentWindow.document}
      if (navigator.vendor == "Apple Computer, Inc." &&
          document.location.protocol == 'file:') {doc = null}
      if (doc) {doc.location.replace(url)} else {this.iframe.src = url}
    },

    /*
     *  Queue items that need to be postponed until jsMath has run
     */
    queue: [],
    Push: function (name,data) {this.queue[this.queue.length] = [name,data]},
    RunStack: function () {
      if (this.tex2math) {jsMath.Autoload.Check2(); return}
      for (var i = 0; i < this.queue.length; i++) {
        var name = this.queue[i][0];
        var data = this.queue[i][1];
        if (data.length == 1) {jsMath[name](data[0])}
          else {jsMath[name](data[0],data[1],data[2],data[3])}
      }
     this.queue = [];
    },

    AfterLoad: function () {jsMath.Autoload.Script.RunStack()},

    /*
     *  Look up the jsMath root directory, if it is not already supplied
     */
    Root: function () {
      if (jsMath.Autoload.root) return;
      var script = document.getElementsByTagName('script');
      if (script) {
        for (var i = 0; i < script.length; i++) {
          var src = script[i].src;
          if (src && src.match('(^|/)plugins/autoload.js$')) {
            jsMath.Autoload.root = src.replace(/plugins\/autoload.js$/,'');
            break;
          }
        }
      }
    }

  },

  /**************************************************************/

  /*
   *  Load tex2math first (so we can call its search functions
   *  to look to see if anything needs to be turned into math)
   *  if it is needed, otherwise go on to the second check.
   */
  Check: function () {
    if (this.checked) return; this.checked = 1;
    if (this.findTeXstrings || this.findLaTeXstrings ||
        this.findCustomStrings || this.findCustomSettings) {
      this.Script.tex2math = 1;
      this.Script.Load('plugins/tex2math.js');
    } else {
      if (!jsMath.tex2math) {jsMath.tex2math = {}}
      this.Check2();
    }
  },
  ReCheck: function () {
    if (jsMath.loaded) return;
    this.InitStubs();
    this.checked = 0;
    this.Script.queue = [];
    this.Check();
  },

  /*
   *  Once tex2math is loaded, use it to check for math that
   *  needs to be tagged for jsMath, and load jsMath if it is needed
   */
  Check2: function () {
    this.Script.tex2math = 0;
    this.needsJsMath = 0; if (this.checkElement == null) {this.checkElement = null}
    if (this.findMathElements) {
      this.needsJsMath = this.areMathElements(this.checkElement);
    }
    jsMath.tex2math.callback = this.tex2mathCallback;
    if (this.findTeXstrings && !this.needsJsMath) {
      jsMath.tex2math.ConvertTeX(this.checkElement);
    }
    if (this.findLaTeXstrings && !this.needsJsMath) {
      jsMath.tex2math.ConvertLaTeX(this.checkElement);
    }
    if (this.findCustomSettings && !this.needsJsMath) {
      jsMath.tex2math.Convert(this.checkElement,this.findCustomSettings);
    }
    if (this.findCustomStrings && !this.needsJsMath) {
      var s = this.findCustomStrings;
      jsMath.tex2math.CustomSearch(s[0],s[1],s[2],s[3]);
      jsMath.tex2math.ConvertCustom(this.checkElement);
    }
    jsMath.tex2math.callback = null;

    if (this.needsJsMath) {
      this.LoadJsMath();
    } else {
      jsMath.Process = function () {};
      jsMath.ProcessBeforeShowing = function () {};
      jsMath.ConvertTeX = function () {};
      jsMath.ConvertTeX2 = function () {};
      jsMath.ConvertLaTeX = function () {};
      jsMath.ConvertCustom = function () {};
      jsMath.CustomSearch = function () {};
      jsMath.Macro = function () {};
      jsMath.Synchronize = function (code,data) {
        if (typeof(code) != 'string') {code(data)} else {eval(code)}
      };
      jsMath.Autoload.Script.RunStack(); // perform pending commands
      jsMath.Autoload.setMessage();
    }
  },

  /*
   *  A callback used in the tex2math searches to signal that
   *  some math has been found.
   */
  tex2mathCallback: function () {
    jsMath.Autoload.needsJsMath = 1;
    return false;
  },

  /*
   *  jsMath.Autoload.Run() can be called to perform the
   *  tex2math calls given by the Autoload parameters.
   */
  Run: function (data) {
    if (jsMath.loaded) {this.Autorun(data)}
                  else {this.Script.Push('Autorun',[data])}
  },

  Autorun: function (element) {
    if (!element) {element = this.checkElement}
    if (this.findTeXstrings) {jsMath.ConvertTeX(element)}
    if (this.findLaTeXstrings) {jsMath.ConvertLaTeX(element)}
    if (this.findCustomSettings) {
      jsMath.Synchronize(function () {
        jsMath.tex2math.Convert(element,jsMath.Autoload.findCustomSettings);
      });
    }
    if (this.findCustomStrings) {
      var s = this.findCustomStrings;
      jsMath.CustomSearch(s[0],s[1],s[2],s[3]);
      jsMath.ConvertCustom(element);
    }
  },

  /*
   *  Look to see if there are SPAN or DIV elements of class "math".
   */
  areMathElements: function (obj) {
    if (!obj) {obj = document}
    if (typeof(obj) == 'string') {obj = document.getElementById(obj)}
    if (!obj.getElementsByTagName) {return false}
    var math = obj.getElementsByTagName('div');
    for (var k = 0; k < math.length; k++)
      {if (math[k].className.match(/(^| )math( |$)/)) {return true}}
    math = obj.getElementsByTagName('span');
    for (var k = 0; k < math.length; k++)
      {if (math[k].className.match(/(^| )math( |$)/)) {return true}}
    return false;
  },

  /*
   *  When math tags are found, load the jsMath.js file,
   *  and afterward, load any auxiliary files or fonts,
   *  and then do any pending commands.
   */
  LoadJsMath: function () {
    if (jsMath.loaded) {this.afterLoad(); return}
    if (this.root) {
      this.setMessage('Loading jsMath...');
      this.Script.AfterLoad = this.afterLoad;
      this.Script.Load('jsMath.js');
    } else {
      alert("Can't determine URL for jsMath.js");
    }
  },
  afterLoad: function () {
    if (jsMath.tex2math.window) {jsMath.tex2math.window.jsMath = jsMath}
    //
    //  Handle MSIE bug where jsMath.window both is and is not the actual window
    //
    if (jsMath.browser == 'MSIE') {window.onscroll = jsMath.window.onscroll};
    var fonts = jsMath.Autoload.loadFonts;
    if (fonts) {
      if (typeof(fonts) != 'object') {fonts = [fonts]}
      for (var i in fonts) {jsMath.Font.Load(fonts[i])}
    }
    var files = jsMath.Autoload.loadFiles;
    if (files) {
      if (typeof(files) != 'object') {files = [files]}
      for (var i in files) {jsMath.Setup.Script(files[i])}
    }
    jsMath.Synchronize(function () {jsMath.Autoload.Script.RunStack()});
    jsMath.Autoload.setMessage();
  },

  /*
   *  Display a message in a small box at the bottom of the screen
   */
  setMessage: function (message) {
    if (message) {
      this.div = document.createElement('div');
      if (!document.body.hasChildNodes) {document.body.appendChild(this.div)}
        else {document.body.insertBefore(this.div,document.body.firstChild)}
      var style = {
        position:'absolute', bottom:'1px', left:'2px',
        backgroundColor:'#E6E6E6', border:'solid 1px #959595',
        margin:'0px', padding:'1px 8px', zIndex:102,
        color:'black', fontSize:'75%', width:'auto'
      };
      for (var id in style) {this.div.style[id] = style[id]}
      this.div.appendChild(jsMath.document.createTextNode(message));
    } else if (this.div) {
      this.div.firstChild.nodeValue = "";
      this.div.style.visibility = 'hidden';
    }
  },

  /*
   *  Queue these so we can do them after jsMath has been loaded
   */
  stubs: {
    Process: function (data) {jsMath.Autoload.Script.Push('Process',[data])},
    ProcessBeforeShowing: function (data) {jsMath.Autoload.Script.Push('ProcessBeforeShowing',[data])},
    ConvertTeX: function (data) {jsMath.Autoload.Script.Push('ConvertTeX',[data])},
    ConvertTeX2: function (data) {jsMath.Autoload.Script.Push('ConvertTeX2',[data])},
    ConvertLaTeX: function (data) {jsMath.Autoload.Script.Push('ConvertLaTeX',[data])},
    ConvertCustom: function (data) {jsMath.Autoload.Script.Push('ConvertCustom',[data])},
    CustomSearch: function (d1,d2,d3,d4) {jsMath.Autoload.Script.Push('CustomSearch',[d1,d2,d3,d4])},
    Synchronize: function (data) {jsMath.Autoload.Script.Push('Synchronize',[data])},
    Macro: function (cs,def,params) {jsMath.Autoload.Script.Push('Macro',[cs,def,params])},
    Autorun: function (data) {jsMath.Autoload.Autorun(data)}
  },

  InitStubs: function () {jsMath.Add(jsMath,jsMath.Autoload.stubs)}

});

/*
 *  Initialize
 */

if (jsMath.Autoload.findMathElements == null) {jsMath.Autoload.findMathElements = 1}
if (jsMath.Autoload.findTeXstrings == null)   {jsMath.Autoload.findTeXstrings = 0}
if (jsMath.Autoload.findLaTeXstrings == null) {jsMath.Autoload.findLaTeXstrings = 0}

jsMath.Autoload.Script.Init();
jsMath.Autoload.InitStubs();
if (document.body) {jsMath.Autoload.Check()}
