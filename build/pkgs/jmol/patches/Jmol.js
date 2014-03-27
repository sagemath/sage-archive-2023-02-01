/* Jmol 12.0 script library Jmol.js 9:48 PM 1/31/2011 Bob Hanson

 checkbox heirarchy -- see http://chemapps.stolaf.edu/jmol/docs/examples-11/check.htm

    based on:
 *
 * Copyright (C) 2004-2005  Miguel, Jmol Development, www.jmol.org
 *
 * Contact: hansonr@stolaf.edu
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This library is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 *  02111-1307  USA.
 */

// for documentation see www.jmol.org/jslibrary

try{if(typeof(_jmol)!="undefined")exit()

// place "?NOAPPLET" on your command line to check applet control action with a textarea
// place "?JMOLJAR=xxxxx" to use a specific jar file

// bob hanson -- jmolResize(w,h) -- resizes absolutely or by percent (w or h 0.5 means 50%)
//    angel herraez -- update of jmolResize(w,h,targetSuffix) so it is not tied to first applet
// bob hanson -- jmolEvaluate -- evaluates molecular math 8:37 AM 2/23/2007
// bob hanson -- jmolScriptMessage -- returns all "scriptStatus" messages 8:37 AM 2/23/2007
// bob hanson -- jmolScriptEcho -- returns all "scriptEcho" messages 8:37 AM 2/23/2007
// bob hanson -- jmolScriptWait -- 11:31 AM 5/2/2006
// bob hanson -- remove trailing separatorHTML in radio groups -- 12:18 PM 5/6/2006
// bob hanson -- adds support for dynamic DOM script nodes 7:04 AM 5/19/2006
// bob hanson -- adds try/catch for wiki - multiple code passes 7:05 AM 5/19/2006
// bob hanson -- auto-initiates to defaultdir/defaultjar -- change as desired.
// bob hanson -- adding save/restore orientation w/ and w/o delay 11:49 AM 5/25/2006
// bob hanson -- adding AjaxJS service 11:16 AM 6/3/2006
// bob hanson -- fix for iframes not available for finding applet
// bob hanson -- added applet fake ?NOAPPLET URL flag
// bob hanson -- added jmolSetCallback(calbackName, funcName) 3:32 PM 6/13/2006
//			used PRIOR to jmolApplet() or jmolAppletInline()
//               added 4th array element in jmolRadioGroup -- title
//               added <span> and id around link, checkbox, radio, menu
//               fixing AJAX loads for MSIE/Opera-Mozilla incompatibility
//            -- renamed Jmol-11.js from Jmol-new.js; JmolApplet.jar from JmolAppletProto.jar
//	 	 renamed Jmol.js for Jmol 11 distribution
//            -- modified jmolRestoreOrientation() to be immediate, no 1-second delay
// bob hanson -- jmolScriptWait always returns a string -- 11:23 AM 9/16/2006
// bh         -- jmolCommandInput()
// bh         -- jmolSetTranslation(TF) -- forces translation even if there might be message callback issues
// bh         -- minor fixes suggested by Angel
// bh         -- adds jmolSetSyncId() and jmolGetSyncId()
// bh 3/2008  -- adds jmolAppendInlineScript() and jmolAppendInlineArray()
// bh 3/2008  -- fixes IE7 bug in relation to jmolLoadInlineArray()
// bh 6/2008  -- adds jmolSetAppletWindow()
// Angel H. 6/2008  -- added html <label> tags to checkboxes and radio buttons [in jmolCheckbox() and _jmolRadio() functions]
// bh 7/2008  -- code fix "for(i..." not "for(var i..."
// bh 12/2008 -- jmolLoadInline, jmolLoadInlineArray, jmolLoadInlineScript, jmolAppendInlineScript, jmolAppendInlineArray all return error message or null (Jmol 11.7.16)
// bh 12/2008 -- jmolScriptWaitOutput() -- waits for script to complete and delivers output normally sent to console

// bh 5/2009  -- Support for XHTML using jmolSetXHTML(id)
// ah & bh 6/2009 -- New jmolResizeApplet() more flexible, similar to jmolApplet() size syntax
// bh 11/2009 -- care in accessing top.document
// bh 12/2009 -- added jmolSetParameter(name, value)
// bh 12/2009 -- added PARAMS=name:value;name:value;name:value... for command line
// bh 12/2009 -- overhaul of target checking
// bh 1/2010  -- all _xxxx() methods ALWAYS have complete argument list
// bh 1/2010  -- adds option to run a JavaScript function from any Jmol control.
//               This is accomplished by passing an array rather than a script:
//               jmolHref([myfunc,"my param 1", "my param 2"], "testing")
//               function myfunc(jmolControlObject, [myfunc,"my param 1", "my param 2"], target){...}
//               and allows much more flexibility with responding to controls
// bh 4/2010  -- added jmolSetMemoryMb(nMb)
// ah 1/2011  -- wider detection of browsers; more browsers now use the object tag instead of the applet tag;
//               fix of object tag (removed classid) accounts for change of behavior in Chrome

var defaultdir = "."
var defaultjar = "JmolApplet.jar"


// Note added 12:41 PM 9/21/2008 by Bob Hanson, hansonr@stolaf.edu:

// JMOLJAR=xxxxx.jar on the URL for this page will override
// the JAR file specified in the jmolInitialize() call.

// The idea is that it can be very useful to test a web page with different JAR files
// Or for an expert user to substitute a signed applet for an unsigned one
// so as to use a broader range of models or to create JPEG files, for example.

// If the JAR file is not in the current directory (has any sort of "/" in its name)
// then the user is presented with a warning and asked whether it is OK to change Jar files.
// The default action, if the user just presses "OK" is to NOT allow the change.
// The user must type the word "yes" in the prompt box for the change to be approved.

// If you don't want people to be able to switch in their own JAR file on your page,
// simply set this next line to read "var allowJMOLJAR = false".


var undefined; // for IE 5 ... wherein undefined is undefined

////////////////////////////////////////////////////////////////
// Basic Scripting infrastruture
////////////////////////////////////////////////////////////////

function jmolInitialize(codebaseDirectory, fileNameOrUseSignedApplet) {
  if (_jmol.initialized)
    return;
  _jmol.initialized = true;
  if(_jmol.jmoljar) {
    var f = _jmol.jmoljar;
    if (f.indexOf("/") >= 0) {
      alert ("This web page URL is requesting that the applet used be " + f + ". This is a possible security risk, particularly if the applet is signed, because signed applets can read and write files on your local machine or network.")
      var ok = prompt("Do you want to use applet " + f + "? ","yes or no")
      if (ok == "yes") {
        codebaseDirectory = f.substring(0, f.lastIndexOf("/"));
        fileNameOrUseSignedApplet = f.substring(f.lastIndexOf("/") + 1);
      } else {
	_jmolGetJarFilename(fileNameOrUseSignedApplet);
        alert("The web page URL was ignored. Continuing using " + _jmol.archivePath + ' in directory "' + codebaseDirectory + '"');
      }
    } else {
      fileNameOrUseSignedApplet = f;
    }
  }
  _jmolSetCodebase(codebaseDirectory);
  _jmolGetJarFilename(fileNameOrUseSignedApplet);
  _jmolOnloadResetForms();
}

function jmolSetTranslation(TF) {
  _jmol.params.doTranslate = ''+TF;
}

function jmolToSigned(){
    _jmolGetJarFilename(true);
}

function _jmolGetJarFilename(fileNameOrFlag) {
  _jmol.archivePath =
    (typeof(fileNameOrFlag) == "string"  ? fileNameOrFlag : (fileNameOrFlag ?  "JmolAppletSigned" : "JmolApplet") + "0.jar");
}

function jmolSetDocument(doc) {
  _jmol.currentDocument = doc;
}

function jmolSetAppletColor(boxbgcolor, boxfgcolor, progresscolor) {
  _jmolInitCheck();
  _jmol.params.boxbgcolor = boxbgcolor;
  if (boxfgcolor)
    _jmol.params.boxfgcolor = boxfgcolor
  else if (boxbgcolor == "white" || boxbgcolor == "#FFFFFF")
    _jmol.params.boxfgcolor = "black";
  else
    _jmol.params.boxfgcolor = "white";
  if (progresscolor)
    _jmol.params.progresscolor = progresscolor;
  if (_jmol.debugAlert)
    alert(" boxbgcolor=" + _jmol.params.boxbgcolor +
          " boxfgcolor=" + _jmol.params.boxfgcolor +
          " progresscolor=" + _jmol.params.progresscolor);
}

function jmolSetAppletWindow(w) {
  _jmol.appletWindow = w;
}

function jmolApplet(size, script, nameSuffix) {
  _jmolInitCheck();
  return _jmolApplet(size, null, script, nameSuffix);
}

////////////////////////////////////////////////////////////////
// Basic controls
////////////////////////////////////////////////////////////////

// undefined means it wasn't there; null means it was explicitly listed as null (so as to skip it)

function jmolButton(script, label, id, title) {
  _jmolInitCheck();
  id != undefined && id != null || (id = "jmolButton" + _jmol.buttonCount);
  label != undefined && label != null || (label = script.substring(0, 32));
  ++_jmol.buttonCount;
  var scriptIndex = _jmolAddScript(script);
  var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input type='button' name='" + id + "' id='" + id +
          "' value='" + label +
          "' onclick='_jmolClick(this," + scriptIndex + _jmol.targetText +
          ")' onmouseover='_jmolMouseOver(" + scriptIndex +
          ");return true' onmouseout='_jmolMouseOut()' " +
          _jmol.buttonCssText + " /></span>";
  if (_jmol.debugAlert)
    alert(t);
  return _jmolDocumentWrite(t);
}

function jmolCheckbox(scriptWhenChecked, scriptWhenUnchecked,
                      labelHtml, isChecked, id, title) {
  _jmolInitCheck();
  id != undefined && id != null || (id = "jmolCheckbox" + _jmol.checkboxCount);
  ++_jmol.checkboxCount;
  if (scriptWhenChecked == undefined || scriptWhenChecked == null ||
      scriptWhenUnchecked == undefined || scriptWhenUnchecked == null) {
    alert("jmolCheckbox requires two scripts");
    return;
  }
  if (labelHtml == undefined || labelHtml == null) {
    alert("jmolCheckbox requires a label");
    return;
  }
  var indexChecked = _jmolAddScript(scriptWhenChecked);
  var indexUnchecked = _jmolAddScript(scriptWhenUnchecked);
  var eospan = "</span>"
  var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input type='checkbox' name='" + id + "' id='" + id +
          "' onclick='_jmolCbClick(this," +
          indexChecked + "," + indexUnchecked + _jmol.targetText +
          ")' onmouseover='_jmolCbOver(this," + indexChecked + "," +
          indexUnchecked +
          ");return true' onmouseout='_jmolMouseOut()' " +
	  (isChecked ? "checked='true' " : "")+ _jmol.checkboxCssText + " />"
  if (labelHtml.toLowerCase().indexOf("<td>")>=0) {
	t += eospan
	eospan = "";
  }
  t += "<label for=\"" + id + "\">" + labelHtml + "</label>" +eospan;
  if (_jmol.debugAlert)
    alert(t);
  return _jmolDocumentWrite(t);
}

function jmolStartNewRadioGroup() {
  ++_jmol.radioGroupCount;
}

function jmolRadioGroup(arrayOfRadioButtons, separatorHtml, groupName, id, title) {
  /*

    array: [radio1,radio2,radio3...]
    where radioN = ["script","label",isSelected,"id","title"]

  */

  _jmolInitCheck();
  var type = typeof arrayOfRadioButtons;
  if (type != "object" || type == null || ! arrayOfRadioButtons.length) {
    alert("invalid arrayOfRadioButtons");
    return;
  }
  separatorHtml != undefined && separatorHtml != null || (separatorHtml = "&nbsp; ");
  var len = arrayOfRadioButtons.length;
  jmolStartNewRadioGroup();
  groupName || (groupName = "jmolRadioGroup" + (_jmol.radioGroupCount - 1));
  var t = "<span id='"+(id ? id : groupName)+"'>";
  for (var i = 0; i < len; ++i) {
    if (i == len - 1)
      separatorHtml = "";
    var radio = arrayOfRadioButtons[i];
    type = typeof radio;
    if (type == "object") {
      t += _jmolRadio(radio[0], radio[1], radio[2], separatorHtml, groupName, (radio.length > 3 ? radio[3]: (id ? id : groupName)+"_"+i), (radio.length > 4 ? radio[4] : 0), title);
    } else {
      t += _jmolRadio(radio, null, null, separatorHtml, groupName, (id ? id : groupName)+"_"+i, title);
    }
  }
  t+="</span>"
  if (_jmol.debugAlert)
    alert(t);
  return _jmolDocumentWrite(t);
}


function jmolRadio(script, labelHtml, isChecked, separatorHtml, groupName, id, title) {
  _jmolInitCheck();
  if (_jmol.radioGroupCount == 0)
    ++_jmol.radioGroupCount;
  var t = _jmolRadio(script, labelHtml, isChecked, separatorHtml, groupName, (id ? id : groupName + "_" + _jmol.radioCount), title ? title : 0);
  if (_jmol.debugAlert)
    alert(t);
  return _jmolDocumentWrite(t);
}

function jmolLink(script, label, id, title) {
  _jmolInitCheck();
  id != undefined && id != null || (id = "jmolLink" + _jmol.linkCount);
  label != undefined && label != null || (label = script.substring(0, 32));
  ++_jmol.linkCount;
  var scriptIndex = _jmolAddScript(script);
  var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><a name='" + id + "' id='" + id +
          "' href='javascript:_jmolClick(this," + scriptIndex + _jmol.targetText + ");' onmouseover='_jmolMouseOver(" + scriptIndex +
          ");return true;' onmouseout='_jmolMouseOut()' " +
          _jmol.linkCssText + ">" + label + "</a></span>";
  if (_jmol.debugAlert)
    alert(t);
  return _jmolDocumentWrite(t);
}

function jmolCommandInput(label, size, id, title) {
  _jmolInitCheck();
  id != undefined && id != null || (id = "jmolCmd" + _jmol.cmdCount);
  label != undefined && label != null || (label = "Execute");
  size != undefined && !isNaN(size) || (size = 60);
  ++_jmol.cmdCount;
  var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input name='" + id + "' id='" + id +
          "' size='"+size+"' onkeypress='_jmolCommandKeyPress(event,\""+id+"\"" + _jmol.targetText + ")'><input type=button value = '"+label+"' onclick='jmolScript(document.getElementById(\""+id+"\").value" + _jmol.targetText + ")' /></span>";
  if (_jmol.debugAlert)
    alert(t);
  return _jmolDocumentWrite(t);
}

function _jmolCommandKeyPress(e, id, target) {
	var keycode = (window.event ? window.event.keyCode : e ? e.which : 0);
	if (keycode == 13) {
		var inputBox = document.getElementById(id)
		_jmolScriptExecute(inputBox, inputBox.value, target)
	}
}

function _jmolScriptExecute(element,script,target) {
	if (typeof(script) == "object")
		script[0](element, script, target)
	else
		jmolScript(script, target)
}

function jmolMenu(arrayOfMenuItems, size, id, title) {
  _jmolInitCheck();
  id != undefined && id != null || (id = "jmolMenu" + _jmol.menuCount);
  ++_jmol.menuCount;
  var type = typeof arrayOfMenuItems;
  if (type != null && type == "object" && arrayOfMenuItems.length) {
    var len = arrayOfMenuItems.length;
    if (typeof size != "number" || size == 1)
      size = null;
    else if (size < 0)
      size = len;
    var sizeText = size ? " size='" + size + "' " : "";
    var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><select name='" + id + "' id='" + id +
            "' onChange='_jmolMenuSelected(this" + _jmol.targetText + ")'" +
            sizeText + _jmol.menuCssText + ">";
    for (var i = 0; i < len; ++i) {
      var menuItem = arrayOfMenuItems[i];
      type = typeof menuItem;
      var script, text;
      var isSelected = undefined;
      if (type == "object" && menuItem != null) {
        script = menuItem[0];
        text = menuItem[1];
        isSelected = menuItem[2];
      } else {
        script = text = menuItem;
      }
      text != undefined && text != null || (text = script);
      if (script=="#optgroup") {
        t += "<optgroup label='" + text + "'>";
	  } else if (script=="#optgroupEnd") {
        t += "</optgroup>";
	  } else {
        var scriptIndex = _jmolAddScript(script);
        var selectedText = isSelected ? "' selected='true'>" : "'>";
        t += "<option value='" + scriptIndex + selectedText + text + "</option>";
      }
    }
    t += "</select></span>";
    if (_jmol.debugAlert)
      alert(t);
    return _jmolDocumentWrite(t);
  }
}

function jmolHtml(html) {
  return _jmolDocumentWrite(html);
}

function jmolBr() {
  return _jmolDocumentWrite("<br />");
}

////////////////////////////////////////////////////////////////
// advanced scripting functions
////////////////////////////////////////////////////////////////

function jmolDebugAlert(enableAlerts) {
  _jmol.debugAlert = (enableAlerts == undefined || enableAlerts)
}

function jmolAppletInline(size, inlineModel, script, nameSuffix) {
  _jmolInitCheck();
  return _jmolApplet(size, _jmolSterilizeInline(inlineModel),
                     script, nameSuffix);
}

function jmolSetTarget(targetSuffix) {
  _jmol.targetSuffix = targetSuffix;
  _jmol.targetText = targetSuffix ? ",\"" + targetSuffix + "\"" : ",0";
}

function jmolScript(script, targetSuffix) {
  if (script) {
    _jmolCheckBrowser();
    if (targetSuffix == "all") {
      with (_jmol) {
	for (var i = 0; i < appletSuffixes.length; ++i) {
	  var applet = _jmolGetApplet(appletSuffixes[i]);
          if (applet) applet.script(script);
        }
      }
    } else {
      var applet=_jmolGetApplet(targetSuffix);
      if (applet) applet.script(script);
    }
  }
}

function jmolLoadInline(model, targetSuffix) {
  if (!model)return "ERROR: NO MODEL"
  var applet=_jmolGetApplet(targetSuffix);
  if (!applet)return "ERROR: NO APPLET"
  if (typeof(model) == "string")
    return applet.loadInlineString(model, "", false);
  else
    return applet.loadInlineArray(model, "", false);
}


function jmolLoadInlineScript(model, script, targetSuffix) {
  if (!model)return "ERROR: NO MODEL"
  var applet=_jmolGetApplet(targetSuffix);
  if (!applet)return "ERROR: NO APPLET"
  return applet.loadInlineString(model, script, false);
}


function jmolLoadInlineArray(ModelArray, script, targetSuffix) {
  if (!model)return "ERROR: NO MODEL"
  script || (script="")
  var applet=_jmolGetApplet(targetSuffix);
  if (!applet)return "ERROR: NO APPLET"
  try {
    return applet.loadInlineArray(ModelArray, script, false);
  } catch (err) {
    //IE 7 bug
    return applet.loadInlineString(ModelArray.join("\n"), script, false);
  }
}

function jmolAppendInlineArray(ModelArray, script, targetSuffix) {
  if (!model)return "ERROR: NO MODEL"
  script || (script="")
  var applet=_jmolGetApplet(targetSuffix);
  if (!applet)return "ERROR: NO APPLET"
  try {
    return applet.loadInlineArray(ModelArray, script, true);
  } catch (err) {
    //IE 7 bug
    return applet.loadInlineString(ModelArray.join("\n"), script, true);
  }
}

function jmolAppendInlineScript(model, script, targetSuffix) {
  if (!model)return "ERROR: NO MODEL"
  var applet=_jmolGetApplet(targetSuffix);
  if (!applet)return "ERROR: NO APPLET"
  return applet.loadInlineString(model, script, true);
}

function jmolCheckBrowser(action, urlOrMessage, nowOrLater) {
  if (typeof action == "string") {
    action = action.toLowerCase();
    action == "alert" || action == "redirect" || action == "popup" || (action = null);
  }
  if (typeof action != "string")
    alert("jmolCheckBrowser(action, urlOrMessage, nowOrLater)\n\n" +
          "action must be 'alert', 'redirect', or 'popup'");
  else {
    if (typeof urlOrMessage != "string")
      alert("jmolCheckBrowser(action, urlOrMessage, nowOrLater)\n\n" +
            "urlOrMessage must be a string");
    else {
      _jmol.checkBrowserAction = action;
      _jmol.checkBrowserUrlOrMessage = urlOrMessage;
    }
  }
  if (typeof nowOrLater == "string" && nowOrLater.toLowerCase() == "now")
    _jmolCheckBrowser();
}

////////////////////////////////////////////////////////////////
// Cascading Style Sheet Class support
////////////////////////////////////////////////////////////////

function jmolSetAppletCssClass(appletCssClass) {
  if (_jmol.hasGetElementById) {
    _jmol.appletCssClass = appletCssClass;
    _jmol.appletCssText = appletCssClass ? "class='" + appletCssClass + "' " : "";
  }
}

function jmolSetButtonCssClass(buttonCssClass) {
  if (_jmol.hasGetElementById) {
    _jmol.buttonCssClass = buttonCssClass;
    _jmol.buttonCssText = buttonCssClass ? "class='" + buttonCssClass + "' " : "";
  }
}

function jmolSetCheckboxCssClass(checkboxCssClass) {
  if (_jmol.hasGetElementById) {
    _jmol.checkboxCssClass = checkboxCssClass;
    _jmol.checkboxCssText = checkboxCssClass ? "class='" + checkboxCssClass + "' " : "";
  }
}

function jmolSetRadioCssClass(radioCssClass) {
  if (_jmol.hasGetElementById) {
    _jmol.radioCssClass = radioCssClass;
    _jmol.radioCssText = radioCssClass ? "class='" + radioCssClass + "' " : "";
  }
}

function jmolSetLinkCssClass(linkCssClass) {
  if (_jmol.hasGetElementById) {
    _jmol.linkCssClass = linkCssClass;
    _jmol.linkCssText = linkCssClass ? "class='" + linkCssClass + "' " : "";
  }
}

function jmolSetMenuCssClass(menuCssClass) {
  if (_jmol.hasGetElementById) {
    _jmol.menuCssClass = menuCssClass;
    _jmol.menuCssText = menuCssClass ? "class='" + menuCssClass + "' " : "";
  }
}

////////////////////////////////////////////////////////////////
// functions for INTERNAL USE ONLY which are subject to change
// use at your own risk ... you have been WARNED!
////////////////////////////////////////////////////////////////
var _jmol = {
  currentDocument: document,

  debugAlert: false,

  codebase: "",
  modelbase: ".",

  appletCount: 0,
  appletSuffixes: [],
  appletWindow: null,
  allowedJmolSize: [25, 2048, 300],   // min, max, default (pixels)
	  /*  By setting the _jmol.allowedJmolSize[] variable in the webpage
	      before calling jmolApplet(), limits for applet size can be overriden.
		    2048 standard for GeoWall (http://geowall.geo.lsa.umich.edu/home.html)
	  */
  buttonCount: 0,
  checkboxCount: 0,
  linkCount: 0,
  cmdCount: 0,
  menuCount: 0,
  radioCount: 0,
  radioGroupCount: 0,

  appletCssClass: null,
  appletCssText: "",
  buttonCssClass: null,
  buttonCssText: "",
  checkboxCssClass: null,
  checkboxCssText: "",
  java_arguments: "-Xmx512m",
  radioCssClass: null,
  radioCssText: "",
  linkCssClass: null,
  linkCssText: "",
  menuCssClass: null,
  menuCssText: "",

  targetSuffix: 0,
  targetText: ",0",
  scripts: [""],
  params: {
	syncId: ("" + Math.random()).substring(3),
	progressbar: "true",
	progresscolor: "blue",
	boxbgcolor: "black",
	boxfgcolor: "white",
	boxmessage: "Downloading JmolApplet ..."
  },
  ua: navigator.userAgent.toLowerCase(),
  // uaVersion: parseFloat(navigator.appVersion),  // not used

  os: "unknown",
  browser: "unknown",
  browserVersion: 0,
  hasGetElementById: !!document.getElementById,
  isJavaEnabled: navigator.javaEnabled(),
  // isNetscape47Win: false,  // not used, N4.7 is no longer supported even for detection
  useIEObject: false,
  useHtml4Object: false,

  windowsClassId: "clsid:8AD9C840-044E-11D1-B3E9-00805F499D93",
  windowsCabUrl:
   "http://java.sun.com/update/1.6.0/jinstall-6u22-windows-i586.cab",

  isBrowserCompliant: false,
  isJavaCompliant: false,
  isFullyCompliant: false,

  initialized: false,
  initChecked: false,

  browserChecked: false,
  checkBrowserAction: "alert",
  checkBrowserUrlOrMessage: null,

  archivePath: null, // JmolApplet0.jar OR JmolAppletSigned0.jar

  previousOnloadHandler: null,

  jmoljar: null,
  useNoApplet: false,

  ready: {}
}

with (_jmol) {
  function _jmolTestUA(candidate) {
    var ua = _jmol.ua;
    var index = ua.indexOf(candidate);
    if (index < 0)
      return false;
    _jmol.browser = candidate;
    _jmol.browserVersion = parseFloat(ua.substring(index+candidate.length+1));
    return true;
  }

  function _jmolTestOS(candidate) {
    if (_jmol.ua.indexOf(candidate) < 0)
      return false;
    _jmol.os = candidate;
    return true;
  }

  _jmolTestUA("konqueror") ||
  _jmolTestUA("webkit") ||
  _jmolTestUA("omniweb") ||
  _jmolTestUA("opera") ||
  _jmolTestUA("webtv") ||
  _jmolTestUA("icab") ||
  _jmolTestUA("msie") ||
  (_jmol.ua.indexOf("compatible") < 0 && _jmolTestUA("mozilla")); //Netscape, Mozilla, Seamonkey, Firefox and anything assimilated

  _jmolTestOS("linux") ||
  _jmolTestOS("unix") ||
  _jmolTestOS("mac") ||
  _jmolTestOS("win");

  isBrowserCompliant = hasGetElementById;
  // known exceptions (old browsers):
  if (browser == "opera" && browserVersion <= 7.54 && os == "mac"
      || browser == "webkit" && browserVersion < 125.12
      || browser == "msie" && os == "mac"
      || browser == "konqueror" && browserVersion <= 3.3
    ) {
    isBrowserCompliant = false;
  }

  // possibly more checks in the future for this
  isJavaCompliant = isJavaEnabled;

  isFullyCompliant = isBrowserCompliant && isJavaCompliant;

  useIEObject = (os == "win" && browser == "msie" && browserVersion >= 5.5);
  useHtml4Object =
   (browser == "mozilla" && browserVersion >= 5) ||
   (browser == "opera" && browserVersion >= 8) ||
   (browser == "webkit" && browserVersion >= 412.2);
 try {
  if (top.location.search.indexOf("JMOLJAR=")>=0)
    jmoljar = top.location.search.split("JMOLJAR=")[1].split("&")[0];
 } catch(e) {
  // can't access top.location
 }
 try {
  useNoApplet = (top.location.search.indexOf("NOAPPLET")>=0);
 } catch(e) {
  // can't access top.document
 }
}

function jmolSetMemoryMb(nMb) {
  _jmol.java_arguments = "-Xmx" + Math.round(nMb) + "m"
}

function jmolSetParameter(name,value) {
  _jmol.params[name] = value
}

function jmolSetCallback(callbackName,funcName) {
  _jmol.params[callbackName] = funcName
}

 try {
// note this is done FIRST, so it cannot override a setting done by the developer
  if (top.location.search.indexOf("PARAMS=")>=0) {
    var pars = unescape(top.location.search.split("PARAMS=")[1].split("&")[0]).split(";");
    for (var i = 0; i < pars.length; i++) {
      var p = pars[i].split(":");
      jmolSetParameter(p[0],p[1]);
    }
  }
 } catch(e) {
  // can't access top.location
 }

function jmolSetSyncId(n) {
  return _jmol.params["syncId"] = n
}

function jmolGetSyncId() {
  return _jmol.params["syncId"]
}

function jmolSetLogLevel(n) {
  _jmol.params.logLevel = ''+n;
}

	/*  AngelH, mar2007:
		By (re)setting these variables in the webpage before calling jmolApplet(),
		a custom message can be provided (e.g. localized for user's language) when no Java is installed.
	*/
if (noJavaMsg==undefined) var noJavaMsg =
        "You do not have Java applets enabled in your web browser, or your browser is blocking this applet.<br />\n" +
        "Check the warning message from your browser and/or enable Java applets in<br />\n" +
        "your web browser preferences, or install the Java Runtime Environment from <a href='http://www.java.com'>www.java.com</a><br />";
if (noJavaMsg2==undefined) var noJavaMsg2 =
        "You do not have the<br />\n" +
        "Java Runtime Environment<br />\n" +
        "installed for applet support.<br />\n" +
        "Visit <a href='http://www.java.com'>www.java.com</a>";
function _jmolApplet(size, inlineModel, script, nameSuffix) {
	/*  AngelH, mar2007
		Fixed percent / pixel business, to avoid browser errors:
		put "px" where needed, avoid where not.

	    Bob Hanson, 1/2010
		Fixed inline escape changing returns to |
	*/
  with (_jmol) {
    nameSuffix == undefined && (nameSuffix = appletCount);
    appletSuffixes.push(nameSuffix);
    ++appletCount;
    script || (script = "select *");
    var sz = _jmolGetAppletSize(size);
    var widthAndHeight = " width='" + sz[0] + "' height='" + sz[1] + "' ";
    var tHeader, tFooter;
    codebase || jmolInitialize(".");
    if (useIEObject || useHtml4Object) {
      params.archive = archivePath;
      params.mayscript = 'true';
      params.codebase = codebase;
      params.code = 'JmolApplet';
      tHeader =
        "<object name='jmolApplet" + nameSuffix +
        "' id='jmolApplet" + nameSuffix + "' " + appletCssText + "\n" +
				widthAndHeight + "\n";
      tFooter = "</object>";
    }
    if (java_arguments)
      params.java_arguments = java_arguments;
    if (useIEObject) { // use MSFT IE6 object tag with .cab file reference
      tHeader += " classid='" + windowsClassId + "'\n" +
      (windowsCabUrl ? " codebase='" + windowsCabUrl + "'\n" : "") + ">\n";
    } else if (useHtml4Object) { // use HTML4 object tag
      tHeader += " type='application/x-java-applet'\n>\n";
				/*	" classid='java:JmolApplet'\n" +	AH removed this
				  Chromium Issue 62076: 	Java Applets using an <object> with a classid paramater don't load.
					http://code.google.com/p/chromium/issues/detail?id=62076
					They say this is the correct behavior according to the spec, and there's no indication at this point
					that WebKit will be changing the handling, so eventually Safari will acquire this behavior too.
					Removing the classid parameter seems to be well tolerated by all browsers (even IE!).
				*/
    } else { // use applet tag
      tHeader =
        "<applet name='jmolApplet" + nameSuffix +
        "' id='jmolApplet" + nameSuffix + "' " + appletCssText + "\n" +
				widthAndHeight + "\n" +
        " code='JmolApplet'" +
        " archive='" + archivePath + "' codebase='" + codebase + "'\n" +
        " mayscript='true'>\n";
      tFooter = "</applet>";
    }
    var visitJava;
    if (useIEObject || useHtml4Object) {
		var szX = "width:" + sz[0]
		if ( szX.indexOf("%")==-1 ) szX+="px"
		var szY = "height:" + sz[1]
		if ( szY.indexOf("%")==-1 ) szY+="px"
      visitJava =
        "<p style='background-color:yellow; color:black; " +
		szX + ";" + szY + ";" +
        // why doesn't this vertical-align work?
	"text-align:center;vertical-align:middle;'>\n" +
		noJavaMsg +
        "</p>";
    } else {
      visitJava =
        "<table bgcolor='yellow'><tr>" +
        "<td align='center' valign='middle' " + widthAndHeight + "><font color='black'>\n" +
		noJavaMsg2 +
        "</font></td></tr></table>";
    }
    params.loadInline = (inlineModel ? inlineModel : "");
    params.script = (script ? _jmolSterilizeScript(script) : "");
    var t = tHeader + _jmolParams() + visitJava + tFooter;
    jmolSetTarget(nameSuffix);
    ready["jmolApplet" + nameSuffix] = false;
    if (_jmol.debugAlert)
      alert(t);
    return _jmolDocumentWrite(t);
  }
}

function _jmolParams() {
 var t = "";
 for (var i in _jmol.params)
	if(_jmol.params[i]!="")
		 t+="  <param name='"+i+"' value='"+_jmol.params[i]+"' />\n";
 return t
}

function _jmolInitCheck() {
  if (_jmol.initChecked)
    return;
  _jmol.initChecked = true;
  jmolInitialize(defaultdir, defaultjar)
}

function _jmolCheckBrowser() {
  with (_jmol) {
    if (browserChecked)
      return;
    browserChecked = true;

    if (isFullyCompliant)
      return true;

    if (checkBrowserAction == "redirect")
      location.href = checkBrowserUrlOrMessage;
    else if (checkBrowserAction == "popup")
      _jmolPopup(checkBrowserUrlOrMessage);
    else {
      var msg = checkBrowserUrlOrMessage;
      if (msg == null)
        msg = "Your web browser is not fully compatible with Jmol\n\n" +
              "browser: " + browser +
              "   version: " + browserVersion +
              "   os: " + os +
              "   isBrowserCompliant: " + isBrowserCompliant +
              "   isJavaCompliant: " + isJavaCompliant +
              "\n\n" + ua;
      alert(msg);
    }
  }
  return false;
}

function jmolSetXHTML(id) {
	_jmol.isXHTML = true
	_jmol.XhtmlElement = null
	_jmol.XhtmlAppendChild = false
	if (id){
		_jmol.XhtmlElement = document.getElementById(id)
		_jmol.XhtmlAppendChild = true
	}
}

function _jmolDocumentWrite(text) {
	if (_jmol.currentDocument) {
		if (_jmol.isXHTML && !_jmol.XhtmlElement) {
			var s = document.getElementsByTagName("script")
			_jmol.XhtmlElement = s.item(s.length - 1)
			_jmol.XhtmlAppendChild = false
		}
		if (_jmol.XhtmlElement) {
			_jmolDomDocumentWrite(text)
		} else {
			_jmol.currentDocument.write(text);
		}
	}
	return text;
}

function _jmolDomDocumentWrite(data) {
	var pt = 0
	var Ptr = []
	Ptr[0] = 0
	while (Ptr[0] < data.length) {
		var child = _jmolGetDomElement(data, Ptr)
		if (!child)break
		if (_jmol.XhtmlAppendChild)
			_jmol.XhtmlElement.appendChild(child)
		else
			_jmol.XhtmlElement.parentNode.insertBefore(child, _jmol.XhtmlElement);
	}
}
function _jmolGetDomElement(data, Ptr, closetag, lvel) {
	var e = document.createElement("span")
	e.innerHTML = data
	Ptr[0] = data.length
	return e

//unnecessary?

	closetag || (closetag = "")
	lvel || (lvel = 0)
	var pt0 = Ptr[0]
	var pt = pt0
	while (pt < data.length && data.charAt(pt) != "<") pt++
	if (pt != pt0) {
		var text = data.substring(pt0, pt)
		Ptr[0] = pt
		return document.createTextNode(text)
	}
	pt0 = ++pt
	var ch
	while (pt < data.length && "\n\r\t >".indexOf(ch = data.charAt(pt)) < 0) pt++
	var tagname = data.substring(pt0, pt)
	var e = (tagname == closetag  || tagname == "/" ? ""
		: document.createElementNS ? document.createElementNS('http://www.w3.org/1999/xhtml', tagname)
		: document.createElement(tagname));
	if (ch == ">") {
		Ptr[0] = ++pt
		return e
	}
	while (pt < data.length && (ch = data.charAt(pt)) != ">") {
		while (pt < data.length && "\n\r\t ".indexOf(ch = data.charAt(pt)) >= 0) pt++
		pt0 = pt
		while (pt < data.length && "\n\r\t =/>".indexOf(ch = data.charAt(pt)) < 0) pt++
		var attrname = data.substring(pt0, pt).toLowerCase()
		if (attrname && ch != "=")
			e.setAttribute(attrname, "true")
		while (pt < data.length && "\n\r\t ".indexOf(ch = data.charAt(pt)) >= 0) pt++
		if (ch == "/") {
			Ptr[0] = pt + 2
			return e
		} else if (ch == "=") {
			var quote = data.charAt(++pt)
			pt0 = ++pt
			while (pt < data.length && (ch = data.charAt(pt)) != quote) pt++
			var attrvalue = data.substring(pt0, pt)
			e.setAttribute(attrname, attrvalue)
			pt++
		}
	}
	Ptr[0] = ++pt
	while (Ptr[0] < data.length) {
		var child = _jmolGetDomElement(data, Ptr, "/" + tagname, lvel+1)
		if (!child)break
		e.appendChild(child)
	}
	return e
}

function _jmolPopup(url) {
  var popup = window.open(url, "JmolPopup",
                          "left=150,top=150,height=400,width=600," +
                          "directories=yes,location=yes,menubar=yes," +
                          "toolbar=yes," +
                          "resizable=yes,scrollbars=yes,status=yes");
  if (popup.focus)
    poup.focus();
}

function _jmolReadyCallback(name) {
  if (_jmol.debugAlert)
    alert(name + " is ready");
  _jmol.ready["" + name] = true;
}

function _jmolSterilizeScript(script) {
  script = script.replace(/'/g, "&#39;");
  if (_jmol.debugAlert)
    alert("script:\n" + script);
  return script;
}

function _jmolSterilizeInline(model) {
  model = model.replace(/\r|\n|\r\n/g, (model.indexOf("|") >= 0 ? "\\/n" : "|")).replace(/'/g, "&#39;");
  if (_jmol.debugAlert)
    alert("inline model:\n" + model);
  return model;
}

function _jmolRadio(script, labelHtml, isChecked, separatorHtml, groupName, id, title) {
  ++_jmol.radioCount;
  groupName != undefined && groupName != null || (groupName = "jmolRadioGroup" + (_jmol.radioGroupCount - 1));
  if (!script)
    return "";
  labelHtml != undefined && labelHtml != null || (labelHtml = script.substring(0, 32));
  separatorHtml || (separatorHtml = "")
  var scriptIndex = _jmolAddScript(script);
  var eospan = "</span>"
  var t = "<span id=\"span_"+id+"\""+(title ? " title=\"" + title + "\"":"")+"><input name='"
	+ groupName + "' id='"+id+"' type='radio' onclick='_jmolClick(this," +
         scriptIndex + _jmol.targetText + ");return true;' onmouseover='_jmolMouseOver(" +
         scriptIndex + ");return true;' onmouseout='_jmolMouseOut()' " +
	 (isChecked ? "checked='true' " : "") + _jmol.radioCssText + " />"
  if (labelHtml.toLowerCase().indexOf("<td>")>=0) {
	t += eospan
	eospan = "";
  }
  t += "<label for=\"" + id + "\">" + labelHtml + "</label>" +eospan + separatorHtml;

  return t;
}

function _jmolFindApplet(target) {
  // first look for the target in the current window
  var applet = _jmolFindAppletInWindow(_jmol.appletWindow != null ? _jmol.appletWindow : window, target);
  // THEN look for the target in child frames
  if (applet == undefined)
    applet = _jmolSearchFrames(window, target);
  // FINALLY look for the target in sibling frames
  if (applet == undefined)
    applet = _jmolSearchFrames(top, target); // look starting in top frame
  return applet;
}

function _jmolGetApplet(targetSuffix){
 var target = "jmolApplet" + (targetSuffix ? targetSuffix : "0");
 var applet = _jmolFindApplet(target);
 if (applet) return applet
 _jmol.alerted || alert("could not find applet " + target);
 _jmol.alerted = true;
 return null
}

function _jmolSearchFrames(win, target) {
  var applet;
  var frames = win.frames;
  if (frames && frames.length) { // look in all the frames below this window
   try{
    for (var i = 0; i < frames.length; ++i) {
      applet = _jmolSearchFrames(frames[i], target);
      if (applet)
        return applet;
    }
   }catch(e) {
	if (_jmol.debugAlert)
		alert("Jmol.js _jmolSearchFrames cannot access " + win.name + ".frame[" + i + "] consider using jmolSetAppletWindow()")
   }
  }
  return applet = _jmolFindAppletInWindow(win, target)
}

function _jmolFindAppletInWindow(win, target) {
    var doc = win.document;
		if (doc.getElementById(target))
      return doc.getElementById(target);
    else if (doc.applets)
      return doc.applets[target];
    else
      return doc[target];
}

function _jmolAddScript(script) {
  if (!script)
    return 0;
  var index = _jmol.scripts.length;
  _jmol.scripts[index] = script;
  return index;
}

function _jmolClick(elementClicked, scriptIndex, targetSuffix) {
  _jmol.element = elementClicked;
  _jmolScriptExecute(elementClicked, _jmol.scripts[scriptIndex], targetSuffix);
}

function _jmolMenuSelected(menuObject, targetSuffix) {
  var scriptIndex = menuObject.value;
  if (scriptIndex != undefined) {
    _jmolScriptExecute(menuObject, _jmol.scripts[scriptIndex], targetSuffix);
    return;
  }
  var len = menuObject.length;
  if (typeof len == "number") {
    for (var i = 0; i < len; ++i) {
      if (menuObject[i].selected) {
        _jmolClick(menuObject[i], menuObject[i].value, targetSuffix);
	return;
      }
    }
  }
  alert("?Que? menu selected bug #8734");
}


_jmol.checkboxMasters = {};
_jmol.checkboxItems = {};

function jmolSetCheckboxGroup(chkMaster,chkBox) {
	var id = chkMaster;
	if(typeof(id)=="number")id = "jmolCheckbox" + id;
	chkMaster = document.getElementById(id);
	if (!chkMaster)alert("jmolSetCheckboxGroup: master checkbox not found: " + id);
	var m = _jmol.checkboxMasters[id] = {};
	m.chkMaster = chkMaster;
	m.chkGroup = {};
	for (var i = 1; i < arguments.length; i++){
		var id = arguments[i];
		if(typeof(id)=="number")id = "jmolCheckbox" + id;
		checkboxItem = document.getElementById(id);
		if (!checkboxItem)alert("jmolSetCheckboxGroup: group checkbox not found: " + id);
		m.chkGroup[id] = checkboxItem;
		_jmol.checkboxItems[id] = m;
	}
}

function _jmolNotifyMaster(m){
	//called when a group item is checked
	var allOn = true;
	var allOff = true;
	for (var chkBox in m.chkGroup){
		if(m.chkGroup[chkBox].checked)
			allOff = false;
		else
			allOn = false;
	}
	if (allOn)m.chkMaster.checked = true;
	if (allOff)m.chkMaster.checked = false;
	if ((allOn || allOff) && _jmol.checkboxItems[m.chkMaster.id])
		_jmolNotifyMaster(_jmol.checkboxItems[m.chkMaster.id])
}

function _jmolNotifyGroup(m, isOn){
	//called when a master item is checked
	for (var chkBox in m.chkGroup){
		var item = m.chkGroup[chkBox]
		item.checked = isOn;
		if (_jmol.checkboxMasters[item.id])
			_jmolNotifyGroup(_jmol.checkboxMasters[item.id], isOn)
	}
}

function _jmolCbClick(ckbox, whenChecked, whenUnchecked, targetSuffix) {
  _jmol.control = ckbox
  _jmolClick(ckbox, ckbox.checked ? whenChecked : whenUnchecked, targetSuffix);
  if(_jmol.checkboxMasters[ckbox.id])
	_jmolNotifyGroup(_jmol.checkboxMasters[ckbox.id], ckbox.checked)
  if(_jmol.checkboxItems[ckbox.id])
	_jmolNotifyMaster(_jmol.checkboxItems[ckbox.id])
}

function _jmolCbOver(ckbox, whenChecked, whenUnchecked) {
  window.status = _jmol.scripts[ckbox.checked ? whenUnchecked : whenChecked];
}

function _jmolMouseOver(scriptIndex) {
  window.status = _jmol.scripts[scriptIndex];
}

function _jmolMouseOut() {
  window.status = " ";
  return true;
}

function _jmolSetCodebase(codebase) {
  _jmol.codebase = codebase ? codebase : ".";
  if (_jmol.debugAlert)
    alert("jmolCodebase=" + _jmol.codebase);
}

function _jmolOnloadResetForms() {
  // must be evaluated ONLY once
  _jmol.previousOnloadHandler = window.onload;
  window.onload =
  function() {
    with (_jmol) {
      if (buttonCount+checkboxCount+menuCount+radioCount+radioGroupCount > 0) {
        var forms = document.forms;
        for (var i = forms.length; --i >= 0; )
          forms[i].reset();
      }
      if (previousOnloadHandler)
        previousOnloadHandler();
    }
  }
}

////////////////////////////////////
/////extensions for getProperty/////
////////////////////////////////////


function _jmolEvalJSON(s,key){
 s=s+""
 if(!s)return []
 if(s.charAt(0)!="{"){
	if(s.indexOf(" | ")>=0)s=s.replace(/\ \|\ /g, "\n")
	return s
 }
 var A = eval("("+s+")")
 if(!A)return
 if(key && A[key])A=A[key]
 return A
}

function _jmolEnumerateObject(A,key){
 var sout=""
 if(typeof(A) == "string" && A!="null"){
	sout+="\n"+key+"=\""+A+"\""
 }else if(!isNaN(A)||A==null){
	sout+="\n"+key+"="+(A+""==""?"null":A)
 }else if(A.length){
    sout+=key+"=[]"
    for(var i=0;i<A.length;i++){
	sout+="\n"
	if(typeof(A[i]) == "object"||typeof(A[i]) == "array"){
		sout+=_jmolEnumerateObject(A[i],key+"["+i+"]")
	}else{
		sout+=key+"["+i+"]="+(typeof(A[i]) == "string" && A[i]!="null"?"\""+A[i].replace(/\"/g,"\\\"")+"\"":A[i])
	}
    }
 }else{
    if(key != ""){
	sout+=key+"={}"
	key+="."
    }

    for(var i in A){
	sout+="\n"
	if(typeof(A[i]) == "object"||typeof(A[i]) == "array"){
		sout+=_jmolEnumerateObject(A[i],key+i)
	}else{
		sout+=key+i+"="+(typeof(A[i]) == "string" && A[i]!="null"?"\""+A[i].replace(/\"/g,"\\\"")+"\"":A[i])
	}
    }
 }
 return sout
}


function _jmolSortKey0(a,b){
 return (a[0]<b[0]?1:a[0]>b[0]?-1:0)
}

function _jmolSortMessages(A){
 if(!A || typeof(A)!="object")return []
 var B = []
 for(var i=A.length-1;i>=0;i--)for(var j=0;j<A[i].length;j++)B[B.length]=A[i][j]
 if(B.length == 0) return
 B=B.sort(_jmolSortKey0)
 return B
}

/////////additional extensions //////////


function _jmolDomScriptLoad(URL){
 //open(URL) //to debug
 _jmol.servercall=URL
 var node = document.getElementById("_jmolScriptNode")
 if (node && _jmol.browser!="msie"){
    document.getElementsByTagName("HEAD")[0].removeChild(node)
    node=null
 }
 if (node) {
   node.setAttribute("src",URL)
 } else {
   node=document.createElement("script")
   node.setAttribute("id","_jmolScriptNode")
   node.setAttribute("type","text/javascript")
   node.setAttribute("src",URL)
   document.getElementsByTagName("HEAD")[0].appendChild(node)
 }
}


function _jmolExtractPostData(url){
 S=url.split("&POST:")
 var s=""
 for(var i=1;i<S.length;i++){
	KV=S[i].split("=")
	s+="&POSTKEY"+i+"="+KV[0]
	s+="&POSTVALUE"+i+"="+KV[1]
 }
 return "&url="+escape(S[0])+s
}

function _jmolLoadModel(targetSuffix,remoteURL,array,isError,errorMessage){
 //called by server, but in client
 //overload this function to customize return
 _jmol.remoteURL=remoteURL
 isError && alert(errorMessage)
 jmolLoadInlineScript(array.join("\n"),_jmol.optionalscript,targetSuffix)
}

//////////user property/status functions/////////

function jmolGetStatus(strStatus,targetSuffix){
 return _jmolSortMessages(jmolGetPropertyAsArray("jmolStatus",strStatus,targetSuffix))
}

function jmolGetPropertyAsArray(sKey,sValue,targetSuffix) {
 return _jmolEvalJSON(jmolGetPropertyAsJSON(sKey,sValue,targetSuffix),sKey)
}

function jmolGetPropertyAsString(sKey,sValue,targetSuffix) {
 var applet = _jmolGetApplet(targetSuffix);
 sValue == undefined && (sValue="");
 return (applet ? applet.getPropertyAsString(sKey,sValue) + "" : "")
}

function jmolGetPropertyAsJSON(sKey,sValue,targetSuffix) {
 sValue == undefined && (sValue = "")
 var applet = _jmolGetApplet(targetSuffix);
 try {
  return (applet ? applet.getPropertyAsJSON(sKey,sValue) + "" : "")
 } catch(e) {
  return ""
 }
}

function jmolGetPropertyAsJavaObject(sKey,sValue,targetSuffix) {
 sValue == undefined && (sValue = "")
 var applet = _jmolGetApplet(targetSuffix);
 return (applet ? applet.getProperty(sKey,sValue) : null)
}


function jmolDecodeJSON(s) {
 return _jmolEnumerateObject(_jmolEvalJSON(s),"")
}


///////// synchronous scripting ////////

function jmolScriptWait(script, targetSuffix) {
  targetSuffix == undefined && (targetSuffix="0")
  var Ret=jmolScriptWaitAsArray(script, targetSuffix)
  var s = ""
  for(var i=Ret.length;--i>=0;)
  for(var j=0;j< Ret[i].length;j++)
	s+=Ret[i][j]+"\n"
  return s
}

function jmolScriptWaitOutput(script, targetSuffix) {
  targetSuffix == undefined && (targetSuffix="0")
  var ret = ""
  try{
   if (script) {
    _jmolCheckBrowser();
    var applet=_jmolGetApplet(targetSuffix);
    if (applet) ret += applet.scriptWaitOutput(script);
   }
  }catch(e){
  }
 return ret;
}

function jmolEvaluate(molecularMath, targetSuffix) {

  //carries out molecular math on a model

  targetSuffix == undefined && (targetSuffix="0")
  var result = "" + jmolGetPropertyAsJavaObject("evaluate", molecularMath, targetSuffix);
  var s = result.replace(/\-*\d+/,"")
  if (s == "" && !isNaN(parseInt(result)))return parseInt(result);
  var s = result.replace(/\-*\d*\.\d*/,"")
  if (s == "" && !isNaN(parseFloat(result)))return parseFloat(result);
  return result;
}

function jmolScriptEcho(script, targetSuffix) {
  // returns a newline-separated list of all echos from a script
  targetSuffix == undefined && (targetSuffix="0")
  var Ret=jmolScriptWaitAsArray(script, targetSuffix)
  var s = ""
  for(var i=Ret.length;--i>=0;)
  for(var j=Ret[i].length;--j>=0;)
        if (Ret[i][j][1] == "scriptEcho")s+=Ret[i][j][3]+"\n"
  return s.replace(/ \| /g, "\n")
}


function jmolScriptMessage(script, targetSuffix) {
  // returns a newline-separated list of all messages from a script, ending with "script completed\n"
  targetSuffix == undefined && (targetSuffix="0")
  var Ret=jmolScriptWaitAsArray(script, targetSuffix)
  var s = ""
  for(var i=Ret.length;--i>=0;)
  for(var j=Ret[i].length;--j>=0;)
        if (Ret[i][j][1] == "scriptStatus")s+=Ret[i][j][3]+"\n"
  return s.replace(/ \| /g, "\n")
}


function jmolScriptWaitAsArray(script, targetSuffix) {
 var ret = ""
 try{
  jmolGetStatus("scriptEcho,scriptMessage,scriptStatus,scriptError",targetSuffix)
  if (script) {
    _jmolCheckBrowser();
    var applet=_jmolGetApplet(targetSuffix);
    if (applet) ret += applet.scriptWait(script);
    ret = _jmolEvalJSON(ret,"jmolStatus")
    if(typeof ret == "object")
	return ret
  }
 }catch(e){
 }
  return [[ret]]
}



////////////   save/restore orientation   /////////////

function jmolSaveOrientation(id, targetSuffix) {
 targetSuffix == undefined && (targetSuffix="0")
 return _jmol["savedOrientation"+id] = jmolGetPropertyAsArray("orientationInfo","info",targetSuffix).moveTo
}

function jmolRestoreOrientation(id, targetSuffix) {
 targetSuffix == undefined && (targetSuffix="0")
 var s=_jmol["savedOrientation"+id]
 if (!s || s == "")return
 s=s.replace(/1\.0/,"0")
 return jmolScriptWait(s,targetSuffix)
}

function jmolRestoreOrientationDelayed(id, delay, targetSuffix) {
 arguments.length < 2 && (delay=1)
 targetSuffix == undefined && (targetSuffix="0")
 var s=_jmol["savedOrientation"+id]
 if (!s || s == "")return
 s=s.replace(/1\.0/,delay)
 return jmolScriptWait(s,targetSuffix)
}

////////////  add parameter /////////////
/*
 * for adding callbacks or other parameters. Use:

   jmolSetDocument(0)
   var s= jmolApplet(....)
   s = jmolAppletAddParam(s,"messageCallback", "myFunctionName")
   document.write(s)
   jmolSetDocument(document) // if you want to then write buttons and such normally

 */

function jmolAppletAddParam(appletCode,name,value){
  return (value == "" ? appletCode : appletCode.replace(/\<param/,"\n<param name='"+name+"' value='"+value+"' />\n<param"))
}

///////////////auto load Research Consortium for Structural Biology (RCSB) data ///////////

function jmolLoadAjax_STOLAF_RCSB(fileformat,pdbid,optionalscript,targetSuffix){

 _jmol.thismodel || (_jmol.thismodel = "1crn")
 _jmol.serverURL || (_jmol.serverURL="http://fusion.stolaf.edu/chemistry/jmol/getajaxjs.cfm")
 _jmol.RCSBserver || (_jmol.RCSBserver="http://www.rcsb.org")
 _jmol.defaultURL_RCSB || (_jmol.defaultURL_RCSB=_jmol.RCSBserver+"/pdb/files/1CRN.CIF")
 fileformat || (fileformat="PDB")
 pdbid || (pdbid=prompt("Enter a 4-digit PDB ID:",_jmol.thismodel))
 if(!pdbid || pdbid.length != 4)return ""
 targetSuffix || (targetSuffix="0")
 optionalscript || (optionalscript="")
 var url=_jmol.defaultURL_RCSB.replace(/1CRN/g,pdbid.toUpperCase())
 fileformat=="CIF" || (url=url.replace(/CIF/,fileformat))
 _jmol.optionalscript=optionalscript
 _jmol.thismodel=pdbid
 _jmol.thistargetsuffix=targetSuffix
 _jmol.thisurl=url
 _jmol.modelArray = []
 url=_jmol.serverURL+"?returnfunction=_jmolLoadModel&returnArray=_jmol.modelArray&id="+targetSuffix+_jmolExtractPostData(url)
 _jmolDomScriptLoad(url)
 return url
}

/////////////// St. Olaf College AJAX server -- ANY URL ///////////

function jmolLoadAjax_STOLAF_ANY(url, userid, optionalscript,targetSuffix){
 _jmol.serverURL="http://fusion.stolaf.edu/chemistry/jmol/getajaxjs.cfm"
 _jmol.thisurlANY || (_jmol.thisurlANY = "http://www.stolaf.edu/depts/chemistry/mo/struc/data/ycp3-1.mol")
 url || (url=prompt("Enter any (uncompressed file) URL:", _jmol.thisurlANY))
 userid || (userid="0")
 targetSuffix || (targetSuffix="0")
 optionalscript || (optionalscript="")
 _jmol.optionalscript=optionalscript
 _jmol.thistargetsuffix=targetSuffix
 _jmol.modelArray = []
 _jmol.thisurl = url
 url=_jmol.serverURL+"?returnfunction=_jmolLoadModel&returnArray=_jmol.modelArray&id="+targetSuffix+_jmolExtractPostData(url)
 _jmolDomScriptLoad(url)
}


/////////////// Mineralogical Society of America (MSA) data /////////

function jmolLoadAjax_MSA(key,value,optionalscript,targetSuffix){

 _jmol.thiskeyMSA || (_jmol.thiskeyMSA = "mineral")
 _jmol.thismodelMSA || (_jmol.thismodelMSA = "quartz")
 _jmol.ajaxURL_MSA || (_jmol.ajaxURL_MSA="http://rruff.geo.arizona.edu/AMS/result.php?mineral=quartz&viewing=ajaxjs")
 key || (key=prompt("Enter a field:", _jmol.thiskeyMSA))
 if(!key)return ""
 value || (value=prompt("Enter a "+key+":", _jmol.thismodelMSA))
 if(!value)return ""
 targetSuffix || (targetSuffix="0")
 optionalscript || (optionalscript="")
 optionalscript == 1 && (optionalscript='load "" {1 1 1}')
 var url=_jmol.ajaxURL_MSA.replace(/mineral/g,key).replace(/quartz/g,value)
 _jmol.optionalscript=optionalscript
 _jmol.thiskeyMSA=key
 _jmol.thismodelMSA=value
 _jmol.thistargetsuffix=targetSuffix
 _jmol.thisurl=url
 _jmol.modelArray = []
 loadModel=_jmolLoadModel
 _jmolDomScriptLoad(url)
 return url
}



function jmolLoadAjaxJS(url, userid, optionalscript,targetSuffix){
 userid || (userid="0")
 targetSuffix || (targetSuffix="0")
 optionalscript || (optionalscript="")
 _jmol.optionalscript=optionalscript
 _jmol.thismodel=userid
 _jmol.thistargetsuffix=targetSuffix
 _jmol.modelArray = []
 _jmol.thisurl = url
 url+="&returnFunction=_jmolLoadModel&returnArray=_jmol.modelArray&id="+targetSuffix
 _jmolDomScriptLoad(url)
}


//// in case Jmol library has already been loaded:

}catch(e){}

///////////////moving atoms //////////////

// HIGHLY experimental!!

function jmolSetAtomCoord(i,x,y,z,targetSuffix){
    _jmolCheckBrowser();
      var applet=_jmolGetApplet(targetSuffix);
      if (applet) applet.getProperty('jmolViewer').setAtomCoord(i,x,y,z)
}

function jmolSetAtomCoordRelative(i,x,y,z,targetSuffix){
    _jmolCheckBrowser();
      var applet=_jmolGetApplet(targetSuffix);
      if (applet) applet.getProperty('jmolViewer').setAtomCoordRelative(i,x,y,z)
}


///////////////applet fake for testing buttons/////////////


if(_jmol.useNoApplet){
	jmolApplet = function(w){
		var s="<table style='background-color:black' width="+w+"><tr height="+w+">"
		+"<td align=center valign=center style='background-color:white'>"
		+"Applet would be here"
		+"<p><textarea id=fakeApplet rows=5 cols=50></textarea>"
		+"</td></tr></table>"
		return _jmolDocumentWrite(s)
	}

	_jmolFindApplet = function(){return jmolApplet0}

	jmolApplet0 = {
	 script: function(script){document.getElementById("fakeApplet").value="\njmolScript:\n"+script}
	,scriptWait: function(script){document.getElementById("fakeApplet").value="\njmolScriptWait:\n"+script}
	,loadInline: function(data,script){document.getElementById("fakeApplet").value="\njmolLoadInline data:\n"+data+"\n\nscript:\n"+script}
	}
}


///////////////////////////////////////////

  //  This should no longer be needed, jmolResizeApplet() is better; kept for backwards compatibility
  /*
	Resizes absolutely (pixels) or by percent of window (w or h 0.5 means 50%).
	targetSuffix is optional and defaults to zero (first applet in page).
	Both w and h are optional, but needed if you want to use targetSuffix.
		h defaults to w
		w defaults to 100% of window
	If either w or h is between 0 and 1, then it is taken as percent/100.
	If either w or h is greater than 1, then it is taken as a size (pixels).
	*/
function jmolResize(w,h,targetSuffix) {
 _jmol.alerted = true;
 var percentW = (!w ? 100 : w <= 1  && w > 0 ? w * 100 : 0);
 var percentH = (!h ? percentW : h <= 1 && h > 0 ? h * 100 : 0);
 if (_jmol.browser=="msie") {
   var width=document.body.clientWidth;
   var height=document.body.clientHeight;
 } else {
   var netscapeScrollWidth=15;
   var width=window.innerWidth - netscapeScrollWidth;
   var height=window.innerHeight-netscapeScrollWidth;
 }
 var applet = _jmolGetApplet(targetSuffix);
 if(!applet)return;
 applet.style.width = (percentW ? width * percentW/100 : w)+"px";
 applet.style.height = (percentH ? height * percentH/100 : (h ? h : w))+"px";
 //title=width +  " " + height + " " + (new Date());
}

// 13 Jun 09 -- makes jmolResize() obsolete  (kept for backwards compatibility)
function jmolResizeApplet(size,targetSuffix) {
 // See _jmolGetAppletSize() for the formats accepted as size [same used by jmolApplet()]
 //  Special case: an empty value for width or height is accepted, meaning no change in that dimension.
 _jmol.alerted = true;
 var applet = _jmolGetApplet(targetSuffix);
 if(!applet)return;
 var sz = _jmolGetAppletSize(size, "px");
 sz[0] && (applet.style.width = sz[0]);
 sz[1] && (applet.style.height = sz[1]);
}

function _jmolGetAppletSize(size, units) {
	/* Accepts single number or 2-value array, each one can be one of:
	   percent (text string ending %), decimal 0 to 1 (percent/100), number, or text string (interpreted as nr.)
	   [width, height] array of strings is returned, with units added if specified.
	   Percent is relative to container div or element (which should have explicitly set size).
	*/
  var width, height;
  if ( (typeof size) == "object" && size != null ) {
    width = size[0]; height = size[1];
  } else {
    width = height = size;
  }
  return [_jmolFixDim(width, units), _jmolFixDim(height, units)];
}

function _jmolFixDim(x, units) {
  var sx = "" + x;
  return (sx.length == 0 ? (units ? "" : _jmol.allowedJmolSize[2])
	: sx.indexOf("%") == sx.length-1 ? sx
  	: (x = parseFloat(x)) <= 1 && x > 0 ? x * 100 + "%"
  	: (isNaN(x = Math.floor(x)) ? _jmol.allowedJmolSize[2]
  		: x < _jmol.allowedJmolSize[0] ? _jmol.allowedJmolSize[0]
  	    : x > _jmol.allowedJmolSize[1] ? _jmol.allowedJmolSize[1]
        : x) + (units ? units : ""));
}




