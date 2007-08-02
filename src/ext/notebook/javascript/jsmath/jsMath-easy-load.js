/*
 *  jsMath-easy-load.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file is used to load jsMath with one easy <SCRIPT>
 *  command in your HTML file.  It is called by the files
 *  in the jsMath/easy/ directory.  It expects that the jsMath.Easy
 *  array has been initialized before it is called.
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

if (!window.jsMath) {window.jsMath = {}}
if (!jsMath.Easy) {jsMath.Easy = {}}
if (!jsMath.tex2math) {jsMath.tex2math = {}}

jsMath.tex2math.doubleDollarsAreInLine = jsMath.Easy.doubleDollarsAreInLine;
jsMath.tex2math.allowDisableTag = jsMath.Easy.allowDisableTag;

if (jsMath.Easy.scale) {
  if (!jsMath.Controls) {jsMath.Controls = {}}
  if (!jsMath.Controls.cookie) {jsMath.Controls.cookie = {}}
  jsMath.Controls.cookie.scale = jsMath.Easy.scale;
}
if (!jsMath.Easy.allowDoubleClicks) {
  if (!jsMath.Click) {jsMath.Click = {}}
  jsMath.Click.CheckDblClick = function () {};
}
if (!jsMath.Easy.showFontWarnings) {
  if (!jsMath.Font) {jsMath.Font = {}}
  jsMath.Font.Message = function () {};
}

if (!jsMath.Easy.allowGlobal) {
  if (!jsMath.Controls) {jsMath.Controls = {}}
  if (!jsMath.Controls.cookie) {jsMath.Controls.cookie = {}}
  jsMath.Controls.cookie.global = 'never';
  jsMath.noGoGlobal = 1;
  jsMath.noChangeGlobal = 1;
  jsMath.noShowGlobal = 1;
}

if (jsMath.Easy.noImageFonts) {
  jsMath.noImgFonts = 1;
  if (!jsMath.Font) {jsMath.Font = {}}
  jsMath.Font.extra_message =
    'Extra TeX fonts not found: <b><span id="jsMath_ExtraFonts"></span></b><br/>'
      + 'Using unicode fonts instead.  This may be slow and might not print well.<br/>\n'
      + 'Use the jsMath control panel to get additional information.';
}

if (jsMath.Easy.processSingleDollars ||
    jsMath.Easy.processDoubleDollars ||
    jsMath.Easy.processSlashParens ||
    jsMath.Easy.processSlashBrackets ||
    jsMath.Easy.fixEscapedDollars) {

  jsMath.Easy.findCustomSettings = {
    processSingleDollars: jsMath.Easy.processSingleDollars,
    processDoubleDollars: jsMath.Easy.processDoubleDollars,
    processSlashParens:   jsMath.Easy.processSlashParens,
    processSlashBrackets: jsMath.Easy.processSlashBrackets,
    fixEscapedDollars:    jsMath.Easy.fixEscapedDollars,
    custom: 0
  }
}

if (!jsMath.Autoload) {jsMath.Autoload = {}}
jsMath.Autoload.root = jsMath.Easy.root+'/';
if (jsMath.Easy.autoload) {
  jsMath.Autoload.findTeXstrings = 0;
  jsMath.Autoload.findLaTeXstrings = 0;
  jsMath.Autoload.findCustomStrings = jsMath.Easy.customDelimiters;
  jsMath.Autoload.findCustomSettings = jsMath.Easy.findCustomSettings;
  jsMath.Autoload.loadFiles = jsMath.Easy.loadFiles;
  jsMath.Autoload.loadFonts = jsMath.Easy.loadFonts;
  jsMath.Autoload.root = jsMath.Easy.root + '/';

  if (!document.body) {jsMath.Easy.autoloadCheck = 1}
  document.write('<script src="'+jsMath.Autoload.root+'plugins/autoload.js"></script>');

} else {
  jsMath.Easy.tex2math =
     (jsMath.Easy.processSingleDollars ||
      jsMath.Easy.processDoubleDollars ||
      jsMath.Easy.processSlashParens ||
      jsMath.Easy.processSlashBrackets ||
      jsMath.Easy.fixEscapedDollars ||
      jsMath.Easy.customDelimiters);

  if (!jsMath.Setup) {jsMath.Setup = {}}
  if (!jsMath.Setup.UserEvent) {jsMath.Setup.UserEvent = {}}
  jsMath.Setup.UserEvent.onload = function () {
    if (jsMath.Easy.tex2math) jsMath.Setup.Script("plugins/tex2math.js");
    var i;
    for (i = 0; i < jsMath.Easy.loadFiles.length; i++)
      jsMath.Setup.Script(jsMath.Easy.loadFiles[i]);
    for (i = 0; i < jsMath.Easy.loadFonts.length; i++)
      jsMath.Font.Load(jsMath.Easy.loadFonts[i]);
  }
  document.write('<script src="'+jsMath.Easy.root+'/jsMath.js"></script>'+"\n");
}

jsMath.Easy.onload = function () {
  if (jsMath.Easy.autoloadCheck) jsMath.Autoload.Check();
  if (jsMath.Easy.tex2math) {
    jsMath.Synchronize(function () {
      if (jsMath.Easy.findCustomSettings)
        jsMath.tex2math.Convert(document,jsMath.Easy.findCustomSettings);
      if (jsMath.Easy.customDelimiters) {
        var s = jsMath.Easy.customDelimiters;
        jsMath.tex2math.CustomSearch(s[0],s[1],s[2],s[3]);
        jsMath.tex2math.ConvertCustom();
      }
    });
  }
  (jsMath[jsMath.Easy.method])();
}

if (window.addEventListener) {window.addEventListener("load",jsMath.Easy.onload,false)}
else if (window.attachEvent) {window.attachEvent("onload",jsMath.Easy.onload)}
else {window.onload = jsMath.Easy.onload}
