/*
 *  CHMmode.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes jsMath work with MicroSoft's HTML Help system
 *  from within .chm files (compiled help archives).
 *
 *  This file should be loaded BEFORE jsMath.js.
 *
 *  ---------------------------------------------------------------------
 *
 *  Copyright 2006-2007 by Davide P. Cervone
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
if (!jsMath.Controls) {jsMath.Controls = {}}
if (!jsMath.Controls.cookie) {jsMath.Controls.cookie = {}}

jsMath.isCHMmode = 1;

jsMath.noChangeGlobal = 1;
jsMath.noShowGlobal = 1;
jsMath.noImgFonts = 1;
jsMath.Controls.cookie.global = 'always';
jsMath.Controls.cookie.hiddenGlobal = 1;

if (window.location.protocol == "mk:") {

  /*
   *  Work around bug in hh.exe that causes it to run at 100% CPU
   *  and not exit if the page is reloaded after an IFRAME is used
   *  to load the controls file, so fake it using XMLHttpRequest.
   *  Load the data into a DIV instead of an IFRAME, and make sure
   *  that the styles are correct for it.  Change the GetPanel()
   *  call to hide the other panel and open the correct one.
   */

  jsMath.Controls.Init = function () {
    this.controlPanels = jsMath.Setup.DIV("controlPanels");
    if (!jsMath.Browser.msieButtonBug) {this.Button()}
    else {setTimeout("jsMath.Controls.Button()",500)}
  }

  jsMath.Controls.Panel = function () {
    jsMath.Translate.Cancel();
    jsMath.Setup.AddStyleSheet({
      '#jsMath_options': jsMath.styles['#jsMath_panel'],
      '#jsMath_options .disabled': jsMath.styles['#jsMath_panel .disabled'],
      '#jsMath_options .infoLink': jsMath.styles['#jsMath_panel .infoLink']
    });
    if (this.loaded) {this.panel = jsMath.Element("panel"); this.Main(); return}
    var html = jsMath.Script.xmlRequest(jsMath.root+"jsMath-controls.html");
    var body = (html.match(/<body>([\s\S]*)<\/body>/))[1];
    this.controlPanels.innerHTML = body;
    var script = (body.match(/<script>([\s\S]*?)<\/script>/))[1];
    jsMath.window.eval(script);
    jsMath.Controls.GetPanel = function (name) {
      if (this.panel) {this.panel.style.display = "none"}
      this.panel = jsMath.Element(name);
    }
    jsMath.Controls.oldClose = jsMath.Controls.Close;
    jsMath.Controls.Close = function () {this.oldClose(); this.panel = null}
    jsMath.Element("options").style.display = "none";
    jsMath.Controls.Main();
    if (!jsMath.Browser.IE7 || jsMath.Browser.quirks) {
      jsMath.window.attachEvent("onscroll",jsMath.Controls.MoveButton);
      if (jsMath.Browser.IE7) jsMath.window.attachEvent("onresize",jsMath.Controls.MoveButton);
    }
  }

}

