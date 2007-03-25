/*
 *  plugins/spriteImageFonts.js
 *
 *  Part of the jsMath package for mathematics on the web.
 *
 *  This file makes jsMath use single files for the image fonts
 *  rather than individual images for each character.
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

if (!window.jsMath) {window.jsMath = {}}
if (!jsMath.styles) {jsMath.styles = []}

jsMath.styles['.typeset .img']     = 'position:relative; display:inline-block; overflow:hidden; background-repeat: no-repeat';
jsMath.styles['.typeset .img .xy'] = 'position:relative; left:0px; top:0px';

// for Mozilla
jsMath.styles['.typeset .mimg']       = 'position:relative';
jsMath.styles['.typeset .mimg .size'] = 'display:-moz-inline-box';
jsMath.styles['.typeset .mimg .wh']   = 'position:absolute; left:0px; top:0px; overflow:hidden';
jsMath.styles['.typeset .mimg .xy']   = 'position:relative; left:0px; top:0px';
// for MSIE
jsMath.styles['.typeset .mimg .h']    = 'position:relative; display:inline-block; width:0px';

/*
 *  Replace the TeXIMG function with one that uses the sprite fonts
 */
if (!jsMath.Box) {jsMath.Box = {}}
  jsMath.Box.TeXIMG = function (font,C,size) {
    var c = jsMath.TeX[font][C];
    if (c.img.reload && jsMath.Img[c.img.reload][font].loaded == 1)
      {delete c.img.reload; c.img.size = null}
    if (c.img.size != null && c.img.size == size &&
        c.img.best != null && c.img.best == jsMath.Img.best) return;
    var mustScale = (jsMath.Img.scale != 1);
    var id = jsMath.Img.best + size - 4;
    if (id < 0) {id = 0; mustScale = 1} else
    if (id >= jsMath.Img.fonts.length) {id = jsMath.Img.fonts.length-1; mustScale = 1}
    var imgFont = jsMath.Img[jsMath.Img.fonts[id]][font];
    if (!imgFont.loaded && jsMath.Browser.waitForImages) {
      // store information so several fonts can be loaded at once
      jsMath.Img.mustLoad[jsMath.Img.mustLoad.length] = [font,jsMath.Img.fonts[id]];
      imgFont.loaded = -1;
    }
    var img = imgFont[C];
    var scale = 1/jsMath.Img.w[jsMath.Img.fonts[id]];
    if (id != jsMath.Img.best + size - 4) {
      if (c.w != null) {scale = c.w/img[0]} else {
        scale *= jsMath.Img.fonts[size]/jsMath.Img.fonts[4]
              *  jsMath.Img.fonts[jsMath.Img.best]/jsMath.Img.fonts[id];
      }
    }

    // get the metrics for the character glyph
    var bScale = jsMath.Browser.imgScale;
    if (img[3] == null) {img[3] = 0}
    var w = (img[0]-img[3])*scale; var h = img[1]*scale; var d = -img[2]*scale;
    var x = img[3]-imgFont.x[C%16]; var y = img[1]-img[2]-imgFont.y[Math.floor(C/16)];
    var wh; var xy; var v;
    var ladjust = ""; var resize = ""; var vadjust; var wadjust;

    if ((mustScale || jsMath.Controls.cookie.scaleImg) && !jsMath.Browser.operaImageFonts) {
      w += 2/jsMath.em; h += 2/jsMath.em; d -= 1/jsMath.em; y += 1; x += 1; // try to adjust for rounding errors
      resize = "width:"+jsMath.HTML.Em(imgFont.wh[0]*scale*bScale)+";";
      wh = "width:"+jsMath.HTML.Em(w*bScale)+";height:"+jsMath.HTML.Em(h*bScale)+";";
      xy = "left:"+jsMath.HTML.Em(x*scale*bScale)+";top:"+jsMath.HTML.Em(y*scale*bScale)+";";
      vadjust = "vertical-align:"+jsMath.HTML.Em(d*bScale)+";";
      v = jsMath.HTML.Em(h+d);
      if (img[3]) {ladjust = "margin-left:"+jsMath.HTML.Em(-img[3]*scale*bScale)+";"}
    } else {
      if (jsMath.Browser.msieAlphaBug && jsMath.Controls.cookie.alpha) {
        resize = "height:"+(imgFont.wh[1]*jsMath.Browser.imgScale)+"px;"
               + "width:"+(imgFont.wh[0]*jsMath.Browser.imgScale)+"px;";
      }
      wh = "width:"+img[0]+"px; height:"+img[1]+"px;";
      xy = "left:"+x+"px; top:"+y+"px;";
      vadjust = "vertical-align:"+(-img[2])+"px;"; v = (img[1]-img[2])+"px";
      if (img[3]) {ladjust = "margin-left:"+(-img[3])+"px;"}
    }

    wadjust = (c.w == null || Math.abs(c.w-w) < .01)? "" : " margin-right:"+jsMath.HTML.Em(c.w-w)+';';
    if (img[2] == 0 || jsMath.Browser.valignBug) {vadjust = ""}

    // get the image
    var URL = jsMath.Img.URL(font,jsMath.Img.fonts[id],C); var IMG;
    if (jsMath.Browser.msieAlphaBug && jsMath.Controls.cookie.alpha) {
      IMG = '<span class="xy" style="'+xy+'">'
          +   '<img src="'+jsMath.blank+'" style="' + resize
          +      'filter:progid:DXImageTransform.Microsoft.AlphaImageLoader(src='
          +      "'" + URL + "', sizingMethod='scale'" + ');" />'
          + '</span>';
    } else {
      IMG = '<img src="'+URL+'" class="xy" style="'+xy+resize+'" />';
    }
    if (imgFont.loaded == -1) {IMG = ""; c.img.reload = jsMath.Img.fonts[id]}

    // get the HTML for cropping the image
    c.c = this.IMG(IMG,wh,vadjust,wadjust,ladjust,v,x,y,URL);
    c.tclass = "normal";
    c.img.bh = h+d; c.img.bd = -d;
    c.img.size = size; c.img.best = jsMath.Img.best;
  };

  /*
   *  Default uses inline-box containing an image
   */
  jsMath.Box.IMG = function (IMG,wh,vadjust,wadjust,ladjust,v,x,y,URL) {
    return '<span class="img" style="'+wh+vadjust+wadjust+ladjust+'">'
             + IMG + '</span>' + jsMath.Browser.msieImgFontBBoxFix;
  };

  /*
   *  Opera bug in inline-block alignment forces use of background image
   */
  jsMath.Box.IMG_opera = function (IMG,wh,vadjust,wadjust,ladjust,v,x,y,URL) {
    var html = '<span class="img" style="background-image: url('+URL+');'
             +    wh + vadjust + 'background-position: '+x+'px '+y+'px"></span>';
    if (wadjust || ladjust)
      {html = '<span style="'+wadjust+ladjust+'">' + html + '</span>'}
    return html;
  };

  /*
   *  Mozilla's -moz-inline-box has top aligned with baseline, so adjust
   */
  jsMath.Box.IMG_mozilla = function (IMG,wh,vadjust,wadjust,ladjust,v,x,y,URL) {
    vadjust = "vertical-align:"+v+";";
    var html = '<span class="mimg" style="'+vadjust+'">'
             +   '<span class="size" style="'+wh+'"></span>'
             +   '<span class="wh" style="'+wh+'">' + IMG + '</span>'
             + '</span>';
    if (wadjust || ladjust)
      {html = '<span style="'+wadjust+ladjust+'">' + html + '</span>'}
    return html;
  };

  /*
   *  MSIE screws up vadjust on inline-box elements, so use absolute
   *  positioning and an extra span to set the width and height
   */
  jsMath.Box.IMG_msie = function (IMG,wh,vadjust,wadjust,ladjust,v,x,y,URL) {
    var html = '<span class="mimg">'
             +   '<span class="h" style="height:'+v+'">'
             +     '<span class="wh" style="'+wh+'">' + IMG + '</span>'
             +   '</span>'
             +   '<span class="img" style="'+wh+vadjust+'"></span>'
             + '</span>';
    if (wadjust || ladjust) {
      html = jsMath.Browser.msieSpaceFix
           +    '<span style="'+wadjust+ladjust+'">' + html + '</span>';
    }
    return html;
  };

if (!jsMath.Img) {jsMath.Img = {}}

  /*
   *  Called by the exta-font definition files to add an image font
   *  into the mix (save offset data and image size)
   */
  jsMath.Img.AddFont = function (size,def) {
    if (!jsMath.Img[size]) {jsMath.Img[size] = {}};
    for (var font in def) {
      def[font].x = def[font][128]; def[font].y = def[font][129];
      def[font].wh = def[font][130];
      delete def[font][128]; delete def[font][129]; delete def[font][130];
    }
    jsMath.Add(jsMath.Img[size],def);
  };

  /*
   *  Get URL to directory for given font and size, based on the
   *  user's alpha/plain setting
   */
  jsMath.Img.URL = function (name,size) {
    if (size == null) {return this.root+name+'/font.js'}
    var type = (jsMath.Controls.cookie.alpha) ? '/alpha/': '/plain/';
    return this.root+name+type+size+'.png';
  };

  /*
   *  Laod the data for an image font
   */
  jsMath.Img.LoadFont = function (name) {
    if (jsMath.browser == 'OmniWeb' && !jsMath.Browser.hasInlineBlock) {
      jsMath.noImgFonts = 1;
      jsMath.Font.Check();
      return;
    }
    if (!this.loaded) this.Init();
    jsMath.Setup.Script(this.URL(name));
  };

  /*
   *  Setup for print mode
   */
  jsMath.Img.Init = function () {
    if ((jsMath.Controls.cookie.print || jsMath.Controls.cookie.stayhires) && !jsMath.Browser.operaImgFonts) {
      jsMath.Controls.cookie.print = jsMath.Controls.cookie.stayhires;
      this.factor *= 3;
      if (!jsMath.Controls.isLocalCookie || !jsMath.Global.isLocal) {jsMath.Controls.SetCookie(0)}
      if (jsMath.Browser.alphaPrintBug) {jsMath.Controls.cookie.alpha = 0}
    }
    this.loaded = 1;
    jsMath.Browser.ImgFontInit();
    jsMath.Img.root = jsMath.root + "fonts-sprite/";
  };


if (!jsMath.Browser) {jsMath.Browser = {}}
  /*
   *  These should be part of the regular browser
   *  test functions
   */
  jsMath.Browser.ImgFontInit = function () {
    this.msieImgFontBBoxFix = '';
    if (jsMath.browser == 'Mozilla') {
      jsMath.Box.IMG = jsMath.Box.IMG_mozilla;
    } else if (jsMath.browser == 'Opera')   {
      this.operaImageFonts = 1;
      jsMath.Box.IMG = jsMath.Box.IMG_opera;
    } else if (jsMath.browser == 'MSIE') {
      if (navigator.platform == 'MacPPC') {
        this.msieImgFontBBoxFix = '<span style="display:none">x</span>'
      } else {
        jsMath.Parser.prototype.oldTypeset = jsMath.Parser.prototype.Typeset;
        jsMath.Parser.prototype.Typeset = jsMath.Parser.prototype.msieTypeset;
        jsMath.Img.mustLoad = [];
        this.msieImageFonts = 1;
        jsMath.Controls.defaults.alpha = 0;
        if (!jsMath.Controls.userSet.alpha) {jsMath.Controls.cookie.alpha = 0}
        jsMath.Box.IMG = jsMath.Box.IMG_msie;
      }
    }
    jsMath.version += "-sp1.0";
  };


if (!jsMath.Parser) {jsMath.Parser = {}}
if (!jsMath.Parser.prototype) {jsMath.Parser.prototype = {}}
  /*
   *  Handle loading of image files needed for this equation
   *  (avoids MSIE bug where it will request the image more
   *  than once if it is used more than once before it is
   *  loaded.)
   */
  jsMath.Parser.prototype.msieTypeset = function () {
    var HTML = this.oldTypeset();
    if (jsMath.Img.mustLoad.length > 0) {
      for (var i = 0; i < jsMath.Img.mustLoad.length; i++) {
        var IMG = jsMath.Img.URL(jsMath.Img.mustLoad[i][0],jsMath.Img.mustLoad[i][1]);
        jsMath.Script.WaitForImage(IMG);
        jsMath.Img[jsMath.Img.mustLoad[i][1]][jsMath.Img.mustLoad[i][0]].loaded = 1;
      }
      jsMath.Img.mustLoad = [];
      jsMath.Translate.restart = 1;
      throw "restart";
    }
    return HTML;
  };

/*
 *  Override the control panel calls in order to
 *  disable scaling in Opera.
 */
if (!jsMath.Controls) {jsMath.Controls = {}}
  /*
   *  Add Opera calls to loading of control panel
   */
  jsMath.Controls.Panel = function () {
    jsMath.Translate.Cancel();
    if (this.loaded) {
      this.Main();
    } else {
      jsMath.Script.delayedLoad(jsMath.root+"jsMath-controls.html");
      jsMath.Script.Push(this,"OperaInit");
    }
  };

  /*
   *  Disable hi-res fonts and image scaling in Opera
   */
  jsMath.Controls.OperaMain = function (init) {
    if (!init) {this.OldMain()}
    jsMath.Element("_resolution").disabled = true;
  };

  jsMath.Controls.OperaOptions = function () {
    this.OldOptions();
    jsMath.Element("_scaleImg").disabled = true;
    jsMath.Element("_scaleImgText").className = "disabled";
  };

  jsMath.Controls.OperaInit = function () {
    if (!jsMath.Browser.operaImageFonts) return;
    this.OldMain = this.Main; this.Main = this.OperaMain;
    this.OldOptions = this.Options; this.Options = this.OperaOptions;
    this.OperaMain(1);
  };
