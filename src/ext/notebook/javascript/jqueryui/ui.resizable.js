(function($) {

  $.fn.resizable = function(options) {
    return this.each(function() {
      var args = Array.prototype.slice.call(arguments, 1);

      if (typeof options == "string") {
        var resize = $.data(this, "ui-resizable");
        resize[options].apply(resize, args);

      } else if(!$(this).is(".ui-resizable"))
        new $.ui.resizable(this, options);

    });
  };

  $.ui.resizable = function(element, options) {
    //Initialize needed constants
    var self = this;

    this.element = $(element);

    $.data(element, "ui-resizable", this);
    this.element.addClass("ui-resizable");

    //Prepare the passed options
    this.options = $.extend({
      preventDefault: true,
      transparent: false,
      minWidth: 10,
      minHeight: 10,
      aspectRatio: false,
      disableSelection: true,
      preserveCursor: true,
      autohide: false
    }, options);

		this.options._aspectRatio = !!(this.options.aspectRatio);

    $(element).bind("setData.resizable", function(event, key, value){
      self.options[key] = value;
    }).bind("getData.resizable", function(event, key){
      return self.options[key];
    });

    var o = this.options;

    //Default Theme
    var aBorder = '1px solid #DEDEDE';

    o.defaultTheme = {
      'ui-resizable': { display: 'block' },
      'ui-resizable-handle': { position: 'absolute', background: '#F5F5F5' },
      'ui-resizable-n': { cursor: 'n-resize', height: '4px', left: '0px', right: '0px', borderTop: aBorder },
      'ui-resizable-s': { cursor: 's-resize', height: '4px', left: '0px', right: '0px', borderBottom: aBorder },
      'ui-resizable-e': { cursor: 'e-resize', width: '4px', top: '0px', bottom: '0px', borderRight: aBorder },
      'ui-resizable-w': { cursor: 'w-resize', width: '4px', top: '0px', bottom: '0px', borderLeft: aBorder },
      'ui-resizable-se': { cursor: 'se-resize', width: '4px', height: '4px', borderRight: aBorder, borderBottom: aBorder },
      'ui-resizable-sw': { cursor: 'sw-resize', width: '4px', height: '4px', borderBottom: aBorder, borderLeft: aBorder },
      'ui-resizable-ne': { cursor: 'ne-resize', width: '4px', height: '4px', borderRight: aBorder, borderTop: aBorder },
      'ui-resizable-nw': { cursor: 'nw-resize', width: '4px', height: '4px', borderLeft: aBorder, borderTop: aBorder }
    };

    //Position the node
    if(!o.proxy && (this.element.css('position') == 'static' || this.element.css('position') == ''))
      this.element.css('position', 'relative');

    var nodeName = element.nodeName;

    //Wrap the element if it cannot hold child nodes
    if(nodeName.match(/textarea|input|select|button|img/i)) {

      //Create a wrapper element and set the wrapper to the new current internal element
      this.element.wrap('<div class="ui-wrapper"  style="overflow: hidden; position: relative; width: '+this.element.outerWidth()+'px; height: '+this.element.outerHeight()+';"></div>');
      var oel = this.element; element = element.parentNode; this.element = $(element);

      //Move margins to the wrapper
      this.element.css({ marginLeft: oel.css("marginLeft"), marginTop: oel.css("marginTop"),
        marginRight: oel.css("marginRight"), marginBottom: oel.css("marginBottom")
      });

      oel.css({ marginLeft: 0, marginTop: 0, marginRight: 0, marginBottom: 0});

      //Prevent Safari textarea resize
      if ($.browser.safari && o.preventDefault) oel.css('resize', 'none');

      o.proportionallyResize = oel.css('position', 'static');

			// fix handlers offset
			this._proportionallyResize();
    }

    if(!o.handles) o.handles = !$('.ui-resizable-handle', element).length ? "e,s,se" : { n: '.ui-resizable-n', e: '.ui-resizable-e', s: '.ui-resizable-s', w: '.ui-resizable-w', se: '.ui-resizable-se', sw: '.ui-resizable-sw', ne: '.ui-resizable-ne', nw: '.ui-resizable-nw' };
    if(o.handles.constructor == String) {

      if(o.handles == 'all')
        o.handles = 'n,e,s,w,se,sw,ne,nw';

      var n = o.handles.split(","); o.handles = {};

      o.zIndex = o.zIndex || 1000;

      var insertions = {
				handle: 'overflow:hidden; position:absolute;',
        n: 'top: 0pt; width:100%;',
        e: 'right: 0pt; height:100%;',
        s: 'bottom: 0pt; width:100%;',
        w: 'left: 0pt; height:100%;',
        se: 'bottom: 0pt; right: 0px;',
        sw: 'bottom: 0pt; left: 0px;',
        ne: 'top: 0pt; right: 0px;',
        nw: 'top: 0pt; left: 0px;'
      };

      for(var i = 0; i < n.length; i++) {
        var d = jQuery.trim(n[i]), t = o.defaultTheme, hname = 'ui-resizable-'+d;

        var rcss = $.extend(t[hname], t['ui-resizable-handle']),
            axis = $(["<div class=\"",hname," ui-resizable-handle\" style=\"",insertions[d], insertions.handle,"\"></div>"].join("")).css(/sw|se|ne|nw/.test(d) ? { zIndex: ++o.zIndex } : {});

        o.handles[d] = '.ui-resizable-'+d;

        this.element.append(
          //Theme detection, if not loaded, load o.defaultTheme
          axis.css( !$.ui.css(hname) ? rcss : {} )
        );
      }
    }

    this._renderAxis = function(target) {
      target = target || this.element;

      for(var i in o.handles) {
        if(o.handles[i].constructor == String)
          o.handles[i] = $(o.handles[i], element).show();

        if (o.transparent)
          o.handles[i].css({opacity:0});

        //Apply pad to wrapper element, needed to fix axis position (textarea, inputs, scrolls)
        if (this.element.is('.ui-wrapper') &&
          nodeName.match(/textarea|input|select|button/i)) {

          var axis = $(o.handles[i], element), padWrapper = 0;

          //Checking the correct pad and border
          padWrapper = /sw|ne|nw|se|n|s/.test(i) ? axis.outerHeight() : axis.outerWidth();

          //The padding type i have to apply...
          var padPos = [ 'padding',
            /ne|nw|n/.test(i) ? 'Top' :
            /se|sw|s/.test(i) ? 'Bottom' :
            /^e$/.test(i) ? 'Right' : 'Left' ].join("");

          if (!o.transparent)
            target.css(padPos, padWrapper);
        }
        if(!$(o.handles[i]).length) continue;
      }
    };

    this._renderAxis(this.element);
    var handlers = $('.ui-resizable-handle', self.element);

    if (o.disableSelection)
      handlers.each(function(i, e) { $.ui.disableSelection(e); });

    //Matching axis name
    handlers.mouseover(function() {
      if (!o.resizing) {
        if (this.className)
          var axis = this.className.match(/ui-resizable-(se|sw|ne|nw|n|e|s|w)/i);
        //Axis, default = se
        o.axis = axis && axis[1] ? axis[1] : 'se';
      }
    });

    //If we want to auto hide the elements
    if (o.autohide) {
      var tLoaded = $.ui.css('ui-resizable-s') || $.ui.css('ui-resizable-e');
      if (!tLoaded) handlers.hide();

      $(self.element).addClass("ui-resizable-autohide").hover(function(){
        if (!tLoaded) handlers.show();
        $(this).removeClass("ui-resizable-autohide");
      }, function(){
        if (!o.resizing) {
          if (!tLoaded) handlers.hide();
          $(this).addClass("ui-resizable-autohide");
        }
      });
    }

    //Initialize mouse events for interaction
    this.element.mouseInteraction({
      executor: this,
      delay: 0,
      distance: 0,
      dragPrevention: ['input','textarea','button','select','option'],
      start: this.start,
      stop: this.stop,
      drag: this.drag,
      condition: function(e) {
        if(this.disabled) return false;
        for(var i in this.options.handles) {
          if($(this.options.handles[i])[0] == e.target) return true;
        }
        return false;
      }
    });
  };

  $.extend($.ui.resizable.prototype, {
    plugins: {},
    ui: function() {
      return {
        instance: this,
        axis: this.options.axis,
        options: this.options
      };
    },
    _proportionallyResize: function() {
			var o = this.options;

			if (!o.proportionallyResize)
				return;

			var prel = o.proportionallyResize;

      var b = [ prel.css('borderTopWidth'), prel.css('borderRightWidth'), prel.css('borderBottomWidth'), prel.css('borderLeftWidth') ];
      var p = [ prel.css('paddingTop'), prel.css('paddingRight'), prel.css('paddingBottom'), prel.css('paddingLeft') ];

      o.borderDif = o.borderDif || $.map(b, function(v, i) {
        var border = parseInt(v,10)||0, padding = parseInt(p[i],10)||0;
        return border + padding;
      });

      prel.css({
        display: 'block', //Needed to fix height autoincrement
        height: (this.element.height() - o.borderDif[0] - o.borderDif[2]) + "px",
        width: (this.element.width() - o.borderDif[1] - o.borderDif[3]) + "px"
      });
    },
    _renderProxy: function() {
      var el = this.element, o = this.options;
      this.offset = el.offset();

      if(o.proxy) {
        this.helper = this.helper || $('<div style="overflow:hidden;"></div>');

				// fix ie6 offset
				var ie6offset = ($.browser.msie && $.browser.version < 7 ? 3 : 0);

        this.helper.addClass(o.proxy).css({
          width: el.outerWidth(),
          height: el.outerHeight(),
          position: 'absolute',
          left: this.offset.left - ie6offset +'px',
          top: this.offset.top - ie6offset +'px',
          zIndex: ++o.zIndex
        });

        this.helper.appendTo("body");

        if (o.disableSelection)
          $.ui.disableSelection(this.helper.get(0));

      } else {
        this.helper = el;
      }
    },
    propagate: function(n,e) {
      $.ui.plugin.call(this, n, [e, this.ui()]);
      this.element.triggerHandler(n == "resize" ? n : "resize"+n, [e, this.ui()], this.options[n]);
    },
    destroy: function() {
      this.element
        .removeClass("ui-resizable ui-resizable-disabled")
        .removeMouseInteraction()
        .removeData("ui-resizable")
        .unbind(".resizable");
    },
    enable: function() {
      this.element.removeClass("ui-resizable-disabled");
      this.disabled = false;
    },
    disable: function() {
      this.element.addClass("ui-resizable-disabled");
      this.disabled = true;
    },
    start: function(e) {
      var o = this.options, iniPos = this.element.position(), el = this.element;
      o.resizing = true;
			o.documentScroll = { top: $(document).scrollTop(),	left: $(document).scrollLeft() };

			// buf fix #1749
      if (el.is('.ui-draggable') || (/absolute/).test(el.css('position'))) {
				// sOffset decides if document scrollOffset will be added to the top/left of the resizable element
				var sOffset = $.browser.msie && !o.containment && (/absolute/).test(el.css('position')) && !(/relative/).test(el.parent().css('position'));
				var dscrollt = sOffset ? o.documentScroll.top : 0, dscrolll = sOffset ? o.documentScroll.left : 0;
	  		el.css({ position: 'absolute', top: (iniPos.top + dscrollt),	left: (iniPos.left + dscrolll)	});
	  	}

      //Opera fixing relative position
      if (/relative/.test(el.css('position')) && $.browser.opera)
        el.css({ position: 'relative', top: 'auto', left: 'auto' });

      this._renderProxy();

			var curleft = parseInt(this.helper.css('left'),10) || 0, curtop = parseInt(this.helper.css('top'),10) || 0;

      //Store needed variables
      $.extend(o, {
        currentSize: { width: el.outerWidth(), height: el.outerHeight() },
        currentSizeDiff: { width: el.outerWidth() - el.width(), height: el.outerHeight() - el.height() },
        startMousePosition: { left: e.pageX, top: e.pageY },
        startPosition: { left: curleft, top: curtop },
        currentPosition: { left: curleft,top: curtop }
      });

			//Aspect Ratio
			var iswlt = o.currentSize.width < o.currentSize.height;
 			o.aspectRatio = (typeof o.aspectRatio == 'number') ? o.aspectRatio : Math.pow(o.currentSize.width / o.currentSize.height, iswlt ? 1 : -1);
			o.aspectRatioTarget = iswlt ? "width" : "height";

      if (o.preserveCursor)
        $('body').css('cursor', o.axis + '-resize');

      if (o.containment) {
        var oc = o.containment,
           ce = (oc instanceof jQuery) ? oc.get(0) :
              (/parent/.test(oc)) ? el.parent().get(0) : null;
        if (ce) {

          var scroll = function(e, a) {
            var scroll = /top/.test(a||"top") ? 'scrollTop' : 'scrollLeft', has = false;
            if (e[scroll] > 0) return true; e[scroll] = 1;
            has = e[scroll] > 0 ? true : false; e[scroll] = 0;
            return has;
          };

          var co = $(ce).offset(), ch = $(ce).innerHeight(), cw = $(ce).innerWidth();
          o.cdata = { e: ce, l: co.left, t: co.top, w: (scroll(ce, "left") ? ce.scrollWidth : cw ), h: (scroll(ce) ? ce.scrollHeight : ch) };
        }
        if (/document/.test(oc) || oc == document) o.cdata = { e: document, l: 0, t: 0, w: $(document).width(), h: $(document).height() };
      }
      this.propagate("start", e);
      return false;

    },
    stop: function(e) {
      this.options.resizing = false;
      var o = this.options;

      if(o.proxy) {
        var style = {
          width: (this.helper.width() - o.currentSizeDiff.width) + "px",
          height: (this.helper.height() - o.currentSizeDiff.height) + "px",
          top: ((parseInt(this.element.css('top'),10) || 0) + ((parseInt(this.helper.css('top'),10) - this.offset.top)||0)),
          left: ((parseInt(this.element.css('left'),10) || 0) + ((parseInt(this.helper.css('left'),10) - this.offset.left)||0))
        };
       	this.element.css(style);
        if (o.proxy) this._proportionallyResize();
        this.helper.remove();
      }

      if (o.preserveCursor)
        $('body').css('cursor', 'auto');

      this.propagate("stop", e);
      return false;

    },
    drag: function(e) {
      //Increase performance, avoid regex
      var el = this.helper, o = this.options, props = {},
						self = this, pRatio = o._aspectRatio || e.shiftKey;

      var change = function(a,b) {
        var isth = (a=="top"||a=="height"), ishw = (a=="width"||a=="height"),
							defAxis = (o.axis=="se"||o.axis=="s"||o.axis=="e");

        var mod = (e[isth ? 'pageY' : 'pageX'] - o.startMousePosition[isth ? 'top' : 'left']) * (b ? -1 : 1);
        var val = o[ishw ? 'currentSize' : 'startPosition'][a] - mod - (!o.proxy && defAxis ? o.currentSizeDiff.width : 0);

				//Preserve ratio
        if (pRatio) {
					var v = val * Math.pow(o.aspectRatio, (isth ? -1 : 1) * (o.aspectRatioTarget == 'height' ? 1 : -1)), locked = false;

					if (isth && v >= o.maxWidth || !isth && v >= o.maxHeight)	locked = true;
					if (isth && v <= o.minWidth || !isth && v <= o.minHeight)	locked = true;

					if (ishw && !locked)	el.css(isth ? "width" : "height", v);

					if (a == "top" && (o.axis == "ne" || o.axis == "nw")) {
						//el.css('top', o.startPosition['top'] - (el.outerHeight() - o.currentSize.height) );
						/*TODO*/ return;
					};
				}
				el.css(a, val);
      };

			var a = o.axis, tminval = 0, tmaxval;

      //Change the height
      if(a=="n"||a=="ne"||a=="nw") change("height");
      if(a=="s"||a=="se"||a=="sw") change("height", 1);

      //Measure the new height and correct against min/maxHeight
			var curheight = parseInt(el.css('height'),10)||0;
      if(o.minHeight && curheight <= o.minHeight) el.css('height', o.minHeight);
      if(o.maxHeight && curheight >= o.maxHeight) el.css('height', o.maxHeight);

      //Change the top position when picking a handle at north
      if(a=="n"||a=="ne"||a=="nw") change("top", 1);

      //Measure the new top position and correct against min/maxHeight
			var curtop = parseInt(el.css('top'),10)||0;

			tminval = (o.startPosition.top + (o.currentSize.height - o.minHeight));
			tmaxval = (o.startPosition.top + (o.currentSize.height - o.maxHeight));
      if(o.minHeight && curtop >= tminval) el.css('top', tminval);
      if(o.maxHeight && curtop <= tmaxval) el.css('top', tmaxval);

      //Change the width
      if(a=="e"||a=="se"||a=="ne") change("width", 1);
      if(a=="sw"||a=="w"||a=="nw") change("width");

      //Measure the new width and correct against min/maxWidth
			var curwidth = parseInt(el.css('width'),10)||0;
      if(o.minWidth && curwidth <= o.minWidth) el.css('width', o.minWidth);
      if(o.maxWidth && curwidth >= o.maxWidth) el.css('width', o.maxWidth);

      //Change the left position when picking a handle at west
      if(a=="sw"||a=="w"||a=="nw") change("left", 1);

      //Measure the new left position and correct against min/maxWidth
			var curleft = parseInt(el.css('left'),10)||0;

			tminval = (o.startPosition.left + (o.currentSize.width - o.minWidth));
			tmaxval = (o.startPosition.left + (o.currentSize.width - o.maxWidth));
      if(o.minWidth && curleft >= tminval) el.css('left', tminval);
      if(o.maxWidth && curleft <= tmaxval) el.css('left', tmaxval);

      if (o.containment && o.cdata.e) {
        if (curleft < 0) {
          el.css('left', 0);
          el.css('width', curwidth + curleft);
        }
        if (curtop < 0) {
          el.css('top', 0);
          el.css('height', curheight + curtop);
        }
        if (curwidth + o.currentSizeDiff.width + curleft >= o.cdata.w)
					el.css('width', o.cdata.w - o.currentSizeDiff.width - (curleft < 0 ? 0 : curleft));
        if (curheight + o.currentSizeDiff.height + curtop >= o.cdata.h)
					el.css('height', o.cdata.h - o.currentSizeDiff.height - (curtop < 0 ? 0 : curtop));
      }
			o.currentPosition = { left: curleft, top: curtop };
      if (!o.proxy) this._proportionallyResize();
      this.propagate("resize", e);
      return false;
    }
  });

})(jQuery);
