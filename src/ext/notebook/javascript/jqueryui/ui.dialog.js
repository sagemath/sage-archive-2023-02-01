;(function($) {

	//If the UI scope is not available, add it
	$.ui = $.ui || {};

	$.fn.extend({
		dialog: function(options, data) {
			var args = Array.prototype.slice.call(arguments, 1);

			return this.each(function() {
				if (typeof options == "string") {
					var dialog = $.data(this, "ui-dialog") ||
						$.data($(this).parents(".ui-dialog:first").find(".ui-dialog-content")[0], "ui-dialog");
					dialog[options].apply(dialog, args);

				// INIT with optional options
				} else if (!$(this).is(".ui-dialog-content"))
					new $.ui.dialog(this, options);
			});
		}
	});

	$.ui.dialog = function(el, options) {

		this.options = options = $.extend({},
			$.ui.dialog.defaults,
			options && options.modal ? {resizable: false} : {},
			options);
		this.element = el;
		var self = this; //Do bindings

		$.data(this.element, "ui-dialog", this);

		$(el).bind("setData.dialog", function(event, key, value){
			options[key] = value;
		}).bind("getData.dialog", function(event, key){
			return options[key];
		});

		var uiDialogContent = $(el).addClass('ui-dialog-content');

		if (!uiDialogContent.parent().length) {
			uiDialogContent.appendTo('body');
		}
		uiDialogContent
			.wrap(document.createElement('div'))
			.wrap(document.createElement('div'));
		var uiDialogContainer = uiDialogContent.parent().addClass('ui-dialog-container').css({position: 'relative'});
		var uiDialog = this.uiDialog = uiDialogContainer.parent().hide()
			.addClass('ui-dialog')
			.css({position: 'absolute', width: options.width, height: options.height, overflow: 'hidden'});

		var classNames = uiDialogContent.attr('className').split(' ');

		// Add content classes to dialog, to inherit theme at top level of element
		$.each(classNames, function(i, className) {
			if (className != 'ui-dialog-content')
				uiDialog.addClass(className);
		});

		if (options.resizable && $.fn.resizable) {
			uiDialog.append('<div class="ui-resizable-n ui-resizable-handle"></div>')
				.append('<div class="ui-resizable-s ui-resizable-handle"></div>')
				.append('<div class="ui-resizable-e ui-resizable-handle"></div>')
				.append('<div class="ui-resizable-w ui-resizable-handle"></div>')
				.append('<div class="ui-resizable-ne ui-resizable-handle"></div>')
				.append('<div class="ui-resizable-se ui-resizable-handle"></div>')
				.append('<div class="ui-resizable-sw ui-resizable-handle"></div>')
				.append('<div class="ui-resizable-nw ui-resizable-handle"></div>');
			uiDialog.resizable({ maxWidth: options.maxWidth, maxHeight: options.maxHeight, minWidth: options.minWidth, minHeight: options.minHeight });
		}

		uiDialogContainer.prepend('<div class="ui-dialog-titlebar"></div>');
		var uiDialogTitlebar = $('.ui-dialog-titlebar', uiDialogContainer);
		var title = (options.title) ? options.title : (uiDialogContent.attr('title')) ? uiDialogContent.attr('title') : '';
		uiDialogTitlebar.append('<span class="ui-dialog-title">' + title + '</span>');
		uiDialogTitlebar.append('<a href="#" class="ui-dialog-titlebar-close"><span>X</span></a>');
		this.uiDialogTitlebarClose = $('.ui-dialog-titlebar-close', uiDialogTitlebar)
			.hover(function() { $(this).addClass('ui-dialog-titlebar-close-hover'); },
			       function() { $(this).removeClass('ui-dialog-titlebar-close-hover'); })
			.mousedown(function(ev) {
				ev.stopPropagation();
			})
			.click(function() {
				self.close();
				return false;
			})
			.keydown(function(ev) {
				var ESC = 27;
				ev.keyCode && ev.keyCode == ESC && self.close();
			});

		var l = 0;
		$.each(options.buttons, function() { l = 1; return false; });
		if (l == 1) {
			uiDialog.append('<div class="ui-dialog-buttonpane"></div>');
			var uiDialogButtonPane = $('.ui-dialog-buttonpane', uiDialog);
			$.each(options.buttons, function(name, value) {
				var btn = $(document.createElement('button')).text(name).click(value);
				uiDialogButtonPane.append(btn);
			});
		}

		if (options.draggable && $.fn.draggable) {
			uiDialog.draggable({
				handle: '.ui-dialog-titlebar',
				start: function() {
					self.activate();
				}
			});
		}
		uiDialog.mousedown(function() {
			self.activate();
		});
		uiDialogTitlebar.click(function() {
			self.activate();
		});

		// TODO: determine if this is necessary for modal dialogs
		options.bgiframe && $.fn.bgiframe && uiDialog.bgiframe();

		this.open = function() {
			options.modal && overlay.show(self, options.overlay);
			uiDialog.appendTo('body');
			var wnd = $(window), doc = $(document), top = doc.scrollTop(), left = doc.scrollLeft();
			if (options.position.constructor == Array) {
				// [x, y]
				top += options.position[1];
				left += options.position[0];
			} else {
				switch (options.position) {
					case 'center':
						top += (wnd.height() / 2) - (uiDialog.height() / 2);
						left += (wnd.width() / 2) - (uiDialog.width() / 2);
						break;
					case 'top':
						top += 0;
						left += (wnd.width() / 2) - (uiDialog.width() / 2);
						break;
					case 'right':
						top += (wnd.height() / 2) - (uiDialog.height() / 2);
						left += (wnd.width()) - (uiDialog.width());
						break;
					case 'bottom':
						top += (wnd.height()) - (uiDialog.height());
						left += (wnd.width() / 2) - (uiDialog.width() / 2);
						break;
					case 'left':
						top += (wnd.height() / 2) - (uiDialog.height() / 2);
						left += 0;
						break;
					default:
						//center
						top += (wnd.height() / 2) - (uiDialog.height() / 2);
						left += (wnd.width() / 2) - (uiDialog.width() / 2);
				}
			}
			top = top < doc.scrollTop() ? doc.scrollTop() : top;
			uiDialog.css({top: top, left: left});
			uiDialog.show();
			self.activate();

			// CALLBACK: open
			var openEV = null;
			var openUI = {
				options: options
			};
			this.uiDialogTitlebarClose.focus();
			$(this.element).triggerHandler("dialogopen", [openEV, openUI], options.open);
		};

		this.activate = function() {
			var maxZ = 0;
			$('.ui-dialog:visible').each(function() {
				maxZ = Math.max(maxZ, parseInt($(this).css("z-index"),10));
			});
			overlay.$el && overlay.$el.css('z-index', ++maxZ);
			uiDialog.css("z-index", ++maxZ);
		};

		this.close = function() {
			options.modal && overlay.hide();
			uiDialog.hide();

			// CALLBACK: close
			var closeEV = null;
			var closeUI = {
				options: options
			};
			$(this.element).triggerHandler("dialogclose", [closeEV, closeUI], options.close);
		};

		if (options.autoOpen)
			this.open();
	};

	$.extend($.ui.dialog, {
		defaults: {
			autoOpen: true,
			bgiframe: false,
			buttons: [],
			draggable: true,
			height: 200,
			minHeight: 100,
			minWidth: 150,
			modal: false,
			overlay: {},
			position: 'center',
			resizable: true,
			width: 300
		}
	});

	// This is a port of relevant pieces of Mike Alsup's blockUI plugin (http://www.malsup.com/jquery/block/)
	// duplicated here for minimal overlay functionality and no dependency on a non-UI plugin
	var overlay = {
		$el: null,
		events: $.map('focus,mousedown,mouseup,keydown,keypress,click'.split(','),
			function(e) { return e + '.ui-dialog-overlay'; }).join(' '),

		show: function(dialog, css) {
			if (this.$el) return;

			this.dialog = dialog;
			this.selects = this.ie6 && $('select:visible').css('visibility', 'hidden');
			var width = this.width();
			var height = this.height();
			this.$el = $('<div/>').appendTo(document.body)
				.addClass('ui-dialog-overlay').css($.extend({
					borderWidth: 0, margin: 0, padding: 0,
					position: 'absolute', top: 0, left: 0,
					width: width,
					height: height
				}, css));

			// prevent use of anchors and inputs
			$('a, :input').bind(this.events, function() {
				if ($(this).parents('.ui-dialog').length == 0) {
					dialog.uiDialogTitlebarClose.focus();
					return false;
				}
			});

			// allow closing by pressing the escape key
			$(document).bind('keydown.ui-dialog-overlay', function(e) {
				var ESC = 27;
				e.keyCode && e.keyCode == ESC && dialog.close();
			});

			// handle window resizing
			$overlay = this.$el;
			function resize() {
				// If the dialog is draggable and the user drags it past the
				// right edge of the window, the document becomes wider so we
				// need to stretch the overlay.  If the user then drags the
				// dialog back to the left, the document will become narrower,
				// so we need to shrink the overlay to the appropriate size.
				// This is handled by resetting the overlay to its original
				// size before setting it to the full document size.
				$overlay.css({
					width: width,
					height: height
				}).css({
					width: overlay.width(),
					height: overlay.height()
				});
			};
			$(window).bind('resize.ui-dialog-overlay', resize);
			dialog.uiDialog.is('.ui-draggable') && dialog.uiDialog.data('stop.draggable', resize);
			dialog.uiDialog.is('.ui-resizable') && dialog.uiDialog.data('stop.resizable', resize);
		},

		hide: function() {
			$('a, :input').add([document, window]).unbind('.ui-dialog-overlay');
			this.ie6 && this.selects.css('visibility', 'visible');
			this.$el = null;
			$('.ui-dialog-overlay').remove();
		},

		height: function() {
			var height;
			if (this.ie6
				// body is smaller than window
				&& ($(document.body).height() < $(window).height())
				// dialog is above the fold
				&& !(document.documentElement.scrollTop
					|| (this.dialog.uiDialog.offset().top
						+ this.dialog.uiDialog.height())
						> $(window).height())) {
				height = $(window).height();
			} else {
				height = $(document).height();
			}
			return height + 'px';
		},

		width: function() {
			var width;
			if (this.ie6
				// body is smaller than window
				&& ($(document.body).width() < $(window).width())
				// dialog is off to the right
				&& !(document.documentElement.scrollLeft
					|| (this.dialog.uiDialog.offset().left
						+ this.dialog.uiDialog.width())
						> $(window).width())) {
				width = $(window).width();
			} else {
				width = $(document).width();
			}
			return width + 'px';
		},

		// IE 6 compatibility
		ie6: $.browser.msie && $.browser.version < 7,
		selects: null
	};

})(jQuery);
