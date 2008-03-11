(function($) {

	$.fn.extend({
		draggable: function(options) {
			var args = Array.prototype.slice.call(arguments, 1);

			return this.each(function() {
				if (typeof options == "string") {
					var drag = $.data(this, "ui-draggable");
					drag[options].apply(drag, args);

				} else if(!$.data(this, "ui-draggable"))
					new $.ui.draggable(this, options);
			});
		}
	});

	$.ui.draggable = function(element, options) {
		//Initialize needed constants
		var self = this;

		this.element = $(element);

		$.data(element, "ui-draggable", this);
		this.element.addClass("ui-draggable");

		//Prepare the passed options
		this.options = $.extend({}, options);
		var o = this.options;
		$.extend(o, {
			helper: o.ghosting == true ? 'clone' : (o.helper || 'original'),
			handle : o.handle ? ($(o.handle, element)[0] ? $(o.handle, element) : this.element) : this.element,
			appendTo: o.appendTo || 'parent'
		});

		$(element).bind("setData.draggable", function(event, key, value){
			self.options[key] = value;
		}).bind("getData.draggable", function(event, key){
			return self.options[key];
		});

		//Initialize mouse events for interaction
		$(o.handle).mouseInteraction({
			executor: this,
			delay: o.delay,
			distance: o.distance || 0,
			dragPrevention: o.prevention ? o.prevention.toLowerCase().split(',') : ['input','textarea','button','select','option'],
			start: this.start,
			stop: this.stop,
			drag: this.drag,
			condition: function(e) { return !(e.target.className.indexOf("ui-resizable-handle") != -1 || this.disabled); }
		});

		//Position the node
		if(o.helper == 'original' && (this.element.css('position') == 'static' || this.element.css('position') == ''))
			this.element.css('position', 'relative');

	};

	$.extend($.ui.draggable.prototype, {
		plugins: {},
		ui: function(e) {
			return {
				helper: this.helper,
				position: this.position,
				absolutePosition: this.positionAbs,
				instance: this,
				options: this.options
			};
		},
		propagate: function(n,e) {
			$.ui.plugin.call(this, n, [e, this.ui()]);
			return this.element.triggerHandler(n == "drag" ? n : "drag"+n, [e, this.ui()], this.options[n]);
		},
		destroy: function() {
			this.handle.removeMouseInteraction();
			this.element
				.removeClass("ui-draggable ui-draggable-disabled")
				.removeData("ui-draggable")
				.unbind(".draggable");
		},
		enable: function() {
			this.element.removeClass("ui-draggable-disabled");
			this.disabled = false;
		},
		disable: function() {
			this.element.addClass("ui-draggable-disabled");
			this.disabled = true;
		},
		recallOffset: function(e) {

			var elementPosition = { left: this.elementOffset.left - this.offsetParentOffset.left, top: this.elementOffset.top - this.offsetParentOffset.top };
			var r = this.helper.css('position') == 'relative';

			//Generate the original position
			this.originalPosition = {
				left: (r ? parseInt(this.helper.css('left'),10) || 0 : elementPosition.left + (this.offsetParent[0] == document.body ? 0 : this.offsetParent[0].scrollLeft)),
				top: (r ? parseInt(this.helper.css('top'),10) || 0 : elementPosition.top + (this.offsetParent[0] == document.body ? 0 : this.offsetParent[0].scrollTop))
			};

			//Generate a flexible offset that will later be subtracted from e.pageX/Y
			this.offset = {left: this._pageX - this.originalPosition.left, top: this._pageY - this.originalPosition.top };

		},
		start: function(e) {

			var o = this.options;
			if($.ui.ddmanager) $.ui.ddmanager.current = this;

			//Create and append the visible helper
			this.helper = typeof o.helper == 'function' ? $(o.helper.apply(this.element[0], [e])) : (o.helper == 'clone' ? this.element.clone().appendTo((o.appendTo == 'parent' ? this.element[0].parentNode : o.appendTo)) : this.element);
			if(this.helper[0] != this.element[0]) this.helper.css('position', 'absolute');
			if(!this.helper.parents('body').length) this.helper.appendTo((o.appendTo == 'parent' ? this.element[0].parentNode : o.appendTo));

			//Find out the next positioned parent
			this.offsetParent = (function(cp) {
				while(cp) {
					if(cp.style && (/(absolute|relative|fixed)/).test($.css(cp,'position'))) return $(cp);
					cp = cp.parentNode ? cp.parentNode : null;
				}; return $("body");
			})(this.helper[0].parentNode);

			//Prepare variables for position generation
			this.elementOffset = this.element.offset();
			this.offsetParentOffset = this.offsetParent.offset();
			var elementPosition = { left: this.elementOffset.left - this.offsetParentOffset.left, top: this.elementOffset.top - this.offsetParentOffset.top };
			this._pageX = e.pageX; this._pageY = e.pageY;
			this.clickOffset = { left: e.pageX - this.elementOffset.left, top: e.pageY - this.elementOffset.top };
			var r = this.helper.css('position') == 'relative';

			//Generate the original position
			this.originalPosition = {
				left: (r ? parseInt(this.helper.css('left'),10) || 0 : elementPosition.left + (this.offsetParent[0] == document.body ? 0 : this.offsetParent[0].scrollLeft)),
				top: (r ? parseInt(this.helper.css('top'),10) || 0 : elementPosition.top + (this.offsetParent[0] == document.body ? 0 : this.offsetParent[0].scrollTop))
			};

			//If we have a fixed element, we must subtract the scroll offset again
			if(this.element.css('position') == 'fixed') {
				this.originalPosition.top -= this.offsetParent[0] == document.body ? $(document).scrollTop() : this.offsetParent[0].scrollTop;
				this.originalPosition.left -= this.offsetParent[0] == document.body ? $(document).scrollLeft() : this.offsetParent[0].scrollLeft;
			}

			//Generate a flexible offset that will later be subtracted from e.pageX/Y
			this.offset = {left: e.pageX - this.originalPosition.left, top: e.pageY - this.originalPosition.top };

			//Call plugins and callbacks
			this.propagate("start", e);

			this.helperProportions = { width: this.helper.outerWidth(), height: this.helper.outerHeight() };
			if ($.ui.ddmanager && !o.dropBehaviour) $.ui.ddmanager.prepareOffsets(this, e);

			//If we have something in cursorAt, we'll use it
			if(o.cursorAt) {
				if(o.cursorAt.top != undefined || o.cursorAt.bottom != undefined) {
					this.offset.top -= this.clickOffset.top - (o.cursorAt.top != undefined ? o.cursorAt.top : (this.helperProportions.height - o.cursorAt.bottom));
					this.clickOffset.top = (o.cursorAt.top != undefined ? o.cursorAt.top : (this.helperProportions.height - o.cursorAt.bottom));
				}
				if(o.cursorAt.left != undefined || o.cursorAt.right != undefined) {
					this.offset.left -= this.clickOffset.left - (o.cursorAt.left != undefined ? o.cursorAt.left : (this.helperProportions.width - o.cursorAt.right));
					this.clickOffset.left = (o.cursorAt.left != undefined ? o.cursorAt.left : (this.helperProportions.width - o.cursorAt.right));
				}
			}

			return false;

		},
		clear: function() {
			if($.ui.ddmanager) $.ui.ddmanager.current = null;
			this.helper = null;
		},
		stop: function(e) {

			//If we are using droppables, inform the manager about the drop
			if ($.ui.ddmanager && !this.options.dropBehaviour)
				$.ui.ddmanager.drop(this, e);

			//Call plugins and trigger callbacks
			this.propagate("stop", e);

			if(this.cancelHelperRemoval) return false;
			if(this.options.helper != 'original') this.helper.remove();
			this.clear();

			return false;
		},
		drag: function(e) {

			//Compute the helpers position
			this.position = { top: e.pageY - this.offset.top, left: e.pageX - this.offset.left };
			this.positionAbs = { left: e.pageX - this.clickOffset.left, top: e.pageY - this.clickOffset.top };

			//Call plugins and callbacks
			this.position = this.propagate("drag", e) || this.position;

			this.helper.css({ left: this.position.left+'px', top: this.position.top+'px' }); // Stick the helper to the cursor
			if($.ui.ddmanager) $.ui.ddmanager.drag(this, e);
			return false;

		}
	});

})(jQuery);
