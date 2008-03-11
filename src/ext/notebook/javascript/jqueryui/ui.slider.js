(function($) {

	$.fn.extend({
		slider: function(options) {
			var args = Array.prototype.slice.call(arguments, 1);

			if ( options == "value" )
				return $.data(this[0], "ui-slider").value(arguments[1]);

			return this.each(function() {
				if (typeof options == "string") {
					var slider = $.data(this, "ui-slider");
					slider[options].apply(slider, args);

				} else if(!$.data(this, "ui-slider"))
					new $.ui.slider(this, options);
			});
		}
	});

	$.ui.slider = function(element, options) {

		//Initialize needed constants
		var self = this;
		this.element = $(element);
		$.data(element, "ui-slider", this);
		this.element.addClass("ui-slider");

		//Prepare the passed options
		this.options = $.extend({}, options);
		var o = this.options;
		$.extend(o, {
			axis: o.axis || (element.offsetWidth < element.offsetHeight ? 'vertical' : 'horizontal'),
			maxValue: !isNaN(parseInt(o.maxValue,10)) ? parseInt(o.maxValue,10) :  100,
			minValue: parseInt(o.minValue,10) || 0,
			startValue: parseInt(o.startValue,10) || 'none'
		});

		//Prepare the real maxValue
		o.realMaxValue = o.maxValue - o.minValue;

		//Calculate stepping based on steps
		o.stepping = parseInt(o.stepping,10) || (o.steps ? o.realMaxValue/o.steps : 0);

		$(element).bind("setData.slider", function(event, key, value){
			self.options[key] = value;
		}).bind("getData.slider", function(event, key){
			return self.options[key];
		});

		//Initialize mouse and key events for interaction
		this.handle = o.handle ? $(o.handle, element) : $('> *', element);
		$(this.handle)
			.mouseInteraction({
				executor: this,
				delay: o.delay,
				distance: o.distance || 0,
				dragPrevention: o.prevention ? o.prevention.toLowerCase().split(',') : ['input','textarea','button','select','option'],
				start: this.start,
				stop: this.stop,
				drag: this.drag,
				condition: function(e, handle) {
					if(!this.disabled) {
						if(this.currentHandle) this.blur(this.currentHandle);
						this.focus(handle,1);
						return !this.disabled;
					}
				}
			})
			.wrap('<a href="javascript:void(0)"></a>')
			.parent()
				.bind('focus', function(e) { self.focus(this.firstChild); })
				.bind('blur', function(e) { self.blur(this.firstChild); })
				.bind('keydown', function(e) {
					if(/(37|39)/.test(e.keyCode))
						self.moveTo((e.keyCode == 37 ? '-' : '+')+'='+(self.options.stepping ? self.options.stepping : (self.options.realMaxValue / self.size)*5),this.firstChild);
				})
		;

		//Position the node
		if(o.helper == 'original' && (this.element.css('position') == 'static' || this.element.css('position') == '')) this.element.css('position', 'relative');

		//Prepare dynamic properties for later use
		if(o.axis == 'horizontal') {
			this.size = this.element.outerWidth();
			this.properties = ['left', 'width'];
		} else {
			this.size = this.element.outerHeight();
			this.properties = ['top', 'height'];
		}

		//Bind the click to the slider itself
		this.element.bind('click', function(e) { self.click.apply(self, [e]); });

		//Move the first handle to the startValue
		if(!isNaN(o.startValue)) this.moveTo(o.startValue, 0);

		//If we only have one handle, set the previous handle to this one to allow clicking before selecting the handle
		if(this.handle.length == 1) this.previousHandle = this.handle;


		if(this.handle.length == 2 && o.range) this.createRange();

	};

	$.extend($.ui.slider.prototype, {
		plugins: {},
		createRange: function() {
			this.rangeElement = $('<div></div>')
				.addClass('ui-slider-range')
				.css({ position: 'absolute' })
				.css(this.properties[0], parseInt($(this.handle[0]).css(this.properties[0]),10) + this.handleSize(0)/2)
				.css(this.properties[1], parseInt($(this.handle[1]).css(this.properties[0]),10) - parseInt($(this.handle[0]).css(this.properties[0]),10))
				.appendTo(this.element);
		},
		updateRange: function() {
				this.rangeElement.css(this.properties[0], parseInt($(this.handle[0]).css(this.properties[0]),10) + this.handleSize(0)/2);
				this.rangeElement.css(this.properties[1], parseInt($(this.handle[1]).css(this.properties[0]),10) - parseInt($(this.handle[0]).css(this.properties[0]),10));
		},
		getRange: function() {
			return this.rangeElement ? this.convertValue(parseInt(this.rangeElement.css(this.properties[1]),10)) : null;
		},
		ui: function(e) {
			return {
				instance: this,
				options: this.options,
				handle: this.currentHandle,
				value: this.value(),
				range: this.getRange()
			};
		},
		propagate: function(n,e) {
			$.ui.plugin.call(this, n, [e, this.ui()]);
			this.element.triggerHandler(n == "slide" ? n : "slide"+n, [e, this.ui()], this.options[n]);
		},
		destroy: function() {
			this.element
				.removeClass("ui-slider ui-slider-disabled")
				.removeData("ul-slider")
				.unbind(".slider");
			this.handles.removeMouseInteraction();
		},
		enable: function() {
			this.element.removeClass("ui-slider-disabled");
			this.disabled = false;
		},
		disable: function() {
			this.element.addClass("ui-slider-disabled");
			this.disabled = true;
		},
		focus: function(handle,hard) {
			this.currentHandle = $(handle).addClass('ui-slider-handle-active');
			if(hard) this.currentHandle.parent()[0].focus();
		},
		blur: function(handle) {
			$(handle).removeClass('ui-slider-handle-active');
			if(this.currentHandle && this.currentHandle[0] == handle) { this.previousHandle = this.currentHandle; this.currentHandle = null; };
		},
		value: function(handle) {
			if(this.handle.length == 1) this.currentHandle = this.handle;
			return ((parseInt($(handle != undefined ? this.handle[handle] || handle : this.currentHandle).css(this.properties[0]),10) / (this.size - this.handleSize())) * this.options.realMaxValue) + this.options.minValue;
		},
		convertValue: function(value) {
			return (value / (this.size - this.handleSize())) * this.options.realMaxValue;
		},
		translateValue: function(value) {
			return ((value - this.options.minValue) / this.options.realMaxValue) * (this.size - this.handleSize());
		},
		handleSize: function(handle) {
			return $(handle != undefined ? this.handle[handle] : this.currentHandle)['outer'+this.properties[1].substr(0,1).toUpperCase()+this.properties[1].substr(1)]();
		},
		click: function(e) {

			// This method is only used if:
			// - The user didn't click a handle
			// - The Slider is not disabled
			// - There is a current, or previous selected handle (otherwise we wouldn't know which one to move)
			var pointer = [e.pageX,e.pageY];
			var clickedHandle = false; this.handle.each(function() { if(this == e.target) clickedHandle = true;  });
			if(clickedHandle || this.disabled || !(this.currentHandle || this.previousHandle)) return;

			//If a previous handle was focussed, focus it again
			if(this.previousHandle) this.focus(this.previousHandle, 1);

			//Move focussed handle to the clicked position
			this.offset = this.element.offset();
			this.moveTo(this.convertValue(e[this.properties[0] == 'top' ? 'pageY' : 'pageX'] - this.offset[this.properties[0]] - this.handleSize()/2));

		},
		start: function(e, handle) {

			var o = this.options;

			this.offset = this.element.offset();
			this.handleOffset = this.currentHandle.offset();
			this.clickOffset = { top: e.pageY - this.handleOffset.top, left: e.pageX - this.handleOffset.left };
			this.firstValue = this.value();

			this.propagate('start', e);
			return false;

		},
		stop: function(e) {
			this.propagate('stop', e);
			if(this.firstValue != this.value()) this.propagate('change', e);
			return false;
		},
		drag: function(e, handle) {

			var o = this.options;
			var position = { top: e.pageY - this.offset.top - this.clickOffset.top, left: e.pageX - this.offset.left - this.clickOffset.left};

			var modifier = position[this.properties[0]];
			if(modifier >= this.size - this.handleSize()) modifier = this.size - this.handleSize();
			if(modifier <= 0) modifier = 0;

			if(o.stepping) {
				var value = this.convertValue(modifier);
				value = Math.round(value / o.stepping) * o.stepping;
				modifier = this.translateValue(value);
			}

			if(this.rangeElement) {
				if(this.currentHandle[0] == this.handle[0] && modifier >= this.translateValue(this.value(1))) modifier = this.translateValue(this.value(1));
				if(this.currentHandle[0] == this.handle[1] && modifier <= this.translateValue(this.value(0))) modifier = this.translateValue(this.value(0));
			}

			this.currentHandle.css(this.properties[0], modifier);
			if(this.rangeElement) this.updateRange();
			this.propagate('slide', e);
			return false;

		},
		moveTo: function(value, handle) {

			var o = this.options;
			if(handle == undefined && !this.currentHandle && this.handle.length != 1) return false; //If no handle has been passed, no current handle is available and we have multiple handles, return false
			if(handle == undefined && !this.currentHandle) handle = 0; //If only one handle is available, use it
			if(handle != undefined) this.currentHandle = this.previousHandle = $(this.handle[handle] || handle);

			if(value.constructor == String) value = /\-\=/.test(value) ? this.value() - parseInt(value.replace('-=', ''),10) : this.value() + parseInt(value.replace('+=', ''),10);
			if(o.stepping) value = Math.round(value / o.stepping) * o.stepping;
			value = this.translateValue(value);

			if(value >= this.size - this.handleSize()) value = this.size - this.handleSize();
			if(value <= 0) value = 0;
			if(this.rangeElement) {
				if(this.currentHandle[0] == this.handle[0] && value >= this.translateValue(this.value(1))) value = this.translateValue(this.value(1));
				if(this.currentHandle[0] == this.handle[1] && value <= this.translateValue(this.value(0))) value = this.translateValue(this.value(0));
			}

			this.currentHandle.css(this.properties[0], value);
			if(this.rangeElement) this.updateRange();

			this.propagate('start', null);
			this.propagate('stop', null);
			this.propagate('change', null);

		}
	});

})(jQuery);