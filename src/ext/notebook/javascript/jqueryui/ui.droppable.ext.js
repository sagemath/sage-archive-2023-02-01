(function($) {

	$.ui.plugin.add("droppable", "activeClass", {
		activate: function(e,ui) {
			$(this).addClass(ui.options.activeClass);
		},
		deactivate: function(e,ui) {
			$(this).removeClass(ui.options.activeClass);
		},
		drop: function(e,ui) {
			$(this).removeClass(ui.options.activeClass);
		}
	});

	$.ui.plugin.add("droppable", "hoverClass", {
		over: function(e,ui) {
			$(this).addClass(ui.options.hoverClass);
		},
		out: function(e,ui) {
			$(this).removeClass(ui.options.hoverClass);
		},
		drop: function(e,ui) {
			$(this).removeClass(ui.options.hoverClass);
		}
	});

})(jQuery);
