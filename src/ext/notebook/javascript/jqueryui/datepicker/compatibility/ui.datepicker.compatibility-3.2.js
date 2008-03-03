/* jQuery UI Date Picker v3.3 compatibility with jQuery UI Date Picker v3.2.
   Written by Marc Grabanski (m@marcgrabanski.com) and Keith Wood (kbwood@virginbroadband.com.au).

   Copyright (c) 2007 Marc Grabanski (http://marcgrabanski.com/code/ui-datepicker)
   Dual licensed under the MIT (MIT-LICENSE.txt)
   and GPL (GPL-LICENSE.txt) licenses.
   Date: 09-03-2007  */

(function($) { // hide the namespace

/* Attach the date picker to a jQuery selection.
   @param  settings  object - the new settings to use for this date picker instance (anonymous)
   @return jQuery object - for chaining further calls */
$.fn.datepicker = function(settings) {
	return this.attachDatepicker(settings);
};

$(document).ready(function() {
	// Add the old functions back again
	$.extend($.datepicker, {
		enableFor: function(control) {
			(control.jquery ? control : $(control)).enableDatepicker();
			return this;
		},

		disableFor: function(control) {
			(control.jquery ? control : $(control)).disableDatepicker();
			return this;
		},

		isDisabled: function(control) {
			return (control.jquery ? control : $(control)).isDisabledDatepicker();
		},

		reconfigureFor: function(control, settings) {
			(control.jquery ? control : $(control)).changeDatepicker(settings);
			return this;
		},

		setDateFor: function(control, date, endDate) {
			(control.jquery ? control : $(control)).setDatepickerDate(date, endDate);
			return this;
		},

		getDateFor: function(control) {
			return (control.jquery ? control : $(control)).getDatepickerDate();
		}
	});
});

})(jQuery);
