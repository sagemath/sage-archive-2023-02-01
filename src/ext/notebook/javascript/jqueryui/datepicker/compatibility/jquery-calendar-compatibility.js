/* jQuery UI Date Picker v3.3 - compatibility with jQuery Calendar v2.7
   Written by Marc Grabanski (m@marcgrabanski.com) and Keith Wood (kbwood@iprimus.com.au).

   Copyright (c) 2007 Marc Grabanski (http://marcgrabanski.com/code/jquery-calendar)
   Dual licensed under the GPL (http://www.gnu.org/licenses/gpl-3.0.txt) and
   CC (http://creativecommons.org/licenses/by/3.0/) licenses. "Share or Remix it but please Attribute the authors."
   Date: 09-03-2007  */

(function($) { // hide the namespace

function PopUpCal() {
	this.regional = $.datepicker.regional;
}

$.extend(PopUpCal.prototype, {

	/* Override the default settings for all instances of the calendar.
	   @param  settings  object - the new settings to use as defaults (anonymous object)
	   @return void */
	setDefaults: function(settings) {
		$.datepicker.setDefaults(convertSettings(settings));
	},

	/* Pop-up the calendar in a "dialog" box.
	   @param  dateText  string - the initial date to display (in the current format)
	   @param  onSelect  function - the function(dateText) to call when a date is selected
	   @param  settings  object - update the dialog calendar instance's settings (anonymous object)
	   @param  pos       int[2] - coordinates for the dialog's position within the screen
			leave empty for default (screen centre)
	   @return void */
	dialogCalendar: function(dateText, onSelect, settings, pos) {
		$.datepicker.dialogDatepicker(dateText, onSelect, convertSettings(settings), pos);
	},

	/* Enable the input field(s) for entry.
	   @param  inputs  element - single input field or
	                   string - the ID or other jQuery selector of the input field(s) or
	                   object - jQuery collection of input fields
	   @return void */
	enableFor: function(inputs) {
		(inputs.jquery ? inputs : $(inputs)).enableDatepicker();
	},

	/* Disable the input field(s) from entry.
	   @param  inputs  element - single input field or
	                   string - the ID or other jQuery selector of the input field(s) or
	                   object - jQuery collection of input fields
	   @return void */
	disableFor: function(inputs) {
		(inputs.jquery ? inputs : $(inputs)).disableDatepicker();
	},

	/* Update the settings for a calendar attached to an input field or division.
	   @param  control   element - the input field or div/span attached to the calendar or
	                     string - the ID or other jQuery selector of the input field
	   @param  settings  object - the new settings to update
	   @return void */
	reconfigureFor: function(control, settings) {
		$(control).changeDatepicker(convertSettings(settings));
	},

	/* Set the date for a calendar attached to an input field or division.
	   @param  control  element - the input field or div/span attached to the calendar or
	                    string - the ID or other jQuery selector of the input field
	   @param  date     Date - the new date
	   @return void */
	setDateFor: function(control, date) {
		$(control).setDatepickerDate(date);
	},

	/* Retrieve the date for a calendar attached to an input field or division.
	   @param  control  element - the input field or div/span attached to the calendar or
	                    string - the ID or other jQuery selector of the input field
	   @return Date - the current date */
	getDateFor: function(control) {
		return $(control).getDatepickerDate();
	},

	/* Pop-up the calendar for a given input field.
	   @param  target  element - the input field attached to the calendar or
	                   string - the ID or other jQuery selector of the input field
	   @return void */
	showFor: function(target) {
		$.datepicker.showFor(target);
	},

	/* Hide the calendar from view.
	   @param  id     string/object - the ID of the current calendar instance,
			or the instance itself
	   @param  speed  string - the speed at which to close the calendar
	   @return void */
	hideCalendar: function(id, speed) {
		$.datepicker.hideDatepicker(speed);
	},

	/* Set as customDate function to prevent selection of weekends.
	   @param  date  Date - the date to customise
	   @return [boolean, string] - is this date selectable?, what is its CSS class? */
	noWeekends: function(date) {
		return $.datepicker.noWeekends(date);
	},

	/* Format a date object into a string value.
	   @param  date  Date - the date to customise */
	formatDate: function(date) {
		return $.datepicker.formatDate($.datepicker._defaults.dateFormat, date);
	}
});

/* Translate the calendar settings. */
function convertSettings(settings) {
	if (settings) {
		if (settings.autoPopUp) {
			settings.showOn = settings.autoPopUp;
			settings.autoPopUp = null;
		}
		if (settings.fieldSettings) {
			settings.beforeShow = settings.fieldSettings;
			settings.fieldSettings = null;
		}
		if (settings.customDate) {
			settings.beforeShowDay = settings.customDate;
			settings.customDate = null;
		}
	}
	return settings;
}

/* Attach the calendar to a jQuery selection.
   Convert to use the new jQuery Date Picker functionality.
   @param  settings  object - the new settings to use for this calendar instance (anonymous)
   @return jQuery object - for chaining further calls */
$.fn.calendar = function(settings) {
	this.each(function() {
		for (attrName in $.datepicker._defaults) {
			var attrValue = this.getAttribute('cal:' + attrName);
			if (attrValue) {
				this.setAttribute('date:' + attrName, attrValue);
				this.removeAttribute('cal:' + attrName);
			}
		}
	});
	return this.attachDatepicker(convertSettings(settings));
};

$.fn.datepicker = $.fn.calendar;

/* Initialise the calendar. */
$(document).ready(function() {
	popUpCal = new PopUpCal(); // singleton instance
});

})(jQuery);
