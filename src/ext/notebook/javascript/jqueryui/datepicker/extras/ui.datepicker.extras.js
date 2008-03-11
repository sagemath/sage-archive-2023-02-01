/* jQuery UI Date Picker Extensions v3.2
   Written by Marc Grabanski (m@marcgrabanski.com) and Keith Wood (kbwood@iprimus.com.au).

   Copyright (c) 2007 Marc Grabanski (http://marcgrabanski.com/code/ui-datepicker)
   Dual licensed under the MIT (MIT-LICENSE.txt)
   and GPL (GPL-LICENSE.txt) licenses.
   Date: 09-03-2007  */

/* Datepicker object providing beforeShowDay function to prevent selection of special days.
   Use it like the following:
   var specials = new SpecialDays(false, {thirteenth: ['13/01', '13/02', '13/03']}, 'dd/mm');
   $('#date').datepicker({beforeShowDay: function(date) { return specials.noSpecials(date); }});
   @param  allowWeekends  boolean - true if weekends are selectable, false if not
   @param  days           Object - object with named fields that contain arrays of
                          strings indicating days and months -
                          field name becomes the CSS class for those dates
   @param  dateFormat     string - format of day and month in days above */
function SpecialDays(allowWeekends, days, dateFormat) {
	this._allowWeekends = allowWeekends;
	this._days = days;
	this._dateFormat = dateFormat;
	for (var name in this._days) {
		var days = this._days[name];
		for (var d in days) {
			days[d] = $.datepicker.parseDate(this._dateFormat, days[d]);
		}
	}
}

$.extend(SpecialDays.prototype, {
	/* Set as datepicker beforeShowDay function to prevent selection
	   of special days and (optionally) weekends.
	   @param  date  Date - the date to customise
	   @return [boolean, string] - is this date selectable?, what is its CSS class? */
	noSpecials: function(date) {
		for (var name in this._days) {
			var days = this._days[name];
			for (var d in days) {
				if (date.getDate() == days[d].getDate() &&
						date.getMonth() == days[d].getMonth()) {
					return [false, name];
				}
			}
		}
		var day = date.getDay();
		return [this._allowWeekends || (day > 0 && day < 6), ''];
	}
});
