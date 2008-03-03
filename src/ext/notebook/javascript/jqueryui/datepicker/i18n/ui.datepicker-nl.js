/* Dutch (UTF-8) initialisation for the jQuery UI date picker plugin. */
jQuery(function($){
	$.datepicker.regional['nl'] = {clearText: 'Wissen', clearStatus: '',
		closeText: 'Sluiten', closeStatus: '',
		prevText: '&lt;Terug', prevStatus: '',
		nextText: 'Volgende&gt;', nextStatus: '',
		currentText: 'Vandaag', currentStatus: '',
		monthNames: ['Januari','Februari','Maart','April','Mei','Juni',
		'Juli','Augustus','September','Oktober','November','December'],
		monthNamesShort: ['Jan','Feb','Maa','Apr','Mei','Jun',
		'Jul','Aug','Sep','Okt','Nov','Dec'],
		monthStatus: '', yearStatus: '',
		weekHeader: 'Wk', weekStatus: '',
		dayNames: ['Zondag','Maandag','Dinsdag','Woensdag','Donderdag','Vrijdag','Zaterdag'],
		dayNamesShort: ['Zon','Maa','Din','Woe','Don','Vri','Zat'],
		dayNamesMin: ['Zo','Ma','Di','Wo','Do','Vr','Za'],
		dayStatus: 'DD', dateStatus: 'D, M d',
		dateFormat: 'dd.mm.yy', firstDay: 0,
		initStatus: '', isRTL: false};
	$.datepicker.setDefaults($.datepicker.regional['nl']);
});