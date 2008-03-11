/* Italian initialisation for the jQuery UI date picker plugin. */
/* Written by Apaella (apaella@gmail.com). */
jQuery(function($){
	$.datepicker.regional['it'] = {clearText: 'Svuota', clearStatus: '',
		closeText: 'Chiudi', closeStatus: '',
		prevText: '&lt;Prec', prevStatus: '',
		nextText: 'Succ&gt;', nextStatus: '',
		currentText: 'Oggi', currentStatus: '',
		monthNames: ['Gennaio','Febbraio','Marzo','Aprile','Maggio','Giugno',
		'Luglio','Agosto','Settembre','Ottobre','Novembre','Dicembre'],
		monthNamesShort: ['Gen','Feb','Mar','Apr','Mag','Giu',
		'Lug','Ago','Set','Ott','Nov','Dic'],
		monthStatus: '', yearStatus: '',
		weekHeader: 'Sm', weekStatus: '',
		dayNames: ['Domenica','Lunedi','Martedi','Mercoledi','Giovedi','Venerdi','Sabato'],
		dayNamesShort: ['Dom','Lun','Mar','Mer','Gio','Ven','Sab'],
		dayNamesMin: ['Do','Lu','Ma','Me','Gio','Ve','Sa'],
		dayStatus: 'DD', dateStatus: 'D, M d',
		dateFormat: 'dd/mm/yy', firstDay: 0,
		initStatus: '', isRTL: false};
	$.datepicker.setDefaults($.datepicker.regional['it']);
});
