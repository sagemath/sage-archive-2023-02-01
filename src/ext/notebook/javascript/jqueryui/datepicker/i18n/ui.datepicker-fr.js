/* French initialisation for the jQuery UI date picker plugin. */
/* Written by Keith Wood (kbwood@iprimus.com.au). */
jQuery(function($){
	$.datepicker.regional['fr'] = {clearText: 'Effacer', clearStatus: '',
		closeText: 'Fermer', closeStatus: '',
		prevText: '&lt;Préc', prevStatus: '',
		nextText: 'Proch&gt;', nextStatus: '',
		currentText: 'En cours', currentStatus: '',
		monthNames: ['Janvier','Février','Mars','Avril','Mai','Juin',
		'Juillet','Août','Septembre','Octobre','Novembre','Décembre'],
		monthNamesShort: ['Jan','Fév','Mar','Avr','Mai','Jun',
		'Jul','Aoû','Sep','Oct','Nov','Déc'],
		monthStatus: '', yearStatus: '',
		weekHeader: 'Sm', weekStatus: '',
		dayNames: ['Dimanche','Lundi','Mardi','Mercredi','Jeudi','Vendredi','Samedi'],
		dayNamesShort: ['Dim','Lun','Mar','Mer','Jeu','Ven','Sam'],
		dayNamesMin: ['Di','Lu','Ma','Me','Je','Ve','Sa'],
		dayStatus: 'DD', dateStatus: 'D, M d',
		dateFormat: 'dd/mm/yy', firstDay: 0,
		initStatus: '', isRTL: false};
	$.datepicker.setDefaults($.datepicker.regional['fr']);
});