var count = 10; // Default count

// PageLoad function
// This function is called when:
// 1. after calling $.historyInit();
// 2. after calling $.historyLoad();
// 3. after pushing "Go Back" button of a browser
function pageload(hash) {
	// hash doesn't contain the first # character.
	if(hash) {
		// restore ajax loaded state
		$("#load").load(hash + ".html");
	} else {
		// start page
		$("#load").empty();
	}
}

function getServerDetails () {
    $('#jobs_table').hide();
    $('#job_details').hide();
    $('#help').hide();

    $('#server_details').load('get_server_details');
    $('#server_details').show();
}

function getJobs (c) {
    count = c
    $('#server_details').hide();
    $('#job_details').hide();
    $('#help').hide();
    $('#jobs_table').load('get_jobs', {'count': count}, function() {
            $("#jobs_table").tablesorter(); } );
    $('#jobs_table').show();
}

function getJobDetails (job_id) {
    $('#jobs_table').hide();
    $('#server_details').hide();
    $('#help').hide();
    $('#job_details').load('get_details', {'job_id': job_id}, function() {
                                        $("#job_details").tablesorter(); } );
    $('#job_details').show();
}

function getHelp () {
    $('#jobs_table').hide();
    $('#server_details').hide();
    $('#job_details').hide();
    $('#help').load('get_help');
    $('#help').show();
}

$(document).ready(function()
    {
        // Initialize history plugin.
		// The callback is called at once by present location.hash.
		$.historyInit(pageload);

		// set onlick event for buttons
		$("a[@rel='history']").click(function(){
			//
			var hash = this.href;
			hash = hash.replace(/^.*#/, '');
			// moves to a new page.
			// pageload is called at once.
			$.historyLoad(hash);
			return false;
		});

        getJobs(count) // Run it the first time
        setInterval('getJobs(count)', 5000); // Update it every 5 seconds.
    }
);