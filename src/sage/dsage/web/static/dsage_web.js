var count = 10; // Default count

function getServerDetails () {
    $('#server_details').load('get_server_details');
    showServer()
}

function showServer () {
    $('#jobs_table').hide();
    $('#job_details').hide();
    $('#help').hide();
    $('#navigator').hide();
    $('#server_details').show();
}
function getJobs (c) {
    count = c // Set global variable count to the selected count.

    // Load jobs table with the given number of jobs.
    $('#jobs_table').load('get_jobs', {'count': c}, function() {
            $("#jobs_table").tablesorter(); } );

    // Change current selected count to red
    $('#navigator > a').css('color', 'black')
    $('#' + c).css('color', 'red')
}

function showJobs () {
    $('#server_details').hide();
    $('#job_details').hide();
    $('#help').hide();
    $('#jobs_table').show();
    $('#navigator').show();
}

function getJobDetails (job_id) {
    $('#job_details').load('get_details', {'job_id': job_id}, function() {
                                        $("#job_details").tablesorter(); } );
    showJobDetails()
}

function showJobDetails () {
    $('#jobs_table').hide();
    $('#navigator').hide();
    $('#server_details').hide();
    $('#help').hide();
    $('#job_details').show();
}

function getHelp () {
    $('#help').load('get_help');
    showHelp()
}

function showHelp () {
    $('#jobs_table').hide();
    $('#navigator').hide();
    $('#server_details').hide();
    $('#job_details').hide();
    $('#help').show();
}

// PageLoad function
// This function is called when:
// 1. after calling $.historyInit();
// 2. after calling $.historyLoad();
// 3. after pushing "Go Back" button of a browser
function pageload(hash) {
	// hash doesn't contain the first # character.
	if(hash) {
		// restore ajax loaded state
        if (hash=='jobs') {
            showJobs()
        }
        else if (hash=='server') {
            getServerDetails()
        }
        else if (hash=='help') {
            getHelp()
        }
	}
	else {
		// start page
		showJobs()
	}
}

function gethash (that) {
    var hash = that.href;
	hash = hash.replace(/^.*#/, '');
	// moves to a new page.
	// pageload is called at once.
    $.historyLoad(hash);
	return false;
}

$(document).ready(function()
    {
        // Initialize history plugin.
		// The callback is called at once by present location.hash.
		$.historyInit(pageload);
		console.log('historyInit')
		// set onlick event for buttons
		$("a[@rel='history']").click(function() {gethash(this); });
		console.log('setonclick')
        getJobs(count) // Run it the first time
        // showJobs() // Show jobs
        // setInterval('getJobs(count)', 5000); // Update it every 5 seconds.
    }
);