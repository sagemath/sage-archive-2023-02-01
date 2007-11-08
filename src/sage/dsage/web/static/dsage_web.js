var count = 10; // Default count

function getServerDetails () {
    $('#server_details').load('get_server_details');
    showServer()
}

function showServer () {
    $('#jobs_table').hide();
    $('#job_details').hide();
    $('#help').hide();
    $('#server_details').show();
    $('#navigator').hide();
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

$(document).ready(function()
    {
        getJobs(count) // Run it the first time
        showJobs() // Show jobs
        // setInterval('getJobs(count)', 5000); // Update it every 5 seconds.
    }
);