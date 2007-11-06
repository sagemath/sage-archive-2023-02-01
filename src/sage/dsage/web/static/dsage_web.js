var count = 10; // Default count

$(document).ready(function()
    {
        getJobs(count) // Run it the first time
        setInterval('getJobs(count)', 5000); // Update it every 5 seconds.
    }
);

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