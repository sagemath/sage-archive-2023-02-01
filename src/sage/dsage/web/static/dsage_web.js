$(document).ready(function()
    {
        getJobs()
        $.ajaxHistory.initialize();
    }
);

function getServerDetails () {
    $('#jobs_table').hide();
    $('#job_details').hide();
    $('#server_details').load('get_server_details');
    $('#server_details').show();
}

function getJobs () {
    $('#jobs_table').load('get_jobs', function() {
            $("#jobs_table").tablesorter(); } );
    $('#jobs_table').show();
    $('#server_details').hide();
    $('#job_details').hide();
}

function getJobDetails (job_id) {
    $('#jobs_table').hide();
    $('#server_details').hide();
    $('#job_details').load('get_details', {'job_id':job_id}, function() {
                                        $("#job_details").tablesorter(); } );
    $('#job_details').show();
}