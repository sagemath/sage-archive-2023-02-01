$(document).ready(function()
    {
        $('#jobs_table').load('get_jobs', function() {
            $("#jobs_table").tablesorter(); } );
    }
);
