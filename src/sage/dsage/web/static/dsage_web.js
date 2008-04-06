var count = 10; // Default count
var page_num = 0;
var pages = 2147483647  ;

function hideAll () {
    // Hides all tables
    $("table").each(function() {$(this).hide();});

    // Hide extra div"s
    $("#help").hide();
    $("#navigator").hide();
    $("#user_management").hide();
    $("#success").hide();
    $("#quick_help").hide();
    $("#page_nav").hide();
}

function showJobs() {
    hideAll()
    $("#jobs_table").show();
    $("#navigator").show();
    $("#page_nav").show();
}

function getServerDetails () {
    $("#server_details").load("get_server_details", function() {
                                    $("#server_details").tablesorter(); });
    $("#workers").load("get_workers",
                       function() { $("#workers").tablesorter();} );
    $("#clients").load("get_clients",
                      function() { $("#clients").tablesorter();} );
    showServer();
}

function userManagement () {
    hideAll()
    $("#user_management").show()
}

function showServer () {
    hideAll()
    $("#server_details").show();
    $("#workers").show()
    $("#clients").show()
}

function colorStatus() {
    // Change color of status
    // completed -> green
    // processing -> blue
    // new -> red
    $(".job_status").each(function() {
        if ($(this).text() == "failed") {
            $(this).css("color", "red");
        }
        if ($(this).text() == "killed") {
            $(this).css("color", "grey");
        }
        else if ($(this).text() == "new") {
            $(this).css("color", "orange")
        }
        else if ($(this).text() == "processing") {
            $(this).css("color", "blue");
        }
        else if ($(this).text() == "completed") {
            $(this).css("color", "green");
        }
    });
}

function check_page() {
    if (page_num == 0) {
        $("#later").unbind('click').click(function() { return false;});
        $("#later").css("color", "grey");
    }
    else {
        $("#later").click(later);
        $("#later").css("color", "blue");
    }
    if (page_num >= pages) {
        $("#earlier").unbind('click').click(function() { return false;});
        $("#earlier").css("color", "grey");
    }
    else {
        $("#earlier").click(earlier);
        $("#earlier").css("color", "blue");
    }
}

function loadJobs() {
    $("#jobs_table tbody").load("get_page", {"count": count, "n": page_num},
        function() {
            colorStatus();
            $("#jobs_table").trigger("update");
            var sorting = [[4,1]];
            $("#jobs_table").trigger("sorton",[sorting]);
        });
}

function setCount(c) {
    check_page();
    count = c;
    $.get('pages', {"count": count}, function(p) {
        pages = parseInt(p);
    });
    $("#navigator > a").css("color", "black");
    $("#" + c).css("color", "red");
}

function updateCurrentPage() {
    loadJobs();
}
function later() {
    if (page_num > 0) {
        page_num = page_num - 1;
    }
    check_page();
    loadJobs();
}

function earlier() {
    page_num = page_num + 1;
    check_page();
    loadJobs();
}

function getJobDetails (job_id) {
    $("#job_details").load("get_details",
                           {"job_id": job_id},
                            function() { $("#job_details").tablesorter();
                           }
                          );
    showJobDetails();
}

function showJobDetails () {
    hideAll()
    $("#job_details").show();
}

function getHelp () {
    $("#help").load("get_help");
    showHelp();
}

function showHelp () {
    hideAll()
    $("#help").show();
    $("#quick_help").show();
}

// PageLoad function
// This function is called when:
// 1. after calling $.historyInit();
// 2. after calling $.historyLoad();
// 3. after pushing "Go Back" button of a browser
function pageload(hash) {
	// hash doesn"t contain the first # character.
	if(hash) {
		// restore ajax loaded state
        if (hash=="jobs") {
            showJobs();
        }
        else if (hash=="server") {
            getServerDetails();
        }
        else if (hash=="help") {
            getHelp();
        }
	}
	else {
		// start page
		showJobs();
	}
}

function gethash (that) {
    var hash = that.href;
	hash = hash.replace(/^.*#/, "");
    $.historyLoad(hash);
	return false;
}

function showAddedUser() {
    hideAll();
    $("#add_client")[0].reset();
    $("#success").show("fast")
}

$(document).ready(function() {
        // Initialize history plugin.
		// The callback is called at once by present location.hash.
		$.historyInit(pageload);
		// set onlick event for buttons
		$("a[@rel='history']").click(function() {gethash(this); });

        // initialize tablesorter
		$("#jobs_table").tablesorter();

		// Initialize jobs table
        setCount(10);
        updateCurrentPage();
        $("#add_client").ajaxForm(function() {
            showAddedUser();
        });
        // Only show jobs if this.href is none or jobs
        if (document.location.hash.length == 11) {
            hash = document.location.hash.replace(/^.*#/, "");
            getJobDetails(hash);
        }
        // Only show jobs if this.href is none or jobs
        if (document.location.hash == "" || document.location.hash == "#jobs")
        {
            showJobs();
        }
        setInterval("updateCurrentPage()", 5000); // Update it every 5 seconds.
    }
);