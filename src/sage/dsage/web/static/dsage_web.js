function add_row() {
    var jobs_table = document.getElementById("jobs_table");
    top_row = jobs_table.insertRow(-1);
    top_row.className = 'tr0';
    first_cell = top_row.insertCell(0);
    first_cell.innerHTML= 'Woot';
}

function init() {

}

function build_jobs_table(transport) {
    var xmldoc = XML.newDocument();
    xmldoc.load(transport.responseText);
}

var get_jobs = new Ajax.PeriodicalUpdater('jobs_table', '/get_jobs',
  {
    method: 'get',
    insertion: Insertion.Bottom,
    frequency: 1,
    decay: 2
  });

// new Ajax.Request("/get_jobs",
//   {
//     method:'post',
//     onSuccess: build_jobs_table(transport),
//     onFailure: function(){ alert('Something went wrong...') }
//   });
// window.onload=init()