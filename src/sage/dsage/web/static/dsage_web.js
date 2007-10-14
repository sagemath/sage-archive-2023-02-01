//enter refresh time in "minutes:seconds" Minutes should range from 0 to inifinity. Seconds should range from 0 to 59
var limit="0:10"

if (document.images){
    var parselimit=limit.split(":")
    parselimit=parselimit[0]*60+parselimit[1]*1
}
function beginrefresh(){
    if (!document.images)
        return
    if (parselimit==1)
        window.location.reload()
    else{
        parselimit-=1
        curmin=Math.floor(parselimit/60)
        cursec=parselimit%60
        if (curmin!=0)
            curtime=curmin+" minutes and "+cursec+" seconds left until page refresh!"
        else
            curtime=cursec+" seconds left until page refresh!"
        window.status=curtime
        setTimeout("beginrefresh()",1000)
    }
}

function add_row() {
    var jobs_table = document.getElementById("jobs_table");
    top_row = jobs_table.insertRow(-1);
    top_row.className = 'tr0'
    first_cell = top_row.insertCell(0)
    first_cell.innerHTML= 'Woot'
}

function init() {
    add_row()
    beginrefresh()
}

new Ajax.Request("/get_details?job_id=fooo",
  {
    method:'post',
    onSuccess: function(transport){
      var response = transport.responseText || "no response text";
      alert("Success! \n\n" + response);
    },
    onFailure: function(){ alert('Something went wrong...') }
  });
window.onload=init()