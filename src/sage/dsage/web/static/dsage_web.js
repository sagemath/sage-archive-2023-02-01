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

window.onload=beginrefresh