"""
Javascript (AJAX) Component of SAGE Notebook

AUTHORS:
    -- William Stein
    -- Tom Boothby
    -- Alex Clemesha


This file is one big raw triple-quoted string that contains a bunch of
javascript.  This javascript is inserted into the head of the notebook
web page.
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################

def javascript():
    return r"""

///////////////////////////////////////////////////////////////////
// An AJAX framework for connections back to the
// SAGE server (written by Tom Boothby).
///////////////////////////////////////////////////////////////////

var asyncObj
function getAsyncObject(handler) {
  var asyncObj=null

  if (navigator.userAgent.indexOf("Opera")>=0) {
    alert("SAGE Notebook requires Firefox or Internet Explorer");
  }

  if (navigator.userAgent.indexOf("MSIE")>=0) {
    var strName = "Msxml2.XMLHTTP";
    if (navigator.appVersion.indexOf("MSIE 5.5")>=0) {
      strName = "Microsoft.XMLHTTP";
    }
    try {
      asyncObj = new ActiveXObject(strName);
      asyncObj.onreadystatechange = handler;
      return objXmlHttp;
    }
    catch(e) {
      alert("Error. Scripting for ActiveX disabled.  Enable ActiveX in your browser controls, or use Firefox.") ;
      return
    }
  }
  if (navigator.userAgent.indexOf("Mozilla")>=0) {
    asyncObj = new XMLHttpRequest();
    asyncObj.onload  = handler;
    asyncObj.onerror = handler;
    return asyncObj;
  }
}

function asyncCallbackHandler(name, callback) {
    function f() {
                 eval('asyncObj = ' + name);
                 try {
                   if ((asyncObj.readyState==4 || asyncObj.readyState=="complete")
                           && asyncObj.status == 200)
                       callback("success", asyncObj.responseText);
                   } catch(e) {
                       callback("failure", e);
                   } finally { }
              };
    return f;
}

function asyncRequest(name, url, callback, postvars) {
  f = asyncCallbackHandler(name, callback);
  asyncObj = getAsyncObject(f);
  eval(name + '=asyncObj;');

  if(postvars != null) {
    asyncObj.open('POST',url,false);
    asyncObj.setRequestHeader('Content-Type','application/x-www-form-urlencoded');
    asyncObj.send(postvars);
  } else {
    asyncObj.open('GET',url,true);
    asyncObj.setRequestHeader('Content-Type',  "text/html");
    asyncObj.send(null);
  }
}

///////////////////////////////////////////////////////////////////
//
// WORKBOOK functions -- for the individual cells
//
///////////////////////////////////////////////////////////////////

function cell_key_event(number, event) {
    var theCode = event.keyCode ? event.keyCode :
                   event.which ? event.which : event.charCode;

    if (theCode == 13 && event.shiftKey) {
       // User pressed shift-enter, so we submit the input text for computation.
       document.forms[number].submit();
       return false;

    }

    else if (theCode == 9) {
       var txt = document.getElementById('in' + number).value;
       var j=0;

       for (i = txt.length-1; i>=0 && (txt.charAt(i) == ' ' || txt.charAt(i) == '\n'); i--) {
       }
       for (i=txt.length-1; i>=0; i--) {
           var c = txt.charAt(i);
           lower_alpha = (c >= 'a' && c <= 'z');
           upper_alpha = (c >= 'A' && c <= 'Z');
           digit = (c >= '0' && c <= '9');
           if (! (c == '.' || c == '_' || c == '\n' || c == '?' || lower_alpha || upper_alpha || digit)) {
               j=i+1;
               break;
           }
       }
       txt = txt.slice(j);
       document.getElementById("search_input").value = txt;
       search_box();
    } else {
       return true;
    }
}

function changeAreaSize(val,id) {
    var el = document.getElementById(id);
    if (val==1)
        el.rows = el.rows + 1;
    else
        el.rows = el.rows - 1;
}

function toggleVisibility(id) {
    var outBox = document.getElementById('out'+id);
    if(outBox.style.display == 'none') {
        outBox.style.display = 'block';
        document.getElementById('tog'+id).innerHTML='H';
    } else {
        outBox.style.display = 'none';
        document.getElementById('tog'+id).innerHTML='S';
    }
}

function showBox(id) {
    var el = document.getElementById(id);
    el.style.display='none';
}


///////////////////////////////////////////////////////////////////
//
// INTROSPECTION functions -- for getting help
//
///////////////////////////////////////////////////////////////////

// SEARCH BOX
function search_box() {
    var s = document.getElementById("search_input").value;
    document.getElementById("search_input").style.color = "#888";
    if (s.indexOf('?') == -1) {
        callback = search_fill_in_completions;
    } else {
        search_fill_in_doc('success','searching...');
        callback = search_fill_in_doc;
    }
    asyncRequest('async_obj_search', '/search', callback, 'query='+s)
}

function search_fill_in_doc(status, response_text) {
    document.getElementById("search_input").style.color = "#000";
    if (status  == "success") {
       expand_doc_box();

       document.getElementById("search_doc_topbar").innerHTML =
           '<table bgcolor="73a6ff" width="100%" height="100%"><tr>  \
           <td align=left class="menubar">Documentation</td><td align=right class="menubar"> \
       <a class="menubar" href="javascript:shrink_doc_box()">&nbsp;&nbsp;X&nbsp;&nbsp</a></td></tr> </table>'

       document.getElementById("search_doc").innerHTML = '<pre>' + response_text + '</pre>';
    }
}

function search_fill_in_completions(status, response_text) {
    document.getElementById("search_input").style.color = "#000";
    if (status  == "success") {
       shrink_doc_box();
       document.getElementById("search_doc").innerHTML = response_text;
       document.getElementById("search_doc_topbar").innerHTML =
           '<table bgcolor="73a6ff" width="100%"><tr><td align=left class="menubar"> \
           Completions \
            </td><td></td></tr></table>';
    }
}

function shrink_doc_box() {
    document.getElementById("search_doc").style.width = '154px';
    document.getElementById("search_doc").style.height = '150px';
    document.getElementById("search_doc").style.font = 'arial';
    document.getElementById("search_doc_topbar").style.width = '158px';
    document.getElementById("search_input").style.width = '160px';
}

function expand_doc_box() {
    document.getElementById("search_doc_topbar").style.width = '700px';
    document.getElementById("search_doc").style.width = '696';
    document.getElementById("search_doc").style.font = 'courier';
    document.getElementById("search_doc").style.height = '80%';
    document.getElementById("search_input").style.width = '702px';
}


///////////////////////////////////////////////////////////////////
//
// COMPUTE functions -- for running computations and querying
// for results
//
///////////////////////////////////////////////////////////////////
"""
