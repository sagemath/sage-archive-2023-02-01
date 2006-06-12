r"""
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
#       Copyright (C) 2006 William Stein <wstein@ucsd.edu>
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

cell_output_delta = 200;

SEP = '___S_A_G_E___';

var asyncObj;
var no_async = false;
var userAgent = navigator.userAgent.toLowerCase();
function getAsyncObject(handler) {
  asyncObj=null
  try {
    if (userAgent.indexOf("msie")!=-1  && userAgent.indexOf("opera")==-1) {
      var s =(userAgent.indexOf("msie 5")!=-1)?"Microsoft.XMLHTTP":"Msxml2.XMLHTTP";
      asyncObj = new ActiveXObject(s);
      asyncObj.onreadystatechange = handler;
      return asyncObj;
    } else {
      asyncObj = new XMLHttpRequest();
      asyncObj.onload  = handler;
      asyncObj.onerror = handler;
      return asyncObj;
    }
  } catch(e) {
    no_async = true;
    return null;
  }
}

function asyncCallbackHandler(name, callback) {
    function f() {
                 eval('asyncObj = ' + name);
                 try {
                   if( (asyncObj.readyState==4 || asyncObj.readyState=="complete")
                       && asyncObj.status == 200 )
                     try {
                       callback('success', asyncObj.responseText);
                     } catch(e) {
                       callback('success', "empty");
                     }
                 } catch(e) {
                   callback("failure", e);
                 } finally { }
              };
    return f;
}

function async_request(name, url, callback, postvars) {
  f = asyncCallbackHandler(name, callback);
  asyncObj = getAsyncObject(f);
  eval(name + '=asyncObj;');

  if(postvars != null) {
    asyncObj.open('POST',url,true);
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
// OBJECT functions -- for managing saved objects
//
///////////////////////////////////////////////////////////////////

function click_on_object(name) {
/*
    o = document.open("/" + name + ".sobj");
    */
}


///////////////////////////////////////////////////////////////////
//
// WORKSHEET functions -- for switching between and managing worksheets
//
///////////////////////////////////////////////////////////////////

worksheet_id=0;   // The current worksheet.

function add_worksheet(name) {
    async_request('async_obj_add_worksheet', '/add_worksheet',
                   add_worksheet_callback, 'name='+name)
}

function add_worksheet_callback(status, response_text) {
    if (status == "success") {
        /* expect response_text to encode a pair consisting of
           the HTML for the updated worksheet list and the
           name of the new worksheet. */
        var X = response_text.split(SEP);
        if (X.length <= 1) {
            alert(X);
        } else {
            set_worksheet_list(X[0]);
            switch_to_worksheet(X[1]);
        }
    } else {
        alert("Possible failure adding workbook: \n" + response_text);
    }
}

function delete_worksheet(name) {
    async_request('async_obj_delete_worksheet', '/delete_worksheet',
                   delete_worksheet_callback, 'name='+name)
}

function delete_worksheet_callback(status, response_text) {
    if (status == "success") {
        /* expect response_text to encode a pair consisting of
           the HTML for the updated worksheet list and the
           id of a worksheet to switch to in case we just
           deleted the current worksheet. */
        var X = response_text.split(SEP);
        if (X.length <= 1) {
            alert(X);
        } else {
            set_worksheet_list(X[0]);
            if (X[1] != -1)
               switch_to_worksheet(X[1]);
        }
    } else {
        alert("Possible failure deleting workbook: \n" + response_text);
    }
}

function set_worksheet_list(worksheets) {
    var wlist = document.getElementById('worksheet_list');
    wlist.innerHTML = worksheets;
}

function show_add_new_worksheet_menu() {
    var add_worksheet_menu = document.getElementById('add_worksheet_menu');
    add_worksheet_menu.style.display = 'block';
    document.getElementById('new_worksheet_box').focus()
}

function hide_add_new_worksheet_menu() {
    var add_worksheet_menu = document.getElementById('add_worksheet_menu');
    add_worksheet_menu.style.display = 'none';
}

function show_delete_worksheet_menu() {
    var delete_worksheet_menu = document.getElementById('delete_worksheet_menu');
    delete_worksheet_menu.style.display = 'block';
    document.getElementById('delete_worksheet_box').focus();
}

function hide_delete_worksheet_menu() {
    var delete_worksheet_menu = document.getElementById('delete_worksheet_menu');
    delete_worksheet_menu.style.display = 'none';
}

function process_new_worksheet_menu_submit() {
    hide_add_new_worksheet_menu();
    var add_worksheet_box = document.getElementById('new_worksheet_box');
    name = add_worksheet_box.value;
    add_worksheet_box.value = '';
    add_worksheet(name);
}

function process_delete_worksheet_menu_submit() {
    hide_delete_worksheet_menu();
    var delete_worksheet_box = document.getElementById('delete_worksheet_box');
    name = delete_worksheet_box.value;
    delete_worksheet_box.value = '';
    delete_worksheet(name);
}

function switch_to_worksheet(id) {
    /* 1. check to see if worksheet is already loaded into the DOM
       2. If not, load it into the dom.
       3. Move it to the front and everything else to the back by changing the css.
    */
  /*  alert('switch to worksheet ' + id);  */
}


///////////////////////////////////////////////////////////////////
//
// CELL functions -- for the individual cells
//
///////////////////////////////////////////////////////////////////

function focus_on(id, id_before) {
       var cell = document.getElementById('cell_input_' + id);
       if (cell) {
          cell.focus();
       }
/*       var cell_before = document.getElementById('cell_div_output_' + id_before);  */
       var cell_before = document.getElementById('cell_input_' + id_before);
       if (cell_before)
          cell_before.scrollIntoView();
}

function cell_input_resize(number) {
   var cell_input = document.getElementById('cell_input_' + number);
   rows = cell_input.value.split('\n').length - 1;
   if (rows <= 1) {
      rows = 2;
   } else {
      /* to avoid bottom chop off */
      rows = rows + 1;
   }
   cell_input.style.height = null;
   /* cell_input.style.height = 1.5*rows + 'em'; */   // this sort of works in konqueror...
   cell_input.rows = rows;
}

function cell_input_minimize_size(number) {
   var cell_input = document.getElementById('cell_input_' + number);
   rows = cell_input.value.split('\n').length ;
   if (rows < 1) {
      rows = 1;
   }
   cell_input.rows = rows;
   if (rows == 1) {
       /* hack because of bug in firefox with 1-row textarea */
       cell_input.style.height = '1.5em';
   }
}

id_to_delete=-1;
function cell_delete_callback(status, response_text) {
    if (status == "failure") {
        cell = document.getElementById('cell_' + id_to_delete);
        var worksheet = document.getElementById('worksheet');
        worksheet.removeChild(cell);
     }
    var X = response_text.split(SEP);
    if (X[0] == 'ignore') {
        return;   /* do not delete, for some reason */
    }
    cell = document.getElementById('cell_' + X[1]);
    var worksheet = document.getElementById('worksheet');
    worksheet.removeChild(cell);
    id_before = X[2];
    var cell_before = document.getElementById('cell_input_' + id_before);
    cell_before.focus();
    cell_before.scrollIntoView();
}

function cell_delete(id) {
   id_to_delete=id;
   async_request('async_obj_cell_delete', '/delete_cell', cell_delete_callback, 'id='+id)
}

function cell_input_key_event(number, event) {

    var the_code = event.keyCode ? event.keyCode :
                   event.which ? event.which : event.charCode;

    var cell_input = document.getElementById('cell_input_' + number);


/*    alert(the_code); */
    if (the_code == 8 && cell_input.value == '') {
        cell_delete(number);
        return false;
    }

    cell_input_resize(number);

    if (the_code == 13) {
       if (event.shiftKey) {
           // User pressed shift-enter
           evaluate_cell(number, 0);
           return false;
       } else if (event.ctrlKey) {
           evaluate_cell(number, 1);
           return false;
       }
    }
/*     else if (the_code == 9 || (the_code==39 && event.ctrlKey)) { */
     else if (the_code == 27 || (the_code==39 && event.ctrlKey)) {
       // command completion: tab or ctrl->
       evaluate_cell(number, 2);
       return false;
    }
    else if (the_code == 99 && event.ctrlKey) {
       interrupt();
    }
    return true;
}

cell_id = 0;
last_action = 0;
function evaluate_cell(id, action) {
    cell_id = id;
    last_action = action;
    cell_input_minimize_size(id);
    evaluated_cell_id = 'cell_input_' + id
    var cell_input = document.getElementById(evaluated_cell_id);
    input = cell_input.value;
    input = escape(input);
    input = input.replace(/\+/g,"__plus__");
    cell_set_running(id);
    document.getElementById('interrupt').className = 'interrupt';
    async_request('async_obj_evaluate', '/eval' + action, evaluate_cell_callback,
            'id=' + id + '&input='+input)
}


updating = 0;
function evaluate_cell_callback(status, response_text) {
    /* update focus and possibly add a new cell to the end */
    var X = response_text.split(SEP);
    if (X[0] == '-1') {
        /* something went wrong -- i.e., the requested cell doesn't exist. */
        alert("You requested to evaluate a cell that, for some reason, the server is unaware of.");
        return;
    }
    if (X[1] != 'no_new_cell') {
        /* add a new cell to the very end */
       var new_cell = document.createElement("div");
       new_cell.innerHTML = X[1];
       new_cell.id = 'cell_' + X[0];
       var worksheet = document.getElementById('worksheet_cell_list');
       worksheet.appendChild(new_cell);
    }
    if (last_action != 2) {
       focus_on(X[0],X[2]);
    }

    /* set check-for-updates process in progress */
    if (!updating) {
        updating = 1;
        check_for_cell_output();
    }
}


function cell_output_click(id, event) {
    var cell_div = document.getElementById('cell_div_output_' + id);
    var cell_output = document.getElementById('cell_output_' + id);
    var cell_output_nowrap = document.getElementById('cell_output_nowrap_' + id);

    if (cell_div.className == 'cell_output_hidden') {
        cell_div.className = 'cell_output';
        cell_output.className = 'cell_output';
    } else if (cell_div.className == 'cell_output' && event.layerX <= 50) {
        if (cell_output_nowrap.className == 'cell_output_nowrap') {
             cell_output_nowrap.className = 'cell_output_nowrap_visible';
             cell_output.className = 'cell_output_nowrap';
        } else {
             cell_output_nowrap.className = 'cell_output_nowrap';
             cell_output.className = 'cell_output';
             cell_div.className = 'cell_output_hidden';
             cell_output.className = 'cell_output_hidden';
        }
    }
}

function cell_set_running(id) {
    var cell_div = document.getElementById('cell_div_output_' + id)
    cell_div.className = 'cell_output_running';
}

function cell_set_done(id) {
    var cell_div = document.getElementById('cell_div_output_' + id)
    cell_div.className = 'cell_output';
}

function check_for_cell_output() {
    async_request('async_obj_check', '/update_cells',
                    update_cell_output, 'worksheet_id='+worksheet_id)
}

function set_output_text(id, text, wrapped_text, status) {
    /* fill in output text got so far */
    var cell_output = document.getElementById('cell_output_' + id);
    var cell_output_nowrap = document.getElementById('cell_output_nowrap_' + id);
    cell_output.className = 'cell_output';
    cell_output.innerHTML = wrapped_text;
    cell_output_nowrap.innerHTML = text;

    if (status == 'd') {
         cell_set_done(id);
    } else {
         cell_set_running(id);
    }
}

function set_variable_list(variables) {
    var varlist = document.getElementById('variable_list');
    varlist.innerHTML = variables;
}

function set_object_list(objects) {
    var objlist = document.getElementById('object_list');
    objlist.innerHTML = objects;
}

function update_cell_output(status, response_text) {
    if (status == "success") {
        if (response_text == 'empty') {
            /* done -- nothing being computed, since queue is empty */
            updating = 0;
            document.getElementById('interrupt').className = 'interrupt_grey';
        } else {
            /* computing output for a cell */
            i = response_text.indexOf(' ');
            id = response_text.substring(1, i);

            D = response_text.slice(i+1).split(SEP);
            output_text = D[0];
            output_text_wrapped = D[1];
            stat = response_text.charAt(0)
            set_output_text(id, output_text, output_text_wrapped, stat);


            if (stat == 'd') {
                variable_list = D[2];
                set_variable_list(variable_list);
                object_list = D[3];
                set_object_list(object_list);
            }

            /* wait for output from next cell */
            setTimeout('check_for_cell_output()', cell_output_delta);
        }
    } else {
        /* error message ? */
    }
}

///////////////////////////////////////////////////////////////////
//  Insert and move cells
///////////////////////////////////////////////////////////////////


function insert_new_cell_before_callback(status, response_text) {
    if (status == "failure") {
        alert(response_text);
        return ;
    }
    /* Insert a new cell _before_ a cell. */
    var X = response_text.split(SEP);
    var new_cell = document.createElement("div");
    new_cell.id = 'cell_' + X[0];
    new_cell.innerHTML = X[1];
    var cell = document.getElementById('cell_' + X[2]);
    var worksheet = document.getElementById('worksheet_cell_list');
    worksheet.insertBefore(new_cell, cell);
    focus_on(X[0]);
}

function insert_new_cell_before(id) {
    async_request('async_obj_add_new_cell', '/new_cell',
                   insert_new_cell_before_callback, 'id='+id);
}


///////////////////////////////////////////////////////////////////
//
// INTROSPECTION functions -- for getting help
//
///////////////////////////////////////////////////////////////////

function inspect_variable(name) {
/*
    alert(name);
*/
}

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
    async_request('async_obj_search', '/search', callback, 'query='+s)
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
// CONTROL functions
//
///////////////////////////////////////////////////////////////////

function interrupt_callback(status, response_text) {
    if (response_text == "restart") {
        alert("The SAGE kernel had to be restarted (your variables are no longer defined).");
    }
    var link = document.getElementById("interrupt");
    link.className = "interrupt";
    link.innerHTML = "Interrupt"
}

function interrupt() {
    var link = document.getElementById("interrupt");
    if (link.className == "interrupt_grey") {
        return;
    }
    link.className = "interrupt_in_progress";
    link.innerHTML = "Interrupt!"
    async_request('async_obj_interrupt', '/interrupt',
                  interrupt_callback, 'worksheet_id='+worksheet_id);
}

function evaluate_all_callback(status, response_text) {
    /* Iterate through every input cell in the document, in order,
       and call the evaluate_cell function on it.
    */
    if (status == "success") {
        var v, i;
        v = response_text.split(' ');
        for(i=0; i<v.length; i++) {
            evaluate_cell(v[i],0);
        }
    } else {
       alert(response_text);
    }
}

function evaluate_all() {
    /* Use async request back to the server to find out the
       *ordered* list of active cells in the current worksheet.
    */
    async_request('async_obj_evaluate_all', '/cell_id_list',
                      evaluate_all_callback, 'worksheet_id='+worksheet_id);
}

function hide_all_callback(status, response_text) {
    /* Iterate through every input cell in the document, in order,
       and hide it.
    */
    if (status == "success") {
        var v = response_text.split(' ');
        var i;
        for(i=0; i<v.length; i++) {
           var id = v[i];
           var cell_div = document.getElementById('cell_div_output_' + id);
           var cell_output = document.getElementById('cell_output_' + id);
           var cell_output_nowrap = document.getElementById('cell_output_nowrap_' + id);
           cell_output_nowrap.className = 'cell_output_nowrap';
           cell_output.className = 'cell_output';
           cell_div.className = 'cell_output_hidden';
           cell_output.className = 'cell_output_hidden';
        }
    }
}

function hide_all() {
    /* Use async request back to the server to find out the
       *ordered* list of active cells in the current worksheet.
    */
    async_request('async_obj_hide_all', '/cell_id_list',
                      hide_all_callback, 'worksheet_id='+worksheet_id);
}

///////////////////////////////////////////////////////////////////
//
// HELP Window
//
///////////////////////////////////////////////////////////////////

function show_help_window() {
    var help = document.getElementById("help_window");
    help.style.display = "block";

}

function hide_help_window() {
    var help = document.getElementById("help_window");
    help.style.display = "none";

}
"""
