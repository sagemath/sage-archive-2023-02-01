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
    s = r"""

///////////////////////////////////////////////////////////////////
//
// Cross-Browser Stuff
//
///////////////////////////////////////////////////////////////////

var userAgent = navigator.userAgent.toLowerCase();
var browser_ie = userAgent.indexOf("msie")!=-1  && userAgent.indexOf("opera")==-1;
var browser_ie5= userAgent.indexOf("msie 5")!=-1

function get_element(id) {
  if(document.getElementById)
    return document.getElementById(id);
  if(document.all)
    return document.all[id];
  if(document.layers)
    return document.layers[id];
}

function get_event(e) {
   return (e==null)?window.event:e;
}

///////////////////////////////////////////////////////////////////
// An AJAX framework for connections back to the
// SAGE server (written by Tom Boothby).
///////////////////////////////////////////////////////////////////

cell_output_delta = 350;

SEP = '___S_A_G_E___';


var current_cell;
var asyncObj;
var no_async = false;

function getAsyncObject(handler) {
  asyncObj=null
  try {
    if (browser_ie) {
      var s =browser_ie5?"Microsoft.XMLHTTP":"Msxml2.XMLHTTP";
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

function generic_callback(status, response_text) {
   /* do nothing */
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

function show_worksheet_menu(worksheet) {
    // not implemented
}

function set_worksheet_list(worksheets) {
    var wlist = get_element('worksheet_list');
    wlist.innerHTML = worksheets;
}

function show_add_new_worksheet_menu() {
    var add_worksheet_menu = get_element('add_worksheet_menu');
    add_worksheet_menu.style.display = 'block';
    get_element('new_worksheet_box').focus()
}

function hide_add_new_worksheet_menu() {
    var add_worksheet_menu = get_element('add_worksheet_menu');
    add_worksheet_menu.style.display = 'none';
}

function show_upload_worksheet_menu() {
    var upload_worksheet_menu = get_element('upload_worksheet_menu');
    upload_worksheet_menu.style.display = 'block';
    get_element('upload_worksheet_filename').focus()
}


function hide_upload_worksheet_menu() {
    var upload_worksheet_menu = get_element('upload_worksheet_menu');
    upload_worksheet_menu.style.display = 'none';
}

function process_upload_worksheet_menu_submit() {
    hide_upload_worksheet_menu();
    var box = get_element('upload_worksheet_filename');
    var filename = box.value;
    box.value = '';
    upload_worksheet(filename);
}

function upload_worksheet(filename) {
   async_request('async_upload', '/upload_worksheet',
       upload_worksheet_callback, 'filename='+filename)
}

function upload_worksheet_callback(status, response_text) {
    if (status == "success") {
        if (response_text.slice(0,5) == "Error") {
            alert(response_text);
        } else {
            set_worksheet_list(response_text);
        }
    } else {
        alert("Possible problem uploading file (file must be local).");
    }
}

function show_delete_worksheet_menu() {
    var delete_worksheet_menu = get_element('delete_worksheet_menu');
    delete_worksheet_menu.style.display = 'block';
    get_element('delete_worksheet_box').focus();
}

function hide_delete_worksheet_menu() {
    var delete_worksheet_menu = get_element('delete_worksheet_menu');
    delete_worksheet_menu.style.display = 'none';
}

function process_new_worksheet_menu_submit() {
    hide_add_new_worksheet_menu();
    var add_worksheet_box = get_element('new_worksheet_box');
    name = add_worksheet_box.value;
    add_worksheet_box.value = '';
    add_worksheet(name);
}

function process_delete_worksheet_menu_submit() {
    hide_delete_worksheet_menu();
    var delete_worksheet_box = get_element('delete_worksheet_box');
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

function focus(id) {
       // make_cell_input_active(id);
       var cell = get_element('cell_input_' + id);
       if (cell && cell.focus) {
          cell.focus();
       }
}

function scroll_view(id) {
       var cell = get_element('cell_input_' + id);
       if (cell && cell.focus) {
          cell.scrollIntoView();
       }
}

function cell_input_resize(cell_input) {
   try {
        var rows = cell_input.value.split('\n').length - 1;
        if (rows <= 1) {
          rows = 2;
        } else {
          /* to avoid bottom chop off */
          rows = rows + 1;
        }
        cell_input.style.height = null;
        /* cell_input.style.height = 1.5*rows + 'em'; */   // this sort of works in konqueror...
        cell_input.rows = rows;
   } catch(e) {
   } finally {}
}

function lstrip(s) {
    var n = s.length;
    var i = 0;
    while (i < n && (s[i] == ' ' || s[i] == '\n' || s[i] == '\t')) {
        i = i + 1;
    }
    return s.slice(i);
}

function cell_input_minimize_size(cell_input) {
    var v = cell_input.value;
    var w = lstrip(v);
    var sl = w.slice(0,5);
    if (sl == '%hide') {
        cell_input.className = 'cell_input_hide';
        return;
    }

    cell_input.className = 'cell_input';
    var rows = v.split('\n').length ;
    if (rows < 1) {
      rows = 1;
    }
    cell_input.rows = rows;
    if (rows == 1) {
       // hack because of bug in firefox with 1-row textarea
       cell_input.style.height = '1.5em';
    }
}

function cell_input_minimize_all() {
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        var cell=get_element('cell_input_' + v[i]);
        cell_input_minimize_size(cell);
    }
}

id_to_delete=-1;
function cell_delete_callback(status, response_text) {
    if (status == "failure") {
        cell = get_element('cell_' + id_to_delete);
        var worksheet = get_element('worksheet_cell_list');
        worksheet.removeChild(cell);
    }
    var X = response_text.split(SEP);
    if (X[0] == 'ignore') {
        return;   /* do not delete, for some reason */
    }
    var cell = get_element('cell_' + X[1]);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.removeChild(cell);
    cell_id_list = eval(X[3]);
    var id_before = X[2];
    focus(eval(id_before), 0);
}

function cell_delete(id) {
   id_to_delete=id;
   async_request('async_obj_cell_delete', '/delete_cell', cell_delete_callback, 'id='+id)
}

function cell_input_key_event(id, event) {
    cell_input = get_element('cell_input_'+id);

    e = get_event(event);

    //alert (e.keyCode);

    if (key_delete_cell(e) && cell_input.value == '') {
        cell_delete(id);
        return false;
    }

    cell_input_resize(cell_input);

    // Will need IE version... if possible.
    if (e.keyCode == 38) {  // up arrow
        var before_cursor = cell_input.value.substr(0,cell_input.selectionEnd);
        var i = before_cursor.indexOf('\n');
        if (i == -1) {
            jump_to_cell(id,-1);
            return false;
        } else {
            return true;
        }
    } else if (e.keyCode == 40) {   // down arrow
        var before_cursor = cell_input.value.substr(cell_input.selectionEnd);
        var i = before_cursor.indexOf('\n');
        if (i == -1) {
            jump_to_cell(id,1);
            return false;
        } else {
            return true;
        }
    } else if (key_send_input(e)) {
           // User pressed shift-enter
       evaluate_cell(id, 0);
       return false;
    } else if (key_send_input_newcell(e)) {
      /* evaluate_cell(id, 1); */
       evaluate_cell(id, 1);
       insert_new_cell_after(id);
       return false;
    } else if (key_request_introspections(e)) {
       // command introspection
       evaluate_cell(id, 2);
       return false;
    } else if (key_interrupt(e)) {
       interrupt();
       return false;
    } else if (key_page_down(e)) {
       jump_to_cell(id, 5);
       return false;
    } else if (key_page_up(e)) {
       jump_to_cell(id, -5);
       return false;
    } else if (key_request_history(e)) {
       history_window();
    } else if (key_request_log(e)) {
       text_log_window(worksheet_filename);
    }
    return true;
}

function id_of_cell_delta(id, delta) {
    if (cell_id_list.length == 0) {
        alert("bug -- no cells.");
        return;
    }
    var i = array_indexOf(cell_id_list, eval(id));
    var new_id;
    if (i == -1) {
        return(cell_id_list[0]);
    } else {
        i = i + delta;
        if (i < 0) {
            i = 0;
        } else if (i >= cell_id_list.length) {
            i = cell_id_list.length - 1;
        }
        return(cell_id_list[i]);
    }
}

/*
old_id = -1;
function make_cell_input_active(id) {
   if (old_id != -1) {
        make_cell_input_inactive(old_id);
   }
   var txt = get_element('cell_input_'+id);
   if (txt.style.display == "inline") {
       return;
   }
   cell_input_resize(txt);
   current_cell = id;
   txt.style.display = "inline";

   var pre = get_element('cell_input_pre_'+id);
   pre.style.display = "none";
   txt.value = pre.innerHTML;
   pre.innerHTML = '';
}

function make_cell_input_inactive(id) {
   var txt = get_element('cell_input_'+id);
   if (txt.style.display != "inline") {
       return;
   }

   txt.style.display = "none";


   var pre = get_element('cell_input_pre_'+id);
   pre.style.display = "inline";
   pre.innerHTML = txt.value;
   txt.value = '';
}
*/

function jump_to_cell(id, delta) {
    var i = id_of_cell_delta(id, delta)
    var j = id_of_cell_delta(i, -2);
    if (j >= i) {
        j = i;
    }
    scroll_view(j);
    focus(i);
}

function escape0(input) {
    input = escape(input);
    input = input.replace(/\+/g,"__plus__");
    return input;
}

var non_word = new RegExp("[^a-zA-Z0-9._]");
var evaluated_cell_id = 0;
var cell_id = 0;
var last_action = 0;
function evaluate_cell(id, action) {
    cell_id = id;
    last_action = action;
    var cell_input = get_element('cell_input_' + id);

    I = cell_input.value

    if(action == 2) { // Introspection
       evaluate_cell_introspection(id);
       return;
    }

    input = escape0(I);

    cell_set_running(id);
    get_element('interrupt').className = 'interrupt';
    async_request('async_obj_evaluate', '/eval' + action, evaluate_cell_callback,
            'id=' + id + '&input='+input)
}

function evaluate_cell_introspection(id) {
    var before_cursor;
    var after_cursor;
    var cell_input = get_element('cell_input_' + id);
    if(browser_ie) {
        //for explorer, we call the
        //   document.selection.createRange().duplicate()
        //to generate a selection range object which does not effect the
        //original input box.
        //Then, we rewind the start point of the range until we encounter
        //a non-word character, or we've rewound past the beginning of
        //the textarea).
        range = document.selection.createRange().duplicate();
        var i = range.text.length;
        while((cell_input.value.match(range.text) || i==0) && !non_word.exec(range.text)) {
            range.moveStart('character', -1);
            i = i + 1;
    }
        before_cursor = range.text;
        /* TODO -- total guess -- stein */
        after_cursor = cell_input.value.substr(i);
    } else if(cell_input.selectionEnd) {
        before_cursor = cell_input.value.substr(0,cell_input.selectionEnd);
        after_cursor = cell_input.value.substr(cell_input.selectionEnd);
    }

    // If the character right before the cursor is blank, we instead
    // send 4 spaces.
    if (is_just_a_tab(before_cursor, id)) {
        cell_input.value = before_cursor + '    ' + after_cursor;
        return;
    }
    before_cursor = escape0(before_cursor);
    after_cursor = escape0(after_cursor);
    async_request('async_obj_evaluate', '/introspect', evaluate_cell_callback,
          'id=' + id + '&before_cursor='+before_cursor + '&after_cursor='+after_cursor)
}

function is_just_a_tab(s, id) {
    var n = s.length
    if (n==0)
       return 1;
    c = s.charAt(n-1);
    if (c == ' ' || c == '\t' || c == '\n' || c == ')') {
       return 1;
    }
    return 0;
}


updating = 0;
function evaluate_cell_callback(status, response_text) {
    /* update focus and possibly add a new cell to the end */
    if (status == "failure") {
        alert(response_text);
        return;
    }
    var X = response_text.split(SEP);
    if (X[0] == '-1') {
        /* something went wrong -- i.e., the requested cell doesn't exist. */
        alert("You requested to evaluate a cell that, for some reason, the server is unaware of.");
        return;
    }
    if (X[1] == 'append_new_cell') {
        /* add a new cell to the very end */
        append_new_cell(X[0],X[2]);
    } else if (X[1] == 'insert_cell') {
        /* insert a new cell after the one with id X[3] */
        do_insert_new_cell_after(X[3], X[2], X[0]);
        jump_to_cell(X[0], 0);
    } else if (last_action != 2) {  /* not a introspection */
       focus(X[0]);
    }

    /* set check-for-updates process in progress */
    if (!updating) {
        updating = 1;
        check_for_cell_output();
    }
}

function cell_output_set_type(id, typ) {
    var cell_div = get_element('cell_div_output_' + id);
    var cell_output = get_element('cell_output_' + id);
    var cell_output_nowrap = get_element('cell_output_nowrap_' + id);
    var cell_output_html = get_element('cell_output_html_' + id);

    cell_div.className = 'cell_output_' + typ;
    cell_output.className = 'cell_output_' + typ;
    cell_output_nowrap.className ='cell_output_nowrap_' + typ;
    cell_output_html.className   ='cell_output_html_' + typ;

    /* Do async request back to the server */
    async_request('async_obj_check', '/cell_output_set',
                    generic_callback, 'id='+id+'&type=' + typ)
}

function cycle_cell_output_type(id) {
    var cell_div = get_element('cell_div_output_' + id);

    if (cell_div.className == 'cell_output_hidden' || cell_div.className=='cell_output_running') {
        cell_output_set_type(id, 'wrap');
        return;
    }

    if (cell_div.className == 'cell_output_wrap') {
        cell_output_set_type(id, 'nowrap');
    } else {
        cell_output_set_type(id, 'hidden');
    }
}

function cell_set_running(id) {
    cell_output_set_type(id, 'wrap');
    var cell_div = get_element('cell_div_output_' + id)
    cell_div.className = 'cell_output_running';
}

function cell_set_done(id) {
    var cell_div = get_element('cell_div_output_' + id)
    cell_div.className = 'cell_output_wrap';
}

function check_for_cell_output() {
    async_request('async_obj_check', '/update_cells',
                    update_cell_output, 'worksheet_id='+worksheet_id)
}

function set_output_text(id, text, wrapped_text, output_html, status) {
    /* fill in output text got so far */
    var cell_output = get_element('cell_output_' + id);
    var cell_output_nowrap = get_element('cell_output_nowrap_' + id);
    var cell_output_html = get_element('cell_output_html_' + id);

    cell_output.innerHTML = wrapped_text;
    cell_output_nowrap.innerHTML = text;
    cell_output_html.innerHTML = output_html;

    if (status == 'd') {
         cell_set_done(id);
    } else {
         cell_set_running(id);
    }
}

function set_input_text(id, text) {
    /* fill in input text */
    var cell_input = get_element('cell_input_' + id);
    cell_input.value = text;
}

function set_variable_list(variables) {
    var varlist = get_element('variable_list');
    varlist.innerHTML = variables;
}

function set_object_list(objects) {
    var objlist = get_element('object_list');
    objlist.innerHTML = objects;
}

function set_attached_files_list(objects) {
    var objlist = get_element('attached_list');
    objlist.innerHTML = objects;
}

function update_cell_output(status, response_text) {
    if (status == "success") {
        if (response_text == 'empty') {
            /* done -- nothing being computed, since queue is empty */
            updating = 0;
            /* get_element('interrupt').className = 'interrupt_grey'; */
        } else {
            /* computing output for a cell */
            i = response_text.indexOf(' ');
            id = response_text.substring(1, i);

            D = response_text.slice(i+1).split(SEP);
            output_text = D[0] + ' ';
            output_text_wrapped = D[1] + ' ';
            output_html = D[2];
            stat = response_text.charAt(0)
            set_output_text(id, output_text, output_text_wrapped, output_html, stat);
            var j = id_of_cell_delta(id,1);

            if (stat == 'd') {
                new_cell_input = D[3];
                if (new_cell_input != '') {
                    set_input_text(id, new_cell_input);
                }
                variable_list = D[4];
                set_variable_list(variable_list);

                object_list = D[5];
                set_object_list(object_list);

                attached_files_list = D[6];
                set_attached_files_list(attached_files_list);
            }

            /* wait for output from next cell */
            setTimeout('check_for_cell_output()', cell_output_delta);
        }
    } else {
        alert("Error updating cell output: " + response_text);
    }
}

///////////////////////////////////////////////////////////////////
//  Insert and move cells
///////////////////////////////////////////////////////////////////

function insertAfter(parent, node, referenceNode) {
	parent.insertBefore(node, referenceNode.nextSibling);
}

function insert_into_array(v, i, x) {
    /* Return a new array with x inserted into position i of v. */
    return (v.slice(0,i).concat([x]).concat(v.slice(i,v.length)));
}

// Array.indexOf( value, begin, strict ) - Return index of the first element that matches value
function array_indexOf(v, x) {
    var len = v.length;
    for(var i=0; i < len; i++) {
         if(v[i] == x) {
             return i;
         }
    }
    return -1;
};

function do_insert_new_cell_before(id, new_html, new_id) {
  /* Insert a new cell with the given new_id and new_html
     before the cell with given id. */
    var new_cell = document.createElement("div");
    new_cell.id = 'cell_' + new_id;
    new_cell.innerHTML = new_html;
    var cell = get_element('cell_' + id);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.insertBefore(new_cell, cell);
    var i = array_indexOf(cell_id_list, eval(id));
    cell_id_list = insert_into_array(cell_id_list, i, eval(new_id));
}

function do_insert_new_cell_after(id, new_html, new_id) {
  /* Insert a new cell with the given new_id and new_html
     after the cell with given id. */
    var new_cell = document.createElement("div");
    new_cell.id = 'cell_' + new_id;
    new_cell.innerHTML = new_html;
    var cell = get_element('cell_' + id);
    var worksheet = get_element('worksheet_cell_list');
    insertAfter(worksheet, new_cell, cell);
    var i = array_indexOf(cell_id_list, eval(id));
    cell_id_list = insert_into_array(cell_id_list, i+1, eval(new_id));
}

function insert_new_cell_before_callback(status, response_text) {
    if (status == "failure") {
        alert(response_text);
        return ;
    }
    /* Insert a new cell _before_ a cell. */
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);
    do_insert_new_cell_before(id, new_html, new_id);
    focus(new_id);
}

function insert_new_cell_before(id) {
    async_request('async_obj_add_new_cell', '/new_cell',
                   insert_new_cell_before_callback, 'id='+id);
}

function insert_new_cell_after_callback(status, response_text) {
    if (status == "failure") {
        alert(response_text);
        return ;
    }
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);
    do_insert_new_cell_after(id, new_html, new_id);
    jump_to_cell(new_id,0);
}

function insert_new_cell_after(id) {
    async_request('async_obj_add_new_cell', '/new_cell_after',
                   insert_new_cell_after_callback, 'id='+id);
}

function append_new_cell(id, html) {
    var new_cell = document.createElement("div");
    new_cell.innerHTML = html;
    new_cell.id = 'cell_' + id;
    var worksheet = get_element('worksheet_cell_list');
    worksheet.appendChild(new_cell);
    cell_id_list = cell_id_list.concat([id]);
    jump_to_cell(id, 0);
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
    var link = get_element("interrupt");
    link.className = "interrupt";
    link.innerHTML = "Interrupt"
}

function interrupt() {
    var link = get_element("interrupt");
    if (link.className == "interrupt_grey") {
        return;
    }
    link.className = "interrupt_in_progress";
    link.innerHTML = "Interrupt!"
    async_request('async_obj_interrupt', '/interrupt',
                  interrupt_callback, 'worksheet_id='+worksheet_id);
}


function evaluate_all() {
    /* Iterate through every input cell in the document, in order,
       and evaluate.
    */
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        evaluate_cell(v[i],0);
    }
}

function hide_all() {
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        cell_output_set_type(v[i],'hidden');
    }
    async_request('async_obj_hide_all', '/hide_all',
                      hide_all_callback, 'worksheet_id='+worksheet_id);
}

function show_all() {
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        cell_output_set_type(v[i],'wrap');
    }
    async_request('async_obj_hide_all', '/show_all',
                      hide_all_callback, 'worksheet_id='+worksheet_id);
}

function login(username,password) {
  document.cookie="username="+username;
  document.cookie="password="+password;
  window.location="/";
}

///////////////////////////////////////////////////////////////////
//
// HELP Window
//
///////////////////////////////////////////////////////////////////

function show_help_window() {
    var help = get_element("help_window");
    help.style.display = "block";

}

function hide_help_window() {
    var help = get_element("help_window");
    help.style.display = "none";

}

///////////////////////////////////////////////////////////////////
//
// LOG windows
//
///////////////////////////////////////////////////////////////////
function history_window() {
    history = window.open ("__history__.html",
      "", "menubar=1,scrollbars=1,width=700,height=600");
}

function worksheet_text_window(worksheet) {
    log = window.open (worksheet+"__plain__.html","",
      "menubar=1,scrollbars=1,width=700,height=600");
}

function doctest_window(worksheet) {
    log = window.open (worksheet+"__doc__.html","",
      "menubar=1,scrollbars=1,width=700,height=600");
}


function print_window(worksheet) {
    log = window.open (worksheet+"__print__.html","",
      "menubar=1,scrollbars=1,width=700,height=600,toolbar=0");
}

//////////////////////////////////
// HELP
/////////////////////////////////
function show_help_window(worksheet) {
    help = window.open ("__help__.html","",
    "menubar=1,scrollbars=1,width=700,height=600");
}

"""

    s += r"""
///////////////////////////////////////////////////////////////////
//
// KeyCodes (auto-generated from config.py and user's sage config
//
///////////////////////////////////////////////////////////////////
%s
"""%build_all_key_codes()

    return s


key_codes = {};
class JSKeyCode:
    def __init__(self, name, key='', alt=False, ctrl=False, shift=False):
        global key_codes
        self.key   = key
        self.alt   = alt
        self.ctrl  = ctrl
        self.shift = shift
        key_codes[name] = self

    def js_test(self, name):
        a = "%s"%self.alt
        c = "%s"%self.ctrl
        s = "%s"%self.shift

        return """
function key_%s(e) {
  return e.keyCode == %d && e.altKey == %s && e.ctrlKey == %s && e.shiftKey == %s;
}
        """%(name, self.key, a.lower(), c.lower(), s.lower())

def build_all_key_codes():
    return ''.join([k.js_test(n) for n, k in key_codes.items()])

