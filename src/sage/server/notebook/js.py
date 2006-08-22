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

from sage.misc.misc import SAGE_URL
import keyboards

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
// GLOBAL VARIABLES
//
// PLEASE define all global variables up here, and if you want to
// set them with anything but simple assignment, 'var' them first,
// and set them later.  Your code might work in your browser, but
// it might break initial setup for other critical pieces in  other
// browsers.  Thanks. (and for the record, I'm guilty of this more
// than anybody else here -- I figure a big block comment might
// help keep me in check)
//
// Exception: keyboard globals are defined at the end
//
///////////////////////////////////////////////////////////////////


// The active cell list.
var active_cell_list = [];

//Browser & OS identification
var browser_op, browser_saf, browser_konq, browser_moz, browser_ie, browser_ie5;
var os_mac, os_lin, os_win;

var update_error_count = 0;
var update_error_threshold = 30;

// in milliseconds
var update_error_delta = 1000;
var update_normal_delta = 256
var cell_output_delta = update_normal_delta;

var SEP = '___S_A_G_E___';   // this had better be the same as in the server
var current_cell = -1;       // gets set on focus / blur
var asyncObj;
var no_async = false; //this isn't really used -- should we think about dealing with this?

// introspection variables
var introspection_loaded = false;
var introspect_id;
var introspection_text = "";
var replacement_text = "";
var replacement_row = 0;
var replacement_col = 0;
var replacing_word = "";
var replacement_word = "";
var replacing = false;
var sub_introspecting = false;

var worksheet_id=0;   // The current worksheet.  Where else does this get set?

var id_to_delete=-1;

//regular expressions used to peek into the cell input for introspection
var non_word = "[^a-zA-Z0-9_]"; //finds any character that doesn't belong in a variable name
var command_pat = "([a-zA-Z_][a-zA-Z._0-9]*)$"; //identifies the command at the end of a string
var function_pat = "([a-zA-Z_][a-zA-Z._0-9]*)\\([^()]*$";
var one_word_pat = "([a-zA-Z_][a-zA-Z._0-9]*)";
try{
  non_word = new RegExp(non_word);
  command_pat = new RegExp(command_pat);
  function_pat = new RegExp(function_pat);
  one_word_pat = new RegExp(one_word_pat);
} catch(e){}

var after_cursor, before_cursor, before_replacing_word;

var update_timeout = -1;

var updating = false;
var update_time = -1;

var jsmath_font_msg = '<a href="SAGE_URL/jsmath">jsMath temporarily disabled while we resolve a windows firefox hang bug</a><br>';
/*var jsmath_font_msg = '<a href="SAGE_URL/jsmath">Click to download and install tex fonts.</a><br>'; */

var cell_id_list; // this gets set in worksheet.py

var input_keypress; //this gets set to a function when we set up the keyboards

var in_slide_mode = false; //whether or not we're in slideshow mode

///////////////////////////////////////////////////////////////////
//
// Cross-Browser Stuff
//
///////////////////////////////////////////////////////////////////
function true_function() {return true;}
input_keypress = true_function;

try{
  var n=navigator;
  var nav=n.appVersion;
  var nan=n.appName;
  var nua=n.userAgent;
  browser_op=(nua.indexOf('Opera')!=-1);
  browser_saf=(nua.indexOf('Safari')!=-1);
  browser_konq=(!browser_saf && (nua.indexOf('Konqueror')!=-1) ) ? true : false;
  browser_moz=( (!browser_saf && !browser_konq ) && ( nua.indexOf('Gecko')!=-1 ) ) ? true : false;
  browser_ie=((nua.indexOf('MSIE')!=-1)&&!browser_op);
  browser_ie5=(browser_ie&&(nua.indexOF('MSIE 5')!=-1));
  os_mac=(nav.indexOf('Mac')!=-1);
  os_win=( ( (nav.indexOf('Win')!=-1) || (nav.indexOf('NT')!=-1) ) && !os_mac)?true:false;
  os_lin=(nua.indexOf('Linux')!=-1);

  get_keyboard();
} catch(e){}

function get_keyboard() {
  var b,o,warn=false;

  input_keypress = cell_input_key_event;
  debug_keypress = debug_input_key_event;

  if(browser_op) {
    b = "o";
  } else if(browser_ie) {
    b = "i";
    alert("listening onkeydown");
    document.onkeydown = key_listen_ie;
    input_keypress = true_function;
    debug_keypress = true_function;
    warn = true;
  } else if(browser_saf) {
    b = "s";
  } else if(browser_konq) {
    b = "k";
    warn = true;
  } else {
    b = "m";
  }

  if(os_mac) {
    o = "m";
  } else if(os_lin) {
    o = "l";
  } else {
    o = "w"
  }

  if(b == null || o == null || warn) {
    alert("Your browser / OS combination is not supported.  \nPlease use Firefox or Opera under linux, windows, or mac OSX, or Safari.")
  }

  async_request('keyboard', '__keyboard_'+b+o+'__.js', get_keyboard_callback, null);
}

function get_keyboard_callback(status, response) {
  if(status == 'success') {
    eval(response);
  }
}

function get_element(id) {
  if(document.getElementById)
    return document.getElementById(id);
  if(document.all)
    return document.all[id];
  if(document.layers)
    return document.layers[id];
}

function set_class(id, cname) {
  e = get_element(id);
  if(e!=null) {
      e.className = cname;
  }
}

function get_event(e) {
   return (e==null)?window.event:e;
}

function key_event(e) {
   if(e==null) e = window.event;
   if(e.modifiers) {
     this.a = e.modifiers | 1;
     this.c = e.modifiers | 2;
     this.s = e.modifiers | 4;
   } else {
     this.a = e.altKey;
     this.c = e.ctrlKey;
     this.s = e.shiftKey;
   }
   this.k = e.keyCode + "," + e.which;
   this.m = this.k + (this.s?'!':'');
   return this;
}

function time_now() {
  return (new Date()).getTime();
}

///////////////////////////////////////////////////////////////////
// An AJAX framework for connections back to the
// SAGE server (written by Tom Boothby).
///////////////////////////////////////////////////////////////////

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
// Misc page functions -- for making the page work nicely
// (this is a crappy descriptor)
///////////////////////////////////////////////////////////////////

function trim(s) {
    m = one_word_pat.exec(s);
    if(m == null)
        return s;
    return m[1];
}

function body_load() {
//  init_menus();
}

function init_menus() {
  for( i = 1; i <= 3; i++) {
    menu = get_element("menu"+i);
    menu.style.display="none";
  }
}

///////////////////////////////////////////////////////////////////
//
// Completions interface stuff
//
///////////////////////////////////////////////////////////////////

function handle_replacement_controls(cell_input, event) {
    deselect_replacement_element();
    if(key_menu_up(event)) {
        if(replacement_row <= 0) {
            halt_introspection();
        } else {
            replacement_row--;
        }
    } else if(key_menu_down(event)) {
        replacement_row++;
        if(!replacement_element_exists())
            replacement_row = 0;
    } else if(key_menu_right(event)) {
        replacement_col++;
        if(!replacement_element_exists())
            replacement_col = 0;
    } else if(key_menu_left(event)) {
        replacement_col--;
        if(!replacement_element_exists()) {
            replacement_col = 1;
            while(replacement_element_exists())
                replacement_col++;
            replacement_col--;
        }
    } else if(key_menu_pick(event)) {
        do_replacement(introspect_id, replacement_word, true);
        return false;
    } else if(key_request_introspections(event)) {
        if(sub_introspecting) {
            introspection_text = replacement_text;
            introspection_loaded = true;
            sub_introspecting = false;
            update_introspection_text();
        } else {
            replacement_text = introspection_text;
            introspection_loaded = false;
            sub_introspecting = true;
        }
    } else {
       halt_introspection();
       return true;
    }
    select_replacement_element();

    if(sub_introspecting) {
        active_cell_list = active_cell_list.concat([introspect_id]);
        evaluate_cell_introspection(introspect_id, before_replacing_word+replacement_word+'?', after_cursor);
    }

    return false;
}

function do_replacement(id, word,do_trim) {
    var cell_input = get_cell(id);

    if(do_trim) //optimization 'cause Opera has a slow regexp engine
        word = trim(word);

    cell_input.value = before_replacing_word + word + after_cursor;
    focus(id);  //reset the cursor (for explorer)

    try{ //firefox, et. al.
        var pos = before_replacing_word.length + word.length;
        cell_input.selectionStart = pos;
        cell_input.selectionEnd = pos;
    } catch(e) {}
    try{ //explorer; anybody else?
        var range = document.selection.createRange();
        range.moveStart('character', -after_cursor.length);
        range.moveEnd('character', -after_cursor.length);
    } catch(e) {}

    if(browser_op || browser_saf)
      focus(id,true);

    halt_introspection();
}

function get_replacement_element() {
    return get_element("completion"+introspect_id + "_" + replacement_row + "_" + replacement_col);
}

function replacement_element_exists() {
    return get_replacement_element() != null;
}

function select_replacement(row, col) {
    deselect_replacement_element();
    replacement_row = row;
    replacement_col = col;
    select_replacement_element();
}

function deselect_replacement_element() {
    e = get_replacement_element();
    if(e==null) return;
    e.className = 'completion_menu_two';
}

function select_replacement_element() {
    var e = get_replacement_element();
    if (e==null) return;
    e.className = 'completion_menu_two completion_menu_selected';
    var l = e.getElementsByTagName('a');
    if(l.length && l[0]!=null) {
        var h = l[0].innerHTML;
        var i = h.indexOf('&nbsp')
        if (i != -1) {
            h = h.substr(0,i);
        }
        replacement_word = h;
    }
}

function update_introspection_text(preserve_cursor) {
  close_introspection_text();
  d = get_element("introspect_div_"+introspect_id);
  if(!d) return;

  if(introspection_loaded) {
    if(introspection_text == "") {
        halt_introspection();
        return;
    }
    d.innerHTML = introspection_text;
    if(replacing)
      select_replacement_element();
  } else {
    d.innerHTML = "loading..."
  }
}

function close_introspection_text() {
  d = get_element("introspect_div_"+introspect_id);
  if(d!=null)
    d.innerHTML = "";
}

function halt_introspection() {
    close_introspection_text();
    introspect_id = null;
    replacing = false;
    sub_introspecting = false;
    introspection_loaded = false;
    replacement_row = replacement_col = 0;
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
            alert("Unable to add worksheet.");
        } else {
            set_worksheet_list(X[0]);
            switch_to_worksheet(X[1]);
        }
    } else {
        alert("Possible failure adding workbook.");
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
            alert("Possible failure deleting workbook.");
        } else {
            set_worksheet_list(X[0]);
            if (X[1] != -1)
               switch_to_worksheet(X[1]);
        }
    } else {
        alert("Possible failure deleting workbook.");
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
    w = window.open("__upload__.html", "upload", "width=500, height=200");
    if(w.focus)
      w.focus();
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
            alert("Error uploading worksheet.");
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

function get_cell(id) {
    return get_element('cell_input_'+ id);
}

function cell_focus(id) {
    e = get_cell(id);
    current_cell = id;
    if(e == null) return;
    e.className="cell_input_active";
    cell_input_resize(e);
    return true;
}
function cell_blur(id) {
    current_cell = -1;
    e = get_cell(id);
    if(e == null) return;
    e.className="cell_input";
    cell_input_minimize_size(e);
    return true;
}

function debug_focus() {
    in_debug_input = true;
    w = get_element('debug_window');
    if(w)
       w.className = 'debug_window_active';
}

function debug_blur() {
    in_debug_input = false;
    w = get_element('debug_window');
    if(w)
        w.className = 'debug_window_inactive';
}

//set and_delay to true if you want to refocus the browser in a keyevent
//which expects a tab -- Opera apparently resists canceling the tab key
//event -- so we can subvert that by breaking out of the call stack with
//a little timeout.
function focus(id, and_delay) {
    // make_cell_input_active(id);
    if(in_slide_mode) {
        slide_hide();
        current_cell = id;
        slide_show();
    }

    var cell = get_cell(id);
    if (cell && cell.focus) {
        cell.focus();
        if(and_delay) {
            setTimeout('focus('+id+',false)', 10);
        }
    }

}

function cell_input_resize(cell_input) {
    var rows = 2;
    //var rows = cell_input.value.split('\n').length - 1;
    var rows = cell_input.value.split('\n').length + 1;
    if (rows <= 1) {
      rows = 2;
    } else {
      /* to avoid bottom chop off */
      rows = rows + 1;
    }
    try {
        cell_input.style.height = rows + 'em';   // this sort of works in konqueror...
    } catch(e) {}
    try{
        cell_input.rows = rows;
    } catch(e) {}
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
        cell_input.style.height = '1.5em';
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
        var cell=get_cell(v[i]);
        cell_input_minimize_size(cell);
    }
}

function cell_delete_callback(status, response_text) {
    if (status == "failure") {
        cell = get_element('cell_outer_' + id_to_delete);
        var worksheet = get_element('worksheet_cell_list');
        worksheet.removeChild(cell);
        jump_to_cell(id_to_delete,-1);
        cell_id_list = delete_from_array(cell_id_list, id_to_delete);
        id_to_delete = -1;
        return;
    }
    var X = response_text.split(SEP);
    if (X[0] == 'ignore') {
        id_to_delete = -1;
        return;   /* do not delete, for some reason */
    }
    var cell = get_element('cell_outer_' + id_to_delete);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.removeChild(cell);
    jump_to_cell(id_to_delete,-1);
    cell_id_list = delete_from_array(cell_id_list, id_to_delete);
    id_to_delete = -1;
}

function cell_delete(id) {
   id_to_delete=id;
   async_request('async_obj_cell_delete', '/delete_cell', cell_delete_callback, 'id='+id)
}

function key_listen_ie() {
    var e = get_event(null);
    if(current_cell != -1) {
        k = new key_event(e);
        if(key_shift(k) || key_ctrl(k) || key_alt(k))
          return true;

        if(!cell_input_key_event(current_cell, e)) {
            void(0);
            e.returnValue=false;
            e.cancelBubble=true;
            return false;
        }
        return true;
    }
    if(in_debug_input) {
        if(!debug_input_key_event(e)) {
            void(0);
            e.returnValue=false;
            e.cancelBubble=true;
            return false;
        }
        return true;
    }
    return true;
}

function debug_input_key_event(e) {
    e = new key_event(e);
    debug_input = get_element('debug_input');

    if (key_down_arrow(e)) {
        var after = text_cursor_split(debug_input)[1];
        var i = after.indexOf('\n');
        if (i == -1 || after == '') {
            focus(cell_id_list[0])
            return false;
        } else {
            return true;
        }
    }
    if (key_send_input(e)) {
        var out = ""
        try {
          out = eval(debug_input.value);
        } catch(err) {
          out = "Error: " + err.description;
        } finally {
          debug_append(out);
          return false;
        }
    }
}

function cell_input_key_event(id, e) {
    cell_input = get_cell(id);
    e = new key_event(e);
    if (e==null) return;

    //alert (e.keyCode);

    if (key_delete_cell(e) && cell_input.value == '') {
        cell_delete(id);
        return false;
    }

    if(introspect_id && introspection_loaded && replacing) {
        if(!handle_replacement_controls(cell_input, e)) {
            if(browser_op) focus(id,true);
            return false;  //otherwise, keep going
        }
        halt_introspection();
    }

    cell_input_resize(cell_input);

    // Will need IE version... if possible.
    if (!in_slide_mode && key_up_arrow(e)) {
        var before = text_cursor_split(cell_input)[0];
        var i = before.indexOf('\n');
        if (i == -1 || before == '') {
            jump_to_cell(id,-1);
            return false;
        } else {
            return true;
        }
    } else if (!in_slide_mode && key_down_arrow(e)) {
        var after = text_cursor_split(cell_input)[1];
        var i = after.indexOf('\n');
        if (i == -1 || after == '') {
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
       return false;
    } else if (key_request_introspections(e)) {
       // command introspection (tab completion, ?, ??)
       evaluate_cell(id, 2);
       return false;
    } else if (key_interrupt(e)) {
       interrupt();
       return false;
    } else if (key_page_down(e)) {
       if(in_slide_mode) {
           jump_to_cell(id, 1);
       } else {
           jump_to_cell(id, 5);
       }
       return false;
    } else if (key_page_up(e)) {
       if(in_slide_mode) {
           jump_to_cell(id, -1);
       } else {
           jump_to_cell(id, -5);
       }
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
        return(id);  /* Better not to move. */
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

function debug_clear() {
    output = get_element("debug_output");
    if(output == null) return;
    output.value = "";
}

function debug_append(txt) {
    output = get_element("debug_output");
    if(output == null) return;
    output.value = txt + "\n" + output.value;
}

/*
old_id = -1;
function make_cell_input_active(id) {
   if (old_id != -1) {
        make_cell_input_inactive(old_id);
   }
   var txt = get_cell(id);
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
   var txt = get_cell(id);
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
    focus(i);
}

function escape0(input) {
    input = escape(input);
    input = input.replace(/\+/g,"__plus__");
    return input;
}

function text_cursor_split(input) {
    if(browser_ie) {
        //for explorer, we call the
        //   document.selection.createRange().duplicate()
        //to generate a selection range object which does not effect the
        //original input box.
        //Then, we rewind the start point of the range until we encounter
        //a non-word character, or we've rewound past the beginning of
        //the textarea).
        var range = document.selection.createRange().duplicate();
        var i = range.text.length;
        while((input.value.match(range.text) || i==0)
               && range.text.length == i) {
            range.moveStart('character', -1);
            i = i + 1;
        }
        if(!input.value.match(range.text))
            range.moveStart('character', 1);
        b = range.text;
    } else {
        b = input.value.substr(0,input.selectionEnd);
    }

    a = input.value.substr(b.length);
    return new Array(b,a);
}

function evaluate_cell(id, action) {
    active_cell_list = active_cell_list.concat([id]);

    if(action == 2) { // Introspection
       evaluate_cell_introspection(id,null,null);
       return;
    }

    if(!in_slide_mode) {
        jump_to_cell(id,1);
    }
    cell_set_running(id);
    var cell_input = get_cell(id);
    var I = cell_input.value;
    var input = escape0(I);

    get_element('interrupt').className = 'interrupt';
    async_request('async_obj_evaluate', '/eval' + action, evaluate_cell_callback,
            'id=' + id + '&input='+input)
}

function evaluate_cell_introspection(id, before, after) {
    var cell_input = get_cell(id);

    replacing = false;
    if(before == null) {
        var in_text = text_cursor_split(cell_input);
        before_cursor = before = in_text[0];
        after_cursor  = after  = in_text[1];
        before_replacing_word = before;

        m = command_pat.exec(before);
        f = function_pat.exec(before);
        if(introspect_id != null)
            halt_introspection();
        introspect_id = id;

        var last_char_before = before.charAt(before.length-1);
        if(last_char_before == "?") {
        } else if(m) {
            replacing = true;
            replacing_word  = m[1];
            before_replacing_word = before.substring(0, before.length-replacing_word.length);
        } else if(f != null) { //we're in an open function paren -- give info on the function
            before = f[1] + "?";
        } else { //just a tab
            do_replacement(id, '    ',false);
            return;
        }
    } else {
        sub_introspecting = true;
    }
    if(!replacing && browser_op)
        focus(id,true);

    update_introspection_text();
    var before_cursor_e = escape0(before);
    var after_cursor_e = escape0(after);
    cell_set_running(id);
    async_request('async_obj_evaluate', '/introspect', evaluate_cell_callback,
          'id=' + id + '&before_cursor='+before_cursor_e + '&after_cursor='+after_cursor_e);
}

function evaluate_cell_callback(status, response_text) {
    /* update focus and possibly add a new cell to the end */
    if (status == "failure") {
        alert("Failure evaluating a cell.");
        return;
    }
    var X = response_text.split(SEP);
    if (X[0] == '-1') {
        /* something went wrong -- i.e., the requested cell doesn't exist. */
        alert("You requested to evaluate a cell that, for some reason, the server is unaware of.");
        return;
    }
    if (X[1] == 'append_new_cell') {
        // add a new cell to the very end
        append_new_cell(X[0],X[2]);
    } else if (X[1] == 'insert_cell') {
        // insert a new cell after the one with id X[3]
        do_insert_new_cell_after(X[3], X[0], X[2]);
        focus(X[0]);
    }
    start_update_check();
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

function cell_set_evaluated(id) {
    var D = get_element('cell_'+id);
    D.className = "cell_evaluated";
}

function cell_set_not_evaluated(id) {
    var D = get_element('cell_'+id);
    D.className = "cell_not_evaluated";
    cell_set_done(id);
}

function cell_set_running(id) {
    set_output_text(id, '', '', '', '', '');
    cell_output_set_type(id, 'wrap');
    var cell_div = get_element('cell_div_output_' + id);
    cell_div.className = 'cell_output_running';
    var cell_number = get_element('cell_number_' + id);
    cell_number.className = 'cell_number_running';
}

function cell_set_done(id) {
    var cell_div = get_element('cell_div_output_' + id)
    cell_div.className = 'cell_output_wrap';
    var cell_number = get_element('cell_number_' + id);
    cell_number.className = 'cell_number';
}

function check_for_cell_update() {
    if (active_cell_list.length == 0) {
        cancel_update_check();
        return;
    }
    var cell_id = active_cell_list[0];
    update_time = time_now();
    async_request('async_obj_cell_update', '/cell_update',
                    check_for_cell_update_callback,
                    'cell_id=' + cell_id + '&worksheet_id='+worksheet_id);
}

function start_update_check() {
    if(updating) return;
    updating = true;
    check_for_cell_update();
}

function cancel_update_check() {
    updating = false;
    clearTimeout(update_timeout);
}

function set_output_text(id, text, wrapped_text, output_html, status, introspect_html) {
    /* fill in output text got so far */
    var cell_output = get_element('cell_output_' + id);
    var cell_output_nowrap = get_element('cell_output_nowrap_' + id);
    var cell_output_html = get_element('cell_output_html_' + id);

    cell_output.innerHTML = wrapped_text;
    cell_output_nowrap.innerHTML = text;
    cell_output_html.innerHTML = output_html;

    if (status == 'd') {
         cell_set_done(id);
         // TODO: should make this not case sensitive!!  how to .lower() in javascript?
         if (text.indexOf('class="math"') != -1 || text.indexOf("class='math'") != -1) {
             try {
                 jsMath.Process(cell_output);
                 /* jsMath.ProcessBeforeShowing(cell_output_nowrap); */
                 /* jsMath.ProcessBeforeShowing(cell_output);
                 jsMath.ProcessBeforeShowing(cell_output_nowrap); */
             } catch(e) {
                 cell_output.innerHTML = jsmath_font_msg + cell_output.innerHTML;
                 cell_output_nowrap.innerHTML = jsmath_font_msg + cell_output_nowrap.innerHTML;
             }
         }
    } else {
    }

    if(introspect_id == id) {
        if (status == 'd') {
            introspection_loaded = true;
            introspection_text = introspect_html;
        }
        update_introspection_text();
    } else if(introspect_html != '') {
        cell_output.innerHTML = '';
        cell_output_nowrap.innerHTML = '';
        cell_output_html.innerHTML = introspect_html;
    }
}

function set_input_text(id, text) {
    /* fill in input text */
    var cell_input = get_cell(id);
    cell_input.value = text;
    focus(id)

    try {
        pos = text.length - after_cursor.length;
        cell_input.selectionStart = pos;
        cell_input.selectionEnd = pos;
    }catch(e){}
    try{
        var range = document.selection.createRange();
        range.moveStart('character', -after_cursor.length);
        range.moveEnd('character', -after_cursor.length);
    }catch(e){}

    return false;
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

function check_for_cell_update_callback(status, response_text) {
    // make sure the update happens again in a few hundred milliseconds,
    // unless a problem occurs below.
    if (status != "success") {
        if(update_error_count>update_error_threshold) {
            cancel_update_check();
            halt_active_cells();
            var elapsed_time = update_cell_error_count*update_error_delta/1000;
            var msg = "Error updating cell output after " + elapsed_time + "s";
            msg += "(canceling further update checks).";
            alert(msg);
            return;
        }
        cell_output_delta = update_error_delta;
        update_error_count++;
        continue_update_check();
        return;
    } else {
        cell_output_delta = update_normal_delta;
    }

    var i = response_text.indexOf(' ');
    var id = response_text.substring(1, i);
    var stat = response_text.substring(0,1)

    if(response_text == 'empty') {
        cancel_update_check();
        return;
    }

    if(stat == 'e') {
        cancel_update_check();
        halt_active_cells();
        return;
    }

    /* compute output for a cell */
    var D = response_text.slice(i+1).split(SEP);
    var output_text         = D[0] + ' ';
    var output_text_wrapped = D[1] + ' ';
    var output_html         = D[2];
    var new_cell_input      = D[3];
    var interrupted         = D[4];
    var variable_list       = D[5];
    var object_list         = D[6];
    var attached_files_list = D[7];
    var introspect_html     = D[8];
    var j = id_of_cell_delta(id,1);

    set_output_text(id, output_text, output_text_wrapped,
                    output_html, stat, introspect_html);
    if (stat == 'd') {
        active_cell_list = delete_from_array(active_cell_list, id);

        if (interrupted == 'false') {
            cell_set_evaluated(id);
        } else {
            halt_active_cells();
            cancel_update_check();
        }

        if(active_cell_list.length == 0)
            cancel_update_check();

        if (new_cell_input != '') {
            set_input_text(id, new_cell_input);
        }

        set_variable_list(variable_list);
        set_object_list(object_list);
        set_attached_files_list(attached_files_list);
    }

    continue_update_check();
}

function continue_update_check() {
    var time_elapsed = time_now() - update_time;
    if(time_elapsed < cell_output_delta) {
        update_timeout = setTimeout('check_for_cell_update()', cell_output_delta-time_elapsed);
    } else {
        check_for_cell_update();
    }
}

///////////////////////////////////////////////////////////////////
//  Slideshow Functions
///////////////////////////////////////////////////////////////////

function slide_mode() {
    in_slide_mode = true;
    set_class('left_pane', 'hidden');
    set_class('cell_controls', 'hidden');
    set_class('slide_controls', 'control_commands');
    set_class('worksheet', 'slideshow');

    for(i = 0; i < cell_id_list.length ; i++) {
        set_class('cell_outer_'+cell_id_list[i], 'hidden');
    }
    if(current_cell == -1 && cell_id_list.length>0)
        current_cell = cell_id_list[0];
    if(current_cell != -1)
        set_class('cell_outer_'+current_cell, 'cell_visible');
}

function cell_mode() {
    in_slide_mode = false;
    set_class('left_pane', 'pane');
    set_class('cell_controls', 'control_commands');
    set_class('slide_controls', 'hidden');
    set_class('worksheet', 'worksheet');

    for(i = 0; i < cell_id_list.length ; i++) {
        set_class('cell_outer_'+cell_id_list[i], 'cell_visible');
    }
}

function slide_hide() {
    set_class('cell_outer_' + current_cell, 'hidden');
}
function slide_show() {
    set_class('cell_outer_' + current_cell, 'cell_visible');
}

function slide_first() {
    slide_hide();
    current_cell = cell_id_list[0];
    slide_show();
}

function slide_last() {
    slide_hide();
    current_cell = cell_id_list[cell_id_list.length-1];
    slide_show();
}


function slide_next() {
    slide_hide();
    current_cell =  id_of_cell_delta(current_cell, 1);
    slide_show();
}

function slide_prev() {
    slide_hide();
    current_cell =  id_of_cell_delta(current_cell, -1);
    slide_show();
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

function delete_from_array(v, x) {
    /* Delete first occurrence of x in v.
       Returns resulting array (creates a new array!).
       No error if x is not in v.
    */
    var i;
    for (i=0; i<v.length; i++)
        if (v[i] == x) {
            return v.slice(0,i).concat(v.slice(i+1,v.length));
        }
    return v;
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
}

function make_new_cell(id, html) {
    var new_cell = document.createElement("div");
    var in_cell = document.createElement("div");
    new_cell.appendChild(in_cell);
    new_cell.id = 'cell_outer_' + id;
    in_cell.id = 'cell_' + id;
    in_cell.innerHTML = html;
    return new_cell;
}

function do_insert_new_cell_before(id, new_id, new_html) {
  /* Insert a new cell with the given new_id and new_html
     before the cell with given id. */
    var new_cell = make_new_cell(new_id, new_html);
    var cell = get_element('cell_outer_' + id);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.insertBefore(new_cell, cell);
    var i = array_indexOf(cell_id_list, eval(id));
    cell_id_list = insert_into_array(cell_id_list, i, eval(new_id));
}

function do_insert_new_cell_after(id, new_id, new_html) {
  /* Insert a new cell with the given new_id and new_html
     after the cell with given id. */
    i = id_of_cell_delta(id,1);
    if(i == id) {
        append_new_cell(new_id,new_html);
    } else {
        do_insert_new_cell_before(i, new_id, new_html);
    }
}

function insert_new_cell_before_callback(status, response_text) {
    if (status == "failure") {
        alert("Problem inserting new input cell before current input cell.");
        return ;
    }
    /* Insert a new cell _before_ a cell. */
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);
    do_insert_new_cell_before(id, new_id, new_html);
    jump_to_cell(new_id,0);
}

function insert_new_cell_before(id) {
    async_request('async_obj_add_new_cell', '/new_cell',
                   insert_new_cell_before_callback, 'id='+id);
}

function insert_new_cell_after_callback(status, response_text) {
    if (status == "failure") {
        alert("Problem inserting new cell after current cell.");
        return ;
    }
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);
    do_insert_new_cell_after(id, new_id, new_html);
    jump_to_cell(new_id,0);
}

function insert_new_cell_after(id) {
    async_request('async_obj_add_new_cell', '/new_cell_after',
                   insert_new_cell_after_callback, 'id='+id);
}

function append_new_cell(id, html) {
    var new_cell = make_new_cell(id, html);
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
    restart_sage_callback("success", "");
    if (response_text == "restart") {
        alert("The SAGE kernel had to be restarted (your variables are no longer defined).");
        restart_sage_callback('success', response_text);
    } else if(status == "success") {
        halt_active_cells()
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
    link.innerHTML = "Interrupt"
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

function hide_all_callback() {
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

function show_all_callback() {
}

function show_all() {
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        cell_output_set_type(v[i],'wrap');
    }
    async_request('async_obj_hide_all', '/show_all',
                      show_all_callback, 'worksheet_id='+worksheet_id);
}

function halt_active_cells() {
    for(i = 0; i < active_cell_list.length; i++)
        cell_set_not_evaluated(active_cell_list[i]);
    active_cell_list = []
}

function restart_sage_callback(status, response_text) {
    for(i = 0; i < cell_id_list.length; i++)
        cell_set_not_evaluated(cell_id_list[i]);
    active_cell_list = []
    var link = get_element("restart_sage");
    link.className = "restart_sage";
    link.innerHTML = "Restart";
    set_variable_list('');
}

function restart_sage() {
    var link = get_element("restart_sage");
    link.className = "restart_sage_in_progress";
    link.innerHTML = "Restart";
    async_request('async_obj_restart_sage', '/restart_sage',
                      restart_sage_callback, 'worksheet_id='+worksheet_id);
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
    "menubar=1,scrollbars=1,width=800,height=600,resizable=1");
}


/********************* js math ***************************/


function jsmath_init() {
    try {
        jsMath.Process();
      /*  jsMath.ProcessBeforeShowing();  */
    } catch(e) {
        font_warning();
    }

}

function font_warning() {
/*    alert(jsmath_font_msg); */
}


"""

    s = s.replace('SAGE_URL',SAGE_URL)

    s += r"""
///////////////////////////////////////////////////////////////////
//
// KeyCodes (auto-generated from config.py and user's sage config
//
///////////////////////////////////////////////////////////////////

//this one isn't auto-generated.  its only purpose is to make
//onKeyPress stuff simpler for text inputs and the whatlike.
function is_submit(e) {
  e = new key_event(e);
  return key_generic_submit(e);
}

%s
"""%keyhandler.all_tests()

    s += keyboards.get_keyboard('')

    return s


class JSKeyHandler:
    """This class is used to make javascript functions to check for
specific keyevents."""

    def __init__(self):
        self.key_codes = {};

    def set(self, name, key='', alt=False, ctrl=False, shift=False):
        """Add a named keycode to the handler.  When built by \code{all_tests()},
it can be called in javascript by \code{key_<key_name>(event_object)}.
The function returns true if the keycode numbered by the \code{key} parameter
was pressed with the appropriate modifier keys, false otherwise."""
        self.key_codes.setdefault(name,[])
        self.key_codes[name] = [JSKeyCode(key, alt, ctrl, shift)]

    def add(self, name, key='', alt=False, ctrl=False, shift=False):
        """Similar to \code{set_key(...)}, but this instead checks if there
is an existing keycode by the specified name, and associates the specified key
combination to that name in addition.  This way, if different browsers don't catch
one keycode, multiple keycodes can be assigned to the same test."""
        try: self.key_codes[name]
        except KeyError: self.key_codes.setdefault(name,[])
        self.key_codes[name].append(JSKeyCode(key,alt,ctrl,shift))

    def all_tests(self):
        """Builds all tests currently in the handler.  Returns a string of
javascript code which defines all functions."""
        tests = ''
        for name, keys in self.key_codes.items():
            tests += """
function key_%s(e) {
  return %s;
}"""%(name, "\n  || ".join([k.js_test() for k in keys]))

        return tests;


class JSKeyCode:
    def __init__(self, key, alt, ctrl, shift):
        global key_codes
        self.key   = key
        self.alt   = alt
        self.ctrl  = ctrl
        self.shift = shift

    def js_test(self):
        t = "(((e.k == %s) || (e.m == %s))"%(self.key, self.key)
        if self.alt:
            t += " && e.a"
        if self.ctrl:
            t += " && e.c"
        if self.shift:
            t += " && e.s"
        t+= ")"
        return t


keyhandler = JSKeyHandler()



