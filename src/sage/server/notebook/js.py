r"""nodoctest
Javascript (AJAX) Component of \sage Notebook

AUTHORS:
    -- William Stein
    -- Tom Boothby
    -- Alex Clemesha


This file is one big raw triple-quoted string that contains a bunch of
javascript.  This javascript is inserted into the head of the notebook
web page.
"""

from sage.misc.misc import SAGE_URL
from compress.JavaScriptCompressor import JavaScriptCompressor
import keyboards

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#                     2006 Tom Boothby <boothby@u.washington.edu>
#
#   Released under the *modified* BSD license.
#     Tom wrote in email to me at wstein@gmail.com on March 2, 2008: "You have my permission
#     to change the license on anything I've contributed to the notebook, to whatever suits you."
#
###########################################################################


def javascript():
    s = async_lib()
    s += notebook_lib()
    s += jmol_lib()

    return s


def jmol_lib():
    s = r"""
function jmol_applet(size, url) {
    jmolSetDocument(cell_writer);
    jmolApplet(size, "script " + url);
}

function jmol_popup(url) {
    win = window.open ("", "jmol viewer", "width=600,height=600,resizable=1,statusbar=0");
    win.document.body.innerHTML = "";
    win.document.title = "Sage 3d Viewer";
    win.document.writeln("<h1 align=center>Sage 3d Viewer</h1>");
    jmolSetDocument(win.document);
    jmolApplet("100%", "script" + url);
    win.focus();
}
    """

    return s

def async_lib():
    s = r"""
///////////////////////////////////////////////////////////////////
// An AJAX framework for connections back to the
// Sage server (written by Tom Boothby and William Stein).
///////////////////////////////////////////////////////////////////


//globals

var async_oblist = [null,null,null,null,null];
var async_idstack= [0,1,2,3,4];

function getAsyncObject(handler) {
  var asyncObj;
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

function asyncCallbackHandler(id) {
    //this was a one-liner, but Opera doesn't like to eval
    // "function() {bla}" -- it needs to be part of an assignment
    //Also, some versions of firefox don't like to see "function()"
    //you need a space between the parentheses.  WTF?
    var f;
    eval("f = function( ) { async_callback("+id+"); }");
    return f;
}

function async_callback(id) {
    var asyncObj = async_oblist[id][0];
    var callback = async_oblist[id][1];
    try {
        if( (asyncObj.readyState==4 || asyncObj.readyState=="complete")
              && asyncObj.status == 200 )
            try {
                callback('success', asyncObj.responseText);
                async_release(id);  //don't release the id until we've tried to capture output
            } catch(e) {
                async_release(id);  //release immediately in case exception was in the callback
                callback('success', "empty");
            }
    } catch(e) {
        if(async_oblist[id] != null) //release immediately
            async_release(id);
        callback("failure", e);
    }
}

function async_request(url, callback, postvars) {
  var id = async_id();
  var f = asyncCallbackHandler(id);
  var asyncObj = getAsyncObject(f);
  async_oblist[id] = [asyncObj,callback];

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

function async_id() {
  if(async_idstack.length == 0) {
    id = async_oblist.length;
    async_oblist.push(null);
  } else {
    id = async_idstack.pop();
  }
  return id
}

function async_release(id) {
  async_oblist[id] = null;
  async_idstack.push(id);
  if(async_idstack.length == async_oblist.length && async_oblist.length > 10) {
    async_oblist = [null,null,null,null,null];
    async_idstack= [0,1,2,3,4];
  }
}

"""
    return s


def notebook_lib():
    s= r"""

///////////////////////////////////////////////////////////////////
//
// GLOBAL VARIABLES
//
// PLEASE define all global variables up here, and if you want to
// set them with anything but simple assignment, 'var' them first,
// and set them later.  Your code might work in your browser, but
// it might break initial setup for other critical pieces in other
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
var update_error_delta = 1024;
//var update_normal_delta = 256;
var update_normal_delta = 512;
var cell_output_delta = update_normal_delta;

var server_ping_time = 30000;  /* Is once very 30 seconds way too fast?  Is it just right?  */


var SEP = '___S_A_G_E___';   // this had better be the same as in the server
var current_cell = -1;       // gets set on focus / blur
var no_async = false; //this isn't really used -- should we think about dealing with this?
var cell_has_changed = false;
var cell_to_focus = -1;

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

// Info about the current worksheet.  These get set in notebook.py
var worksheet_id=0;
var worksheet_filename='';
var worksheet_name='';
var user_name='';

//regular expressions used to peek into the cell input for introspection
var non_word = "[^a-zA-Z0-9_]"; //finds any character that doesn't belong in a variable name
var command_pat = "([a-zA-Z_][a-zA-Z._0-9]*)$"; //identifies the command at the end of a string
var function_pat = "([a-zA-Z_][a-zA-Z._0-9]*)\\([^()]*$";
var one_word_pat = "([a-zA-Z_][a-zA-Z._0-9]*)";
var unindent_pat = "^\\s{0,4}(.*)$";
var uncomment_pat = "^([^\\#]*)\\#{0,1}(.*)$";  //the # doesn't need a slash for now... but let's give it one anyway
var whitespace_pat = "(\\s*)";

try{
  non_word = new RegExp(non_word);
  command_pat = new RegExp(command_pat);
  function_pat = new RegExp(function_pat);
  one_word_pat = new RegExp(one_word_pat);
  whitespace_pat = new RegExp(whitespace_pat);
  unindent_pat = new RegExp(unindent_pat);
  uncomment_pat = new RegExp(uncomment_pat);
} catch(e){}

var after_cursor, before_cursor, before_replacing_word;

var update_timeout = -1;

var updating = false; var update_time = -1;

var jsmath_font_msg = '<a href="SAGE_URL/jsmath">Click to download and install tex fonts.</a><br>';

jsMath = {Font: {Message: function () {}}}

var cell_id_list; // this gets set in worksheet.py

var input_keypress; //this gets set to a function when we set up the keyboards
var input_keydown; //this gets set to a function when we set up the keyboards
var debug_keypress; //this gets set to a function when we set up the keyboards

var in_debug_input = false;
var in_slide_mode = false; //whether or not we're in slideshow mode
var slide_hidden = false; //whether the current slide has the hidden input class

var worksheet_locked;

var original_title;

var line_height=1.2;

//var title_spinner = ['    ', '.   ', '..  ', '... '];
//var title_spinner = ['[ ] ', '[.] ', '[:] ', '[.] '];
//var title_spinner = ['S ', 'SA ', 'SAG ', 'SAGE '];
var title_spinner = ['/ ', '\\ '];
//var title_spinner = ['[   ] ', '[.  ] ', '[.. ] ', '[...] '];
//var title_spinner = ['[-] ','[/] ','[|] ','[\\] '];
var title_spinner_i = 0;


///////////////////////////////////////////////////////////////////
//
// Cross-Browser Stuff
//
///////////////////////////////////////////////////////////////////


function initialize_the_notebook(){
    try{
        original_title = document.title;
    } catch(e) {}

    try{
      [].indexOf || (Array.prototype.indexOf = function(v,n){
        n = (n==null)?0:n; m = this.length;
        for(var i = n; i < m; i++)
          if(this[i] == v)
             return i;
        return -1;
      });
    } catch(e){}


    try{
      var n=navigator;
      var nav=n.appVersion;
      var nap=n.appName;
      var nua=n.userAgent;
      browser_op=(nua.indexOf('Opera')!=-1);
      browser_saf=(nua.indexOf('Safari')!=-1);
      browser_konq=(!browser_saf && (nua.indexOf('Konqueror')!=-1) ) ? true : false;
      browser_moz=( (!browser_saf && !browser_konq ) && ( nua.indexOf('Gecko')!=-1 ) ) ? true : false;
      browser_ie=((nap.indexOf('Internet Explorer') != -1)&&!browser_op);
      browser_ie5=(browser_ie&&(nua.indexOf('MSIE 5')!=-1));
      os_mac=(nav.indexOf('Mac')!=-1);
      os_win=( ( (nav.indexOf('Win')!=-1) || (nav.indexOf('NT')!=-1) ) && !os_mac)?true:false;
      os_lin=(nua.indexOf('Linux')!=-1);
//      show_exception(n);
    } catch(e){
      show_exception(e);
    }

    get_keyboard();

    jsmath_init();
}


function true_function() {return true;}
input_keypress = true_function;

function get_keyboard() {
  var b,o,warn=false;

  input_keypress = cell_input_key_event;
  input_keydown = true_function;
  debug_keypress = debug_input_key_event;

  if(browser_op) {
    b = "o";
  } else if(browser_ie) {
    b = "i";
    input_keypress = true_function;
    input_keydown = cell_input_key_event;
    debug_keypress = true_function;
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

  async_request('/javascript/keyboard/'+b+o, get_keyboard_callback, null);
}


function get_keyboard_callback(status, response_text) {
  if(status == 'success') {
    eval(response_text);
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

function get_class(id) {
  e = get_element(id);
  if(e!=null) {
      return e.className;
  }
  return null
}

function set_html(id, html) {
  e = get_element(id);
  if(e!=null) {
      e.innerHTML = html;
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

function current_selection(input) {
    if(browser_ie) {
        var range = document.selection.createRange();
        return range.text;
    } else {
        return input.value.substring(input.selectionStart,input.selectionEnd);
    }
}

function get_selection_range(input) {
    if(browser_ie) {
        var start, end;
        var range = document.selection.createRange();

        var tmprange = range.duplicate();
        tmprange.moveToElementText(input);
        tmprange.setEndPoint("endToStart", range);
        start = tmprange.text.length;

        tmprange = range.duplicate();
        tmprange.moveToElementText(input);
        tmprange.setEndPoint("endToEnd", range);
        end = tmprange.text.length;

        return Array(start, end);
    } else {
        return Array(input.selectionStart, input.selectionEnd);
    }
}

function set_selection_range(input, start, end) {
    if(browser_ie) {
        input.value = input.value.replaceAll("\r\n", "\n");
        var range = document.selection.createRange();
        range.moveToElementText(input);
        range.moveStart('character', start);
        range.setEndPoint("endToStart", range)
        range.moveEnd('character', end-start);
        range.select();
    } else {
        input.selectionStart = start;
        input.selectionEnd = end;
    }
}

///////////////////////////////////////////////////////////////////
//
// Misc page functions -- for making the page work nicely
// (this is a crappy descriptor)
///////////////////////////////////////////////////////////////////


// Replaces all instances of the given substring.
// From http://www.bennadel.com/blog/142-Ask-Ben-Javascript-String-Replace-Method.htm

String.prototype.replaceAll = function(strTarget, strSubString ) {
	var strText = this;
	var intIndexOfMatch = strText.indexOf( strTarget );
	// Keep looping while an instance of the target string
	// still exists in the string.
	while (intIndexOfMatch != -1) {
		// Replace out the current instance.
		strText = strText.replace( strTarget, strSubString )
		// Get the index of any next matching substring.
		intIndexOfMatch = strText.indexOf( strTarget );
	}
	return( strText );
}

function is_whitespace(s) {
    m = whitespace_pat.exec(s);
    return (m[1] == s);
}

function trim(s) {
    m = one_word_pat.exec(s);
    if(m == null)
        return s;
    return m[1];
}

function body_load() {
// init_menus();
}

function init_menus() {
  for( i = 1; i <= 3; i++) {
    menu = get_element("menu"+i);
    menu.style.display="none";
  }
}

function toggle_menu(name) {
  if(get_class(name) == "hidden") {
    set_class(name, name);
    set_html(name+'_hider', '[-]');
  } else {
    set_class(name, 'hidden');
    set_html(name+'_hider', '[+]');
  }
}

function toggle_top() {
  toggle('topbar')
}

function toggle(el) {
  var el = get_element(el)
  if ( el.style.display != 'none' ) {
    el.style.display = 'none';
  } else {
    el.style.display = '';
  }
}

function toggle_top() {
  toggle('topbar')
}

function toggle(el) {
  var el = get_element(el)
  if ( el.style.display != 'none' ) {
    el.style.display = 'none';
  } else {
    el.style.display = '';
  }
}

function toggle_left_pane() {
  if(get_class('left_pane') == "hidden") {
    set_class('left_pane', 'pane');
    set_class('worksheet', 'worksheet');
  } else {
    set_class('left_pane', 'hidden');
  }
}





function show_exception(e) {
/*    var mess = "";
    for(var i = 0; i < e.length; i++) {
        mess+= i+"\t"+e[i]+"\n";
    }
    for(var i in e) {
        mess+= i+"\t"+e[i]+"\n";
    }
    alert(mess);
/*/
    var h = "<table>\n";
    for(var i = 0; i < e.length; i++) {
        h+= "<tr><td>"+i+"</td><td>"+e[i]+"</td></tr>\n";
    }
    for(var i in e) {
        h+= "<tr><td>"+i+"</td><td>"+e[i]+"</td></tr>\n";
    }
    h+="</table>";
    var d = document.createElement("div");
    var n = Math.floor(Math.random()*32000);
    d.setAttribute("id", "exception_div_"+n);
    d.setAttribute("onClick", "this.parentNode.removeChild(this);");
    d.innerHTML = h;
    d.setAttribute("style", "position:absolute; display:block; z-index:1000; top:0px; left:0px;");
    document.body.insertBefore(d, document.body.childNodes[0]);
// */
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
    cell_focus(id, false);

    if(do_trim) //optimization 'cause Opera has a slow regexp engine
        word = trim(word);

    cell_input.value = before_replacing_word + word + after_cursor;

    var pos = before_replacing_word.length + word.length;

    //note for explorer:  may need to focus cell first.
    set_selection_range(cell_input,pos,pos);

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
   // o = document.open("/" + name + ".sobj");
}


///////////////////////////////////////////////////////////////////
//
// WORKSHEET functions -- for switching between and managing worksheets
//
///////////////////////////////////////////////////////////////////

function new_worksheet() {
    open("/new_worksheet")
}

function set_worksheet_list_checks() {
    /* Go through and set all check boxes the same as they are in the control box */
    var C, i, id, X;
    C = get_element("controlbox");
    for(i=0; i<worksheet_filenames.length; i++) {
        id = worksheet_filenames[i];
        X  = get_element(id);
        X.checked = C.checked;
    }
}

function worksheet_list_button(action, desc) {
    /* For each filename listed in worksheet_filenames, look up the corresponding
       input check box, see if it is checked, and if so, do the corresponding
       action.
     */
    var i, id, X, filenames;
    filenames = "";
    for(i=0; i<worksheet_filenames.length; i++) {
        id = worksheet_filenames[i];
        X  = get_element(id);
        if (X.checked) {
            filenames = filenames + worksheet_filenames[i] + SEP;
            X.checked = 0;
        }
    }
    async_request(action, worksheet_list_button_callback, 'filenames='+filenames + '&sep='+SEP);
}

function worksheet_list_button_callback(status, response_text) {
   if (status == 'success') {
      if (response_text != '') {
          alert(response_text);
      }
   } else {
      alert("Failure deleting worksheet." + response_text);
   }
  window.location.reload(true);
}

function delete_button() {
    worksheet_list_button("/send_to_trash", "--> trash");
}

function make_active_button() {
    worksheet_list_button("/send_to_active", "(--> active");
}

function archive_button() {
    worksheet_list_button("/send_to_archive", "--> archived");
}


function history_window() {
    window.open ("/history",
      "", "menubar=1,scrollbars=1,width=800,height=600, toolbar=1,resizable=1");
}

function upload_worksheet_button() {
    window.location.replace("/upload");
}

function copy_worksheet() {
    window.location.replace(worksheet_command("copy"));
}

function rate_worksheet(rating) {
    comment = get_element("rating_comment").value;
    window.location.replace(worksheet_command("rate?rating="+rating + "&comment="+escape0(comment)));
}

function download_worksheet(base_filename) {
    open(worksheet_command("download/" + base_filename + '.sws'));
}

function worksheet_settings() {
    window.location.replace(worksheet_command("settings"));
}

function share_worksheet() {
    window.location.replace(worksheet_command("share"));
}

function publish_worksheet() {
    window.open(worksheet_command("publish"), "",
      "menubar=1,location=1,scrollbars=1,width=800,height=600,toolbar=1,  resizable=1");
}

function save_as(typ) {
    open(worksheet_command('save_as') + '?typ=' +typ);
}

function edit_worksheet() {
    window.location.replace(worksheet_command(""));
}

function save_worksheet() {
    async_request(worksheet_command('save_snapshot'), save_worksheet_callback, null);
}

function save_worksheet_callback(status, response_text) {
   if (status != 'success') {
       alert("Failed to save worksheet.");
       return;
   }
}

function close_callback(status, response_text) {
   if (status != 'success') {
       alert(response_text);
       return;
   }
    window.location.replace('/');
}

function save_worksheet_and_close() {
    async_request(worksheet_command('save_snapshot'), close_callback, null);
}

function worksheet_discard() {
    async_request(worksheet_command('revert_to_last_saved_state'), close_callback, null);
}

function rename_worksheet() {
   var new_worksheet_name = prompt('Enter new worksheet name:',worksheet_name);
   if (new_worksheet_name == null) return;
   var T = get_element("worksheet_title");
   var set_name;
   if (new_worksheet_name.length >= 30) {
       set_name = new_worksheet_name.slice(0,30) + ' ...';
   } else {
       set_name = new_worksheet_name;
   }
   T.innerHTML = set_name;
   worksheet_name = new_worksheet_name;
   async_request(worksheet_command('rename'), null, 'name='+escape0(new_worksheet_name));
}

function entsub_ws(event, typ) {
  if (event && event.which == 13)
     search_worksheets(typ);
  else
     return true;
}


function search_worksheets(typ) {
    X = get_element('search_worksheets');
    url = '?typ=' + typ + '&search=' + escape0(X.value);
    window.location.replace(url);
}

function go_system_select(theform, original_system) {
   with(theform) {
      var system = options[selectedIndex].value;
      if (confirm("All cells will be evaluted using " + system + " until you change the system back.")) {
          system_select(system);
      } else {
          options[original_system].selected = 1;
      }
   }
}

function system_select(s) {
    async_request(worksheet_command('system/'+s), null, null);
}

function go_pretty_print_check(theform) {
   with(theform) {
          pretty_print_check(checked);
   }
}

function pretty_print_check(s) {
    async_request(worksheet_command('pretty_print/'+s), null, null);
}



function go_data(theform) {
   var value;
   with(theform) {
      value = options[selectedIndex].value;
      if(value == "__upload_data_file__") {
          window.location.replace(worksheet_command("upload_data"));
      } else {
          window.location.replace("/home/" + worksheet_filename + "/" + value);
      }
      options[0].selected = 1;
   }
}

function add_worksheet(name) {
    open("/home/" + user_name + "/" + name)
}

function add_worksheet_callback(status,response_text) {
    if (status == "success") {
        /* expect response_text to encode a pair consisting of
           the HTML for the updated worksheet list and the
           name of the new worksheet. */
        var X = response_text.split(SEP);
        if (X.length <= 1) {
            alert("Unable to add worksheet.");
        } else {
            set_worksheet_list(X[0]);
        }
    } else {
        alert("Possible failure adding worksheet.");
    }
}

function delete_worksheet(name) {
    async_request('/send_to_trash', delete_worksheet_callback, 'filename='+escape0(name))
}

function delete_worksheet_callback(status, response_text) {
    if (status == "success") {
        window.location.replace("/?typ=trash");

    } else {
        alert("Possible failure deleting worksheet.");
    }
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
    window.open("__upload__.html","","location=1,menubar=1,scrollbars=0,width=800,height=700,toolbar=1,resizable=1");
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
   async_request('/upload_worksheet', upload_worksheet_callback, 'filename='+filename)
}

function upload_worksheet_callback(status, response_text) {
    if (status == "success") {
        if (response_text.slice(0,5) == "Error") {
            alert("Error uploading worksheet.");
        } else {
            set_worksheet_list(response_text);
        }
    } else {
        alert("Possible problem uploading file.");
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
   /* hide_add_new_worksheet_menu(); */
    var add_worksheet_box = get_element('new_worksheet_box');
    name = add_worksheet_box.value;
    if (name == '') {
       alert("Enter a worksheet name in the box and click new to create a new worksheet.");
       return;
    }
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



function unlock_worksheet() {
    lock = get_element("worksheet_lock");
    lock.innerHTML = 'Enter Passcode: <input onKeyPress="return unlock_worksheet_submit(event,value);" id="lock_input" type="password">';
    lock.innerHTML+= '<span id="unlock_error" class="red"></span>';
    lock_input = get_element("lock_input");
    lock_input.focus();
}

function unlock_worksheet_submit(e,passcode) {
    if(is_submit(e)) {
        document.cookie = "ws_"+worksheet_filename+"_passcode="+passcode;
        async_request('/unlock_worksheet', unlock_worksheet_callback, 'worksheet_id='+worksheet_id);
        return false;
    }
    return true;
}

function unlock_worksheet_callback(status, response_text) {
    if(status == 'success' && response_text == 'ok') {
        lock = get_element("worksheet_lock");
        lock.parentNode.removeChild(lock);
        worksheet_locked = false;
    } else {
        lock_input = get_element("lock_input");
        lock_input.value = "";
        lock_input.focus();
        txt = get_element('unlock_error');
        if(txt)
            txt.innerHTML = 'incorrect';
    }
}

function sync_active_cell_list() {
    async_request('/get_queue', sync_active_cell_list_callback, 'worksheet_id='+worksheet_id);
}

function sync_active_cell_list_callback(status, response_text) {
    if(status == 'success') {
        if(response_text == "")
            return;
        active_cell_list = response_text.split(",");
        for(var i = 0; i < active_cell_list.length; i++)
            cell_set_running(active_cell_list[i]);
        start_update_check();
    }
}

///////////////////////////////////////////////////////////////////
//
// WORKSHEET list functions -- i.e., functions on a specific
// worksheet in the list of worksheets display.
//
///////////////////////////////////////////////////////////////////

function refresh() {
    window.location.replace(location.href);
}

function go_option(theform) {
   with(theform) {
      eval(options[selectedIndex].value);
      options[0].selected = 1;
   }
}

function link_datafile(target_worksheet_filename, filename) {
   open(worksheet_command("link_datafile?filename=" + escape0(filename) +
         "&target="+escape0(target_worksheet_filename)));
}


function list_rename_worksheet(filename, curname) {
   var new_name = prompt('Enter new worksheet name:', curname);
   async_request('/home/' + filename + '/' + 'rename',
            list_rename_worksheet_callback, 'name='+ escape0(new_name));
}

function list_rename_worksheet_callback(status, response_text) {
   refresh();
}


function list_edit_worksheet(filename) {
    window.location.replace('/home/' + filename);
}

function list_copy_worksheet(filename) {
    async_request('/home/' + filename + '/copy', list_copy_worksheet_callback, null);
}

function list_copy_worksheet_callback(status, response_text) {
    window.location.replace('/');
}

function list_share_worksheet(filename) {
   window.location.replace('/home/' + filename + '/share');
}

function list_publish_worksheet(filename) {
   window.open('/home/' + filename + '/publish', "",
      "menubar=1,scrollbars=1,width=800,height=600,toolbar=1,  resizable=1");
}

function list_revisions_of_worksheet(filename) {
   window.location.replace('/home/' + filename + '/revisions');
}

function list_preview_worksheet(filename) {
   window.location.replace('/home/' + filename + '/preview');
}

///////////////////////////////////////////////////////////////////
//
// Tell server I am alive
//
///////////////////////////////////////////////////////////////////

/* This pings the server every 30 seconds to announce that we are
   still viewing this page.   If it fails, it should probably perform
   some action to indicate that there is no server running.
*/

function server_ping_while_alive() {
    async_request(worksheet_command('alive'), server_ping_while_alive_callback, null);
    setTimeout("server_ping_while_alive();", server_ping_time);
}

function server_ping_while_alive_callback(status, response_text) {
    if (status == "failure") {
        server_down();
    } else {
        server_up();
    }
}

function server_down() {
   set_class("ping", "pingdown");
}

function server_up() {
   set_class("ping", "ping");
}

///////////////////////////////////////////////////////////////////
//
// CELL functions -- for the individual cells
//
///////////////////////////////////////////////////////////////////

function get_cell(id) {
    return get_element('cell_input_'+ id);
}

function cell_blur(id) {
    var cell = get_cell(id);
    if(cell == null) return;

    setTimeout("set_class('eval_button"+id+"','eval_button')", 100); //this is unclickable if we don't add a little delay.

    /* Disable coloring and change to div for now */
    cell.className="cell_input";
    cell_input_minimize_size(cell);
    return true;  /* disable for now */

    cell.className="hidden";

   /* if(!in_slide_mode)
        current_cell = -1; */

    var t = cell.value.replaceAll("<","&lt;");

    if(cell_has_changed)
        send_cell_input(id);
    return true;
}

function send_cell_input(id) {
    cell = get_cell(id)
    if(cell == null) return;

    async_request("/set_cell_input", generic_callback, "cell_id="+id+"&input="+cell.value);
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

function refocus_cell() {
    if(cell_to_focus < 0) return;
    var c = cell_to_focus;  //make a temp variable so we don't trigger another body focus event
    cell_to_focus = -1;     //and cause an infinite loop.
    cell_focus(c);
}

//set and_delay to true if you want to refocus the browser in a keyevent
//which expects a tab -- Opera apparently resists canceling the tab key
//event -- so we can subvert that by breaking out of the call stack with
//a little timeout.  Safari also has this problem.
function cell_focus(id, bottom) {
    var cell = get_cell(id);
    if (cell) {
        cell.focus();
//        cell.className="cell_input_active";
        cell_input_resize(cell);
        if (!bottom)
            move_cursor_to_top_of_cell(cell);
        current_cell = id;
        cell.focus();
    }
    current_cell = id;
    cell_has_changed = false;

    return true;
}

function move_cursor_to_top_of_cell(cell) {
    set_selection_range(cell, 0,0);
}

function focus_delay(id,bottom) {
    if(!bottom)
         setTimeout('cell_focus('+id+',false)', 10);
    else
         setTimeout('cell_focus('+id+',true)', 10);
}

function number_of_rows(txt, ncols) {
    var r;
    r = txt.split('\n');
    var e, i, k, nrows;
    nrows = r.length;
    for(i=0; i < nrows; i++) {
        try {
            nrows += Math.floor(r[i].length/ncols);
        } catch(e) {

        };
    }
    return (nrows);
}

function cell_input_resize(cell_input) {
    var rows = number_of_rows(cell_input.value, cell_input.cols);
    if (rows <= 1) {
      rows = 1;
    }
    if (browser_saf) {
       rows += 1;
    }
    try {
        cell_input.style.height = 0.5 + rows*line_height + 'em';
    } catch(e) {}
    try{
        cell_input.rows = rows;
    } catch(e) {}

    if(slide_hidden) {
        cell_input.className="cell_input_active";
        slide_hidden = false;
    }
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
        cell_input.style.height = '1em';
        return;
    }

    cell_input.className = 'cell_input';
    var rows = number_of_rows(v, cell_input.cols);
    if (rows < 1) {
      rows = 1;
    }
    if (rows >= 25) {
      rows = 25;
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
     // cell = get_element('cell_outer_' + id_to_delete);
     // var worksheet = get_element('worksheet_cell_list');
     // worksheet.removeChild(cell);
     // jump_to_cell(id_to_delete,-1);
     // cell_id_list = delete_from_array(cell_id_list, id_to_delete);
     // id_to_delete = -1;
        return;
    }
    var X = response_text.split(SEP);
    if (X[0] == 'ignore') {
        return; /* do not delete, for some reason */
    }
    var cell = get_element('cell_outer_' + X[1]);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.removeChild(cell);
    jump_to_cell(X[1],-1);
    cell_id_list = delete_from_array(cell_id_list, X[1]);
}


function cell_delete(id) {
   async_request(worksheet_command('delete_cell'), cell_delete_callback, 'id='+id)
}

function debug_input_key_event(e) {
    e = new key_event(e);
    debug_input = get_element('debug_input');

    if (key_down_arrow(e)) {
        var after = text_cursor_split(debug_input)[1];
        var i = after.indexOf('\n');
        if (i == -1 || after == '') {
            jump_to_cell(cell_id_list[0],0)
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

    if (key_delete_cell(e) && is_whitespace(cell_input.value)) {
        cell_delete(id);
        return false;
    }

    if((introspect_id == id) && introspection_loaded && replacing) {
        if(!handle_replacement_controls(cell_input, e)) {
            if(browser_op) focus_delay(id,true);
            return false; //otherwise, keep going
        }
        halt_introspection();
    }

//    cell_input_resize_(cell_input); // nix this for now, onInput *should* do the trick
//    cell_input.focus();  //why is this here?  cell shouldn't lose focus unless we've gotten a tab, and that gets checked later
/*    if (browser_saf)   {
        cell_input.scrollIntoView();
    }
    */

    var selection_range = get_selection_range(cell_input);
    var selection_is_empty = (selection_range[0] == selection_range[1]);

    // Will need IE version... if possible.
    if (!in_slide_mode && key_up_arrow(e) && selection_is_empty) {
        var before = cell_input.value.substring(0,selection_range[0])
        var i = before.indexOf('\n');
        if (i == -1 || before == '') {
            jump_to_cell(id,-1, true);
            return false;
        } else {
            return true;
        }
    } else if (!in_slide_mode && key_down_arrow(e) && selection_is_empty) {
        var after = cell_input.value.substring(selection_range[0])
        var i = after.indexOf('\n');
        if (i == -1 || after == '') {
            jump_to_cell(id,1);
            return false;
        } else {
            return true;
        }
    } else if (key_send_input(e)) {
       // User pressed shift-enter (or whatever the submit key is)
       evaluate_cell(id, 0);
       return false;
    } else if (key_send_input_newcell(e)) {
       evaluate_cell(id, 1);
       return false;
    } else if (key_comment(e) && !selection_is_empty) {
       return comment_cell(cell_input);
    } else if (key_uncomment(e) && !selection_is_empty) {
       return uncomment_cell(cell_input);
    } else if (key_unindent(e) && !selection_is_empty) { //unfortunately, shift-tab needs to get caught before not-shift tab
       unindent_cell(cell_input);
       return false;
    } else if (key_request_introspections(e) && selection_is_empty) {
       // command introspection (tab completion, ?, ??)
       evaluate_cell(id, 2);
       focus_delay(id,true);
       return false;
    } else if (key_indent(e) && !selection_is_empty) {
       indent_cell(cell_input);
       return false;
    } else if (key_interrupt(e)) {
       interrupt();
       return false;
    } else if (key_page_down(e)) {
       if(in_slide_mode) {
           slide_next();
       } else {
           jump_to_cell(id, 5);
       }
       return false;
    } else if (key_page_up(e)) {
       if(in_slide_mode) {
           slide_prev();
       } else {
           jump_to_cell(id, -5);
       }
       return false;
    } else if (key_request_history(e)) {
       history_window();
    } else if (key_request_log(e)) {
       text_log_window(worksheet_filename);
    }

    cell_has_changed = true;
    return true;
}

function id_of_cell_delta(id, delta) {
    if (cell_id_list.length == 0) {
        alert("bug -- no cells.");
        return;
    }
    var i = cell_id_list.indexOf(eval(id));
    var new_id;
    if (i == -1) {
        return(id); /* Better not to move. */
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
    output.innerHTML = "";
}

function debug_append(txt) {
    output = get_element("debug_output");
    if(output == null) return;
    output.innerHTML = txt + "\n" + output.innerHTML;
}

function jump_to_cell(id, delta, bottom) {
     if(delta != 0)
        id = id_of_cell_delta(id, delta)
    if(in_slide_mode) {
        jump_to_slide(id);
    } else {
        cell_focus(id, bottom);
    }
}

function escape0(input) {
    input = escape(input);
    input = input.replace(/\+/g,"%2B");
    return input;
}

function text_cursor_split(input) {
    var a,b;
    var R = get_selection_range(input);
    b = input.value.substr(0,R[1]);
    a = input.value.substr(b.length);
    return new Array(b,a);
}

function indent_cell(input) {
    var R = get_selection_range(input);
    var start = 1+input.value.lastIndexOf("\n", R[0]);
    var a = input.value.substring(0, start);
    var b = input.value.substring(start, R[1]);
    var c = input.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++)
        lines[i] = "    "+lines[i];
    b = lines.join("\n");
    input.value = a+b+c;
    set_selection_range(input, a.length, a.length+b.length);
}

function unindent_cell(input) {
    var R = get_selection_range(input);
    var start = 1+input.value.lastIndexOf("\n", R[0]);
    var a = input.value.substring(0, start);
    var b = input.value.substring(start, R[1]);
    var c = input.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++)
        lines[i] = unindent_pat.exec(lines[i])[1];  //square brackets pull the captured pattern
    b = lines.join("\n");
    input.value = a+b+c;
    set_selection_range(input, a.length, a.length+b.length);
}

function comment_cell(input) {
    var R = get_selection_range(input);
    if(R[0] == R[1]) return true;
    var start = 1+input.value.lastIndexOf("\n", R[0]);
    var a = input.value.substring(0, start);
    var b = input.value.substring(start, R[1]);
    var c = input.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++)
        lines[i] = "#"+lines[i];
    b = lines.join("\n");
    input.value = a+b+c;
    set_selection_range(input, a.length, a.length+b.length);
}

function uncomment_cell(input) {
    var R = get_selection_range(input);
    if(R[0] == R[1]) return true;
    var start = 1+input.value.lastIndexOf("\n", R[0]);
    var a = input.value.substring(0, start);
    var b = input.value.substring(start, R[1]);
    var c = input.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++){
        m = uncomment_pat.exec(lines[i]);
        lines[i] = m[1]+m[2];
    }
    b = lines.join("\n");
    input.value = a+b+c;
    set_selection_range(input, a.length, a.length+b.length);
}




function worksheet_command(cmd) {
   return ('/home/' + worksheet_filename + '/' + cmd);
}

function evaluate_cell(id, action) {
    if(worksheet_locked) {
        alert("This worksheet is read only.  Please make a copy or contact the owner to change it.")
        return;
    }

    active_cell_list = active_cell_list.concat([id]);

    if(action == 2) { // Introspection
       evaluate_cell_introspection(id,null,null);
       return;
    }

    cell_has_changed = false; //stop from sending the input twice.
    if(!in_slide_mode) {
       jump_to_cell(id,1);
    }
    cell_set_running(id);

    var cell_input = get_cell(id);
    var I = cell_input.value;
    var input = escape0(I);
    async_request(worksheet_command('eval'), evaluate_cell_callback,
            'newcell=' + action + '&id=' + id + '&input='+input);
}

function evaluate_cell_introspection(id, before, after) {
    var cell_input = get_cell(id);

    replacing = false;
    if(before == null) {
        var in_text = text_cursor_split(cell_input);
        before_cursor = before = in_text[0];
        after_cursor = after = in_text[1];
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
            replacing_word = m[1];
            before_replacing_word = before.substring(0, before.length-replacing_word.length);
        } else if(f != null) { //we're in an open function paren -- give info on the function
            before = f[1] + "?";
        } else { //just a tab
            cell_has_changed = true;
            do_replacement(id, '    ',false);
            return;
        }
    } else {
        sub_introspecting = true;
    }
    if(!replacing && (browser_op || browser_saf))
        focus_delay(id);

    update_introspection_text();
    var before_cursor_e = escape0(before);
    var after_cursor_e = escape0(after);
    cell_set_running(id);
    async_request(worksheet_command('introspect'), evaluate_cell_callback,
          'id=' + id + '&before_cursor='+before_cursor_e + '&after_cursor='+after_cursor_e);
}

function evaluate_cell_callback(status, response_text) {
    /* update focus and possibly add a new cell to the end */
    if (status == "failure") {
       /* alert("Failure evaluating a cell."); */
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
        jump_to_cell(X[0],0);
    }
    start_update_check();
}

function cell_output_set_type(id, typ, do_async) {
    /* We do this specifically because manipulate cells do not work at all when
       displayed in nowrap mode, which is VERY BAD.  So instead for manipulates
       one gets a toggle to and from hidden.
    */
    if (typ=="nowrap" && get_element("cell-manipulate-" + id)) {
        /* if the type is nowrap and the cell-manipulate-[id] div exists (i.e., we are manipulating)
           then just make the thing hidden. */
        typ = "hidden";
    }

    /* OK, now set the sell output type.  */

    set_class('cell_div_output_' + id,    'cell_div_output_' + typ)
    set_class('cell_output_' + id,        'cell_output_' + typ)
    set_class('cell_output_nowrap_' + id, 'cell_output_nowrap_' + typ)
    set_class('cell_output_html_' + id,   'cell_output_html_' + typ)

    /* Do async request back to the server */
    if(do_async != false)
        async_request(worksheet_command('set_cell_output_type'), generic_callback, 'id='+id+'&type=' + typ)
}

function cycle_cell_output_type(id) {
    var cell_div = get_element('cell_div_output_' + id);

    if (cell_div.className == 'cell_div_output_hidden' || cell_div.className=='cell_div_output_running') {
        cell_output_set_type(id, 'wrap');
        return;
    }

    if (cell_div.className == 'cell_div_output_wrap') {
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
    set_output_text(id, '', '', '', '', '', 1);   // the 1 means no manipulation dynamics
    cell_output_set_type(id, 'wrap');
    var cell_div = get_element('cell_div_output_' + id);
    cell_div.className = 'cell_output_running';
    var cell_number = get_element('cell_number_' + id);
    cell_number.className = 'cell_number_running';
}

function cell_set_done(id) {
    var cell_div = get_element('cell_div_output_' + id)
    cell_div.className = 'cell_div_output_wrap';
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
    async_request(worksheet_command('cell_update'),
                    check_for_cell_update_callback,
                    'id=' + cell_id);
    try{
        title_spinner_i = (title_spinner_i+1)%title_spinner.length;
        document.title = title_spinner[title_spinner_i] + original_title;
    } catch(e){}
}

function start_update_check() {
    if(updating) return;
    updating = true;
    check_for_cell_update();
    set_class('interrupt', 'interrupt')
}

function cancel_update_check() {
    updating = false;
    clearTimeout(update_timeout);
    set_class('interrupt', 'interrupt_grey')
    document.title = original_title;
}

function contains_jsmath(text) {
    // TODO: should make this not case sensitive!!  how to .lower() in javascript?
    return (text.indexOf('class="math"') != -1 || text.indexOf("class='math'") != -1);
}

function set_output_text(id, text, wrapped_text, output_html, status, introspect_html, no_manip) {
    var cell_manip = get_element("cell-manipulate-" + id);
    if (!no_manip && cell_manip) {
        if (status  != 'd') return;
        var i = wrapped_text.indexOf('<?START>');
        var j = wrapped_text.indexOf('<?END>');
        if (i == -1 || j == -1) {
            alert("bug -- manipulate wrapped text is invalid" + wrapped_text);
            return;
        }
        var new_manip_output = wrapped_text.slice(i+8,j);

        /* An error occured accessing the data for this cell.  Just force reload
           of the cell, which will certainly define that data. */
        if (new_manip_output.indexOf('__SAGE_MANIPULATE_RESTART__') != -1) {
            evaluate_cell(id, 0);
        } else {
            cell_manip.innerHTML = new_manip_output;
            if (contains_jsmath(new_manip_output)) {
               jsMath.ProcessBeforeShowing(cell_manip);
            }
        }
    } else {
        /* fill in output text got so far */
        var cell_output = get_element('cell_output_' + id);
        var cell_output_nowrap = get_element('cell_output_nowrap_' + id);
        var cell_output_html = get_element('cell_output_html_' + id);

        cell_output.innerHTML = wrapped_text;
        cell_output_nowrap.innerHTML = text;
        cell_output_html.innerHTML = output_html;

        /* Did we just create or evaluate a new manipulate cell? */
        var cell_manip = get_element("cell-manipulate-" + id);
        /* If so, trigger it so that we see the evaluated version
           of the manipulate cell. */
        if (cell_manip) {
            /*****************************************************************
              This is the first time that the underlying Python manipulate function is
               actually called!
             *****************************************************************/
            manipulate(id, 'sage.server.notebook.manipulate.state[' + id + ']["function"]()');
        }
    }

    if (status == 'd') {
         cell_set_done(id);
         if (contains_jsmath(text)) {
             try {
                 /* jsMath.Process(cell_output); */
                 /* jsMath.ProcessBeforeShowing(cell_output_nowrap); */
                 jsMath.ProcessBeforeShowing(cell_output);
                 /* jsMath.ProcessBeforeShowing(cell_output_nowrap); */
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
    jump_to_cell(id,0)

    pos = text.length - after_cursor.length;
    set_selection_range(cell_input, pos, pos);

    return false;
}

function set_object_list(objects) {
    var objlist = get_element('object_list');
    objlist.innerHTML = objects;
}

function set_attached_files_list(objects) {
    var objlist = get_element('attached_list');
    objlist.innerHTML = objects;
}

/* When the page is loaded, let javascript write
 * directly to the document. After that, make sure
 * javascript writes to a CellWriter object. */

function CellWriter() {
    function write(s) {
        this.buffer += s;
    }
    this.write = write;
    this.buffer = "";
}

cell_writer = document;

function eval_script_tags(text) {
   var s = text; //text.replaceAll('\n','');
   var i = s.indexOf('<'+'script>');
   while (i != -1) {
       var j = s.indexOf('<'+'/script>');
       var code = s.slice(8+i,j);
       try {
           cell_writer = new CellWriter();
           window.eval(code);
       } catch(e) {
           alert(e);
       }
       s = s.slice(0,i) + cell_writer.buffer + s.slice(j+9);
       i = s.indexOf('<'+'script>');
   }
   return s;
}

function check_for_cell_update_callback(status, response_text) {
    // make sure the update happens again in a few hundred milliseconds,
    // unless a problem occurs below.
    if (status != "success") {
        if(update_error_count>update_error_threshold) {
            cancel_update_check();
            halt_active_cells();
            var elapsed_time = update_error_count*update_error_delta/1000;
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
        update_error_count = 0;
        cell_output_delta = update_normal_delta;
    }

    var i = response_text.indexOf(' ');
    var id = response_text.substring(1, i);
    var stat = response_text.substring(0,1)

    if(response_text == 'empty') {
    /* TODO  -- hack -- we are sometimes getting something nothing back from
       twisted for some reason.  Ignoring it seems to work....
       */
       continue_update_check();
       return;
  /*        cancel_update_check();
        return; */
    }

    if(stat == 'e') {
        cancel_update_check();
        halt_active_cells();
        return;
    }

    /* compute output for a cell */
    var D = response_text.slice(i+1).split(SEP);
    var output_text = D[0] + ' ';
    var output_text_wrapped = D[1] + ' ';
    var output_html = D[2];
    var new_cell_input = D[3];
    var interrupted = D[4];
    var introspect_html = D[5];
    var j = id_of_cell_delta(id,1);

    /* Evaluate javascript */
    output_text = output_text.replace(/<script.*?>(.|\n|\r)*?<\/script>/gim, '&lt;script&gt;');
    output_text_wrapped = eval_script_tags(output_text_wrapped);
    output_html = eval_script_tags(output_html);

    set_output_text(id, output_text, output_text_wrapped,
                    output_html, stat, introspect_html);

    if (stat == 'd') {
        active_cell_list = delete_from_array(active_cell_list, id);

        if (interrupted == 'restart') {
            restart_sage();
        } else if (interrupted == 'false') {
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
// Slideshow Functions
///////////////////////////////////////////////////////////////////

/* Switch into slide mode. */
function slide_mode() {
    in_slide_mode = true;
    set_class('left_pane', 'hidden');
    set_class('cell_controls', 'hidden');
    set_class('slide_controls', 'slide_control_commands');
    set_class('left_pane_bar', 'hidden');

    for(i = 0; i < cell_id_list.length ; i++) {
        set_class('cell_outer_'+cell_id_list[i], 'hidden');
    }
    slide_show();
}


function cell_mode() {
    in_slide_mode = false;
    set_class('left_pane', 'pane');
    set_class('cell_controls', 'control_commands');
    set_class('slide_controls', 'hidden');
    set_class('worksheet', 'worksheet');
    set_class('left_pane_bar', 'left_pane_bar');

    for(i = 0; i < cell_id_list.length ; i++) {
        set_class('cell_outer_'+cell_id_list[i], 'cell_visible');
    }
}

function slide_hide() {
    set_class('cell_outer_' + current_cell, 'hidden');
}

function slide_show() {
    if(current_cell != -1) {
        set_class('cell_outer_' + current_cell, 'cell_visible');
    } else {
        if(cell_id_list.length>0)
            current_cell = cell_id_list[0];
        set_class('cell_outer_' + current_cell, 'cell_visible');
    }
    if(current_cell != -1) {
        input = get_cell(current_cell);
        if(input != null) {
            s = lstrip(input.value).slice(0,5)
            cell_focus(current_cell, false);
            if (s == '%hide') {
                slide_hidden = true;
                input.className = 'cell_input_hide';
                input.style.height = '1.5em';
            }
        }
    }
    update_slideshow_progress();
}

function slide_first() {
    jump_to_slide(cell_id_list[0]);
}

function slide_last() {
    jump_to_slide(cell_id_list[cell_id_list.length-1]);
}


function slide_next() {
    jump_to_slide(id_of_cell_delta(current_cell, 1));
}

function slide_prev() {
    jump_to_slide(id_of_cell_delta(current_cell, -1));
}

function jump_to_slide(id) {
    slide_hide();
    current_cell = id;
    slide_show();
}


function update_slideshow_progress() {
    var i = cell_id_list.indexOf(current_cell) + 1;
    var n = cell_id_list.length;
    var bar = get_element("slideshow_progress_bar")
    if(bar != null)
        bar.style.width = "" + 100*i/n + "%";
    text = get_element("slideshow_progress_text")
    if(text != null)
        text.innerHTML = i + " / " + n;
}


///////////////////////////////////////////////////////////////////
// Insert and move cells
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
    var i = cell_id_list.indexOf(eval(id));
    cell_id_list = insert_into_array(cell_id_list, i, eval(new_id));
}

function insert_new_cell_after(id) {
    async_request(worksheet_command('new_cell_after'), insert_new_cell_after_callback, 'id='+id);
}

function insert_new_cell_after_callback(status, response_text) {
    if (status == "failure") {
        alert("Problem inserting new input cell after current input cell.\n" + response_text);
        return ;
    }
    if (response_text == "locked") {
        alert("Worksheet is locked.  Cannot insert cells.");
        return;
    }
    /* Insert a new cell _after_ a cell. */
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);
    do_insert_new_cell_after(id, new_id, new_html);
    jump_to_cell(new_id,0);
}


function do_insert_new_cell_after(id, new_id, new_html) {
  /* Insert a new cell with the given new_id and new_html
     after the cell with given id.      */

    i = id_of_cell_delta(id,1);
    if(i == id) {
        append_new_cell(new_id,new_html);
    } else {
        do_insert_new_cell_before(i, new_id, new_html);
    }
}




function insert_new_cell_before(id) {
    async_request(worksheet_command('new_cell_before'), insert_new_cell_before_callback, 'id='+id);
}

function insert_new_cell_before_callback(status, response_text) {
    if (status == "failure") {
        alert("Problem inserting new input cell before current input cell.");
        return ;
    }
    if (response_text == "locked") {
        alert("Worksheet is locked.  Cannot insert cells.");
        return;
    }
    /* Insert a new cell _before_ a cell. */
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);
    do_insert_new_cell_before(id, new_id, new_html);
    jump_to_cell(new_id,0);
}


function append_new_cell(id, html) {
    var new_cell = make_new_cell(id, html);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.appendChild(new_cell);
    cell_id_list = cell_id_list.concat([eval(id)]);
    if(in_slide_mode) {
        set_class('cell_outer_'+id, 'hidden');
        update_slideshow_progress();
    }else {
        jump_to_cell(id, 0);
    }
}


///////////////////////////////////////////////////////////////////
//
// CONTROL functions
//
///////////////////////////////////////////////////////////////////

function interrupt_callback(status, response_text) {
    if (response_text == "failed") {
       alert('Unable to immediately interrupt calculation.');
       return;
    } else if(status == "success") {
        halt_active_cells()
    }
    set_class("interrupt", "interrupt");
}

function interrupt() {
/*    var link = get_element("interrupt");
    if (link.className == "interrupt_grey") {
        return;
    }
    set_class("interrupt", "interrupt_in_progress");
*/
    async_request(worksheet_command('interrupt'), interrupt_callback);
}


function evaluate_all() {
    /* Iterate through every input cell in the document, in order,
       and evaluate.
    */
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        var cell_input = get_cell(v[i]);
        var I = cell_input.value;
        if (trim(I).length > 0) {
            evaluate_cell(v[i],0);
        }
    }
}

function hide_all_callback() {
}

function hide_all() {
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        cell_output_set_type(v[i],'hidden', false);
    }
    async_request(worksheet_command('hide_all'), hide_all_callback);
}

function show_all_callback() {
}

function show_all() {
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        cell_output_set_type(v[i],'wrap', false);
    }
    async_request(worksheet_command('show_all'), show_all_callback);
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
/*    var link = get_element("restart_sage");
    link.className = "restart_sage";
    link.innerHTML = "Restart";
*/
    sync_active_cell_list();
}

function restart_sage() {
    async_request(worksheet_command('restart_sage'), restart_sage_callback);
}

function quit_sage() {
    async_request(worksheet_command('quit_sage'), restart_sage_callback);
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

////////////////////////////////////////////
//
// doc-browser related
//
////////////////////////////////////////////

function show_doc_browser() {
    window.open("/doc_browser?/?index.html","","location=1,menubar=1,scrollbars=0,width=700,height=600,toolbar=1,resizable=1");
}

////////////////////////////////////////////


////////////////////////////////////////////
//
// wiki-window related stuff
//
///////////////////////////////////////////

function insert_new_doc_html_after(id,doc_htmlin) {
    async_request('/new_doc_html_after', insert_doc_html_callback, 'id=' + id + '&doc_htmlin=' + doc_htmlin);
}

//-------------------------------------------------
/* This is a early hack attempt at putting html into the worksheet cell area from the wiki window.*/
//
function insert_doc_html_callback(status, response_text) {
    if (status == "failure") {
        alert("Problem inserting new doc html after current cell.");
        return ;
    }
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);
    do_insert_doc_html(id, new_id, new_html);
    jump_to_cell(new_id,0);
}

function do_insert_doc_html(id,new_id,html) {
    var new_doc_html = make_new_doc_html(id, html);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.appendChild(new_doc_html);
    cell_id_list = cell_id_list.concat([new_id]);
}

function make_new_doc_html(id, html) {
    var new_doc_html = document.createElement("div");
    new_doc_html.id = 'doc_html_'+id;
    new_doc_html.innerHTML = html;
    return new_doc_html;
}

function upload_doc_html(lastid,doc_html) {
    insert_new_doc_html_after(lastid,doc_html);
}
//
//------------------------------------------------------

function get_cell_list() {
    return cell_id_list;
}

function hide_wiki_window() {
    var wiki = get_element("wiki_window");
    wiki.style.display = "none";
}

function show_wiki_window(worksheet) {
    window.open (worksheet+"__wiki__.html","", "location=1,menubar=1,scrollbars=1,width=750,height=700,toolbar=1,resizable=1");
}

function insert_cells_from_wiki(text,do_eval) {
    var eval_param = "&eval=0";
    if(do_eval)
        eval_param = "&eval=1";
    async_request("/insert_wiki_cells", insert_cells_from_wiki_callback,
                  "worksheet_id="+worksheet_id+"&text="+escape0(text)+eval_param);
}

function get_worksheet_id(){
    return worksheet_id;
}

function insert_cells_from_wiki_callback(status, response_text) {
    if(status == "success") {
        var X = response_text.split(SEP);
        var do_eval = eval(X[0]);
        var new_cell_id_list = eval(X[1]);
        var old_first_cell = cell_id_list[0];
        for(var i = 2; i < X.length; i++)
            do_insert_new_cell_before(old_first_cell, new_cell_id_list[i-2], X[i]);
        if(do_eval)
            evaluate_all();
    }
}


///////////////////////////////////////////////////////////////////
//
// LOG windows
//
///////////////////////////////////////////////////////////////////

function history_window() {
    history = window.open ("/history",
      "", "menubar=1,scrollbars=1,width=800,height=600, toolbar=1,resizable=1");
}


function doctest_window(worksheet) {
    log = window.open ("/home/" + worksheet+"/plain","",
    "menubar=1,scrollbars=1,width=800,height=600,toolbar=1, resizable=1");
}


function print_worksheet(worksheet) {
    log = window.open (worksheet_command("print"),"",
      "menubar=1,scrollbars=1,width=800,height=600,toolbar=1,  resizable=1");
}

function help(worksheet) {
    log = window.open ("/help","",
    "menubar=1,location=1,scrollbars=1,width=800,height=650,toolbar=1,  resizable=1");
}


//////////////////////////////////
// HELP
/////////////////////////////////

function show_help_window(worksheet) {
    help = window.open ("/help","",
    "menubar=1,scrollbars=1,width=800,height=600,resizable=1, toolbar=1");
}


///////////////////////////////////////////////////////////////////
// Manipulate
///////////////////////////////////////////////////////////////////

function manipulate(id, input) {
    active_cell_list = active_cell_list.concat([id]);
    async_request(worksheet_command('eval'), evaluate_cell_callback,
            'newcell=0' + '&id=' + id + '&input='+escape0('%manipulate\n' + input));
}

/********************* js math ***************************/


function jsmath_init() {
    try {
    jsMath.Process();
    /*   jsMath.ProcessBeforeShowing(); */
    } catch(e) {
/*        font_warning(); */
    }

}

function font_warning() { /* alert(jsmath_font_msg); */
}

"""

    s = s.replace('SAGE_URL',SAGE_URL)

    s += r"""
///////////////////////////////////////////////////////////////////
//
// KeyCodes (auto-generated from config.py and user's sage config
//
///////////////////////////////////////////////////////////////////

//this one isn't auto-generated.  its only purpose is to make //onKeyPress stuff simpler for text inputs and the whatlike.
function is_submit(e) {
  e = new key_event(e);
  return key_generic_submit(e);
}

%s
"""%keyhandler.all_tests()
    s += keyboards.get_keyboard('')
    return s






class JSKeyHandler:
    """This class is used to make javascript functions to check for specific keyevents."""

    def __init__(self):
        self.key_codes = {};

    def set(self, name, key='', alt=False, ctrl=False, shift=False):
        """Add a named keycode to the handler.  When built by \code{all_tests()}, it can be called in javascript by
\code{key_<key_name>(event_object)}.  The function returns true if the keycode numbered by the \code{key} parameter
was pressed with the appropriate modifier keys, false otherwise."""
        self.key_codes.setdefault(name,[])
        self.key_codes[name] = [JSKeyCode(key, alt, ctrl, shift)]

    def add(self, name, key='', alt=False, ctrl=False, shift=False):
        """Similar to \code{set_key(...)}, but this instead checks if there is an existing keycode by the specified
name, and associates the specified key combination to that name in addition.  This way, if different browsers don't
catch one keycode, multiple keycodes can be assigned to the same test."""
        try: self.key_codes[name]
        except KeyError: self.key_codes.setdefault(name,[])
        self.key_codes[name].append(JSKeyCode(key,alt,ctrl,shift))

    def all_tests(self):
        """Builds all tests currently in the handler.  Returns a string of javascript code which defines all
functions."""
        tests = ''
        for name, keys in self.key_codes.items():
            tests += """ function key_%s(e) {
  return %s;
}"""%(name, "\n || ".join([k.js_test() for k in keys]))

        return tests;


class JSKeyCode:
    def __init__(self, key, alt, ctrl, shift):
        global key_codes
        self.key = key
        self.alt = alt
        self.ctrl = ctrl
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




