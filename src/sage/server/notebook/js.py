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
        try {
            callback("failure", e);
        } catch(e) {
            /* In some cases the failure report can't be done as above because
               callback itself is not a function. */
        }
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

/* DOCSTRINGS:

All functions below must be documented using the following format,
with the docstring in C-style comments.

Short description.

Longer description.

INPUT:
    each input variable name -- description of it.
GLOBAL INPUT:
    each global variable that significantly impacts
    the behavior of this function and how
OUTPUT:
    decription of output or side effects.
*/

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

var update_count = 0;
var update_falloff_threshold = 20;
var update_falloff_level = 0;
var update_falloff_deltas = [250, 500, 1000, 5000];

var update_error_count = 0;
var update_error_threshold = 30;

// in milliseconds
var update_error_delta = 1024;
var update_normal_delta = update_falloff_deltas[0];
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


/* If this flag is set, then the next call to jump_to_cell is ignored and
   this flag is cleared.  We use this in some places to avoid changing the
   focus. */
var ignore_next_jump = false;

var control_key_pressed = 0;

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

original_title = document.title;

function initialize_the_notebook() {
    /*
    Do the following:
        1. Make sure that arrays have an indexOf method, since we use that a lot.
           We do this since not all browsers support this method, so we insert
           it in.
        2. Determine the browser OS, type e.g., opera, safari, etc.; we set global
           variables for each type.
        3. Figure out which keyboard the user has.
        4. Initialize jsmath.
    */
    try{
      [].indexOf || (Array.prototype.indexOf = function(v,n){
        n = (n==null)?0:n; m = this.length;
        for(var i = n; i < m; i++)
          if(this[i] == v)
             return i;
        return -1;
      });
    } catch(e){}


    // Determine the browser, OS and set global variables.
    try {
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
    } catch(e){
      alert(e);
    }

    // Get the keyboard codes for our browser/os combination
    get_keyboard();

    // Attempt to render any jsmath in this page.
    jsmath_init();
}


// The function that always returns true.
function true_function() { return true; }
input_keypress = true_function;

function get_keyboard() {
  /*
  Determine which keycodes we want, then make a request back to the
  server for those keycodes.  When the server returns the javascript
  with exactly those keycodes, we eval that javascript.

  OUTPUT:
      set some global variables that record platform specific key codes
  */
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
    input_keypress = true_function;
    input_keydown = cell_input_key_event;
    debug_keypress = true_function;
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
    /*
    This is called right after we get the list of keycodes for our
    browser back from the server. We eval them hence setting a bunch of
    (global) variables that define platform-specific keycodes.
    */
    if(status == 'success') {
        eval(response_text);
    }
}


function get_element(id) {
    /*
    Return the DOM element with the given id.

    INPUT:
        id -- a string
    OUTPUT:
        a DOM element.
    */
    if(document.getElementById)
        return document.getElementById(id);
    if(document.all)
        return document.all[id];
    if(document.layers)
        return document.layers[id];
}

function set_class(id, cname) {
    /*
    Set the class of the DOM element with given id to cname.

    INPUT:
        id -- a string
        cname -- a string
    OUTPUT:
        Sets the class of the DOM element with the
        given id to be class.
    */
    e = get_element(id);

    if(e!=null) {
        e.className = cname;
    }
}


function get_class(id) {
    /*
    Get the clas of the DOM element with the given id.

    INPUT:
        id -- a string
    OUTPUT:
        string, or null if there is no such DOM element.
    */
    e = get_element(id);
    if(e!=null) {
        return e.className;
    }
    return null
}

function set_html(id, html) {
    /*
    Set the inner HTML of the DOM element with given id, if there is
    such an element.

    INPUT:
        id -- an integer
        html -- string
    OUTPUT:
        changes the DOM
    */
    e = get_element(id);
    if(e!=null) {
        e.innerHTML = html;
    }
}

function get_event(e) {
    /*
    Just returns e unless the browser is IE (or maybe some old
    Netscape browser), in which case we have get what should be e but
    from window.event.

    INPUT:
        e -- a javascript event.
    OUTPUT:
        either e or window.event
    */
    return (e==null)?window.event:e;
}

function key_event(e) {
    /*
    Normalizes the different possible keyboard even structures for different browsers.

    INPUT:
        e -- a javascript event
    OUTPUT:
        Sets properties of the DOM object in a uniform way.
        The properties set are a, c, s, k, m.

    NOTE: We use key_event as an object.  Also, we use only 1-letter variables
    names here specifically to decrease file size.
   */

    // This is exactly as in the get_event function; see the docs there.
    if(e==null) { e = window.event; }

    // Here we set a, c, s, which tell whether the alt, control, or shift
    // keys have been pressed.
    if(e.modifiers) {
        this.a = e.modifiers | 1;
        this.c = e.modifiers | 2;
        this.s = e.modifiers | 4;
    } else {
        this.a = e.altKey;
        this.c = e.ctrlKey;
        this.s = e.shiftKey;
    }

    // we set the specific key that was pressed (no modifier), which is
    // string as a string pair n,m
    this.k = e.keyCode + "," + e.which;

    // Finally we set m, which the key but with '!' at the end if shift is pressed.
    // We do this because that's how we differentiate certain keys for browsers.
    // Look in keycodes.py for more.
    this.m = this.k + (this.s?'!':'');
    return this;
}

function time_now() {
    /*
    Return the time right now as an integer since Unix epoch in
    milliseconds.

    OUTPUT:
        an integer
    */
    return (new Date()).getTime();
}

function current_selection(input) {
    /*
    Return the text that is currently selected in a given text area.

    INPUT:
        input -- a DOM object (a textarea)
    OUTPUT:
        a string
    */
    if(browser_ie) {
        var range = document.selection.createRange();
        return range.text;
    } else {
        return input.value.substring(input.selectionStart,input.selectionEnd);
    }
}

function get_selection_range(input) {
    /*
    Return the start and end positions of the currently selected text
    in the input text area (a DOM object).

    INPUT:
        input -- a DOM object (a textarea)

    OUTPUT:
        an array of two nonnegative integers
    */
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
    /*
    Select a range of text in a given textarea.
    INPUT:
        input -- a DOM input text area
        start -- an integer
        end -- an integer
    OUTPUT:
        changes the state of the input textarea.
    */
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

function get_cursor_position(cell) {
    /* Return an integer that gives the position of the text cursor
       in the cells input field.
    INPUT:
        cell -- an input cell (not the id but the actual DOM element)
    OUTPUT:
        a single integer
    */
    return get_selection_range(cell)[1];
}

function set_cursor_position(cell, n) {
    /*
    Move the cursor position in the cell to position n.

    WARNING: Does nothing when n is 0 on Opera at present.

    INPUT:
        cell -- an actual cell in the DOM, returned by get_cell
        n -- a non-negative integer
    OUTPUT:
        changes the position of the cursor.
    */
    if (browser_op && !n) {
        // program around a "bug" in opera where using this
        // hack to position the cursor selects the entire
        // text area (which is very painful, since then the
        // user accidently deletes all their code).
        return;
    }
    // TODO: note for explorer:  may need to focus cell first.
    set_selection_range(cell, n, n);
}


///////////////////////////////////////////////////////////////////
//
// Misc page functions -- for making the page work nicely
//
///////////////////////////////////////////////////////////////////

String.prototype.replaceAll = function(strTarget, strSubString ) {
    /*
    Replace all instances of the given substring by another string>

    From http://www.bennadel.com/blog/142-Ask-Ben-Javascript-String-Replace-Method.htm

    INPUT:
        this -- the string having part of itself replaced
        strTarget -- the string that will be replaced
        strSubString -- the string that replaces strTarget
    */
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
    return(strText);
}

function is_whitespace(s) {
    /*
    Return true precisely if the input string s consists only of whitespace,
    e.g., spaces, tabs, etc.

    INPUT:
        s -- a string
    OUTPUT:
        true or false
    */

    // We check using the whitespace_pat regular expression defined at the top of
    // this file.
    m = whitespace_pat.exec(s);
    return (m[1] == s);
}

function first_variable_name_in_string(s) {
    /*
    This function returns the first valid variable name in a string.

    INPUT:
        s -- a string
    OUTPUT:
        a string
    */
    m = one_word_pat.exec(s);
    if(m == null)
        return s;
    return m[1];
}

function toggle_top() {
    /*
    Called when one clicks to toggle the top control bar in the worksheet view.
    */
    toggle_displayed('topbar');
}

function toggle_displayed(id) {
    /*
    Toggle whether or not the DOM element with the given
    id is displayed.

    INPUT:
        id -- a string DOM identifier
    OUTPUT:
        changes the given element's display style to/from none.
    */
    var el = get_element(id)
    if ( el.style.display != 'none' ) {
        el.style.display = 'none';
    } else {
        el.style.display = '';
    }
}


///////////////////////////////////////////////////////////////////
//
// Completions interface stuff
//
///////////////////////////////////////////////////////////////////

function handle_replacement_controls(cell_input, event) {
    /*
    This function handles the keyboard controls for the tab-completion
    pop-up menu.

    It's really just a not-so-good attempt to modularize the keyboard
    handling code somewhat.

    INPUT:
        cell_input -- the input textarea where the completion is happening.
        event -- the keypress event
    */

    // First change the currently selected element so it isn't highlighted.
    deselect_replacement_element();



    if(key_menu_up(event)) {            // Press the up arrow
        replacement_row--;
        // Wrap around vertically.
        // Ugly code since we don't know the size.
        if(!replacement_element_exists()) {
            replacement_row = 1;
            while(replacement_element_exists())
                replacement_row++;
            replacement_row--;
        }
    } else if(key_menu_down(event)) {  // press the down key
        replacement_row++;
        if(!replacement_element_exists())   // going down past the
            replacement_row = 0;
    } else if(key_menu_right(event)) { // press the right arrow key
        replacement_col++;
        if(!replacement_element_exists())
            replacement_col = 0;
    } else if(key_menu_left(event)) {   // left arrow key
        replacement_col--;
        // Check if we have to wrap around horizontally.
        // Ugly code since we don't know the size.
        if(!replacement_element_exists()) {
            replacement_col = 1;
            while(replacement_element_exists())
                replacement_col++;
            replacement_col--;
        }
    } else if(key_menu_pick(event)) {  // press the enter key, so we do the replacement.
        do_replacement(introspect_id, replacement_word, true);
        return false;
    } else if(key_request_introspections(event)) {
        // instead of browsing through a list of options, here we are viewing
        // the docstring on a funtion.
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

    // highlight the correct word.
    select_replacement_element();

    if(sub_introspecting) { // do the actual replacement.
        active_cell_list = active_cell_list.concat([introspect_id]);
        evaluate_cell_introspection(introspect_id, before_replacing_word+replacement_word+'?', after_cursor);
    }

    return false;
}

function do_replacement(id, word, do_trim) {
    /*

    INPUT:
        id -- an integer, the id of an input cell
        word -- a string
        do_trim -- true or false
    */

    // Get the input cell and focus on it.
    var cell_input = get_cell(id);
    cell_focus(id, false);

    // If necessary get only the first word out of the input word string.
    if(do_trim) {
         word = first_variable_name_in_string(word);
    }

    // Do the actual completion
    cell_input.value = before_replacing_word + word + after_cursor;

    // Put the cursor right after the word we just put in.
    var pos = before_replacing_word.length + word.length;
    set_cursor_position(cell_input,pos);

    // Done completing, so get rid of the completions menu.
    halt_introspection();
}

function get_replacement_element() {
    /*
    Return the DOM element that is currently highlighted in the tab
    completion popup window.

    GLOBAL INPUT:
        introspect_id -- integer; id of the input cell in which we're doing the introspection
        replacement_row -- integer; the row position of the cell we are currently selected
        replacement_col -- integer; the column position of the cell we are currently selected
    OUTPUT:
        DOM element -- that is currently highlighted
    */
    return get_element("completion"+introspect_id + "_" + replacement_row + "_" + replacement_col);
}

function replacement_element_exists() {
    /*
    Return whether or not the global variables that define the current
    row/column of the tab completion menu actually define an entry
    in the menu.  This is used to implement, e.g., wrapping around
    the sides.
    */
    return get_replacement_element() != null;
}

function select_replacement(row, col) {
    /*
    INPUT:
        row, col -- integers
    OUTPUT:
        -- move the popup highlighted menu item in the completions menu to position (row,col)
        -- set the global variables replacement_row and replacement_col.
    */
    deselect_replacement_element();
    replacement_row = row;
    replacement_col = col;
    select_replacement_element();
}

function deselect_replacement_element() {
    /*
    Change the currently selected highlighted word in the
    tab-completion popup menu so that it is no longer selected.

    This is done simply by setting the CSS className of the currently
    selected element.
    */
    e = get_replacement_element();
    if(e==null) return;
    e.className = 'completion_menu_two';
}

function select_replacement_element() {
    /*
    Highlight the currently selected completions item, and set the global
    variable replacement_word equal to this item (so it can be used if
    a selection is actually mode).

    OUTPUT:
        modifies the DOM and global variable replacement_word
    */
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

function update_introspection_text() {
    /*
    Set the contexts of the introspection (tab completion or help)
    window, or display "loading..." if we are waiting for data from
    the server.

    GLOBAL INPUTS:
        introspection_text -- string; this is what gets put in the
                              introspection window
        introspect_id -- integer; id of the input cell where we
                         are doing introspection.
    */

    // Delete the current introspection text window contexts.
    close_introspection_text();

    // Get the DOM object corresponding to this introspection window.
    d = get_element("introspect_div_"+introspect_id);
    if(!d) return;

    // Set the new introspection text.
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
    /*
    Delete the text in the introspect window.

    GLOBAL INPUT:
        introspect_id -- integer; id of the input cell where we
                         are doing introspection.
    */
    d = get_element("introspect_div_"+introspect_id);
    if(d!=null) {
        d.innerHTML = "";
    }
}

function halt_introspection() {
    /*
    We are done doing the introspection, so we close the completions
    or documentation popup window, and set several global variables to
    indicate the state of this window is closed.

    OUTPUT:
    Each of these global variables is changed:
        introspect_id, replacing, sub_introspecting, introspection_loaded,
        replacement_row, replacement_col
    */
    close_introspection_text();
    introspect_id = null;
    replacing = false;
    sub_introspecting = false;
    introspection_loaded = false;
    replacement_row = replacement_col = 0;
}

///////////////////////////////////////////////////////////////////
//
// WORKSHEET functions -- for switching between and managing worksheets
//
///////////////////////////////////////////////////////////////////

function new_worksheet() {
    /*
    Ask the server to create a new worksheet, which is then opened
    replacing the current worksheet.
    */
    open("/new_worksheet")
}

function set_worksheet_list_checks() {
    /*
    Go through and set all check boxes the same as they are in the
    control box.

    This is called when the user clicks the checkbox in the top left
    of the list of worksheets, which toggles all the checkboxes below
    it to be either on or off (select all or none).

    GLOBAL INPUT:
        worksheet_filenames -- list of strings
    */
    var C, i, id, X;
    C = get_element("controlbox");
    for(i=0; i<worksheet_filenames.length; i++) {
        id = worksheet_filenames[i];
        X  = get_element(id);
        X.checked = C.checked;
    }
}

function worksheet_list_button(action) {
    /*
    For each filename listed in worksheet_filenames, look up the
    corresponding input check box, see if it is checked, and if so, do
    the corresponding action.

    INPUT:
        action -- url that defines a message to send to the server
    GLOBAL INPUT:
        worksheet_filenames -- list of strings
        SEP -- separator string used when encoding tuples of data to send
               back to the server.
    OUTPUT:
        calls the server and requests an action be performened on all the
        listed worksheets
    */
    var i, id, X, filenames;
    filenames = "";

    // Concatenate the list of all worksheet filenames that are checked
    // togethers separated by the separator string.
    for(i=0; i<worksheet_filenames.length; i++) {
        id = worksheet_filenames[i];
        X  = get_element(id);
        if (X.checked) {
            filenames = filenames + worksheet_filenames[i] + SEP;
            X.checked = 0;
        }
    }
    // Send the list of worksheet names and requested action back to
    // the server.
    async_request(action, worksheet_list_button_callback,
                  'filenames='+filenames + '&sep='+SEP);
}

function worksheet_list_button_callback(status, response_text) {
    /*
    Handle result of performing some action on a list of worksheets.

    INPUT:
        status, response_text -- standard AJAX return values
    OUTPUT:
        display an alert if something goes wrong; refresh this
        browser window no matter what.
    */
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
    /*
    This javascript function is called when the worksheet list delete
    button is pressed.  Each worksheet whose box is checked gets sent
    to the trash.
    */
    worksheet_list_button("/send_to_trash");
}

function make_active_button() {
    /*
    Sends each checked worksheet to the active worksheets folder.
    */
    worksheet_list_button("/send_to_active");
}

function archive_button() {
    /*
    Sends each checked worksheet to the archived worksheets folder.
    */
    worksheet_list_button("/send_to_archive");
}


function history_window() {
    /*
    Display the history popup window, which displays the last few hundred
    commands typed into any worksheet.
    */
    window.open ("/history",
      "", "menubar=1,scrollbars=1,width=800,height=600, toolbar=1,resizable=1");
}

function upload_worksheet_button() {
    /*
    Replace the current display window with the upload entry box.
    */
    window.location.replace("/upload");
}

function copy_worksheet() {
    /*
    Make a copy of the current worksheet then load the copy into the
    current frame.
    */
    window.location.replace(worksheet_command("copy"));
}

function rate_worksheet(rating) {
    /*
    Save the comment and rating that the uses chooses for a public worksheet.

    INPUT:
        rating -- integer
    */
    comment = get_element("rating_comment").value;
    window.location.replace(worksheet_command("rate?rating="+rating + "&comment="+escape0(comment)));
}

function download_worksheet(base_filename) {
    /*
    Download the current worksheet to the file with given name.

    INPUT:
        base_filename
    */
    open(worksheet_command("download/" + base_filename + '.sws'));
}

function worksheet_settings() {
    /*
    Bring up the worksheet settings menu.
    */
    window.location.replace(worksheet_command("settings"));
}

function share_worksheet() {
    /*
    Display the worksheet sharing window.
    */
    window.location.replace(worksheet_command("share"));
}

function publish_worksheet() {
    /*
    Public the current worksheet.
    */
    window.open(worksheet_command("publish"), "",
      "menubar=1,location=1,scrollbars=1,width=800,height=600,toolbar=1,  resizable=1");
}

function save_as(typ) {
    /*
    Save the current worksheet to a file.
    */
    open(worksheet_command('save_as') + '?typ=' +typ);
}

function edit_worksheet() {
    /*
    Edit the current worksheet as a plain text file.
    */
    window.location.replace(worksheet_command(""));
}

function save_worksheet() {
    /*
    Save a snapshot of the current worksheet.
    */
    async_request(worksheet_command('save_snapshot'), save_worksheet_callback, null);
}

function save_worksheet_callback(status, response_text) {
    /*
    Verify that saving the current worksheet worked.
    */
    if (status != 'success') {
        alert("Failed to save worksheet.");
        return;
    }
}

function close_callback(status, response_text) {
    /*
    Called when we successfully close the current worksheet and
    want to display the user home screen (i.e., worksheet list).
    */
    if (status != 'success') {
        alert(response_text);
        return;
    }
    window.location.replace('/');
}

function save_worksheet_and_close() {
    /*
    Send message back to the server saving the current
    worksheet and quitting the Sage process, then
    close the current window returning to the home screen.
    */
    async_request(worksheet_command('save_and_quit'), close_callback, null);
}

function worksheet_discard() {
    /*
    Discard the current worksheet and quit the currently
    running Sage process, then close the current window and
    replace it by the home screen .
    */
    async_request(worksheet_command('discard_and_quit'), close_callback, null);
}

function rename_worksheet() {
    /*
    Rename the current worksheet.  This pops up a dialog that asks for
    the new worksheet name, then sets it in the browser, and finally
    sends a message back to the server stating that the worksheet has
    been renamed.
    */
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
    original_title = worksheet_name + ' (Sage)';
    document.title = original_title;
    async_request(worksheet_command('rename'), null, 'name='+escape0(new_worksheet_name));
}

function entsub_ws(event, typ) {
    /*
    Full text search through the worksheets after pressing return in
    the search input box.

    INPUT:
        event -- keyboard event; checked if it is return
        typ -- the type of worksheets to search through (active, archived, etc.)
    */
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
    /* Return the input cell as a DOM element with the given integer id.
    INPUT:
        id -- integer
    OUTPUT:
        a DOM element
    */
    return get_element('cell_input_'+ id);
}

function cell_blur(id) {
    /* This function is called when the cell with the given
       id is blurred.  It removes whitespace around the input,
       and if the cell has changed sends the changed input
       back to the server.
    INPUT:
        id -- integer
    OUTPUT:
       true -- to avoid infinite recursion.
    */
    var cell = get_cell(id);
    if(cell == null) { return true; }
    cell_input_minimize_size(cell);
    if(cell_has_changed)
        send_cell_input(id);

    /* It is very important to return true here, or one gets an
       infinite javascript recursion. */
    return true;
}

function send_cell_input(id) {
    /* Send the input text of the current cell back to the server.  This is
       typically called when the cursor leaves the current cell.
    INPUT:
       id -- an integer, id of the current cell
    OUTPUT:
       makes an async call back to the server sending the input text.
    */
    cell = get_cell(id)
    if(cell == null) return;

    async_request(worksheet_command('eval'), generic_callback, "save_only=1&id="+id+"&input="+escape0(cell.value));
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

// Set and_delay to true if you want to refocus the browser in a keyeven
// which expects a tab -- Opera apparently resists canceling the tab key
// event -- so we can subvert that by breaking out of the call stack with
// a little timeout.  Safari also has this problem.
function cell_focus(id, bottom) {
    var cell = get_cell(id);
    if (cell) {
        cell.focus();
//        cell.className="cell_input_active";
        cell_input_resize(cell);
        if (!bottom)
            move_cursor_to_top_of_cell(cell);
        cell.focus();
    }
    cell_has_changed = false;

    return true;
}

//this gets called when the cell object has gotten focus
function cell_focused(cell, id) {
    cell.className = "cell_input_active";
    if(current_cell == id) return;
    if (current_cell != -1) {
        set_class("eval_button"+current_cell,"eval_button");
    }
    current_cell = id;
    set_class("eval_button"+id,"eval_button_active");
}

function move_cursor_to_top_of_cell(cell) {
    /* Move the cursor to the first position in the given input cell.
    INPUT:
        cell -- an input cell as a DOM element
    */
    set_cursor_position(cell, 0);
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
        return;
    }
    var X = response_text.split(SEP);
    if (X[0] == 'ignore') {
        return; /* do not delete, for some reason */
    }
    var cell = get_element('cell_outer_' + X[1]);
    var worksheet = get_element('worksheet_cell_list');
    worksheet.removeChild(cell);
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
    /*
    This function is called each time a key is pressed when the cursor is inside
    an input cell.

    INPUT:
        id -- the id of the input cell
        e -- determines the keyboard event, i.e., the key press
    GLOBAL_INPUT:
        control_key_pressed -- used to detect if the control key was pressed; this is
                  really only relevant to handling Opera's quirky even model.
    OUTPUT:
        All kinds of interesting things can happen:
            - introspection
            - cell join
            - cell split
            - cell delete
            - a cell may be evaluated
    */
    cell_input = get_cell(id);
    e = new key_event(e);
    if (e==null) return;

    /*********** SPLIT AND JOIN HANDLING ********/

    /* Record that just the control key was pressed.  We do this since on Opera
       it is the only way to catch control + key. */
    if (key_control(e)) {
        control_key_pressed = 1;
        return;
    }
    /* Check for the split and join keystrokes. */
    if (key_split_cell(e) || (key_enter(e) && control_key_pressed)) {
        control_key_pressed = 0;
        split_cell(id);
        return false;
    } else if (key_join_cell(e) || (key_delete_cell(e) && control_key_pressed)) {
        control_key_pressed = 0;
        join_cell(id);
        return false;
    }

// IT IS SAFE TO DELETE THIS...
// Deleting an empty cell (a special case of join).
//    if ((key_join_cell(e)&& is_whitespace(cell_input.value) ) {
        // If we press backspace at the beginning of a cell,
        // join with the previous cell.
//        if (get_cursor_position(cell_input) == 0) {
//            join_cell(id);
//            return false;
//        }
//    }

    /* Turn off recording that the control key may have pressed last, since
       we *only* would use that in the above if statement.  NOTE: This
       is only needed on Opera.  */
    control_key_pressed = 0;

    /*********** END of SPLIT AND JOIN HANDLING ********/

    if((introspect_id == id) && introspection_loaded && replacing) {
        if(!handle_replacement_controls(cell_input, e)) {
            if(browser_op) { focus_delay(id,true); }
            return false; //otherwise, keep going
        }
        halt_introspection();
    }

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
       evaluate_cell(id, false);
       return false;
    } else if (key_send_input_newcell(e)) {
       evaluate_cell(id, true);
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
       evaluate_cell_introspection(id,null,null);
       if (browser_op) { focus_delay(id,true); }
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
     /* Put the focus and cursor in the cell that is positioned delta
     spots above or below the cell with given id.  If bottom is true
     the cursor is positioned at the bottom of the cell that is put
     in focus.

     INPUT:
         id -- an integer
         delta -- an integer (default or 0: just focus on the cell
                  with the given id).
         bottom -- if true, puts the cursor at the end of the cell
                   rather than the beginnning
     GLOBAL INPUT:
         ignore_next_jump -- if this variable
            is set globally to true, then this function immediately
            returns after setting it to false.  This is used because
            several functions create new cells with unknown id's then
            jump to them (example, when inserting a new cell after the
            current one).  In some cases, e.g., when splitting or
            joining cells, it is necessary to temporarily disable this
            behavior, even though we do call those functions.   This
            is done by simply setting ignore_next_jump to true.
     OUTPUT:
         Changes the focused cell.  Does not send any information
         back to the server.
     */
     if (ignore_next_jump) {
          ignore_next_jump = false;
          return;
     }
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

function join_cell(id) {
    /* Join the cell with given id to the cell before it.

       The output of the resulting joined cells is the output of the
       second cell, *unless* the input of the second cell is only
       whitespace, in which case the output is the output of the first
       cell.  We do this since a common way to delete a cell is to
       empty its input, then hit backspace.   It would be very confusing
       if the output of the second cell were retained.

       INPUT:
          id -- integer cell id.
       OUTPUT:
          change the state of the worksheet in the DOM, global variables,
          etc., and updates the server on this change.
    */


    var id_prev = id_of_cell_delta(id, -1);
    if(id_prev == id) return;

    var cell = get_cell(id);
    var cell_prev = get_cell(id_prev);

    // We delete the cell above the cell with given id except in the
    // one case when the cell with id has empty input, in which case
    // we just delete that cell.
    if (is_whitespace(cell.value)) {
        cell_delete(id);
        cell_prev.focus();
        return;
    }


    // The lower cell in the join is now not empty.  So we
    // delete the previous cell and put its contents in the
    // bottom cell.
    var val_prev = cell_prev.value;

    cell_delete(id_prev);
    cell.focus();

    // The following is so that joining two cells keeps a newline
    // between the input contents.
    var n = val_prev.length;
    if(val_prev[n-1] != '\n') {
        val_prev += '\n';
        n += 1;
    }
    cell.value = val_prev + cell.value;

    // Send a message back to the server reporting that the cell
    // has changed (as a result of joining).
    send_cell_input(id);

    // Set the cursor position in the joined cell to about
    // where it was before the join.
    set_cursor_position(cell, n);

    // Finally resize the joined cell to account for its new text.
    cell_input_resize(cell);
}

function split_cell(id) {
    /*
    Split the cell with the given id into two cells, inserting a new
    cell after the current one and placing the cursor at the beginning
    of the new cell.

    INPUT:
        id -- an integer

    OUTPUT:
        changes the state of the worksheet, DOM, and sends a message
        back to the server.
    */
    var cell = get_cell(id);
    var txt = text_cursor_split(cell)
    if (txt[1].length > 0 && txt[1][0] == '\n') {
        txt[1] = txt[1].slice(1);
    }

    cell.value = txt[1];
    cell_input_resize(cell);
    send_cell_input(id);  /* tell the server about how the input just got split in half. */

    set_cursor_position(cell,0);

    /* Make sure that the cursor doesn't move to the new cell. */
    ignore_next_jump = true;
    insert_new_cell_before(id,txt[0]);
}


function worksheet_command(cmd) {
    /* Create a string formatted as a URL to send back to the server
       and execute the given cmd on the current worksheet.

    INPUT:
        cmd -- string
    OUTPUT:
        a string
    */
    return ('/home/' + worksheet_filename + '/' + cmd);
}

function evaluate_cell(id, newcell) {
    /*
    Evaluate the given cell, and if newcell is true (the default),
    insert a new cell after the current one.

    INPUT:
        id -- an integer that identifies a cell
        newcell -- false -- do not insert a new cell after the current one
                     1 -- do insert a new cell
    GLOBAL INPUT:
        worksheet_locked -- if true, pop up an alert and return immediately
    OUTPUT:
        a message is sent to the server and the "check for updates"
        loop is started if it isn't already going; typically this
        will result in output being generated that we get later
    */
    if(worksheet_locked) {
        alert("This worksheet is read only.  Please make a copy or contact the owner to change it.")
        return;
    }

    // append that cell id is currently having some sort of computation
    // possibly occuring.  Note that active_cell_list is a global variable.
    active_cell_list = active_cell_list.concat([id]);

    // Stop from sending the input again to the server when we leave focus and the
    // send_cell_input function is called.
    cell_has_changed = false;

    // We only advance to the next cell if the worksheet is not in
    // single cell mode.
    if(!in_slide_mode) {
       jump_to_cell(id,1);
    }

    // Clear the output text and set the css to indicate that this
    // is a running cell.
    cell_set_running(id);

    // Finally make the request back to the server to do the actual calculation.
    var cell_input = get_cell(id);
    var I = cell_input.value;
    var input = escape0(I);
    if (newcell) { newcell = 1; } else { newcell = 0; }
    async_request(worksheet_command('eval'), evaluate_cell_callback,
            'newcell=' + newcell + '&id=' + id + '&input='+input);
}

function evaluate_cell_introspection(id, before, after) {
    /*
    Do an introspection in the cell with given id.

    INPUT:
        id -- integer; the id a cell
        before -- null, or if given all the text before the cursor
        after -- null, or if given, all the text after the cursor
    OUTPUT:
        sends a message back to the server to do an introspection on
        this cell; also set the cell running.
    */
    var cell_input = get_cell(id);

    active_cell_list = active_cell_list.concat([id]);
    cell_set_running(id);

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
    if(!replacing && browser_op)
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
    /* We do this specifically because interact cells do not work at all when
       displayed in nowrap mode, which is VERY BAD.  So instead for interacts
       one gets a toggle to and from hidden.
    */
    if (typ=="nowrap" && get_element("cell-interact-" + id)) {
        /* if the type is nowrap and the cell-interact-[id] div exists (i.e., we are interacting)
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
    set_output_text(id, '', '', '', '', '', 1);   // the 1 means no interact dynamics
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
    update_count = 0;
    update_falloff_level = 0;
    cell_output_delta = update_falloff_deltas[0];
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

function set_output_text(id, text, wrapped_text, output_html, status, introspect_html, no_interact) {
    if (id < 0) {
        /* negative id's come up for special internal usage. */
        return;
    }
    var cell_interact = get_element("cell-interact-" + id);
    if (!no_interact && cell_interact) {
        if (status  != 'd') return;
        var i = wrapped_text.indexOf('<?__SAGE__START>');
        var j = wrapped_text.indexOf('<?__SAGE__END>');
        if (i == -1 || j == -1) {
            alert("Bug in notebook -- interact wrapped text is invalid" + wrapped_text);
            return;
        }
        var new_interact_output = wrapped_text.slice(i+16,j);

        /* An error occured accessing the data for this cell.  Just force reload
           of the cell, which will certainly define that data. */
        if (new_interact_output.indexOf('__SAGE_INTERACT_RESTART__') != -1) {
            evaluate_cell(id, 0);
        } else {
            cell_interact.innerHTML = new_interact_output;
            if (contains_jsmath(new_interact_output)) {
               jsMath.ProcessBeforeShowing(cell_interact);
            }
        }
    } else {
        /* fill in output text got so far */
        var cell_output = get_element('cell_output_' + id);
        if (!cell_output) {
            alert("Bug in notebook -- missing output for cell with id "+id);
            return;
        }
        var cell_output_nowrap = get_element('cell_output_nowrap_' + id);
        var cell_output_html = get_element('cell_output_html_' + id);

        cell_output.innerHTML = wrapped_text;
        cell_output_nowrap.innerHTML = text;
        cell_output_html.innerHTML = output_html;

        if (status == 'd' && introspect_html=="") {
            /* Did we just create or evaluate a new interact cell? */
            var cell_interact = get_element("cell-interact-" + id);
            /* If so, trigger it so that we see the evaluated version
               of the interact cell. */
            if (cell_interact) {
                /*****************************************************************
                  This is the first time that the underlying Python interact function is
                   actually called!
                 *****************************************************************/
                interact(id, 'sage.server.notebook.interact.state[' + id + ']["function"]()');
            }
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
    set_cursor_position(cell_input, pos);

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
        if(update_error_count > 0) {
            update_error_count = 0;
            update_count = 0;
            update_falloff_level = 1;
            cell_output_delta = update_falloff_deltas[1];
        }
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

        update_count = 0;
        update_falloff_level = 0;
        cell_output_delta = update_falloff_deltas[0];
    } else {
        if(  update_count > update_falloff_threshold &&
             update_falloff_level+1 < update_falloff_deltas.length) {
            update_falloff_level+= 1;
            update_count = 0;
            cell_output_delta = update_falloff_deltas[update_falloff_level];
        } else {
            update_count += 1;
        }
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

function insert_new_cell_after(id, input) {
    if(input == null) input = "";
    async_request(worksheet_command('new_cell_after'), insert_new_cell_after_callback, 'id='+id+'&input='+escape0(input));
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




function insert_new_cell_before(id, input) {
    if(input == null) input = "";
    async_request(worksheet_command('new_cell_before'), insert_new_cell_before_callback, 'id='+id+'&input='+escape0(input));
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
        if (first_variable_name_in_string(I).length > 0) {
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
// Interact
///////////////////////////////////////////////////////////////////

function interact(id, input) {
    active_cell_list = active_cell_list.concat([id]);
    /* the __sage_interact__ string appears also in cell.py */
    async_request(worksheet_command('eval'), evaluate_cell_callback,
            'newcell=0' + '&id=' + id + '&input='+escape0('%__sage_interact__\n' + input));
}

///////////////////////////////////////////////////////////////////
// Base 64 encoding and decoding (mainly used for manipulate).
///////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
// The following header applies to the encode64 and decode64 functions
// This code was written by Tyler Akins and has been placed in the
// public domain.  It would be nice if you left this header intact.
// Base64 code from Tyler Akins -- http://rumkin.com

var keyStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=";

function encode64(input) {
   /* I had to add this, since otherwise when input is numeric there are
      errors below. */
   try {
       input = input.toString();
   } catch(e) {
       return input;
   }
   var output = "";
   var chr1, chr2, chr3;
   var enc1, enc2, enc3, enc4;
   var i = 0;

   do {
      chr1 = input.charCodeAt(i++);
      chr2 = input.charCodeAt(i++);
      chr3 = input.charCodeAt(i++);

      enc1 = chr1 >> 2;
      enc2 = ((chr1 & 3) << 4) | (chr2 >> 4);
      enc3 = ((chr2 & 15) << 2) | (chr3 >> 6);
      enc4 = chr3 & 63;

      if (isNaN(chr2)) {
         enc3 = enc4 = 64;
      } else if (isNaN(chr3)) {
         enc4 = 64;
      }

      output = output + keyStr.charAt(enc1) + keyStr.charAt(enc2) +
         keyStr.charAt(enc3) + keyStr.charAt(enc4);
   } while (i < input.length);

   return output;
}

function decode64(input) {
   var output = "";
   var chr1, chr2, chr3;
   var enc1, enc2, enc3, enc4;
   var i = 0;

   // remove all characters that are not A-Z, a-z, 0-9, +, /, or =
   input = input.replace(/[^A-Za-z0-9\+\/\=]/g, "");

   do {
      enc1 = keyStr.indexOf(input.charAt(i++));
      enc2 = keyStr.indexOf(input.charAt(i++));
      enc3 = keyStr.indexOf(input.charAt(i++));
      enc4 = keyStr.indexOf(input.charAt(i++));

      chr1 = (enc1 << 2) | (enc2 >> 4);
      chr2 = ((enc2 & 15) << 4) | (enc3 >> 2);
      chr3 = ((enc3 & 3) << 6) | enc4;

      output = output + String.fromCharCode(chr1);

      if (enc3 != 64) {
         output = output + String.fromCharCode(chr2);
      }
      if (enc4 != 64) {
         output = output + String.fromCharCode(chr3);
      }
   } while (i < input.length);

   return output;
}

///////////////////////////////////////////////////////////////////
// Trash
///////////////////////////////////////////////////////////////////

function empty_trash() {
    /* This asks for confirmation from the user then sends a request back to the
       server asking that the trash be emptied for this user. The request to the
       server goes by accessing the url /emptytrash.  After that finishes, the
       empty trash folder is displayed.  */
    if(confirm('Emptying the trash will permanently delete all items in the trash. Continue?')) {
        window.location.replace("/emptytrash");
        window.location.replace("/?typ=trash");
    }
}


/********************* js math ***************************/

function jsmath_init() {
    try {
         jsMath.Process();
    } catch(e) {
    }

}

"""

    s = s.replace('SAGE_URL',SAGE_URL)

    s += r"""

///////////////////////////////////////////////////////////////////
//
// KeyCodes (auto-generated from config.py and user's sage config
//
///////////////////////////////////////////////////////////////////

%s

"""%keyhandler.all_tests()

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




