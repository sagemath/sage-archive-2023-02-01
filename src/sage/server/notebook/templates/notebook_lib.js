/*
INDENTATION:
All code should have 4-space indentation, exactly like in our Python code.

DOCSTRINGS:
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
    description of output or side effects.
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
var browser_op, browser_saf, browser_konq, browser_moz, browser_ie, browser_ie5, browser_iphone;
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

var keypress_resize_delay = 250; // don't resize the cell too often! without this, typing in Safari is worse than typing on a 286!
var last_keypress_resize = 0;
var will_resize_soon = false;

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

var updating = false
var update_time = -1;

var jsmath_font_msg = '<a href="{{ SAGE_URL }}/jsmath">Click to download and install tex fonts.</a><br>';

var cell_id_list; // this gets set in worksheet.py

var input_keypress; //this gets set to a function when we set up the keyboards
var input_keydown; //this gets set to a function when we set up the keyboards
var debug_keypress; //this gets set to a function when we set up the keyboards
var skip_keyup = false; //this gets set to work around a bug

var in_debug_input = false;
var in_slide_mode = false; //whether or not we're in slideshow mode
var slide_hidden = false; //whether the current slide has the hidden input class

var doing_split_eval = false; // whether or not we're splitting a cell and evaluating it.

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

var evaluating_all = false;
var evaluating_all_cursor = 0;

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
    } catch(e) {}

    //TODO: Maybe use jquery browser plugin (http://plugins.jquery.com/project/Browser)?
    // Determine the browser, OS and set global variables.
    try {
        var n=navigator;
        var nav=n.appVersion;
        var nap=n.appName;
        var nua=n.userAgent;
        browser_op=(nua.indexOf('Opera')!=-1);
        browser_saf=(nua.indexOf('Safari')!=-1);
        browser_iphone=(nua.indexOf('iPhone')!=-1);
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

    // Trigger cell resize when the notebook resizes
    window.onresize = resize_all_cells;
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
        alert("Your browser / OS combination is not supported.  \nPlease use Firefox or Opera under Linux, Windows, or Mac OS X, or Safari.")
    }

    async_request('/javascript/keyboard/'+b+o, get_keyboard_callback);
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
    } else { /* alert("Error in set_class: no element " + id); */ }
}


function get_class(id) {
    /*
    Get the class of the DOM element with the given id.

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
    Normalizes the different possible keyboard event structures for different browsers.

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

    // Here we set this.v which tell whether the alt, control, or shift
    // keys have been pressed.
    if(e.modifiers) {
        this.v = e.modifiers;
    } else {
        this.v = 0;
        if(e.altKey) this.v+=1
        if(e.ctrlKey) this.v+=2
        if(e.shiftKey) this.v+=4
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
        // user accidentally deletes all their code).
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


function refresh() {
    /*
    This function simply refreshes the current HTML page, thus completely
    reloading the DOM from the server.
    */
    window.location.replace(location.href);
}

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


function lstrip(s) {
    /*
    Given a string s, strip leading whitespace from s and return the resulting string.

    INPUT:
        s -- a string
    OUTPUT:
        a string
    */
    var n = s.length;
    var i = 0;
    while (i < n && (s[i] == ' ' || s[i] == '\n' || s[i] == '\t')) {
        i = i + 1;
    }
    return s.slice(i);
}


function resize_all_cells() {
    /*
    Resizes all cells that do not start with %hide;
    called whenever the window gets resized.

    GLOBAL INPUT:
        cell_id_list -- a list of integers
    */
    var i,id;
    for(i=0;i<cell_id_list.length;i++) {
        // Get the id of the cell to resize
        id = cell_id_list[i];
        // Make sure it is not hidden, and if not resize it.
        if (get_cell(id).className != "cell_input_hide") {
           cell_input_resize(id);
        }
    }
}

function input_keyup(id, event) {
    /*
    Resize the cell once and a while.  Not too often.
    INPUT:
        id -- an integer
        event -- a keyup event

    GLOBAL INPUT:
        keypress_resize_delay -- amount of time to wait between resizes
        last_keypress_resize -- last time we resized
        will_resize_soon -- if a keyup event is ignored for the purpose of resizing,
                            then we queue one up.  Avoid a timeout-flood with this lock.
    */
    if(skip_keyup) { skip_keyup = false; return false; }
    var t = time_now()

    if((t - last_keypress_resize) > keypress_resize_delay) {
        last_keypress_resize = t;
        cell_input_resize(id);
    } else if(!will_resize_soon) {
        will_resize_soon = true;
        setTimeout("cell_input_resize("+id+"); will_resize_soon=false;", keypress_resize_delay)
    }

    // automatic indentation
    if (browser_iphone) return;
    var e = new key_event(event);
    if (e==null) return;
    if (key_enter(e)) {
        var cell = get_cell(id)
        //warning!  very subtle regular expression (for nonjaphs)
        // (?:\n|^)        -- starting from the last line ending (or beginning of the cell), (don't capture contents)
        // ( *)            -- as many spaces as we can find (capture this, we'll find it in RegExp.$1)
        // (?:.*?)         -- everything else in the string, but save room for the following terms (don't capture contents)
        // (:?)            -- optionally, capture a colon before the following term (find it in RegExp.$2)
        // [ \t\r\v\f]*\n$ -- ignore whitespace at the end of the line
        re = /(?:\n|^)( *)(?:.*?)(:?)[ \t\r\v\f]*\n$/;

        var position = get_cursor_position(cell);
        var text = text_cursor_split(cell);
        re.test(text[0])
        var indent = RegExp.$1
        var colon = RegExp.$2
        if (colon == ':') { indent = indent + "    " }
        get_cell(id).value = text[0] + indent + text[1];
        set_cursor_position(cell, position + indent.length)
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
        skip_keyup = true;
        return false;
    } else if(key_request_introspections(event)) {
        // instead of browsing through a list of options, here we are viewing
        // the docstring on a function.
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

        if (contains_jsmath(introspection_text)) {
            try {
                jsMath.ProcessBeforeShowing(d);
            } catch(e) {
                text_cell.innerHTML = jsmath_font_msg + d.innerHTML;
            }
        }

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

function checked_worksheet_filenames() {
    /*
    For each filename listed in worksheet_filenames, look up the
    corresponding input check box, see if it is checked, and if so,
    add it to the list.

    GLOBAL INPUT:
        worksheet_filenames -- list of strings
        SEP -- separator string used when encoding tuples of data to send
               back to the server.
    OUTPUT:
        string of worksheet filenames that are checked, separated by SEP
    */
    var i, id, X, filenames;
    filenames = "";

    // Concatenate the list of all worksheet filenames that are checked
    // together separated by the separator string.
    for(i=0; i<worksheet_filenames.length; i++) {
        id = worksheet_filenames[i];
        X  = get_element(id);
        if (X.checked) {
            filenames = filenames + worksheet_filenames[i] + SEP;
            X.checked = 0;
        }
    }
    return filenames;
}

function worksheet_list_button(action) {
    /*
    For each filename listed in worksheet_filenames, look up the
    corresponding input check box, see if it is checked, and if so, do
    the corresponding action.

    INPUT:
        action -- URL that defines a message to send to the server
    GLOBAL INPUT:
        worksheet_filenames -- list of strings
        SEP -- separator string used when encoding tuples of data to send
               back to the server.
    OUTPUT:
        calls the server and requests an action be performed on all the
        listed worksheets
    */
    // Send the list of worksheet names and requested action back to
    // the server.
    async_request(action, worksheet_list_button_callback,
                  {filenames: checked_worksheet_filenames(), sep: SEP});
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
        alert("Error applying function to worksheet(s)." + response_text);
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

function stop_worksheets_button() {
    /*
    Saves and then quits sage process for each checked worksheet.
    */
    worksheet_list_button("/send_to_stop");
}

function download_worksheets_button() {
    /*
    Downloads the set of checked worksheets as a zip file.
    */
    window.location.replace("/download_worksheets?filenames=" + checked_worksheet_filenames() + "&sep=" + SEP);
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
    async_request(worksheet_command('save_snapshot'), save_worksheet_callback);
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
    async_request(worksheet_command('save_and_quit'), close_callback);
}

function worksheet_discard() {
    /*
    Discard the current worksheet and quit the currently
    running Sage process, then close the current window and
    replace it by the home screen .
    */
    async_request(worksheet_command('discard_and_quit'), close_callback);
}

function rename_worksheet() {
    /*
    Rename the current worksheet.  This pops up a dialog that asks for
    the new worksheet name, then sets it in the browser, and finally
    sends a message back to the server stating that the worksheet has
    been renamed.
    */
    var new_worksheet_name = prompt('Enter new worksheet name:',worksheet_name);
    if (new_worksheet_name == null || new_worksheet_name == "") return;
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
    async_request(worksheet_command('rename'), null,
                  {name: new_worksheet_name});
}

function search_worksheets_enter_pressed(event, typ) {
    /*
    Full text search through the worksheets after pressing enter in
    the search input box.

    INPUT:
        event -- keyboard event; checked if it is return
        typ -- the type of worksheets to search through (active, archived, trashed, etc.)
    */
    if (event && event.which == 13)
        search_worksheets(typ);
    else
        return true;
}


function search_worksheets(typ) {
    /*
    Send a request back to the server and open in the current window the results
    of doing a full-text search through all worksheets of a given type.

    INPUT:
        typ -- the type of worksheets (active, archived, etc.)
    */
    X = get_element('search_worksheets');
    url = '?typ=' + typ + '&search=' + escape0(X.value);
    window.location.replace(url);
}

function go_system_select(theform, original_system) {
    /*
    Switch the current input system from one system to another (e.g., form Sage to Pari or Python).
    A confirmation box is displayed.

    INPUT:
        theform -- the drop down with the list of systems
        original_system -- the system we're switching *from*, so that we can
                           change back to it if the user decides not to do
                           the switch.
    */
    with(theform) {
        var system = options[selectedIndex].value;
        system_select(system);
/*        if (confirm("All cells will be evaluated using " + system + " until you change the system back.")) {
            system_select(system);
        } else {
            options[original_system].selected = 1;
        }
*/
    }
}

function system_select(s) {
    /*
    Send a message back to the server stating that we're switching to evaluating
    all cells using the new system s.
    INPUT:
        s -- a string
    */
    async_request(worksheet_command('system/'+s));
}

function pretty_print_check(s) {
    /*
    Send a message back to the server either turn pretty typeset printing on or off.
    INPUT:
        s -- true or false; true if the pretty print selection is checked.
    */
    async_request(worksheet_command('pretty_print/'+s));
}



function handle_data_menu(theform) {
    /*
    Handle what happens when the user clicks on the worksheet data
    menu and selects an option.

    INPUT:
        theform -- the form in the worksheet with the data drop down menu.
    */
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

function delete_worksheet(name) {
    /*
    Send the worksheet with the given name to the trash.
    INPUT:
        name -- string
    */
    async_request('/send_to_trash', delete_worksheet_callback,
                  {filename: name});
}

function delete_worksheet_callback(status, response_text) {
    /*
    Replace the current page by a page that shows the worksheet in the trash,
    or if the delete worksheet function failed display an error.
    */
    if (status == "success") {
        window.location.replace("/?typ=trash");

    } else {
        alert("Possible failure deleting worksheet.");
    }
}

///////////////////////////////////////////////////////////////////
//
// WORKSHEET list functions -- i.e., functions on a specific
// worksheet in the list of worksheets display.
//
///////////////////////////////////////////////////////////////////

function go_option(theform) {
    /*
    This is called when the user selects a menu item.  This just
    evaluates the corresponding value, which results in running some
    javascript code to do the action.
    */
    with(theform) {
        eval(options[selectedIndex].value);
        options[0].selected = 1;
    }
}

function link_datafile(target_worksheet_filename, filename) {
    /*
    Tell the server to create a symbolic link from the given data file
    to the target worksheet.  This is used to share data
    between multiple worksheets.

    INPUT:
        target_worksheet_filename -- string; the name of the worksheet to link this file to
        filename -- string; the name of this file
    */
    open(worksheet_command("link_datafile?filename=" + escape0(filename) +
         "&target="+escape0(target_worksheet_filename)), process=false);
}


function list_rename_worksheet(filename, curname) {
    /*
    Prompt for a new worksheet name, then send the requested new
    name back to the server, thus changing the worksheet name.

    INPUT:
        filename -- string; the filename of this worksheet to rename
        curname -- string; the current name of this worksheet
    */
    var new_name = prompt('Enter new worksheet name:', curname);
    async_request('/home/' + filename + '/' + 'rename',
        refresh, {name: new_name});
}


function list_edit_worksheet(filename) {
    /*
    In the list of all worksheets, when the user selects "Edit" from
    the menu to edit a worksheet, this function is called, which
    simply loads that worksheet.

    INPUT:
        filename -- string
    */
    window.location.replace('/home/' + filename);
}

function list_copy_worksheet(filename) {
    /*
    When the user selects "Copy" from the list of worksheets, this
    function is called.  It simply sends a message back to the server
    asking that a copy of the worksheet is made.  The worksheet
    list is then refreshed.

    INPUT:
        filename -- string; filename of the worksheet to share
    */
    async_request('/home/' + filename + '/copy?no_load', refresh);
}

function list_share_worksheet(filename) {
    /*
    Bring up menu that allows one to share the selected worksheet from
    the worksheet list with other users.

    INPUT:
        filename -- string; filename of the worksheet to share
    */
    window.location.replace('/home/' + filename + '/share');
}

function list_publish_worksheet(filename) {
    /*
    Publish the given worksheet, when this is selected from the
    worksheet list, and popup the published worksheet.

    INPUT:
        filename -- string; filename of the worksheet to share
    */
    window.open('/home/' + filename + '/publish', "",
             "menubar=1,scrollbars=1,width=800,height=600,toolbar=1,  resizable=1");
}

function list_revisions_of_worksheet(filename) {
    /*
    Display all revisions of the selected worksheet.  This brings
    up the revisions browser.

    INPUT:
        filename -- string; filename of the worksheet to share
    */
    window.location.replace('/home/' + filename + '/revisions');
}


///////////////////////////////////////////////////////////////////
//
// Server pinging support, so server knows page is being viewed.
//
///////////////////////////////////////////////////////////////////

function server_ping_while_alive() {
    /*
    Ping the server every server_ping_time milliseconds to announce
    that we are still viewing this page.
    */
    async_request(worksheet_command('alive'), server_ping_while_alive_callback);
    setTimeout("server_ping_while_alive();", server_ping_time);
}

function server_ping_while_alive_callback(status, response_text) {
    /*
    Whenever the server ping callback occurs, this function runs, and
    if the server didn't respond it calls server_down(); otherwise it
    calls server_up().
    */
    if (status == "failure") {
        server_down();
    } else {
        server_up();
    }
}

function server_down() {
    /*
    Display a warning that indicates that the server is down (this is
    done via CSS).
    */
    set_class("ping", "pingdown");
}

function server_up() {
    /*
    Set the warning message so it's in server up state (this is done via CSS).
    */
    set_class("ping", "ping");
}

///////////////////////////////////////////////////////////////////
//
// CELL functions -- for the individual cells
//
///////////////////////////////////////////////////////////////////
var cell_element_cache = [];
function get_cell(id) {
    /* Return the input cell as a DOM element with the given integer id.
    INPUT:
        id -- integer
    OUTPUT:
        a DOM element
    GLOBAL INPUT:
        cell_element_cache -- an associative array that maps ids to elements
    */
   var v = cell_element_cache[id];
   if(v == undefined)
       v = cell_element_cache[id] = get_element('cell_input_'+ id);
   return v;
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

    // Set the style back to the non-active input cell
    cell.className = 'cell_input';

    // If cell starts with %hide, hide the input.
    var v = lstrip(cell.value).slice(0,5);
    if (v == '%hide') {
        cell.className = 'cell_input_hide';
        cell.style.height = '1em';
    }

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

    // When the input changes we set the CSS to indicate that
    // the cell with this new text has not been evaluated.
    cell_set_not_evaluated(id);

    async_request(worksheet_command('eval'), null,
                  {save_only: 1, id: id, input: cell.value});
}

function evaluate_text_cell_input(id,value, settings) {
    /* Send the input text of the current cell back to the server.
    INPUT:
       id -- The id of the cell
       value -- The new text for the cell
    OUTPUT:
       makes an async call back to the server sending the input text.
    */
    async_request(worksheet_command('eval'), evaluate_text_cell_callback,
                  {text_only: 1, id: id, input: value});
}

function evaluate_text_cell_callback(status, response_text) {
    /*
    Display the new content of a text cell, parsing for math if needed.

    INPUT:
        response_text -- string that is of the form
             [id][cell_html]
             id -- string (integer) of the current text cell

             cell_html -- the html to put in the cell
    */
    if (status == "failure") {
        // Failure evaluating a cell.
        return;
    }
    var X = response_text.split(SEP);
    if (X[0] == '-1') {
        // something went wrong -- i.e., the requested cell doesn't exist.
        alert("You requested to evaluate a cell that, for some reason, the server is unaware of.");
        return;
    }
    id = X[0];
    text = X[1];
    text_cell = get_element('cell_text_'+id);
    var new_html =separate_script_tags(text);
    $(text_cell).replaceWith(new_html[0]);
    // Need to get the new text cell.
    text_cell = get_element('cell_text_'+id);
    setTimeout(new_html[1], 50);

    if (contains_jsmath(text)) {
        try {
            jsMath.ProcessBeforeShowing(text_cell);
        } catch(e) {
            text_cell.innerHTML = jsmath_font_msg + text_cell.innerHTML;
        }
    }

}

function debug_focus() {
    /*
    Called when the Javascript debugging window gets focus.  This window
    is displayed when the notebook server is run with the show_debug option.
    */
    in_debug_input = true;
    w = get_element('debug_window');
    if(w)
       w.className = 'debug_window_active';
}

function debug_blur() {
    /*
    Called when the Javascript debugging window looses focus.
    */
    in_debug_input = false;
    w = get_element('debug_window');
    if(w)
        w.className = 'debug_window_inactive';
}

function cell_focus(id, leave_cursor) {
    /*
    Set the focus on the cell with the given id.

    INPUT:
        id -- integer; id of a cell
        leave_cursor -- if true, do not move the cursor to the top
                        left of the input cell
    OUTPUT:
        focuses on a given cell, possibly moves the cursor, and sets
        the global variable cell_has_changed to false, since we have
        just entered this cell and it hasn't been changed yet.
    */
    var cell = get_cell(id);
    if (cell) {

        // focus on the cell with the given id and resize it
        cell_input_resize(id);
        cell.focus();

        // Possibly also move the cursor to the top left in this cell.
        if (!leave_cursor)
            move_cursor_to_top_of_cell(cell);
    }
    // set since we're now in a new cell, whose state hasn't changed yet.
    cell_has_changed = false;
    return true;
}

function cell_focused(cell, id) {
    /*
    This function is called when the cell gets focus.  It sets the
    CSS so that the input part of the cell is highlighted.

    INPUT:
        cell -- DOM element
        id -- integer

    OUTPUT:
        sets the global variable current_cell and update display of the evaluate link.
    */
    cell.className = "cell_input_active";

    // This makes sure the input textarea is resized right when it is
    // clicked on.
    cell_input_resize(id);

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

function focus_delay(id, leave_cursor) {
    /*
    Set the focus on the cell with given input after a
    10 milliseconds delay.

    NOTE: I (Stein) think this sort of use of Javascript is really bad, and code
    that uses this should be rewritten to not use it.

    INPUT:
        id -- an integer
        leave_cursor -- see docs for the cell_focus function
    */
    setTimeout('cell_focus('+id+','+leave_cursor+')', 10);
}


function cell_input_resize(id) {
    /*
    Resize the given input cell so that it has the right
    number of rows for the number of rows currently typed into
    the input cell.

    INPUT:
        id -- a cell id
    OUTPUT:
        changes the height of the corresponding DOM object to fit the input

    ALGORITHM:
    Create a hidden div with the same style as the textarea, then copy
    all the text into it, set the height of the textarea in pixels
    based on the height of the div, then delete the div.
    */

    var resizer = get_element('cell_resizer');
    var cell_input = get_cell(id);
    resizer.style.width = cell_input.offsetWidth + 'px';
    resizer.innerHTML = cell_input.value.replace(/&/g,'&amp;').replace(/</g,'&lt;').replace(/\r?\n/g,'<br>').replace(/\s\s/g,' &nbsp;') + '&nbsp;';
    cell_input.style.height = resizer.offsetHeight + 'px';

    if(slide_hidden) {
        cell_input.className="cell_input_active";
        slide_hidden = false;
    }
    return;
}

function cell_delete(id) {
    /*
    Send a request back to the server that we would like to delete the
    cell with given id.

    INPUT:
        id -- an integer
    */
    if (active_cell_list.indexOf(id) != -1) {
        // Deleting a running cell causes evaluation to be interrupted.
        // In most cases this avoids potentially tons of confusion.
        async_request(worksheet_command('interrupt'));
    }
    async_request(worksheet_command('delete_cell'), cell_delete_callback,
                  {id: id});
}

function cell_delete_callback(status, response_text) {
    /*
    When a cell is deleted this callback is called after the
    server hopefully does the deletion.  This function then
    removes the cell from the DOM and cell_id_list.

    INPUT:
        response_text -- [command]SEP[id]
               command -- string (empty or 'ignore')
               id -- id of cell being deleted.
    */
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

    /*
    If we are in slide mode, we call slide_mode() again
    to recalculate the slides.
    */
    if (in_slide_mode) {
        current_cell = -1;
        slide_mode();
    }
}

function debug_input_key_event(e) {
    /*
    Handle an input key even when we're in debug mode.
    INPUT:
        e -- the key event
    */
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
    if (browser_iphone) return;
    cell_input = get_cell(id);
    e = new key_event(e);
    if (e==null) return;

    /*********** SPLIT AND JOIN HANDLING ********/

    // Record that just the control key was pressed.  We do this since on Opera
    // it is the only way to catch control + key.
    if (key_control(e)) {
        control_key_pressed = 1;
        return;
    }
    // Check for the split and join keystrokes.
    // The extra control_key_pressed cases are needed for Safari.
    if (key_split_cell(e) || (key_split_cell_noctrl(e) && control_key_pressed)) {
        doing_split_eval = false;
        split_cell(id);
        return false;
    } else if (key_spliteval_cell(e) || (key_enter(e) && control_key_pressed)) {
        doing_split_eval = true;
        jump_to_cell(id, 1);
        control_key_pressed = 0;
        split_cell(id);
        return false;
    } else if (key_join_cell(e) || (key_delete_cell(e) && control_key_pressed) ||
                                   (key_delete_cell(e) && is_whitespace(get_cell(id).value))) {
        control_key_pressed = 0;
        join_cell(id);
        return false;
    }

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
       doing_split_eval = false;
       evaluate_cell(id, false);
       return false;
    } else if (key_send_input_newcell(e)) {
       doing_split_eval = false;
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



    // An actual non-controlling character was sent, which means this cell has changed.
    // When the cursor leaves the cell, we'll use this to know to send the changed
    // version back to the server.
    // We do still have to account for the arrow keys which don't change the text.
    if (! (key_up_arrow(e) || key_down_arrow(e) || key_menu_right(e) || key_menu_left(e)) )
        cell_has_changed = true;
    return true;
}

function id_of_cell_delta(id, delta) {
    /*
    Return the id of the cell that is delta positions from the cell
    with given id, where delta is an integer, either positive,
    negative, or 0.

    INPUT:
        id -- integer
        delta -- integer
    */
    if (cell_id_list.length == 0) {
        /* alert("bug -- no cells."); */
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
    /*
    Clear the debug window.
    */
    output = get_element("debug_output");
    if(output == null) return;
    output.innerHTML = "";
}

function debug_append(txt) {
    /*
    Append output to the debug window.

    INPUT:
        txt -- a string.
    */
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
                   rather than the beginning
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
    /*
    Escape the string for sending via a URL; also replace all +'s by %2B.

    INPUT:
        input -- a string
    OUTPUT:
        a string
    */
    input = escape(input);
    input = input.replace(/\+/g,"%2B");
    return input;
}

function text_cursor_split(cell) {
    /*
    Returns a pair of substrings, the first from the start of the cell
    to the cursor position, and the second from the cursor position to
    the end of the cell.

    INPUT:
        cell -- an input cell (DOM textarea element)
    OUTPUT:
        Javascript Array with two strings in it.
    */
    var a,b;
    var R = get_selection_range(cell);
    b = cell.value.substr(0,R[1]);
    a = cell.value.substr(b.length);
    return new Array(b,a);
}

function indent_cell(cell) {
    /*
    Indent all the highlighted text in the given input cell by 4 spaces.

    INPUT:
        cell -- an input cell (DOM textarea element)
    */
    var R = get_selection_range(cell);
    var start = 1+cell.value.lastIndexOf("\n", R[0]);
    var a = cell.value.substring(0, start);
    var b = cell.value.substring(start, R[1]);
    var c = cell.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++)
        lines[i] = "    "+lines[i];
    b = lines.join("\n");
    cell.value = a+b+c;
    set_selection_range(cell, a.length, a.length+b.length);
}

function unindent_cell(cell) {
    /*
    Unindent all the highlighted text in the given input cell by 4 spaces.

    INPUT:
        cell -- an input cell (DOM textarea element)
    */
    var R = get_selection_range(cell);
    var start = 1+cell.value.lastIndexOf("\n", R[0]);
    var a = cell.value.substring(0, start);
    var b = cell.value.substring(start, R[1]);
    var c = cell.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++)
        lines[i] = unindent_pat.exec(lines[i])[1];  //square brackets pull the captured pattern
    b = lines.join("\n");
    cell.value = a+b+c;
    set_selection_range(cell, a.length, a.length+b.length);
}

function comment_cell(cell) {
    /*
    Comment out all the highlighted (selected) text in the given input cell.

    INPUT:
        cell -- an input cell (DOM textarea element)
    */
    var R = get_selection_range(cell);
    if(R[0] == R[1]) return true;
    var start = 1+cell.value.lastIndexOf("\n", R[0]);
    var a = cell.value.substring(0, start);
    var b = cell.value.substring(start, R[1]);
    var c = cell.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++)
        lines[i] = "#"+lines[i];
    b = lines.join("\n");
    cell.value = a+b+c;
    set_selection_range(cell, a.length, a.length+b.length);
}

function uncomment_cell(cell) {
    /*
    Uncomment the highlighted (selected) text in the given input cell.

    INPUT:
        cell -- an input cell (DOM textarea element)
    */
    var R = get_selection_range(cell);
    if(R[0] == R[1]) return true;
    var start = 1+cell.value.lastIndexOf("\n", R[0]);
    var a = cell.value.substring(0, start);
    var b = cell.value.substring(start, R[1]);
    var c = cell.value.substring(R[1]);
    var lines = b.split("\n");
    for(var i = 0; i < lines.length; i++){
        m = uncomment_pat.exec(lines[i]);
        lines[i] = m[1]+m[2];
    }
    b = lines.join("\n");
    cell.value = a+b+c;
    set_selection_range(cell, a.length, a.length+b.length);
}

function join_cell(id) {
    /*
    Join the cell with given id to the cell before it.

    The output of the resulting joined cells is the output of the
    second cell, *unless* the input of the second cell is only
    whitespace, in which case the output is the output of the first
    cell.  We do this since a common way to delete a cell is to empty
    its input, then hit backspace.  It would be very confusing if the
    output of the second cell were retained.  WARNING: Backspace
    on the first cell if empty deletes it.

    INPUT:
        id -- integer cell id.
    OUTPUT:
        change the state of the worksheet in the DOM, global variables,
        etc., and updates the server on this change.
    */
    var id_prev = id_of_cell_delta(id, -1);
    var cell = get_cell(id);

    // The top cell is a special case.  Here we delete the top cell
    // if it is empty.  Otherwise, we simply return doing nothing.
    if(id_prev == id) {
        // yes, top cell
        if (is_whitespace(cell.value)) {
            // Special case -- deleting the first cell in a worksheet and its whitespace
            // get next cell
            var cell_next = get_cell(id_of_cell_delta(id,1));
            // put cursor on next one
            cell_next.focus();
            // delete this cell
            cell_delete(id);
            return;
        } else {
            return;
        }
    }

    var cell_prev = get_cell(id_prev);

    // We delete the cell above the cell with given id except in the
    // one case when the cell with id has empty input, in which case
    // we just delete that cell.
    if (is_whitespace(cell.value)) {
        cell_prev.focus();
        cell_delete(id);
        return;
    }


    // The lower cell in the join is now not empty.  So we
    // delete the previous cell and put its contents in the
    // bottom cell.
    var val_prev = cell_prev.value;

    cell.focus();
    cell_delete(id_prev);

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
    cell_input_resize(id);
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
    cell_input_resize(id);
    send_cell_input(id);  /* tell the server about how the input just got split in half. */

    set_cursor_position(cell,0);

    /* Make sure that the cursor doesn't move to the new cell. */
    ignore_next_jump = true;
    insert_new_cell_before(id,txt[0]);
}


function worksheet_command(cmd) {
    /*
    Create a string formatted as a URL to send back to the server
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
                   true -- do insert a new cell
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
    // possibly occurring.  Note that active_cell_list is a global variable.
    active_cell_list = active_cell_list.concat([id]);

    // Stop from sending the input again to the server when we leave focus and the
    // send_cell_input function is called.
    cell_has_changed = false;

    // Clear the output text and set the CSS to indicate that this
    // is a running cell.
    cell_set_running(id);

    // Finally make the request back to the server to do the actual calculation.
    var cell_input = get_cell(id);
    if (newcell) { newcell = 1; } else { newcell = 0; }
    async_request(worksheet_command('eval'), evaluate_cell_callback,
            {newcell: newcell, id: id, input: cell_input.value});
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
    active_cell_list = active_cell_list.concat([id]);
    cell_set_running(id);
    async_request(worksheet_command('introspect'), evaluate_cell_callback,
          {id: id, before_cursor: before, after_cursor: after});
}

function evaluate_cell_callback(status, response_text) {
    /*
    Update the focus and possibly add a new cell.  If evaluate all
    has been clicked, start evaluating the next cell (and don't
    add a new cell).

    INPUT:
        response_text -- string that is of the form
             [id][command][new_html][new_cell_id]
             id -- string (integer) current cell
             command -- string 'append_new_cell' or 'insert_cell' or 'no_new_cell' or 'introspect'
             new_html -- string
             new_cell_id -- optional (if command is 'insert_cell')
                            string (integer); id of new cell to create
    */
    if (status == "failure") {
        // Failure evaluating a cell.
        return;
    }
    var X = response_text.split(SEP);
    if (X[0] == '-1') {
        // something went wrong -- i.e., the requested cell doesn't exist.
        alert("You requested to evaluate a cell that, for some reason, the server is unaware of.");
        return;
    }

    if (evaluating_all) {
        if(evaluating_all_cursor >= cell_id_list.length) {
            evaluating_all = false;
        } else {
            evaluate_cell(cell_id_list[evaluating_all_cursor], false);
            evaluating_all_cursor++;
        }
    } else if (X[1] == 'append_new_cell') {
        // add a new cell to the very end
        append_new_cell(X[0],X[2]);
    } else if (X[1] == 'insert_cell') {
        // insert a new cell after the one with id X[3]
        do_insert_new_cell_after(X[3], X[0], X[2]);
        jump_to_cell(X[0],0);
    } else if (X[1] != 'introspect' && !in_slide_mode && !doing_split_eval) {
        // move to the next cell after the one that we just evaluated.
        if (is_interacting_cell(current_cell)) {
            jump_to_cell(current_cell);
        } else {
            jump_to_cell(current_cell, 1);
        }
    }
    start_update_check();
}

function is_interacting_cell(id) {
    /*
    Return true if the cell with given id is currently an @interact cell.

    INPUT:
        id -- an integer
    */
    return (get_element("cell-interact-" + id) != null);
}

function cell_output_set_type(id, typ, do_async) {
    /*
    Set the output type of the cell with given id.

    INPUT:
        id -- integer
        typ -- 'wrap', 'nowrap', 'hidden'
        do_async -- true or false; if true tell the server about the change.
    */

    // We do the following specifically because interact cells do not work
    // at all when displayed in nowrap mode, which is VERY BAD.  So instead
    // for interacts one gets a toggle to and from hidden.

    if (typ=="nowrap" && is_interacting_cell(id)) {
        /* if the type is nowrap and the cell-interact-[id] div exists (i.e., we are interacting)
           then just make the thing hidden. */
        typ = "hidden";
    }

    /* OK, now set the sell output type.  */

    set_class('cell_div_output_' + id,    'cell_div_output_' + typ)
    set_class('cell_output_' + id,        'cell_output_' + typ)
    set_class('cell_output_nowrap_' + id, 'cell_output_nowrap_' + typ)
    set_class('cell_output_html_' + id,   'cell_output_html_' + typ)

    // Do async request back to the server
    if(do_async)
        async_request(worksheet_command('set_cell_output_type'), null, {id: id, type: typ})
}

function cycle_cell_output_type(id) {
    /*
    When called the cell with given id has its output cycled from one type to the next.
    There are three types: word wrap, no word wrap, hidden.

    INPUT:
        id -- an integer
    */
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
    /*
    Set the cell with given id to be evaluated.  This is purely a CSS style setting.

    INPUT:
        id -- an integer
    */
    var D = get_element('cell_'+id);
    D.className = "cell_evaluated";
}

function cell_set_not_evaluated(id) {
    /*
    Set the cell with given id to be not evaluated.  This is purely a CSS style setting.

    INPUT:
        id -- an integer
    */
    var D = get_element('cell_'+id);
    D.className = "cell_not_evaluated";
    cell_set_done(id);
}

function cell_set_running(id) {
    /*
    Start the cell with given id running -- this is purely a style and
    content; the server is not contacted by this function.

    INPUT:
        id -- an integer
    */
    // Blank the output text
    set_output_text(id, '', '', '', '', '', 1);   // the 1 means not @interact

    // If the output type is hidden, toggle it to be visible.  Otherwise
    // we leave it alone.
    if (get_element('cell_div_output_' + id).className == 'cell_div_output_hidden') {
        cycle_cell_output_type(id);
    }

    // Set the CSS
    var cell_div = get_element('cell_div_output_' + id);
    cell_div.className = 'cell_output_running';
    var cell_number = get_element('cell_number_' + id);
    cell_number.className = 'cell_number_running';
}

function cell_set_done(id) {
    /*
    Change the CSS for the cell with given id to indicate
    that it is no longer computing.

    INPUT:
        id -- integer
    */
    var cell_div = get_element('cell_div_output_' + id)
    cell_div.className = 'cell_div_output_wrap';
    var cell_number = get_element('cell_number_' + id);
    cell_number.className = 'cell_number';
}

function check_for_cell_update() {
    /*
    Ask the server if there is any new output that should be placed
    in an output cell.

    OUTPUT:
        * if the active cell list is empty, cancel update checking.
        * makes an async request
        * causes the title bar compute spinner to spin
    */

    // cancel update checks if no cells are doing computations.
    if (active_cell_list.length == 0) {
        cancel_update_check();
        return;
    }

    // record in a global variable when the last update occurred.
    update_time = time_now();

    // check on the cell currently computing to see what's up.
    var cell_id = active_cell_list[0];
    async_request(worksheet_command('cell_update'),
                    check_for_cell_update_callback,
                    {id: cell_id});

    // spin the little title spinner in the title bar.
    try{
        title_spinner_i = (title_spinner_i+1)%title_spinner.length;
        document.title = title_spinner[title_spinner_i] + original_title;
    } catch(e){}
}

function check_for_cell_update_callback(status, response_text) {
    /*
    Callback after the server responds for our request for updates.

    INPUT:
        status -- string
        responese_test -- string that encodes three variables, with this format (no []'s):
[status (1-letter)][id] [output_text]SEP[output_text_wrapped]SEP[output_html]SEP[new_cell_input]SEP[interrupted]SEP[introspect_html]
             status --    'e' -- empty; no more cells in the queue
                          'd' -- done; actively computing cell just finished
                          'w' -- still working
             id -- string; an integer (encoded as a decimal string)
             output_text -- string; the output text so far
             output_text_wrapped -- string; word wrapped version of output text
             output_html -- string; html output
             new_cell_input -- string; if the input to the cell should be changed (e.g., when
                               doing a tab completion), this gives the new input
             interrupted -- string ('restart' or 'false');
                            whether the computation of this cell was interrupted and
                            if so why.
             introspect_html -- new introspection html to be placed in the introspection window
    */
    // make sure the update happens again in a few hundred milliseconds,
    // unless a problem occurs below.

    if (status != "success") {
        // a problem occurs -- stop trying to evaluate.
        if(update_error_count>update_error_threshold) {
            cancel_update_check();
            halt_active_cells();
            var elapsed_time = update_error_count*update_error_delta/1000;
            var msg = "Error updating cell output after " + elapsed_time + "s";
            msg += "(canceling further update checks).";
            /* alert(msg); */
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

    if(response_text == 'empty') {
        // if the server returns nothing, we just ignore that
        // response and try again later.
       continue_update_check();
       return;
    }

    var i = response_text.indexOf(' ');
    var id = response_text.substring(1, i);
    var stat = response_text.substring(0,1)

    if(stat == 'e') {
        cancel_update_check();
        halt_active_cells();
        return;
    }

    // compute output for a cell
    var D = response_text.slice(i+1).split(SEP);
    var output_text = D[0] + ' ';
    var output_text_wrapped = D[1] + ' ';
    var output_html = D[2];
    var new_cell_input = D[3];
    var interrupted = D[4];
    var introspect_html = D[5];
    var j = id_of_cell_delta(id,1);

    // Evaluate javascript, but *only* after the entire
    // cell output has been loaded (hence the stat == 'd') below.
    var cell_is_not_an_interact_update = ! get_element("cell-interact-" + id);
    if (stat == 'd' && cell_is_not_an_interact_update) {
        output_text_wrapped = eval_script_tags(output_text_wrapped);
        output_html = eval_script_tags(output_html);
    }

    // Set the latest output text got from the server.
    set_output_text(id, output_text, output_text_wrapped,
                    output_html, stat, introspect_html);

    if (stat == 'd') {
        active_cell_list = delete_from_array(active_cell_list, id);

        if (interrupted == 'restart') {
            restart_sage();
        } else if (interrupted == 'false') {
            cell_set_evaluated(id);
        } else {
            cancel_update_check();
            halt_active_cells();
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
        if(update_count > update_falloff_threshold &&
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
    /*
    If enough time has elapsed, check for more output from the server.
    If not, wait longer and try again later.

    GLOBAL INPUT:
        update_time -- global variable that records when last update check occurred.
    */
    var time_elapsed = time_now() - update_time;
    if(time_elapsed < cell_output_delta) {
        update_timeout = setTimeout('check_for_cell_update()', cell_output_delta-time_elapsed);
    } else {
        check_for_cell_update();
    }
}

function start_update_check() {
    /*
    Start the updating check system.  This system checks for update
    from the server with an exponential backup strategy.
    */
    if(updating) return;

    // set several global variables that cells are computing so we have to
    // check for output.
    updating = true;
    update_count = 0;
    update_falloff_level = 0;

    // The starting value for how long we wait between checks for new updates.
    cell_output_delta = update_falloff_deltas[0];

    // Do one initial check without waiting, since some calculations are very
    // fast and doing this feels snappy.
    check_for_cell_update();
}

function cancel_update_check() {
    /*
    Turn off the loop that checks for now output updates in the worksheet.

    This just cancels the updating timer and gets rid of the spinning
    "active" indicator in the title bar.
    */
    updating = false;
    clearTimeout(update_timeout);
    document.title = original_title;
}

function contains_jsmath(text) {
    /*
    Returns true if text contains some jsmath text.  This function
    sucks, since it really just looks for class="math" and is easy
    to throw off.  Fix this!

    INPUT:
        text -- a string
    */
    // TODO: should make this not case sensitive!!  how to .lower() in javascript?
    // TODO: Or write this using a regexp.
    return (text.indexOf('class="math"') != -1 || text.indexOf("class='math'") != -1);
}

function set_output_text(id, text, wrapped_text, output_html,
                         status, introspect_html, no_interact) {
    /*
    Set the output text for the cell with given id to the given text.

    INPUT:
        id -- an integer; the id of a cell
        text -- string
        wrapped_text -- string; word wrapped version of text
        output_html -- string; html formatted output
        status -- letter (length 1 string); 'd' -- done; anything else -- working
        introspect_html -- when user is introspecting this html will go in
                           the introspection dialog box
        no_interact -- true or false; if true then this is not an @interact output.
    */
    if (id < 0) {
        // negative id's come up for special internal usage, and should be ignored.
        return;
    }
    var cell_interact = get_element("cell-interact-" + id);
    if (!no_interact && cell_interact) {
        // Uncomment to change so that only show output at the end.
        if (status  != 'd') return;
        var i = wrapped_text.indexOf('<?__SAGE__START>');
        var j = wrapped_text.indexOf('<?__SAGE__END>');
        if (i == -1 || j == -1) {
            /* alert("Bug in notebook -- interact wrapped text is invalid" + wrapped_text); */
            return;
        }

        var new_interact_output = wrapped_text.slice(i+16,j);
        new_interact_output = eval_script_tags(new_interact_output);

        // An error occurred accessing the data for this cell.  Just force reload
        // of the cell, which will certainly define that data.
        if (new_interact_output.indexOf('__SAGE_INTERACT_RESTART__') != -1) {
            evaluate_cell(id, 0);
        } else {
            cell_interact.innerHTML = new_interact_output;
            if (contains_jsmath(new_interact_output)) {
               jsMath.ProcessBeforeShowing(cell_interact);
            }
        }
    } else {
        // fill in output text got so far
        var cell_output = get_element('cell_output_' + id);
        if (!cell_output) {
            // This can happen, e.g., if a cell is deleted from the DOM, but
            // the server has some output it still wants to put in the cell.
            // This happens because once a cell is running there is no stopping
            // it beyond an explicit interrupt (since interrupt may or may not
            // succeed -- this is the real world with hard to kill C code, etc.).
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
                 // This is the first time that the underlying Python interact function is
                 // actually called!
                interact(id, 'sage.server.notebook.interact.recompute(' + id + ')');
            }
        }
    }

    if (status == 'd') {
         cell_set_done(id);
         if (contains_jsmath(text)) {
             try {
                 jsMath.ProcessBeforeShowing(cell_output);
             } catch(e) {
                 cell_output.innerHTML = jsmath_font_msg + cell_output.innerHTML;
                 cell_output_nowrap.innerHTML = jsmath_font_msg + cell_output_nowrap.innerHTML;
             }
         }
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
        if (contains_jsmath(introspect_html)) {
            try {
                jsMath.ProcessBeforeShowing(cell_output_html);
            } catch(e) {
                cell_output.innerHTML = jsmath_font_msg + cell_output_html.innerHTML;
            }
        }
    }
}

function set_input_text(id, text) {
    /*
    Fill in input text for the cell with given id.  This is used by
    the tab completion system, so it also sets the cell with given id
    to be in focus and positions the cursor in exactly the right spot.

    INPUT:
        id -- an integer
        text -- a string
    */
    var cell_input = get_cell(id);
    cell_input.value = text;

    jump_to_cell(id,0)
    pos = text.length - after_cursor.length;
    set_cursor_position(cell_input, pos);

    return false;
}


///////////////////////////////////////////////////////////////////
// Dynamic evaluation of javascript related in cell output.
///////////////////////////////////////////////////////////////////

function CellWriter() {
    /*
    When a new cell is loaded, this class is used to let javascript
    write directly to the document. After that, make sure javascript
    writes to a CellWriter object.  This is used in order to get jmol
    to work.
    */
    function write(s) {
        this.buffer += s;
    }
    this.write = write;
    this.buffer = "";
}

// The global cell_writer target.
cell_writer = document;

function eval_script_tags(text) {
    /*
    Find all the tags in the given script and eval them, where
    tags are javascript code in <script>...</script> tags.
    This allows us put javascript in the output of computations
    and have it evaluated.

    INPUT:
        text -- a string
    OUTPUT
        string -- like text, but with all script tags removed.
    */
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

function separate_script_tags(text) {
    /*
    Find all the tags in the given script and return a list of two strings.  The first
    string is the html in text, the second is the contents of any <script>...</script> tags.

    INPUT:
        text -- a string
    OUTPUT
        a list of two strings.  The first is the text without script tags,
        the second is the contents of the script tags.
    */

    var script = '';
    var s = text; //text.replaceAll('\n','');
    var i = s.indexOf('<'+'script>');
    while (i != -1) {
        var j = s.indexOf('<'+'/script>');
        script += s.slice(8+i,j);

        s = s.slice(0,i) + s.slice(j+9);
        i = s.indexOf('<'+'script>');
    }
    return [s, script];
}


///////////////////////////////////////////////////////////////////
// Single Cell Functions
///////////////////////////////////////////////////////////////////

function slide_mode() {
    /*
    Switch into single cell mode.
    This involves changing a bunch of CSS and some global variables.
    */
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
    /*
    Switch from single cell mode back to normal worksheet mode.
    This involves changing CSS and global variables.
    */
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
    /*
    Hide the currently displayed slide.
    GLOBAL INPUT:
        current_cell -- integer
    */
    set_class('cell_outer_' + current_cell, 'hidden');
}

function slide_show() {
    /*
    Switch into slide show mode.
    This involves changing a lot of CSS in the DOM.
    */
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
    /*
    Move to the first input cell in single cell mode.
    GLOBAL INPUT:
        the first cell is the first entry in the cell_id_list.
    */
    jump_to_slide(cell_id_list[0]);
}

function slide_last() {
    /*
    Move to the last input cell in single cell mode.
    GLOBAL INPUT:
        the last cell is the last entry in the cell_id_list.
    */
    jump_to_slide(cell_id_list[cell_id_list.length-1]);
}


function slide_next() {
    /*
    Move to the next cell in single cell mode.
    */
    jump_to_slide(id_of_cell_delta(current_cell, 1));
}

function slide_prev() {
    /*
    Move to the previous cell in single cell mode.
    */
    jump_to_slide(id_of_cell_delta(current_cell, -1));
}

function jump_to_slide(id) {
    /*
    Move to the display only the cell with the given id.
    INPUT:
        id -- an integer
    OUTPUT:
        sets the global variable current_cell to the id.
    */
    slide_hide();
    current_cell = id;
    slide_show();
}


function update_slideshow_progress() {
    /*
    There is a bar at the top of the screen that shows how far through
    the worksheet we currently are in the list of cells (in single
    cell mode).  This function updates the CSS of that "progress
    meter" to have the right percentage filled in.
    */
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

function insert_into_array(v, i, x) {
    /*
    This is a generic function that return a new array with x inserted
    into position i of v.

    INPUT:
        v -- an array
        i -- an integer
        x -- object
    OUTPUT:
        an array
    */
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
    /*
    Create a new cell in the DOM with given id and html defining it.    This does
    not send a message back to the server.
    INPUT:
        id -- integer
        html -- string
    */
    var new_html =separate_script_tags(html);
    var new_cell = document.createElement("div");
    var in_cell = document.createElement("div");
    new_cell.appendChild(in_cell);
    new_cell.id = 'cell_outer_' + id;
    in_cell.id = 'cell_' + id;
    in_cell.innerHTML = new_html[0];
    setTimeout(new_html[1], 50);
    return new_cell;
}

function make_new_text_cell(id, html) {
    /*
    Create a new cell in the DOM with given id and html defining it.    This does
    not send a message back to the server.
    INPUT:
        id -- integer
        html -- string
    */
    var new_html =separate_script_tags(html);
    var new_cell = $(new_html[0])
    setTimeout(new_html[1], 50);
    return new_cell;
}


function do_insert_new_cell_before(id, new_id, new_html) {
    /*
    Insert a new cell with the given new_id and new_html
    before the cell with given id.

    INPUT:
        id -- integer; id of a cell
        new_id -- integer
        new_html -- text to put in new cell
    GLOBAL INPUT:
        doing_split_eval -- if true, evaluate both cells id and new_id.
    */
    var new_cell = make_new_cell(new_id, new_html);
    $('#cell_outer_'+id).before(new_cell)

    var i = cell_id_list.indexOf(eval(id));
    cell_id_list = insert_into_array(cell_id_list, i, eval(new_id));

    // Deal with one special case when this function is called, where
    // we evaluate both of the cells that result from a split.
    // It is annoying to put this code here, but it is much simpler
    // than coding it in any other way, because we don't know the new_id
    // until this point.
    if (doing_split_eval) {
        evaluate_cell(id, false);
        evaluate_cell(new_id, false);
    }
}


function do_insert_new_text_cell_before(id, new_id, new_html) {
    /*
    Insert a new cell with the given new_id and new_html
    before the cell with given id.

    INPUT:
        id -- integer; id of a cell
        new_id -- integer
        new_html -- text to put in new cell
    */
    var new_cell = make_new_text_cell(new_id, new_html);
    $('#cell_outer_'+id).before(new_cell)
}



function insert_new_cell_after(id, input) {
    /*
    Send a message back to the server requesting that a new cell with
    the given input be inserted before the cell with given id.

    INPUT:
        id -- an integer
        input -- text (default: "")
    */
    if(input == null) input = "";
    async_request(worksheet_command('new_cell_after'), insert_new_cell_after_callback, {id: id, input: input});
}

function insert_new_cell_after_callback(status, response_text) {
    /*
    Callback that is called when the server inserts a new cell
    after a given cell.

    INPUT:
        response_text -- 'locked': means that the user is not allowed to
                         insert new cells into this worksheet
                      -- or a string that encodes several variables:
                             [new_id]SEP[new_html]SEP[id]
                         where
                            new_id -- string representation of an integer
                            new_html -- new HTML output for new cell
                            id -- id of cell before the new one being inserted
    */

    if (status == "failure") {
        alert("Problem inserting new input cell after current input cell.\n" + response_text);
        return ;
    }
    if (response_text == "locked") {
        alert("Worksheet is locked.  Cannot insert cells.");
        return;
    }

    // Extract the input variables that are encoded in the response_text.
    var X = response_text.split(SEP);
    var new_id = eval(X[0]);
    var new_html = X[1];
    var id = eval(X[2]);

    // Insert a new cell _after_ a cell.
    do_insert_new_cell_after(id, new_id, new_html);
    jump_to_cell(new_id,0);
}

function insert_new_text_cell_after(id, input) {
     /*
    Insert a new text cell after the cell with given id.

    This sends a message to the server requesting that a new cell be
    inserted, then via a callback modifies the DOM.

    INPUT:
        id -- integer
        input -- string
    */
    if(input == null) input = "";
    async_request(worksheet_command('new_text_cell_after'), insert_new_text_cell_after_callback, {id: id, input: input});
}

function insert_new_text_cell_after_callback(status, response_text) {
    /*
    Callback that is called when the server inserts a new cell
    after a given cell.

    INPUT:
        response_text -- 'locked': means that the user is not allowed to
                         insert new cells into this worksheet
                      -- or a string that encodes several variables:
                             [new_id]SEP[new_html]SEP[id]
                         where
                            new_id -- string representation of an integer
                            new_html -- new HTML output for new cell
                            id -- id of cell before the new one being inserted
    */
    if (status == "failure") {
        alert("Problem inserting new text cell before current input cell.");
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
    do_insert_new_text_cell_after(id, new_id, new_html);
}

function do_insert_new_cell_after(id, new_id, new_html) {
    /*
    Insert a new cell with the given new_id and new_html after the
    cell with given id.

    INPUT:
        id -- integer
        new_id -- integer
        new_html -- string
    */
    // Find the cell id of the cell right after the cell with id.
    i = id_of_cell_delta(id,1);

    if(i == id) {
        // i is the last cell.
        append_new_cell(new_id,new_html);
    } else {
        do_insert_new_cell_before(i, new_id, new_html);
    }
}

function do_insert_new_text_cell_after(id, new_id, new_html) {
    /*
    Insert a new text cell with the given new_id and new_html after the
    cell with given id.

    INPUT:
        id -- integer
        new_id -- integer
        new_html -- string
    */
    // Find the cell id of the cell right after the cell with id.
    i = id_of_cell_delta(id,1);

    if(i == id) {
        // i is the last cell.
        append_new_text_cell(new_id,new_html);
    } else {
        do_insert_new_text_cell_before(i, new_id, new_html);
    }
}


function insert_new_cell_before(id, input) {
    /*
    Insert a new cell before the cell with given id.

    This sends a message to the server requesting that a new cell be inserted, then
    via a callback modifies the DOM.

    INPUT:
        id -- integer
        input -- string
    */
    if(input == null) input = "";
    async_request(worksheet_command('new_cell_before'), insert_new_cell_before_callback, {id: id, input: input});
}

function insert_new_cell_before_callback(status, response_text) {
    /*
    See the documentation for insert_new_cell_after_callback, since
    response_text is encoded in exactly the same way there.
    */
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

function insert_new_text_cell_before(id, input) {
     /*
    Insert a new text cell before the cell with given id.

    This sends a message to the server requesting that a new cell be
    inserted, then via a callback modifies the DOM.

    INPUT:
        id -- integer
        input -- string
    */
    if(input == null) input = "";
    async_request(worksheet_command('new_text_cell_before'), insert_new_text_cell_before_callback, {id: id, input: input});
}

function insert_new_text_cell_before_callback(status, response_text) {
    /*
    See the documentation for insert_new_text_cell_after_callback, since
    response_text is encoded in exactly the same way there.
    */
    if (status == "failure") {
        alert("Problem inserting new text cell before current input cell.");
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
    do_insert_new_text_cell_before(id, new_id, new_html);
}


function append_new_cell(id, html) {
    /*
    Append a new cell with given id to the end of the list of cells, then
    position the cursor in that cell.

    This modifies the DOM and nothing else.

    INPUT:
        id -- an integer
        html -- html text that goes in the new cell.
    */
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


function append_new_text_cell(id, html) {
    /*
    Append a new text cell with given id to the end of the list of cells, then
    position the cursor in that cell.

    This modifies the DOM and nothing else.

    INPUT:
        id -- an integer
        html -- html text that goes in the new cell.
    */
    var new_cell = make_new_text_cell(id, html);
    $('#worksheet_cell_list').append(new_cell)
}


///////////////////////////////////////////////////////////////////
//
// CONTROL functions
//
///////////////////////////////////////////////////////////////////

function interrupt() {
    /*
    Send a message to the server that we would like to interrupt
    all running calculations in the worksheet.
    */
    async_request(worksheet_command('interrupt'), interrupt_callback);
}

function interrupt_callback(status, response_text) {
    /*
    Callback called after we send the interrupt signal to the
    server.  If the interrupt succeeds, we change the CSS/DOM
    to indicate that no cells are currently computing.  If
    it fails, we display an annoying alert (it might be a good
    idea to change this, e.g., to blink red or something instead
    of an alert).
    */
    if (response_text == "failed") {
       alert('Unable to immediately interrupt calculation.');
       return;
    } else if(status == "success") {
        halt_active_cells()
    }
}


function evaluate_all() {
    /*
    Iterate through every input cell in the document, in order, and
    evaluate it.  Previously, we just called evaluate on everything
    all at once.  This is undesirable, since packets often arrive
    out-of-order, so the cells get evaluated out-of-order.

    Set the global variable evaluating_all = true.  Then, we kick off
    evaluations by evaluating the first cell.  In cell_evaluate_callback,
    we check to see if evaluating_all is set, and proceed from there.
    This way, each cell is evaluated immediately after the server
    acknowledges that it has received the previous request.

    */

    evaluating_all = true;
    evaluating_all_cursor = 1; //start at 1 since we kick-off with zero
    evaluate_cell(cell_id_list[0],false);
}

function hide_all() {
    /*
    Hide every output cell in the worksheet (via CSS) then send a
    message back to the server recording that we hid all cells, so if
    we refresh the browser or visit the page with another browser,
    etc., the cells are still hidden.
    */
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        cell_output_set_type(v[i],'hidden', false);
    }
    async_request(worksheet_command('hide_all'));
}

function show_all() {
    /*
    Show ever output cell in the worksheet, and send a message to the
    server reporting that we showed every cell.
    */
    var v = cell_id_list;
    var n = v.length;
    var i;
    for(i=0; i<n; i++) {
        cell_output_set_type(v[i],'wrap', false);
    }
    async_request(worksheet_command('show_all'));
}

function delete_all_output() {
    /*
    Delete the contents of every output cell in the worksheet (in the DOM) then
    send a message back to the server recording that we deleted all cells, so
    if we refresh the browser or visit the page with another browser,
    etc., the cells are still deleted.

    Things that could go wrong:
         1. Message to server to actually do the delete is not received or fails.
            Not so bad, since no data is lost; a mild inconvenience.
         2. User accidentally clicks on delete all.  There is no confirm dialog.
            Not so bad, since we save a revision right before the delete all, so
            they can easily go back to the previous version.
    */
    var v = cell_id_list;
    var n = v.length;
    var i, id;
    /* Iterate over each cell in the worksheet. */
    for(i=0; i<n; i++) {
        id = v[i];
        /* First delete the actual test from the output of each cell. */
        get_element('cell_output_' + id).innerHTML = "";
        get_element('cell_output_nowrap_' + id).innerHTML = "";
        get_element('cell_output_html_' + id).innerHTML = "";
        /* Then record that the cell hasn't been evaluated and produced that output. */
        cell_set_not_evaluated(id);
    }
    /* Finally tell the server to do the actual delete.
       We first delete from DOM then contact the server for maximum
       snappiness of the user interface. */
    async_request(worksheet_command('delete_all_output'));
}

function halt_active_cells() {
    /*
    Set all cells so they do not look like they are being evaluates or
    queued up for evaluation, and empty the list of active cells from
    the global active_cell_list variable.
    */
    for(i = 0; i < active_cell_list.length; i++)
        cell_set_not_evaluated(active_cell_list[i]);
    active_cell_list = []
}

function set_all_cells_to_be_not_evaluated() {
    /*
    Change the CSS so that all cells are displayed as
    having not been evaluated.
    */
    for(i = 0; i < cell_id_list.length; i++)
        cell_set_not_evaluated(cell_id_list[i]);
}

function restart_sage() {
    /*
    Restart the running Sage process that supports calculations in this
    worksheet.

    This function immediately changes the DOM so it looks like no cells
    are running and none have been evaluated, then it sends a message
    back to the server requesting that the worksheet Sage process
    actually be stopped.
    */
    halt_active_cells();
    set_all_cells_to_be_not_evaluated();
    async_request(worksheet_command('restart_sage'));
}

function quit_sage() {
    /*
    Called when the worksheet process is terminated.  All actively
    computing cells are stopped, and a request is sent to the server
    to quit the worksheet process.
    */
    halt_active_cells();
    set_all_cells_to_be_not_evaluated();
    async_request(worksheet_command('quit_sage'), restart_sage_callback);
}

function login(username,password) {
    /*
    Set the username and password for this user in a cookie.

    This is called when the user logs in from the login screen to set
    the user's credentials.
    INPUT:
        username -- string
        password -- string
    */
    document.cookie="username="+username;
    document.cookie="password="+password;
    window.location="/";
}

///////////////////////////////////////////////////////////////////
//
// Various POPUP WINDOWS
//
///////////////////////////////////////////////////////////////////

function history_window() {
    /*
    Popup the history window.
    */
    history = window.open ("/history",
      "", "menubar=1,scrollbars=1,width=800,height=600, toolbar=1,resizable=1");
}

function print_worksheet() {
    /*
    Display a version of this worksheet that is suitable for printing.
    */
    log = window.open (worksheet_command("print"),"",
      "menubar=1,scrollbars=1,width=800,height=600,toolbar=1,  resizable=1");
}

function help() {
    /*
    Popup the help window.
    */
    log = window.open ("/help","",
    "menubar=1,location=1,scrollbars=1,width=800,height=650,toolbar=1,  resizable=1");
}

function bugreport() {
    /*
    Popup the bug report window.
    */
    log = window.open ("http://spreadsheets.google.com/viewform?key=pCwvGVwSMxTzT6E2xNdo5fA","",
    "menubar=1,location=1,scrollbars=1,width=800,height=650,toolbar=1,  resizable=1");
}



///////////////////////////////////////////////////////////////////
// Interact
///////////////////////////////////////////////////////////////////

function interact(id, input) {
    /*
    Cancels any current computations, then sends an interact request back to the
    server.  This is called by individual interact controls.

    INPUT:
        id -- an integer
        input -- string
    */
//    async_request(worksheet_command('interrupt'));

    active_cell_list = active_cell_list.concat([id]);
    cell_has_changed = false;
    current_cell = id;

    // Delete the old images, etc., that might be sitting
    // in the output from the previous evaluation of this cell.
    get_element('cell_output_html_' + id).innerHTML = "";

    var cell_number = get_element('cell_number_' + id);
    cell_number.className = 'cell_number_running';

    // the __sage_interact__ string appears also in cell.py
    async_request(worksheet_command('eval'), evaluate_cell_callback,
            {newcell: 0, id: id, input: '%__sage_interact__\n' + input});
}

///////////////////////////////////////////////////////////////////
// Base 64 encoding and decoding (mainly used for @interact).
///////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////
// The following header applies to the encode64 and decode64 functions
// This code was written by Tyler Akins and has been placed in the
// public domain.  It would be nice if you left this header intact.
// Base64 code from Tyler Akins -- http://rumkin.com
//////////////////////////////////////////////////////////////////
var keyStr = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/=";

function encode64(input) {
    /*
    Base 64 encode the given input.
    INPUT:
        input -- string
    OUTPUT:
        string
    */
    // I had to add this, since otherwise when input is numeric there are
    // errors below.
    try {
        input = input.toString();
    } catch(e) {
        return input;
    }
    var output = "";
    var chr1, chr2, chr3;
    var enc1, enc2, enc3, enc4;
    var i = 0;

     while (i < input.length) {
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
    }

    return output;
}

function decode64(input) {
    /*
    Base 64 decode the given input.
    INPUT:
        input -- string
    OUTPUT:
        string
    */
    var output = "";
    var chr1, chr2, chr3;
    var enc1, enc2, enc3, enc4;
    var i = 0;

    // remove all characters that are not A-Z, a-z, 0-9, +, /, or =
    input = input.replace(/[^A-Za-z0-9\+\/\=]/g, "");

    while (i < input.length) {
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
    }

    return output;
}

///////////////////////////////////////////////////////////////////
// Trash
///////////////////////////////////////////////////////////////////

function empty_trash() {
    /*
       This asks for confirmation from the user then sends a request back to the
       server asking that the trash be emptied for this user. The request to the
       server goes by accessing the url /emptytrash.  After that finishes, the
       empty trash folder is displayed.
    */
    if(confirm('Emptying the trash will permanently delete all items in the trash. Continue?')) {
        window.location.replace("/emptytrash");
        window.location.replace("/?typ=trash");
    }
}


/********************* js math ***************************/

function jsmath_init() {
    /*
    Process all the jsmath in this page.
    */
    try {
         jsMath.Process();
    } catch(e) {
    }
}

///////////////////////////////////////////////////////////////////
//
// KeyCodes (auto-generated from config.py and user's sage config
//
///////////////////////////////////////////////////////////////////



{{ KEY_CODES }}


{% include "jmol_lib.js" %}

{% include "canvas3d_lib.js" %}

{% include "async_lib.js" %}
