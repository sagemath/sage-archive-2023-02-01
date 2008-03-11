"""nodoctest
Sage Notebook CSS
"""


#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

import os

from sage.misc.misc import DOT_SAGE

def css(color='default'):
    r"""
    Return the CSS header used by the \sage Notebook.

    INPUT:
        color -- string or pair of html colors, e.g.,
                    'gmail'
                    'grey'
                    \code{('\#ff0000', '\#0000ff')}

    EXAMPLES:
        sage: import sage.server.notebook.css as c
        sage: type(c.css())
        <type 'str'>
    """
    s = r"""
body {
  background-color: white;
}

div.hidden {
  display:none;
}
span.hidden{
  display:none;
}

div.fivepix {
  height:5px;
}

textarea.hidden {
   display:none;
}

pre.hidden {
   display:none;
}

/**** TOP CONTROL BAR ************************/

div.top_control_bar {
   z-index: 0;
   background-color: white;
   position: fixed;
   left: 0px;
   width: 100%;
   top: 0px;
   padding-left: 2ex;
}

span.control_commands {
   position: fixed;
   top:1ex;
   right:1ex;
   text-align:right;
   color:blue;
   font-weight:normal;
   font-family:arial;
   font-size:12px;
/*    text-decoration:underline; */
}

span.worksheet_control_commands {
   position: relative;
   top:0px;
   right:0px;
   text-align:right;
   color:blue;
   font-weight:normal;
   font-family:arial;
   font-size:12px;
   padding:5px;
}

div.slide_control_commands {
   float:right;
   position: fixed;
   width:500px;
   top:1ex;
   right:45%;
   text-align:right;
   color:blue;
   font-weight:normal;
   font-family:arial;
   font-size:12px;
/*    text-decoration:underline; */
}

span.vbar {
   height:1.5ex;
   border-left:1px solid black;
   width:1px;
}

div.top_control_bar a {
   color: #0000BB;
   text-decoration: none;
   padding:5px;
}

div.top_control_bar a:hover {
   cursor:pointer;
}

div.top_control_bar a.upload_worksheet {
}

div.top_control_bar a.worksheets_button {
}

div.top_control_bar a.upload_worksheet:hover {
}

div.top_control_bar a.restart_sage {
}

div.top_control_bar a.restart_sage:hover {
   background-color:#ff0000;
}

div.top_control_bar  a.restart_sage_in_progress {
   background-color:#ff0000;
   text-decoration:blink;
}

div.top_control_bar  a.interrupt {
/*   text-decoration: underline;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   color:#000000;
   */
}

div.top_control_bar a.interrupt:hover {
   background-color:#bb0000;
}

div.top_control_bar  a.interrupt_grey {
   color:#888888;
}

div.top_control_bar  a.interrupt_in_progress {
   color:#FFFFFF;
   background-color:#bb0000;
   text-decoration:blink;
}

div.top_control_bar a.help {
}

div.top_control_bar a.help:hover{
   cursor:pointer;
}

a.slide_mode {
}

a.slide_mode:hover{
   cursor:pointer;
}

a.cell_mode {
}

a.cell_mode:hover{
   cursor:pointer;
}

a.slide_arrow {
}

a.slide_arrow:hover{
   cursor:pointer;
}


div.top_control_bar  a.history_link {
}

div.top_control_bar a.history_link:hover {
}


/* links darker! no underlines! */

span.worksheet_control_commands a {
   color: #0000BB;
   text-decoration: none;
   padding:5px;
}

span.worksheet_control_commands a:hover {
   cursor:pointer;
}

span.worksheet_control_commands  a.plain_text {
}

span.worksheet_control_commands a.plain_text:hover {
}

span.worksheet_control_commands  a.doctest_text {
}

span.worksheet_control_commands a.doctest_text:hover {
}

span.worksheet_control_commands  a.download_sws {
}

span.worksheet_control_commands a.download_sws:hover {
}

span.worksheet_control_commands a.evaluate {
}

span.worksheet_control_commands a.evaluate:hover {
}

span.worksheet_control_commands  a.hide{
}

span.worksheet_control_commands a.hide:hover {
}





div.slideshow_control {
  float:right;
}

div.slideshow_control:hover {
   cursor:pointer;
}


div.slideshow_progress {
  float:right;
  background-color:white;
  padding:1px;
  border:1px solid <color2>;
  width:20%;
}
div.slideshow_progress_bar{
  z-index:1;
  position:relative;
  background-color: <color1>;
}
div.slideshow_progress_bar:hover {
   cursor:pointer;
}
div.slideshow_progress_text{
  position:absolute;
  z-index:2;
  top:2px;
  text-align:center;
  color:black;
  width:20%;
}

/************** Introspection ****************/

.completion_menu_selected {
  background-color: #8888ff;
}

div.introspection {
}

pre.introspection {
  font-family: monospace;
  font-size:15px;
  background-color: #efefef;
  border: solid 1px black;
  padding:8px;
  margin:8px;
}


ul.completion_menu_one {
  list-style: none;
  position: absolute;
  z-index: 2;
  background-color: #efefef;
  border: solid 1px black;
  display:inline;
  margin: 5px;
  font-family: monospace;
  font-size:15px;
  padding:5px;
}
li.completion_menu_one {
  display:inline;
  position: relative;
  float:left;
  margin: 0px;
}

ul.completion_menu_two {
  display:inline;
  position: relative;
  list-style: none;
  margin: 0px;
}

li.completion_menu_two{
  display:block;
  position:relative;
  margin: 3px;
  padding-left: 3px;
  padding-right: 3px;
}

li.completion_menu_two:hover{
   background-color: #8888bb;
   cursor:pointer;
}


/***** SEARCH / HELP AREA *********************************/

span.search_doc_topbar {
   z-index: 12;
   height: 24px;
   font-family:monospace;
   font-size: 12px;
   width:158px;
   top: 40px;
   left: 5px;
   position: fixed;
   border:1px solid #387CAF;
   background-color: #73a6ff;
}

td.menubar{
   text-decoration: none;
   font-family:arial;
   font-size:15px;
   font-weight:bold;
   color:#FFFFFF;
}

a.menubar{
   text-decoration: none;
   font-family:arial;
   font-size:15px;
   font-weight:bold;
   color:#FFFFFF;
   background-color:#73a6ff;
}

input.search_input {
   position: fixed;
   left: 5px;
   top: 65px;
   height: 32px;
   width: 160px;
   padding: 4px;
   z-index: 12;
   font-family:monospace;
   font-size:14px;
   color: #222222;
   color: #808080;
   border: 3px solid #387CAF;
   background: #FFF;
}

span.search_doc {
   z-index: 12;
   font-family:arial;
   font-size:12px;
   overflow:auto;
   position: fixed;
   top: 96px;
   left: 5px;
   width: 154px;
   height: 150px;
   margin: 0px;
   border:1px solid #387CAF;
   background-color: white;
   padding: 2px;
}


/************ INFO PANES **************************/

span.pane {
   z-index:30;
   font-family: monospace;
   font-size:12px;
   position: fixed;
   left: 5px;
   top: 33px;
   width: 180px;
   height:100%;
   margin: 0px;
   padding-right: 2px;
   padding-left: 0px;
   padding-top: 0px;
   bottom: 0ex;
}

span.plusminus {
  color:black;
  font-size:8pt;
  font-family:monospace;
}
span.plusminus:hover {
    cursor:pointer;
}

span.controltoggle {
  color:#0000ff;
  font-size:10pt;
  text-decoration:underline;
}
span.controltoggle:hover {
    cursor:pointer;
}

div.left_pane_bar {
  position:fixed;
  left: 0px;
  top:36px;
  background-color: white;
  width:8px;
  height:100%;
  z-index:100;
}
div.left_pane_bar:hover {
  background-color: #8888fe;  /* #000000; */
}


/************ VARIABLES **************************/

span.pane div.variables_topbar {
   color:black;
   background: url('corner.png') no-repeat top left;
   background-color: <color1>;
   font-family:arial;
   text-decoration: none;
   font-size:13px;
   height: 2ex;
   padding-left: 10px;
   padding-bottom:1px;
   width: 174px;
}

span.pane div.variable_list {
   font-size:11px;
   top:0ex;
   height:20ex;
   border:2px solid <color1>;
   width: 180px;
   overflow:auto;
}

div.variable_name {
   padding-left:1ex;
   border-top:1px solid #d3e9ff;
}

/*div.variable_name:hover {
   background-color:<color1>;
   cursor:pointer;
}*/

span.varname {
}

span.vartype {
  /* color:#888888; */
  color:#657d6c;
}

/************ ATTACHED **************************/

span.pane div.attached_topbar {
   color: black;
   height: 2ex;
   top: 0ex;
   background: url('corner.png') no-repeat top left;
   background-color: <color1>;
   text-decoration: none;
   font-size:13px;
   font-family:arial;
   padding-left: 10px;
   padding-bottom:1px;
   width: 174px;
}

span.pane div.attached_list {
   font-size:11px;
   top:0ex;
   height:20ex;
   border:2px solid  <color1>;
   width: 180px;
   overflow:auto;
}

div.attached_filename {
   padding-left:1ex;
   border-top:1px solid #d3e9ff;
}
/*
div.attached_filename:hover {
   background-color:<color1>;
   cursor:pointer;
}
*/

/************ WORKSHEETS **************************/

span.pane div.worksheets_topbar {
   color:black;
   height: 3ex;
   top: 0ex;
   background: url('/images/corner.png') no-repeat top left;
   background-color: <color2>;
   text-decoration: none;
   font-size:15px;
   font-family:arial;
   padding-left: 10px;
   padding-top:10px;
   width: 174px;
}

a.left_panel_hide {
   position: relative;
   top:0px;
   right:-1px;
   text-align:right;
   color:blue;
   font-weight:normal;
   font-family:arial;
   font-size:12px;
}

a.left_panel_hide:hover {
   cursor:pointer;
}

span.X {
   color:white;
   font-family:arial monospace;
   font-weight:bold;
   cursor:pointer;
}

span.pane div.add_new_worksheet_menu {
   position:relative;
   color:black;
   padding-top: 0.5ex;
   padding-bottom: 0.5ex;
   left: 0ex;
   background-color: white;
   text-decoration: none;
   font-size:11px;
   font-family:arial;
   padding-left: 0px;
   width: 174px;
}

input.add_new_worksheet_menu {
   width:100%
}

button.add_new_worksheet_menu {
   font-size:14px;
   font-family:sans-serif;
}

span.pane div.upload_worksheet_menu {
   color:black;
   top: 0ex;
   background-color: <color2>;
   text-decoration: none;
   font-size:11px;
   font-family:arial;
   padding-left: 10px;
   width: 174px;
   display:none;
}

button.upload_worksheet_menu {
   font-size:11px;
   font-family:arial;
}

input.upload_worksheet_menu {

}

span.pane div.delete_worksheet_menu {
   color:black;
   top: 0ex;
   background-color: <color2>;
   text-decoration: none;
   font-size:11px;
   font-family:arial;
   padding-left: 10px;
   width: 174px;
   display:none;
}

input.delete_worksheet_menu {
   width:50%
}

button.delete_worksheet_menu {
   font-size:11px;
   font-family:arial;
   background-color: #ffcccc;
}

span.pane div.worksheet_list {
   position:fixed;
   overflow:scroll;
   font-size:12px;
   top:25ex;
   bottom:2ex;
   left:1ex;
   border:2px solid <color2>;
   width: 180px;
}

a.new_worksheet {
   font-family: arial, monospace;
   font-size:12pt;
   text-align:right;
   color: #0000aa
}

a.new_worksheet:hover {
  cursor:pointer;
}

div.worksheet_bottom_padding {
   height:50%;
}

div.worksheet_top_padding {
   height:5%;
}

div.worksheet_menu {
   top:50px;
}

a.worksheet_title {
   text-decoration:none;
   font-size:20px;
   font-family:arial;
   font-weight:bold;
   color:#000000;
}

a.worksheet_title:hover {
   background-color:#ffffcc;
   cursor:pointer;
}

div.worksheet_title {
   padding-left: 1em;
   background-color: #ffffff;
   color:black;
}


div.worksheet_print_title {
   text-decoration:none;
   font-size:24px;
   font-family:arial;
   font-weight:bold;
   color:#000000;
   text-align:center;
}

div.worksheet_title_under {
   z-index:0;
   padding-top: 3px;
   padding-left: 1em;
   background-color: <color1>;
   /*width: 100%;*/
   font-family:arial;
   font-size: 22px;
   font-weight:bold;
   color:black;
}

div.worksheet_cell_list {
   padding-left:0.5ex;
}

a.delete_worksheet {
   font-family: arial, monospace;
   font-size:8pt;
   text-decoration:underline;
   text-align:right;
   color: #0000aa
}

a.delete_worksheet:hover {
  cursor:pointer;
}


a.upload_worksheet:hover {
  cursor:pointer;
}

span.pane a.worksheet_current {
   font-size:14px;
   padding-left:1ex;
   background-color:<color2>;
   text-decoration:none;
   color:black;
 }

span.pane a.worksheet_current_computing {
   font-size:14px;
   padding-left:1ex;
   background-color:#ffd1d1;
   text-decoration:none;
   color:black;
 }

span.pane a.worksheet_other {
   font-size:14px;
   padding-left:1ex;
   background-color:white;
   text-decoration:none;
   color:black;
}

span.pane a.worksheet_other:hover {
   background-color:<color2>;
   text-decoration:none;
   cursor:pointer;
}

span.pane a.worksheet_other_computing {
   font-size:14px;
   padding-left:1ex;
   background-color:#ffd1d1;
   text-decoration:none;
   color:black;
}

/*********** DOC-BROWSER************************/


.verbatim {
    background-color: #fafad2;
    border-style: solid;
    border-width: 1px 1px;
    border-color: black;
}


/************ OBJECTS **************************/

span.pane div.objects_topbar {
   color:black;
   height: 2ex;
   top: 0ex;
   background: url('corner.png') no-repeat top left;
   background-color: <color2>;
   text-decoration: none;
   font-size:13px;
   font-family:arial;
   padding-left: 10px;
   padding-bottom:1px;
   width: 174px;
}

span.pane div.object_list {
   font-size:11px;
   height:20ex;
   border:2px solid <color2>;
   width: 180px;
   overflow:auto;
}

a.object_name {
   padding-left:1ex;
   border-top:1px solid <color2>;
   background-color:white;
   text-decoration:none;
   color:black;
}

a.object_name:hover {
   background-color:<color2>;
   text-decoration:none;
   color:black;
   cursor:pointer;
}



/************ CONTROLS **************************/

div.control_area{
    vertical-align: top;
}

span.control {
    border:1px solid white;
    font-family: monospace;
    font-size:14pt;
    font-weight:bold;
}

span.control a.cs {
    color:#777777;
    text-decoration:none;
    border:0px solid white;
}

span.control:hover a.cs, span.control a:hover.cs {
    color:black;
    border:1px solid #333333;
}

/************ WORKSHEET **************************/

div.worksheet {
  overflow:auto;
  background-color: white;
  border:1px solid #aaa;
}


div.banner{
  background-color:white;
  font-family: sans-serif;
  font-size:18px;
  text-decoration: none;
  color: #1950c8;
}

div.banner a.banner{
    text-decoration:none;
    border:none;
    margin-top:2px;
}

a.banner:visited {
  color: #1950c8;
}

div.banner a.banner img{
    text-decoration:none;
    border:none;
    margin-top:2px;
}

div.banner a.banner:hover {
}

input.btn {
  font-family: monospace;
  font-size:13pt;
  font-weight:bold;
  color:#808080;
  text-decoration:none;
  background: white;
  padding:0px;
  margin:0px;
  border:1px solid white;
}
input.btn:hover {
  color: black;
  text-decoration: none;
  background: white;
  padding: 0px;
  margin: 0px;
  border: 1px solid #333333;
}

/************ CELL INPUT **************************/

div.cell_visible {
    display:block;
}

div.cell_evaluated {
    border-left: 1px solid white;
    padding-left:3px;
}

div.cell_not_evaluated {
    border-left: 1px solid #ff8888;
    padding-left:3px;
}

td.cell_number {
   font-size:12pt;
   font-family:arial, monospace;
   color:#bbbbbb;
   text-align:left;
}

td.cell_number:hover {
   color:#555555;
   cursor:pointer;
}

td.cell_number_running {
   font-size:12pt;
   font-family:arial, monospace;
   color:#bbbbbb;
   border-left:4px solid #aaffaa;
   /* background-color: #eeeeee; */
   text-align:left;
}

td.cell_number_running:hover {
  cursor:wait;
}

td.output_cell {
   /* width: 100% */
   height:3px;
}

div.cellbox {
  z-index:2;
  background-color: white;
  padding-left: .5em;
  padding-top: 4em;
}


textarea.cell_input {
  color:#000000;
  background-color: white;
  border-left: 1px solid  #a8a8a8;
  border-bottom: 1px solid  #a8a8a8;
  border-top: 1px solid  #a8a8a8;
  border-right: 1px solid  #a8a8a8;
  font-family: monospace;
  font-size:12pt;
  overflow:auto;
  padding-left:5px;
  padding-top:3px;
  padding-bottom:0px;
  width: 95%;
  margin-bottom:0px;
  margin-top:0px;
  line-height:1.2em;
}

pre.cell_input {
  color:#000000;
  background-color: white;
  border-left: 1px solid  #a8a8a8;
  border-bottom: 1px solid  #a8a8a8;
  border-top: 1px solid  #a8a8a8;
  border-right: 1px solid  #a8a8a8;
  font-family: monospace;
  font-size:12pt;
  padding-left:5px;
  padding-top:3px;
  padding-bottom:0px;
  width: 95%;
  margin-bottom:0px;
  margin-top:0px;

}
pre.cell_input:hover {
  cursor:text;
}

textarea.cell_input_hide {
  background-color: white;
  border: 0px solid white;
  font-family: monospace;
  font-size:12pt;
  overflow:hidden;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  /* width: 100%; */
  height:0.5em;
  margin:0px;
}


pre.cell_input_hide {
  background-color: white;
  border: 2px solid #e8e8e8;
  font-family: monospace;
  font-size:12pt;
  overflow:hidden;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  /* width: 100%; */
  height:1em;
  margin:0px;
}

pre.cell_input_hide:hover {
  cursor:text;
}

textarea.cell_input_active {
  background-color: white;
  border-left: 1px solid  #8888fe;
  border-bottom: 1px solid  #8888fe;
  border-top: 1px solid  #8888fe;
  border-right: 1px solid  #8888fe;
  font-family: monospace;
  font-size:12pt;
  overflow:auto;
  padding-left:5px;
  padding-top:3px;
  padding-bottom:0px;
  margin-top:0px;
  margin-bottom:0px;
  line-height:1.2em;
  width: 95%;
}

textarea.cell_input:hover{
  cursor:text;
}

a.eval_button {
  display:none;
}
a.eval_button_active {
  display: block;
  position: relative;
  top: 2px;
  margin:0px;
  padding:0px;
  font-size:10pt;
}


/************ CELL OUTPUT **************************/

div.cell_div_output {
  font-family: monospace;
  font-size:12pt;
  /*width: 95%;*/
  margin-top:-5px;
  margin-bottom:5px;
  padding-bottom:5px;

 /* border-left: 1px solid #aaaaff;  */
}

table.cell_output_box {
  margin:0px;
  padding:0px;
}

/*table.cell_output_box:hover {
  background-color: #fafafa;
}
*/

div.cell_div_output_wrap {
  font-size:12pt;
  margin:0px;
  padding-left:0px;
  color:#0000aa;
}

div.cell_output_wrap pre.cell_output_print_wrap {
  font-size:12pt;
  margin:0px;
  padding:0px;
  color:#0000aa;
}

div.cell_output_print_wrap {
   font-size:10pt;
}

div.cell_output_nowrap {
  display:none;
}
div.cell_output_print_nowrap {
  display:none;
}

div.cell_output_hidden {
  display:none;
}


div.cell_output_nowrap_wrap {
 display:none;
}
div.cell_output_print_nowrap_wrap {
 display:none;
}
div.cell_output_nowrap_nowrap {
  font-size:12pt;
  margin:0px;
  padding:0px;
  color:#0000aa;
}
div.cell_output_nowrap_hidden {
  display:none;
}

div.cell_output_html_wrap {
  font-family: monospace;
  font-size:12pt;
}
div.cell_output_html_nowrap {
  font-family: monospace;
  font-size:12pt;
}
div.cell_output_html_hidden {
   display:none;
}

div.cell_div_output_running {
  font-family: monospace;
  font-size:12pt;
  /* width: 100%; */
  margin:0px;
  background-color:#ffffff;
  padding:0px;
}

div.cell_div_output_running:hover {
  cursor:wait;
}


div.cell_div_output_hidden {
  width: 100%;
  height: 3px;
  margin:0px;
  border-left: 4em solid #aaaaaa;
/*   border-top: 1px solid <color1>;
  border-bottom: 1px solid <color1>;
  */
}

pre.shrunk {
  font-size:12pt;
  margin:0px;
}

pre.cell_output_hidden {
  display: none;
}

pre.cell_output_hide {
  display:none;
}

a.file_link {
  text-decoration:underline;
}


/************ INSERTING NEW CELLS **************************/

div.insert_new_cell {
  height:6px;
  /* width:100%; /*
  /* border-top: 4px solid white; */
  display:block;
  margin:3px;
}

div.insert_new_cell:hover {
  /* border-top: 4px solid #000000; */
  background-color: #8888fe;  /* #000000; */
  margin:3px;
  /* background-color:#eeeeee; */
}


/************ DEBUG WINDOW **************************/

div.debug_window_active {
  background-color: white;
  border: 1px solid #fe8888;
  overflow:auto;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  /* width: 100%; */
}
div.debug_window_inactive {
  background-color: white;
  border: 1px solid #888888;
  overflow:auto;
  padding-left:0px;
  padding-top:0px;
  padding-bottom:0px;
  /* width: 100%; */
}


div.debug_output {
  background-color: white;
  border: 0px;
  font-family: monospace;
  font-size:10pt;
  overflow:scroll;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  height: 10em;
  /* width: 100%; */
}


textarea.debug_input {
  background-color: white;
  border: 1px solid #8888fe;
  font-family: monospace;
  font-size:12pt;
  overflow:scroll;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  /* width: 100%; */
}

span.red{
  color:red;
}

/****************** Other **************/

a.worksheetname{
   text-decoration: none;
}

a.worksheetname_moved{
   color:#888888;
   text-decoration: none;
   font-weight:normal;
}

span.worksheet_buttons {
    position:relative;
    top: -20ex;
    right: 0ex;
}

.flush-right {
    position:absolute;      /* All browsers */
    top: auto;              /* Standards  browsers */
    right: 0;               /* All except IE */
}

.thin-right {
    position:absolute;      /* All browsers */
    top: auto;              /* Standards  browsers */
    right: 0;               /* All except IE */
    width: 70%;
}

/************ User Home **************************/

span.username {
  font-family: sans-serif;
  font-weight:bold;
  font-size:14px;
  padding:1ex;
}

span.ratingmsg {
  color: #112abb;
  padding:0.6ex;
  font-size:14px;
}

span.pubmsg {
  font-family: sans-serif;
  color: #112abb;
  padding:0.6ex;
  font-size:12px;
}

a.usercontrol {
  color: #112abb;
  padding:0.6ex;
  font-size:14px;
  text-decoration:underline;
}

a.usercontrol:hover {
  cursor:pointer;
}

span.usercontrol {
  color: #112abb;
  padding:0.6ex;
  font-size:14px;
}


a.boldusercontrol {
  color: #112abb;
  padding:1ex;
  font-weight:bold;
  font-size:14px;
}

a.control, a.control-select {
  background-color:#7799bb;
  font-family: sans-serif;
  color: #ffffff;
  padding-top:0.5ex;
  padding-bottom:0.5ex;
  padding-left:1ex;
  padding-right:1ex;
  font-size:15px;
  font-weight:bold;
  text-decoration: none;
}

a.control:hover {
   cursor:pointer;
}

a.control-select {
 background-color:#4477aa;
}

a.control-select:hover {
   cursor:pointer;
}

span.sharebar {
  background-color:#4477aa;
  font-family: sans-serif;
  color: #ffffff;
  position:absolute;
  left:1ex;
  right:0ex;
  padding-top:1ex;
  padding-bottom:1ex;
  padding-left:4ex;
  font-size:18px;
  font-weight:bold;
}

textarea.edit {
    font-family: courier, monospace;
    font-size:10pt;
    border: 1px solid #8cacbb;
    color: black;
    background-color: white;
    padding: 3px;
    overflow:auto;
    /* width: 100%; */
    margin-top: 0.5em;
}

a.listcontrol {
  padding:1ex;
  color: #112abb;
  font-weight:bold;
  font-size:14px;
  text-decoration:none;
}

hr.usercontrol {
   border: 0;
   width: 99%;
   color: #c9d7f1;
   background-color: #c9d7f1;
   height: 1px;
}

hr.greybar hr.negative_greybar {
   border: 0;
   width: 99%;
   color: #aaa;
   background-color: #aaa;
   height: 1px;
}

hr.negative_greybar {
   top:-1em;
   position:relative;
}

span.checkcol {
  position:relative;
  left:0%;
  width:10%;
}

span.leftcol {
  position:relative;
  left:10%;
  width:20%;
}

span.middlecol {
  position:relative;
  left:30%;
  width:20%;
}

span.rightcol {
  position:relative;
  left:50%;
  width:20%;
}

tr.greybox {
   background-color: #e8eef7;
}

td.entry {
   padding:4px;
}

div.thinspace {
   border: 0;
   /*width: 100%;*/
   height: 2px;
}

tr.thingreybox {
   background-color: #aaa;
}

div.ultrathinspace {
   border: 0;
   /*width: 100%;*/
   height: 0px;
}

span.lastedit {
  font-family: sans-serif;
  font-size:10px;
  color: #717171;
}

span.revs {
  font-family: sans-serif;
  font-size:12px;
  font-weight:bold;
  color: #333333;
}

span.users {
  font-family: sans-serif;
  font-size:13px;
  color: #222222;
}

a.share {
  font-family: sans-serif;
  font-size:10px;
  color: #7777cc;
}

select.worksheet {
  width:6em;
  border: #aaaaaa;
  border-style: solid;
  border-top-width: 1px;
  border-right-width: 1px;
  border-bottom-width: 1px;
  border-left-width: 1px
}

select.worksheet_list {
  width:5em;
    border: #aaaaaa;
  border-style: solid;
  border-top-width: 1px;
  border-right-width: 1px;
  border-bottom-width: 1px;
  border-left-width: 1px
}

select.worksheet_edit {
  width:5em;
    border: #aaaaaa;
  border-style: solid;
  border-top-width: 1px;
  border-right-width: 1px;
  border-bottom-width: 1px;
  border-left-width: 1px
}

td.worksheet_link {
  font-family: sans-serif;
  font-size:12px;
  font-weight:bold;
  color: #000000;
}

td.archived_worksheet_link {
  font-family: sans-serif;
  font-size:12px;
  color: #000000;
}

td.owner_collab {
  font-family: sans-serif;
  font-size:12px;
  color: #000000;
}

td.last_edited {
  font-family: sans-serif;
  font-size:12px;
  color: #000000;
}

span.addtext {
  font-family: sans-serif;
  font-size:13px;
  color: #222;
}

textarea.plaintextedit {
    font-family: courier, monospace;
    font-size:10pt;
    border: 1px solid #8cacbb;
    color: black;
    background-color: white;
    overflow: auto;
    width: 99%;
    height: 60%;
}

pre.plaintext {
    overflow:auto;
    font-family: courier, monospace;
    font-size:10pt;
    border: 1px solid #8cacbb;
    color: black;
    background-color: white;
    margin-top: 0.5em;
}

div.docidx {
  text-align:center;
  font-family: sans-serif;
  font-size:16px;
  color: #222;
  font-weight:bold;
}

span.ping {
   display:none;
}

span.pingdown {
   font-family:sans-serif;
   font-size:15px;
   font-weight:bold;
   color:white;
   background-color: #990000;
}

"""
    if color == 'gmail':
        color1 = '#c3d9ff'
        color2 = '#b5edbc'
    elif color == 'grey':
        color1 = '#aaaaaa'
        color2 = '#888888'
    elif color == 'default' or color == None:
        color1 = '#dcdcdc'
        color2 = '#cccccc'
        #color1 = '#aaaaff'
        #color2 = '#b5edbc'
        #color2 = '#6cc755'
    elif isinstance(color, (tuple,list)):
        color1, color2 = color
    else:
        raise ValueError, "unknown color scheme %s"%color

    s = s.replace('<color1>',color1).replace('<color2>',color2)
    user_css = DOT_SAGE + '/notebook.css'
    if os.path.exists(user_css):
        s += '\n' + open(user_css).read()

    return s
