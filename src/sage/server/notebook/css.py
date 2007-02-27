"""
SAGE Notebook CSS
"""

import os

from sage.misc.misc import DOT_SAGE

def css(color='default'):
    r"""
    Return the CSS header used by the SAGE Notebook.

    INPUT:
        color -- string or pair of html colors, e.g.,
                    'gmail'
                    'grey'
                    \code{('#ff0000', '#0000ff')}

    EXAMPLES:
        sage: import sage.server.notebook.css as c
        sage: type(c.css())
        <type 'str'>
    """
    s = r"""
div.hidden {
  display:none;
}
span.hidden{
  display:none;
}

div.fivepix {
  height:5px;
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
   position: fixed;
   top:36px;
   right:0px;
   text-align:right;
   color:blue;
   font-weight:normal;
   font-family:arial;
   font-size:12px;
   padding:5;
}

div.slide_control_commands {
   float:right;
   position: fixed;
   width:500px;
   top:1ex;
   right:1ex;
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
   padding:5;
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
}

div.top_control_bar a.slide_mode {
}

div.top_control_bar a.slide_mode:hover{
}

div.top_control_bar a.cell_mode {
}

div.top_control_bar a.cell_mode:hover{
}

div.top_control_bar a.slide_arrow {
}

div.top_control_bar a.slide_arrow:hover{
}


div.top_control_bar  a.history_link {
}

div.top_control_bar a.history_link:hover {
}


/* links darker! no underlines! */

span.worksheet_control_commands a {
   color: #0000BB;
   text-decoration: none;
   padding:5;
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
  font-family:courier, monospace;
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
  font-family:courier, monospace;
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
   font-family:courier;
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
   font-family:courier;
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
   font-family:courier, monospace;
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
  font-family:courier, monospace;
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
   height: 2ex;
   top: 0ex;
   background: url('corner.png') no-repeat top left;
   background-color: <color2>;
   text-decoration: none;
   font-size:12px;
   font-family:arial;
   padding-left: 10px;
   width: 174px;
}

span.X {
   color:white;
   font-family:arial monospace;
   font-weight:bold;
   cursor:pointer;
}

span.pane div.add_new_worksheet_menu {
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

input.add_new_worksheet_menu {
   width:60%
}

button.add_new_worksheet_menu {
   font-size:11px;
   font-family:arial;
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
   font-size:12px;
   top:0ex;
   height:25ex;
   border:2px solid <color2>;
   overflow:auto;
   width: 180px;
}

a.new_worksheet {
   font-family: arial, monospace;
   font-size:8pt;
   text-decoration:underline;
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

div.worksheet_title {
   z-index:2;
   top:36px;
   height:31px;
   padding-top: 3px;
   padding-left: 1em;
   background: url('corner.png') no-repeat top left;
   background-color: <color1>;
   width: 100%;
   font-family:arial;
   font-size: 16px;
   font-weight:bold;
   color:black;
   position: fixed;
}

div.worksheet_title_under {
   z-index:0;
   padding-top: 3px;
   padding-left: 1em;
   background-color: <color1>;
   width: 100%;
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
   font-size:12px;
   padding-left:1ex;
   border-top:1px solid <color2>;
   background-color:<color2>;
   text-decoration:none;
   color:black;
 }

span.pane a.worksheet_current_computing {
   font-size:12px;
   padding-left:1ex;
   border-top:1px solid <color1>;
   background-color:#ffd1d1;
   text-decoration:none;
   color:black;
 }

span.pane a.worksheet_other {
   font-size:12px;
   padding-left:1ex;
   border-top:1px solid <color2>;
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
   font-size:12px;
   padding-left:1ex;
   border-top:1px solid <color1>;
   background-color:ffd1d1;
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
    font-family: courier, monospace;
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
  position:fixed;
  overflow:auto;
  z-index:1;
  background-color: white;
  border-top: 0px;
  border-left: 2px solid <color1>;
  top: 70px;
  bottom: 0ex;
  right: 0ex;
  left: 198px;
  padding-left: 0ex;
  float: right;
  padding-top: 0ex;
}

div.slideshow {
  position:fixed;
  overflow:auto;
  z-index:1;
  background-color: white;
  border-top: 0px;
  border-left: 2px solid <color1>;
  top: 70px;
  bottom: 0ex;
  right: 0ex;
  left: 5px;
  padding-left: 0ex;
  float: right;
  padding-top: 0ex;
}


span.banner{
  background-color:white;
  /* font-family:arial;
  font-size:30px;
  text-decoration: none;
  font-weight: bold;
  color: #387CAF; */
}

span.banner a.banner img{
    text-decoration:none;
    border:none;
    margin-top:2px;
}

span.banner a.banner:hover {
}

input.btn {
  font-family: courier;
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
    border-left: 3px solid white;
    padding-left:3px;
}

div.cell_not_evaluated {
    border-left: 2px dotted #ff8888;
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
   background-color: 99ff99;
   text-align:left;
}

td.cell_number_running:hover {
  cursor:wait;
}

td.output_cell {
   width:100%;
   height:3px;
}

div.cellbox {
  z-index:2;
  background-color: white;
  padding-left: .5em;
  padding-top: 4em;
}

pre.cell_input_pre {
  background-color: white;
  border: 0px solid white;
  font-family: courier, monospace;
  font-size:12pt;
  overflow:hidden;
  padding-left:0px;
  padding-top:0px;
  padding-bottom:0px;
  margin:0px;
  display:inline;
  width: 100%;
}

textarea.cell_input {
  color:#000000;
  background-color: #e8e8e8;
  border: 2px solid white;
  font-family: courier, monospace;
  font-size:12pt;
  overflow:auto;
  padding-left:5px;
  padding-top:3px;
  padding-bottom:0px;
  width: 100%;
  margin-bottom:0px;
  margin-top:0px;
}


textarea.cell_input_hide {
  background-color: white;
  border: 0px solid white;
  font-family: courier, monospace;
  font-size:12pt;
  overflow:hidden;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  width: 100%;
  height:0.5em;
  margin:0px;
}


textarea.cell_input_active {
  background-color: white;
  border: 2px solid  #8888fe;
  font-family: courier, monospace;
  font-size:12pt;
  overflow:auto;
  padding-left:5px;
  padding-top:3px;
  padding-bottom:0px;
  margin-top:0px;
  margin-bottom:0px;
  width: 100%;
}


span.cell_evaluate {
  position: relative;
  top: 2px;
  cursor:pointer;
}


/************ CELL OUTPUT **************************/

div.cell_output {
  font-family: courier, monospace;
  font-size:12pt;
  width: 95%;
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

div.cell_output_wrap {
  font-size:12pt;
  margin:0px;
  padding-left:0px;
  color:#0000aa;
}

span.cell_output_wrap {
  font-size:12pt;
  margin:0px;
  padding:0px;
  color:#0000aa;
}
span.cell_output_nowrap {
  display:none;
}
span.cell_output_hidden {
  display:none;
}


span.cell_output_nowrap_wrap {
 display:none;
}
span.cell_output_nowrap_nowrap {
  font-size:12pt;
  margin:0px;
  padding:0px;
  color:#0000aa;
}
span.cell_output_nowrap_hidden {
  display:none;
}

span.cell_output_html_wrap {
  font-family: courier, monospace;
  font-size:12pt;
}
span.cell_output_html_nowrap {
  font-family: courier, monospace;
  font-size:12pt;
}
span.cell_output_html_hidden {
   display:none;
}

div.cell_output_running {
  font-family: courier, monospace;
  font-size:12pt;
  width: 100%;
  margin:0px;
  background-color:#ffffff;
  padding:0px;
}

div.cell_output_running:hover {
  cursor:wait;
}


div.cell_output_hidden {
  width: 100%;
  height: 3px;
  margin:0px;
  border-left: 4em solid #aaaaaa;
/*   border-top: 1px solid <color1>;
  border-bottom: 1px solid <color1>;
  */
}

pre.shrunk {
/*   height:0px; */
   display:inline;
}

pre.cell_output_hidden {
  display: none;
}

pre.cell_output_hide {
  display:none;
}

a.file_link {
  text_decoration:underline;
}


/************ INSERTING NEW CELLS **************************/

div.insert_new_cell {
  height:6px;
  width:100%;
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
  width: 100%;
}
div.debug_window_inactive {
  background-color: white;
  border: 1px solid #888888;
  overflow:auto;
  padding-left:0px;
  padding-top:0px;
  padding-bottom:0px;
  width: 100%;
}


div.debug_output {
  background-color: white;
  border: 0px;
  font-family: courier, monospace;
  font-size:10pt;
  overflow:scroll;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  height: 10em;
  width: 100%;
}


textarea.debug_input {
  background-color: white;
  border: 1px solid #8888fe;
  font-family: courier, monospace;
  font-size:12pt;
  overflow:scroll;
  padding-left:3px;
  padding-top:0px;
  padding-bottom:0px;
  width: 100%;
}

span.red{
  color:red;
}


/***********************************************************/
/*                     wiki css styling                    */
/***********************************************************/



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
