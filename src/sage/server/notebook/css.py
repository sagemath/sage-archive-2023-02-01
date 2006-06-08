def css():
    return r"""


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
   text-decoration:underline;
}

span.vbar {
   height:1.5ex;
   border-left:1px solid black;
}

div.top_control_bar a.evaluate {
   background-color:white;
   padding:5;
}

div.top_control_bar a.evaluate:hover {
   background-color:#00bb00;
   color:#FFFFFF;
   cursor:pointer;
}

div.top_control_bar  a.interrupt {
   text-decoration: underline;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   color:#000000;
   padding:5;
   background-color:white;
}

div.top_control_bar a.interrupt:hover {
   background-color:#bb0000;
   color:#FFFFFF;
   cursor:pointer;
}

div.top_control_bar  a.interrupt_grey {
   color:#888888;
   padding:5;
   background-color:white;
}

div.top_control_bar  a.interrupt_in_progress {
   color:#FFFFFF;
   padding:5;
   background-color:#bb0000;
   text-decoration:blink;
}

div.top_control_bar  a.hide{
   background-color:white;
   padding:5;
}

div.top_control_bar a.hide:hover {
   background-color:#0000bb;
   color:#FFFFFF;
   cursor:pointer;
}

div.top_control_bar  a.help {
   padding:5;
   background-color:white;
}

div.top_control_bar a.help:hover {
   background-color:#00bb00;
   color:#FFFFFF;
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
   left: 1em;
   top: 33px;
   width: 180px;
   height:100%;
   margin: 0px;
   padding-right: 2px;
   padding-left: 2px;
   padding-top: 0px;
   bottom: 0ex;
}



/************ VARIABLES **************************/

span.pane div.variables_topbar {
   color:black;
   background-color: #c3d9ff;
   font-family:arial;
   text-decoration: none;
   font-size:13px;
   height: 2ex;
   padding-left: 10px;
   margin:0;
   width: 174px;
}

span.pane div.variables_list {
   font-size:11px;
   top:0ex;
   height:25ex;
   border:2px solid #c3d9ff;
   width: 180px;
   overflow:auto;
}

div.variable_name {
   padding-left:1ex;
   border-top:1px solid #d3e9ff;
}

div.variable_name:hover {
   background-color:#c3d9ff;
   cursor:pointer;
}

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
   background-color: #c3d9ff;
   text-decoration: none;
   font-size:13px;
   font-family:arial;
   padding-left: 10px;
   width: 174px;
}

span.pane div.attached_list {
   font-size:11px;
   top:0ex;
   height:25ex;
   border:2px solid  #c3d9ff;
   width: 180px;
   overflow:auto;
}

div.attached_filename {
   padding-left:1ex;
   border-top:1px solid #d3e9ff;
}

div.attached_filename:hover {
   background-color:#c3d9ff;
   cursor:pointer;
}

/************ WORKSHEETS **************************/

span.pane div.worksheets_topbar {
   color:black;
   height: 2ex;
   top: 0ex;
   background-color: #b5edbc;
   text-decoration: none;
   font-size:13px;
   font-family:arial;
   padding-left: 10px;
   width: 174px;
}

span.pane div.worksheet_list {
   font-size:11px;
   top:0ex;
   height:25ex;
   border:2px solid #b5edbc;
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

a.upload_worksheet {
   font-family: arial, monospace;
   font-size:8pt;
   text-decoration:underline;
   text-align:right;
   color: #0000aa
}

a.upload_worksheet:hover {
  cursor:pointer;
}

span.pane a.worksheet_current {
   font-size:11px;
   padding-left:1ex;
   border-top:1px solid #b5edbc;
   background-color:#b5edbc;
   text-decoration:none;
   color:black;
 }

span.pane a.worksheet_current_computing {
   font-size:11px;
   padding-left:1ex;
   border-top:1px solid #c3d9ff;
   background-color:#ffd1d1;
   text-decoration:none;
   color:black;
 }

span.pane a.worksheet_other {
   font-size:11px;
   padding-left:1ex;
   border-top:1px solid #b5edbc;
   background-color:white;
   text-decoration:none;
   color:black;
}

span.pane a.worksheet_other:hover {
   background-color:#b5edbc;
   text-decoration:none;
   cursor:pointer;
}

span.pane a.worksheet_other_computing {
   font-size:11px;
   padding-left:1ex;
   border-top:1px solid #c3d9ff;
   background-color:ffd1d1;
   text-decoration:none;
   color:black;
}

/************ OBJECTS **************************/

span.pane div.objects_topbar {
   color:black;
   height: 2ex;
   top: 0ex;
   background-color: #b5edbc;
   text-decoration: none;
   font-size:13px;
   font-family:arial;
   padding-left: 10px;
   width: 174px;
}

span.pane div.object_list {
   font-size:11px;
   height:25ex;
   border:2px solid #b5edbc;
   width: 180px;
   overflow:auto;
}

a.object_name {
   padding-left:1ex;
   border-top:1px solid #b5edbc;
   background-color:white;
   text-decoration:none;
   color:black;
}

a.object_name:hover {
   background-color:#b5edbc;
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
  border-top: 32px solid #c3d9ff;
  border-left: 10px solid #c3d9ff;
  top: 36px;
  bottom: 0ex;
  right: 0ex;
  left: 205px;
  padding-top: 0ex;
  padding-left: 1ex;
  float: right;
}

span.worksheet_title {
   padding-top: 3px;
   font-family:arial;
   font-size:22px;
   font-weight:bold;
   color:black;
   display:inline;
   position: fixed;
}


span.banner{
  background-color:white;
  font-family:arial;
  font-size:30px;
  text-decoration: none;
  font-weight: bold;
  color: #387CAF;
  margin: 0px;
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

/************ CELLS **************************/

textarea.cell_input {
  background-color: white;

  border: 2px solid #ffffff;

  font-family: courier, monospace;
  font-size:12pt;

  overflow:hidden;

  padding-left:0px;
  padding-top:0px;
  padding-bottom:0px;

  width: 100%;
}

textarea.cell_input_active {
  background-color: white;

  border: 2px solid #73a6ff;

  font-family: courier, monospace;
  font-size:12pt;

  overflow:hidden;

  padding-left:0px;
  padding-top:0px;
  padding-bottom:0px;

  width: 100%;
}

/*textarea.cell_input:hover {
  border: 2px solid #73a6ff;
}
*/

div.cell_output {
  width: 100%;
  margin: 0px;
  padding: 2px;
  border-left: 2px solid #000088;
/*  background: #eff2ff */
}

div.cell_output:hover {
  /* border-right: 1px solid #000088;
  border-top: 1px solid #000088;
  border-bottom: 1px solid #000088;
  border-left: 2px solid #000088; */
}

table.cell_output {
  margin-left:2em;
}

pre.cell_output {
  font-family: courier, monospace;
  font-size:12pt;
  color:#000088;
  display:inline;
}

pre.cell_output_nowrap {
  font-family: courier, monospace;
  font-size:10pt;
  color:#000088;
  display:none;
}

pre.cell_output_nowrap_visible {
  font-family: courier, monospace;
  font-size:10pt;
  color:#000066;
  display:inline;
}

div.cell_output_running {
  width: 100%;
  margin: 0px;
  padding: 5px;
  border-left: 5px solid #880000;
  background-color: #ffeeee;
  /*
  border-left: 5px solid #008800;
  background-color: #d1ffd1;
  */
}

div.cell_output_running:hover {
  cursor:wait;
}

div.cell_output_hidden {
  width: 100%;
  height: 10px;
  margin: 0px;
  /* padding: 5px;*/
  border-left: 4em solid #c3d9ff;
}

div.cell_output_hidden:hover {
  /* border-right: 1px solid #000088;
  border-top: 1px solid #000088;
  border-bottom: 1px solid #000088;
  border-left: 8px solid #000088;
  */
}

pre.cell_output_hidden {
  display: none;
}


/************ INSERTING NEW CELLS **************************/

div.insert_new_cell {
  height:8px;
  width:100%;
  border-top: 1px solid white;
  /* border: 2px solid white; */
  display:block;
  font-size:10pt;
  text-align:left;
  /* padding: 2px; */
}

div.insert_new_cell:hover {
  /* border: 2px solid #dddddd; */
  border-top: 1px solid black;
}


/****************** HELP WINDOW ***********************/
div.help_window {
    z-index:60;
    position:fixed;
    overflow:auto;
    background-color:white;
    border: 3px solid #3d86d0;
    top: 10ex;
    bottom:10%;
    left:25%;
    right:15%;
    padding:2ex;
    display:none;
/*    opacity:0.9;   enable this when computers are faster in a few years?...*/
}


div.help_window_title {
    z-index:60;
    position:fixed;
    overflow:auto;
    background-color: #3d86d0;
    font-size:15px;
    font-weight:bold;
    font-family:arial monospace;
    padding: 2px;
    color:#FFFFFF;
    border: 3px solid #3d86d0;
    top: 7ex;
    height: 2ex;
    left:25%;
    right:15%;
}

div.help_window_close {
    z-index:60;
    position:fixed;
    overflow:auto;
    background-color: #5d99d5;
    font-size:18px;
    font-family:courier monospace;
    padding: 1px;
    color:#FFFFFF;
    border: 3px solid #3d86d0;
    top: 6.3ex;
    height: 2ex;
    width:3ex;
    right:15%;
    text-align:right;
}

div.help_window_close:hover {
    z-index:60;
    background-color: #6d79c5;
    cursor:pointer;
}

table.help_window {
  /*  border: 1px solid #000000; */
    width:100%;
}

td.help_window_cmd {
    background-color: #f5e0aa;
    width:30%;
    padding:1ex;
}

td.help_window_how {
    padding:1ex;
    width:70%;
}

"""
