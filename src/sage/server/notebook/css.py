def css():
    return r"""


/**** TOP CONTROL BAR ************************/

span.top_control_bar {
   z-index: 5;
   background-color: white;
   position: fixed;
   left: 12.6em;
   width: 100%;
   top: 0px;
   height: 5ex;
  /* border: 2px solid #c3d9ff; */
   border-bottom: 3px solid #c3d9ff;
}


span.top_control_bar  input.search_fulltext {
   position: fixed;
   left: 14em;
   top: 5px;
   width: 15em;
   background-color:white;
   font-family:arial;
   font-size:14px;
   color:#808080;
}

div.workbook_title {
   z-index: 6;
   position:fixed;
   left: 15em;
   top: 5px;
   font-family:arial;
   font-size:20px;
   color:black;
   padding:5;
   background:white;
   display:inline;
}

span.top_control_bar  a.evaluate{
   text-decoration: underline;
   position: fixed;
   right: 10em;
   top: 1px;
   background-color:white;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   padding:5;
}

span.top_control_bar a.evaluate:hover {
   background-color:#00bb00;
   color:#FFFFFF;
}

span.top_control_bar  a.interrupt {
   text-decoration: underline;
   position: fixed;
   right: 17em;
   top: 1px;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   color:#000000;
   padding:5;
   background-color:white;
}

span.top_control_bar a.interrupt:hover {
   background-color:#bb0000;
   color:#FFFFFF;
}

span.top_control_bar  a.interrupt_grey {
   text-decoration: underline;
   position: fixed;
   right: 17em;
   top: 1px;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   color:#888888;
   padding:5;
   background-color:white;
}

span.top_control_bar  a.interrupt_in_progress {
   text-decoration: underline;
   position: fixed;
   right: 17em;
   top: 1px;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   color:#FFFFFF;
   padding:5;
   background-color:#bb0000;
   text-decoration:blink;
}

span.top_control_bar  a.hide{
   text-decoration: underline;
   position: fixed;
   right: 3em;
   top: 1px;
   background-color:white;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   padding:5;
}

span.top_control_bar a.hide:hover {
   background-color:#0000bb;
   color:#FFFFFF;
}

span.top_control_bar  a.help {
   text-decoration: underline;
   position: fixed;
   right: 23em;
   top: 1px;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   color:#000000;
   padding:5;
   background-color:white;
}

span.top_control_bar a.help:hover {
   background-color:#00bb00;
   color:#FFFFFF;
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
/*   overflow: auto; */
   overflow:none;
/*   position: fixed; */
   position: absolute;
   left: 1em;
   top: -3ex;
   width: 12em;
   height: 100%;
   margin: 1px;
   /* border:1px solid #387CAF;  */
   padding: 2px;
}

/* span.pane:hover {
   width:30em;
} */


/************ VARIABLES **************************/

span.pane div.variables_topbar {
   color:#FFFFFF;
   height: 2ex;
   top: 0ex;
   background-color: #73a6ff;
   font-family:arial;
   text-decoration: none;
   font-size:15px;
   font-weight:bold;
   padding: 2px;
}

span.pane div.variables_list {
   top:0ex;
   height:20ex;
   border:1px solid #387CAF;
   padding:1ex;
   overflow:auto;
}

div.variable_name {
}

div.variable_name:hover {
   background-color:#dfdfdf;
   font-size:18px;
}

span.varname {
}

span.vartype {
   color:blue;
}

/************ ATTACHED **************************/

span.pane div.attached_topbar {
   color:#FFFFFF;
   height: 2ex;
   top: 0ex;
   background-color: #73a6bb;
   text-decoration: none;
   font-size:15px;
   font-weight:bold;
   font-family:arial;
   padding: 2px;
}

span.pane div.attached_list {
   top:0ex;
   height:20ex;
   border:1px solid #387CAF;
   padding:1ex;
   overflow:auto;
}

/************ WORKBOOKS **************************/

span.pane div.workbooks_topbar {
   color:#FFFFFF;
   height: 2ex;
   background-color: #73aaa6;
   text-decoration: none;
   font-size:15px;
   font-weight:bold;
   font-family:arial;
   padding: 2px;
}

span.pane div.workbooks_list {
   height:20ex;
   border:1px solid #387CAF;
   padding:1ex;
}

div.workbook_name {

}

div.workbooks_name:hover {
   background-color:#cfcfcf;
}

/************ OBJECTS **************************/

span.pane div.objects_topbar {
   color:#FFFFFF;
   height: 2ex;
   top: 0ex;
   background-color: #aa73a6;
   text-decoration: none;
   font-size:15px;
   font-weight:bold;
   font-family:arial;
   padding: 2px;
}

span.pane div.objects_list {
   height:20ex;
   border:1px solid #387CAF;
   padding:1ex;
}

div.object_name {

}

div.object_name:hover {
   background-color:#cfcfcf;
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

/************ WORKBOOK **************************/

div.workbook {
  position:fixed;
  overflow:auto;
  z-index:1;
  background-color: white;
  border: 1px solid #387CAF;
  top: 6ex;
  bottom: 0ex;
  right: 0ex;
  left: 12em;
  padding: 2ex;
  float: right;
  /* width:65em; */
/*  width: 80%;
  height:85%;*/
}

div.workbook0 {
  z-index:1;
  position:absolute;
  text-align:center;
  /* background-color: #c3d9ff; */
  background-color: white;
  /* border-left: 1px solid #387CAF;  */
  top: 4ex;
  left: 12em;
  padding: 2ex;
  float: right;
  width:100%;
}

span.banner{
  z-index:5;
  background-color:white;
  font-family:arial;
  font-size:30px;
  text-decoration: none;
  font-weight: bold;
  color: #387CAF;
  margin: 0px;
  position: fixed;
  top: -5px;
  left: 6em;
  padding: 5px;
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
  /* background-color: #efefef; */
  background-color: white;

  border-right: 1px #ffffff;
  border-top: 1px solid #ffffff;
  border-bottom: 1px solid #ffffff;
  border-left: 1px solid #ffffff;

  font-family: courier, monospace;
  font-size:12pt;

  overflow:hidden;

  /* margin-left: 1em;*/
  padding-left: 1em;
  padding: 0px;

  width: 100%;
}

textarea.cell_input:hover {
  border-right: 1px solid #c3d9ff;
  border-top: 1px solid #c3d9ff;
  border-bottom: 1px solid #c3d9ff;
  border-left: 1px solid #c3d9ff;
}

div.cell_output {
  width: 100%;
  margin: 0px;
  padding: 2px;
  border-left: 2px solid #000088;
  background: #efefff
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
  /* background-color: #fff1f1; */
  background-color: #ffd1d1;
}

div.cell_output_hidden {
  width: 100%;
  height: 5px;
  margin: 0px;
  padding: 5px;
  /* border-left: 2em solid #000088; */
  border-left: 4em solid #aaaaaa;
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
    background-color: #6d79c5;
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
