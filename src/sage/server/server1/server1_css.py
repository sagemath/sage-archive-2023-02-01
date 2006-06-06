def css():
    return r"""


/**** TOP CONTROL BAR ************************/

span.top_control_bar {
   z-index: 5;
   background-color: white;
   position: fixed;
   left: 170px;
   width: 100%;
   top: 0px;
   height:40px;
   border: 2px solid #c3d9ff;
}


span.top_control_bar  input.search_fulltext {
   position: fixed;
   left: 175px;
   top: 5px;
   width: 150px;
   background-color:white;
   font-family:arial;
   font-size:14px;
   color:#808080;
}

div.workbook_title {
   z-index: 6;
   position:fixed;
   left: 350px;
   top: 1px;
   font-family:arial;
   font-size:20px;
   color:black;
   padding:5;
   background:white;
}

span.top_control_bar  a.evaluate{
   position: fixed;
   right: 90px;
   top: 5px;
   background-color:white;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   padding:5;
}

span.top_control_bar  a.interrupt{
   position: fixed;
   right: 20px;
   top: 5px;
   background-color:#bb0000;
   font-family:arial;
   font-size:12px;
   font-weight:bold;
   color:#FFFFFF;
   padding:5;
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
   overflow: auto;
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


/************ VARIABLES **************************/

span.variables_topbar {
   z-index: 4;
   height: 24px;
   color:#FFFFFF;
   width:158px;
   top:260px;
   left: 5px;
   position: fixed;
   border:1px solid #387CAF;
   background-color: #73a6ff;
   text-decoration: none;
   font-family:arial;
   font-size:15px;
   font-weight:bold;
}

span.variable_list {
   z-index: 4;
   font-family:arial;
   font-size:12px;
   overflow: auto;
   position: fixed;
   top: 284px;
   left: 5px;
   width: 154px;
   height: 150px;
   margin: 0px;
   border:1px solid #387CAF;
   background-color: white;
   padding: 2px;
}

/************ WORKBOOKS **************************/

span.workbooks_topbar {
   z-index: 4;
   height: 24px;
   color:#FFFFFF;
   width:158px;
   top:450px;
   left: 5px;
   position: fixed;
   border:1px solid #387CAF;
   background-color: #73a6ff;
   text-decoration: none;
   font-family:arial;
   font-size:15px;
   font-weight:bold;
}

span.workbooks_list {
   z-index: 4;
   font-family:arial;
   font-size:12px;
   overflow: auto;
   position: fixed;
   top: 474px;
   left: 5px;
   width: 154px;
   height: 150px;
   margin: 0px;
   border:1px solid #387CAF;
   background-color: white;
   padding: 2px;
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

span.workbook {
  position:absolute;
  text-align:center;
  background-color: #c3d9ff;
  border-left: 1px solid #387CAF;
  top: 40px;
  left: 170px;
  margin: 0px;
  padding: 8px;
  float: right;
  horizontal-align: center;
}


span.banner{
  margin: 0px;
  position: fixed;
  top: 0px;
  left: 10px;
  padding: 5px;
}

a.banner {
  background-color:white;
  font-family:arial;
  font-size:30px;
  text-decoration: none;
  weight: bold;
  color: #387CAF;
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

input.txtarea {
  color: black;
  text-decoration: none;
  background: white;
  padding: 0px;
  margin: 0px;
  border: 1px solid #387CAF;
  width: 70%;
}

div.output_box {
  overflow: auto;
  height: 60%;
}

"""
