var jmol_count = 0
function jmol_applet(size, url) {
    jmolSetDocument(cell_writer);
    jmolApplet(size, "script " + url, jmol_count);
    var s = " <a href='#' onclick='jmol_image("+jmol_count+");return false;'>Get Image</a>";
    cell_writer.write(s);
    jmol_count += 1;
    return s
}

function jmol_image(jmol_count) {
  var myImage = jmolGetPropertyAsString("image","",jmol_count)
  mywindow = window.open("","Jmol Image","menubar=no,width=600,height=600,toolbar=no");
  s = '<HTML><TITLE>Jmol Image</TITLE><BODY>';
  s += '<img src="data:image/jpeg;base64,' + myImage + '">'
  s += '<p>To save this image, you can try right-clicking on the image to copy it or save it to a file, or you may be able to just drag the image to your desktop.</p>'
  s += '</BODY></HTML>';
  mywindow.document.write(s);
}

function jmol_popup(url) {
    win = window.open ("", "jmol viewer", "width=600,height=600,resizable=1,statusbar=0");
    win.document.body.innerHTML = "";
    win.document.title = "Sage 3d Viewer";
    win.document.writeln("<h1 align=center>Sage 3d Viewer</h1>");
    jmolSetDocument(win.document);
    jmolApplet("100%", "script" + url, jmol_count);
    win.focus();
}
