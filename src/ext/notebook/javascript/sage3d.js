/*
 *
 * This manages the 3d viewer applet.
 * This is done centrally because we only want one instance of
 * the actual applet to ever get run.
 *
 * AUTHOR:
 *     -- Robert Bradshaw
 *
 */


var applet_tag = 'Loading Java 3d Libraries... <applet id="sage3d" code="org.jdesktop.applet.util.JNLPAppletLauncher"  \
width="500" height="22"                                                       \
codebase="/java/java3d"                                                          \
archive="sage3d.jar,                                                      \
         sun-libs/applet-launcher.jar,sun-libs/j3dcore.jar,sun-libs/j3dutils.jar,sun-libs/vecmath.jar,sun-libs/jogl.jar,sun-libs/gluegen-rt.jar">      \
<param name="codebase_lookup" value="false">                                  \
<param name="subapplet.displayname" value="Java 3D Viewer">                   \
<param name="jnlpNumExtensions" value="1">                                    \
<param name="progressbar" value="true">                                       \
<param name="noddraw.check" value="true">                                     \
                                                                              \
<param name="subapplet.classname" value="org.sagemath.sage3d.ObjectViewerApplet">       \
<param name="jnlpExtension1" value="'+document.location.protocol+'//'+document.location.host+'/java/sun-libs/java3d-latest.jnlp">              \
</applet>'

var test_applet_tag = '<applet id="sage3d" code="org.sagemath.TestApplet"     \
width="100" height="20"                                                       \
codebase="/java/java3d"                                                          \
archive="sage3d.jar"                                                      \
</applet>xxxx'


var sage3d_div;
var launcher_applet;
var sage3d_applet;

function sage3d_init() {
  sage3d_div = document.createElement("div");
  sage3d_div.innerHTML = applet_tag;
  //document.body.appendChild(div);
  document.body.insertBefore(sage3d_div, document.body.firstChild)
  launcher_applet = document.getElementById("sage3d");
}

function sage3d_show(url, cell, name) {
  //if (!confirm("sage3d_show")) return;
  if (launcher_applet == undefined) {
    sage3d_init();
  }
  sage3d_applet = launcher_applet.getSubApplet();
  if (sage3d_applet == undefined) {
    setTimeout(function() {sage3d_show(url, cell, name); }, 500);
  }
  else if (sage3d_div.style.visibility != "hidden") {
    sage3d_div.style.visibility = "hidden";
    sage3d_div.style.position = "absolute";
    setTimeout(function() {sage3d_show(url, cell, name); }, 500);
  }
  else {
    sage3d_applet.showView(url, cell, name);
  }
}

