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

var applet_tag = '<applet id="sage3d" code="org.jdesktop.applet.util.JNLPAppletLauncher"  \
width="100" height="20"                                                       \
codebase="/java/3d/"                                                          \
archive="lib/sage3d.jar,                                                      \
         sun-libs/applet-launcher.jar,sun-libs/j3dcore.jar,sun-libs/j3dutils.jar,sun-libs/vecmath.jar,sun-libs/jogl.jar,sun-libs/gluegen-rt.jar">      \
<param name="codebase_lookup" value="false">                                  \
<param name="subapplet.displayname" value="Java 3D Viewer">                   \
<param name="jnlpNumExtensions" value="1">                                    \
<param name="jnlpExtension1" value="extlibs/java3d-latest.jnlp">              \
<param name="progressbar" value="true">                                       \
<param name="noddraw.check" value="true">                                     \
                                                                              \
<param name="subapplet.classname" value="org.sagemath.sage3d.ObjectViewerApplet">                 \
</applet>'

var test_applet_tag = '<applet id="sage3d" code="org.sagemath.TestApplet"     \
width="100" height="20"                                                       \
codebase="/java/3d/"                                                          \
archive="lib/sage3d.jar"                                                      \
</applet>'


var sage3d_applet;

function sage3d_init() {
  div = document.createElement("div");
  div.style.visibility = "hidden";
  div.innerHTML = test_applet_tag;
  document.body.appendChild(div);
  sage3d_applet = document.getElementById("sage3d");
}

function sage3d_show(url, cell, name) {
  if (sage3d_applet == undefined) {
    sage3d_init();
  }
  sage3d_applet.showView(url, cell, name);
}
