from sage.server.notebook.applet import Applet

sun_libs = ["libs/objtest.jar", "applet-launcher.jar", "j3dcore.jar", "j3dutils.jar", "vecmath.jar", "jogl.jar", "gluegen-rt.jar"]
sun_params = { "codebase_lookup": "false",
               "subapplet.classname": "org.sagemath.sage3d.ObjectViewerApplet",
               "subapplet.displayname": "SAGE Java 3D Viewer",
               "jnlpNumExtensions": 1,
               "jnlpExtension1": "/java/3d/sun-libs/java3d-latest.jnlp",
               "progressbar": "true",
               "noddraw.check":"true"
               }


class Java3DApplet(Applet):
    def __init__(self, width=400, height=400, inline=None):
        Applet.__init__(self,
                        code = "org.jdesktop.applet.util.JNLPAppletLauncher",
                        width = width,
                        height = height,
                        archive = ["lib/sage3d.jar"] + ["sun-libs/" + lib for lib in sun_libs],
                        codebase = "3d",
                        params = sun_params)
