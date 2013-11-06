r"""
Interactive OpenGL based 3D

TESTS:

Make sure this module can be imported properly::

    sage: import sage.plot.java3d
"""

#*****************************************************************************
#       Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************


from sagenb.notebook.applet import Applet

sun_libs = ["libs/objtest.jar", "applet-launcher.jar", "j3dcore.jar", "j3dutils.jar", "vecmath.jar", "jogl.jar", "gluegen-rt.jar"]
sun_params = { "codebase_lookup": "false",
               "subapplet.classname": "org.sagemath.sage3d.ObjectViewerApplet",
               "subapplet.displayname": "Sage Java 3-D Viewer",
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

