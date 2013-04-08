# This file is part of the OLD Sage notebook and is NOT actively developed,
# maintained, or supported.  As of Sage v4.1.2, all notebook development has
# moved to the separate Sage Notebook project:
#
# http://nb.sagemath.org/
#
# The new notebook is installed in Sage as an spkg (e.g., sagenb-0.3.spkg).
#
# Please visit the project's home page for more information, including directions on
# obtaining the latest source code.  For notebook-related development and support,
# please consult the sage-notebook discussion group:
#
# http://groups.google.com/group/sage-notebook


class Applet:

    def __init__(self, id, code, archive, codebase="", width=400, height=400, params={}):
        self.id = id
        self.code = code
        self.archive = archive
        self.width = width
        self.height = height
        self.params = params
        self.codebase = "/java/" + codebase

    def html_tag(self):
        params_text = "\n".join(["""<param name="%s" value="%s"/>""" % x for x in self.params.iteritems()])
        tag = """
        <applet id="%s", code="%s" width="%s" height="%s" codebase="%s" archive="%s" MAYSCRIPT>
          %s
        </applet>
        """ % (self.id,
               self.code,
               self.width,
               self.height,
               self.codebase,
               ",".join(self.archive),
               params_text)
        return tag

