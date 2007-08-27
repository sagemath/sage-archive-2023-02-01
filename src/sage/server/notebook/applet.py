

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
