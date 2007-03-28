# -*- coding: utf-8 -*-
"""
    Colorize - Python source formatter that outputs Python code in XHTML.
    This script is based on MoinMoin - The Python Source Parser.

    FROM: Modified version of
      http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/442482
"""

from sage.misc.misc import SAGE_EXTCODE


# Imports
import cgi
import string
import sys
import cStringIO
import keyword
import token
import tokenize
import re
import os

#Set up basic values.
_KEYWORD = token.NT_OFFSET + 1
_TEXT    = token.NT_OFFSET + 2

_classes = {
    token.NUMBER:       'token_number',
    token.OP:           'token_op',
    token.STRING:       'token_string',
    tokenize.COMMENT:   'token_comment',
    token.NAME:         'token_name',
    token.ERRORTOKEN:   'token_error',
    _KEYWORD:           'keyword',
    _TEXT:              'text',
}

class Parser:
    """ Send colored python source.
    """

    def __init__(self, raw, out = sys.stdout):
        """ Store the source text.
        """
        self.raw = string.strip(string.expandtabs(raw))
        self.out = out

    def format(self, formatter, form):
        """ Parse and send the colored source.
        """
        # store line offsets in self.lines
        self.lines = [0, 0]
        pos = 0
        while 1:
            pos = string.find(self.raw, '\n', pos) + 1
            if not pos: break
            self.lines.append(pos)
        self.lines.append(len(self.raw))

        # parse the source and write it
        self.pos = 0
        text = cStringIO.StringIO(self.raw)
        try:
            tokenize.tokenize(text.readline, self)
        except tokenize.TokenError, ex:
            msg = ex[0]
            line = ex[1][0]
            self.out.write("<h3>ERROR: %s</h3>%s\n" % (
                msg, self.raw[self.lines[line]:]))


    def __call__(self, toktype, toktext, (srow,scol), (erow,ecol), line):
        """ Token handler.
        """
        if 0:
            print "type", toktype, token.tok_name[toktype], "text", toktext,
            print "start", srow, scol, "end", erow, ecol, "<br />"

        # calculate new positions
        oldpos = self.pos
        newpos = self.lines[srow] + scol
        self.pos = newpos + len(toktext)

        # handle newlines
        if toktype in [token.NEWLINE, tokenize.NL]:
            self.out.write('\n')
            return

        # send the original whitespace, if needed
        if newpos > oldpos:
            self.out.write(self.raw[oldpos:newpos])

        # skip indenting tokens
        if toktype in [token.INDENT, token.DEDENT]:
            self.pos = newpos
            return

        # map token type to a color/class group
        if token.LPAR <= toktype and toktype <= token.OP:
            toktype = token.OP
        elif toktype == token.NAME and keyword.iskeyword(toktext):
            toktype = _KEYWORD
        classval = _classes.get(toktype, _classes[_TEXT])

        style = ''
        if toktype == token.ERRORTOKEN:
            style = ' style="border: solid 1.5pt #FF0000;"'

        # send text
        self.out.write('<span class="%s"%s>' % (classval, style))
        self.out.write(cgi.escape(toktext))
        self.out.write('</span>')


def colorize(source):
    import os, sys
    # write colorized version to "[filename].py.html"
    html = cStringIO.StringIO()
    Parser(source, html).format(None, None)
    html.flush()
    html.seek(0)
    return html.read()
