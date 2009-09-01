"""
Live Documentation in the Notebook

Processes Sage documentation into notebook worksheet format with
evaluable examples.

This takes in any HTML document, i.e., Sage documentation, and returns
it in the editable format (like the notebook edit window). It also
returns a string representing the CSS link for the document.  The SGML
parser is setup to return only the body of the HTML documentation page
and to re-format Sage examples and type-setting.

Note: This extension of sgmllib.SGMLParser was partly inspired by Mark
Pilgrim's 'Dive Into Python' examples.

Author:

- Dorian Raymer (2006): first version

- William Stein (2007-06-10): rewrite to work with twisted Sage notebook

- Mike Hansen (2008-09-27): Rewrite to work with Sphinx HTML documentation
"""
#############################################################################
#       Copyright (C) 2007 William Stein <wstein@gmail.com> and Dorian Raimer
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
#############################################################################

from sgmllib import SGMLParser
from urllib import splittag
from htmlentitydefs import entitydefs

class SphinxHTMLProcessor(SGMLParser):
    def reset(self):
        """
        Initialize necessary variables.  Called by
        :meth:`SGMLParser.__init__`.

        EXAMPLES::

            sage: from sage.server.notebook.docHTMLProcessor import SphinxHTMLProcessor
            sage: d = SphinxHTMLProcessor()
            sage: d.bodyQ
            False
            sage: d.in_highlight_div
            False
            sage: d.temp_pieces
            []
            sage: d.all_pieces
            ''
            sage: d.cellcount
            0
        """
        # flags
        self.bodyQ = False #don't keep anything before the <body> tag
        self.in_highlight_div = False

        # lists of what the parser keeps
        self.temp_pieces = []
        self.all_pieces = ''

        # counters
        self.cellcount = 0

        SGMLParser.reset(self)

    def process_doc_html(self, doc_in):
        """
        Returns processed HTML input as HTML output.  This is the only
        method that needs to be called externally.

        INPUT:

        - ``doc_in`` - a string containing properly formed HTML

        OUTPUT:

        - a string; the processed HTML
        """
        # self.feed() is a SGMLParser method and starts everything
        # off; Most of the functions here are extensions to
        # SGMLParser, and may never actually be visibly called here.
        self.feed(doc_in) #SGMLParser call
        self.close()     #SGMLParser call
        self.hand_off_temp_pieces('to_doc_pieces')
        self.all_pieces = self.all_pieces[:-16]  # drop </body></html>
        return self.all_pieces


    def hand_off_temp_pieces(self, piece_type):
        """
        To separate the documentation's content from the Sage
        examples, everything is split into one of two cell types.
        This method puts the current ``self.temp_pieces`` into
        ``self.all_pieces``.

        INPUT:

        - ``piece_type`` - a string; indicates the type of and how to
          process the current ``self.temp_pieces``
        """
        pieces = "".join(self.temp_pieces)
        pieces = pieces.lstrip()
        if piece_type == 'to_doc_pieces':
            self.all_pieces += pieces
            self.temp_pieces = []
        elif piece_type == 'ignore':
            self.temp_pieces = []
        else:
            pieces = self.process_cell_input_output(pieces)
            self.all_pieces += pieces
            self.temp_pieces = []

    def get_cellcount(self):
        """
        Return the current cell count and increment it by one.

        OUTPUT:

        - an int

        EXAMPLES::

            sage: from sage.server.notebook.docHTMLProcessor import SphinxHTMLProcessor
            sage: d = SphinxHTMLProcessor()
            sage: d.get_cellcount()
            0
            sage: d.get_cellcount()
            1
        """
        self.cellcount += 1
        return self.cellcount - 1

    def process_cell_input_output(self, cell_piece):
        """
        Process and return a ``cell_piece``.

        All divs with CSS class="highlight" contain code examples.
        They include

        - Models of how the function works.  These begin with, e.g.,
          'INPUT:' and are re-styled as divs with
          class="usage_model".

        - Actual Sage input and output.  These begin with 'sage:'.
          The input and output are separated according to the
          Notebook edit format.

        INPUT:

        - ``cell_piece`` - a string; a cell piece

        OUTPUT:

        - a string; the processed cell piece
        """
        if cell_piece[:5] != 'sage:' and cell_piece[:12] != '&gt;'*3:
            piece = '<div class="highlight"><pre>'
            piece += cell_piece
            piece = piece.replace('{','{&nbsp;')
            piece = piece.replace('}','}&nbsp;')
            piece += '</pre></div>'
        else:
            # group and format inputs and outputs
            pieces = cell_piece.split('\n')
            output_flag = False
            piece = '{{{id=%s|\n'%self.get_cellcount()
            for p in pieces:

                if p[:5] == 'sage:' and not output_flag:
                    piece += p[5:].lstrip() + '\n'
                elif p[:5] == 'sage:' and output_flag:
                    piece += '}}}\n{{{id=%s|\n'%self.get_cellcount() + p[5:].lstrip() + '\n'
                    output_flag = False
                elif p[:12] == '&gt;'*3 and not output_flag:
                    piece += p[12:].lstrip() + '\n'
                elif p[:12] == '&gt;'*3 and output_flag:
                    piece += '}}}\n{{{id=%s|\n'%self.get_cellcount() + p[12:].lstrip() + '\n'
                    output_flag = False
                elif p[:3] == '...':
                    piece += p[3:] + '\n'
                else:
                    # in an output string. replace escaped html
                    # strings so they don't get converted twice.
                    p = p.replace('&lt;', '<')
                    p = p.replace('&gt;', '>')
                    p = p.replace('&amp;', '&')
                    p = p.replace('&#39;', "'")
                    # first occurrence of an output string
                    # write /// denoting output
                    if output_flag == False:
                        piece += '///\n'
                        piece += p + '\n'
                        output_flag = True
                    # multiple output lines exist, don't need /// repeated
                    else:
                        piece += p + '\n'
            piece += '}}}\n'
        return piece

    ##############################################
    ## General tag handlers
    ## These just append their HTML to
    ## self.temp_pieces.

    def unknown_starttag(self, tag, attrs):
        if self.bodyQ:
            strattrs = "".join([' %s="%s"' % (key, value) for key, value in attrs])
            self.temp_pieces.append("<%(tag)s%(strattrs)s>" % locals())

    def unknown_endtag(self, tag):
        if self.bodyQ:
            self.temp_pieces.append("</%(tag)s>" % locals())

    def handle_data(self, data):
        if self.bodyQ:
            self.temp_pieces.append(data)

    def handle_charref(self, ref):
        if self.bodyQ:
            self.temp_pieces.append("&#%(ref)s;" % locals())

    def handle_entityref(self, ref):
        if self.bodyQ:
            self.temp_pieces.append("&%(ref)s" % locals())
            if entitydefs.has_key(ref):
                self.temp_pieces.append(';')

    def handle_comment(self, data):
        if self.bodyQ:
            self.temp_pieces.append("<!--%(data)s-->" % locals())

    def handle_pi(self, text):
        if self.bodyQ:
            self.temp_pieces.append("<?%(text)s>" % locals())

    def handle_decl(self, text):
        if self.bodyQ:
            self.temp_pieces.append("<!%(text)s>" % locals())


    #############################################
    ## Specific tag handlers
    ##
    def start_body(self, attrs):
        """
        Set ``self.bodyQ`` to True upon finding the opening body tag.

        INPUT:

        - ``attrs`` - a string:string dictionary containing the
          element's attributes

        EXAMPLES::

            sage: from sage.server.notebook.docHTMLProcessor import SphinxHTMLProcessor
            sage: d = SphinxHTMLProcessor()
            sage: d.bodyQ
            False
            sage: d.start_body(None)
            sage: d.bodyQ
            True
        """
        self.bodyQ = True

    def start_div(self, attrs):
        #Find out if we are starting a highlighted div
        for name, value in attrs:
            if name.lower()=='class' and value.lower()=='highlight':
                self.in_highlight_div = True
                return
        self.unknown_starttag('div', attrs)

    def end_div(self):
        #Once we end the highlighted div, convert all of the pieces
        #to cells
        if self.in_highlight_div:
            self.in_highlight_div = False
            self.hand_off_temp_pieces('to_cell_pieces')
            return
        self.temp_pieces.append("</div>")

    def start_pre(self, attrs):
        #Once we hit the <pre> tag in a highlighted block,
        #hand of all of the pieces we've encountered so far
        #and ignore the tag.
        if self.in_highlight_div:
            self.hand_off_temp_pieces('to_doc_pieces')
            return
        self.unknown_starttag('pre',attrs)

    def end_pre(self):
        #Ignore the pre tags in highlighted blocks
        if self.in_highlight_div:
            return
        self.unknown_endtag('pre')

    #Ignore forms
    def start_form(self, attrs):
        #Hand of everything we've accumulated so far
        self.hand_off_temp_pieces('to_doc_pieces')
        return

    def end_form(self):
        #Ignore all of the pieces since we started
        #the form.
        self.hand_off_temp_pieces('ignore')
        return

    def start_span(self, attrs):
        #Ignore all spans that occur within highlighted blocks
        if self.in_highlight_div:
            return
        self.unknown_starttag('span', attrs)

    def end_span(self):
        #Ignore all spans that occur within highlighted blocks
        if self.in_highlight_div:
            return
        self.unknown_endtag('span')
