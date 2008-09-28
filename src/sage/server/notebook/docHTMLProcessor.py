"""
Live Documentation in the Notebook

Processes Sage documentation into notebook worksheet format with
evaluatable examples.

This takes in any HTML document, i.e. Sage documentation, and returns
it in the editable format (like the notebook edit window). It also
returns a string representing the css link for the document.  The SGML
parser is setup to return only the body of the html documentation page
and to re-format Sage examples and type-setting.

Note: This extension of sgmllib.SGMLParser was partly inspired by Mark
Pilgrim's 'Dive Into Python' examples.

Author:
    -- Dorian Raymer (2006): first version
    -- William Stein (2007-06-10): rewrite to work with twisted SAGE notebook
    -- Mike Hansen (2008-09-27): Rewrite to work with Sphinx HTML documentation
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
        This function is called by SGMLParser.__init__ so all necessary things
        are initiallized here.

        EXAMPLES:
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
        process_doc_html is the only function that needs to be called
        externally.  docin should be a properly marked up html file.

        self.feed() is a
        SGMLParser method and starts everything off; Most of the
        functions here are extensions to SGMLParser, and may never
        actually be visibly called here.
        """
        self.feed(doc_in) #SGMLParser call
        self.close()     #SGMLParser call
        self.hand_off_temp_pieces('to_doc_pieces')
        self.all_pieces = self.all_pieces[:-16]  # drop </body></html>
        return self.all_pieces


    def hand_off_temp_pieces(self, piece_type):
        """
        To seperate documentation content from sage examples,
        everything is split into one of two cell types.  This function
        is called to put the current self.temp_pieces into
        self.all_pieces.
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
        Returns the current cell count and increments it
        by one.

        EXAMPLES:
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
        All class='highlight' div's contain code examples.
        Some examples are models of how the function works;
        those begin with INPUT: or something.
        The rest of the examples should have sage:input and
        output. If the example is a model, it is made into a
        div class='usage_model' so it can be stylized.
        If it is actuall input/output, the input is seperated
        from the output according to the Notebook edit format.
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
                p = p.lstrip()

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
                    # first occurrence of an output string
                    # write /// denoting output
                    if output_flag == False:
                        piece += '///\n'
                        piece += p.lstrip() + '\n'
                        output_flag = True
                    # multiple output lines exist, don't need /// repeated
                    else:
                        piece += p.lstrip() + '\n'
            piece += '}}}\n'
        return piece

    ##############################################
    ## General tag handlers
    ## These just append their HTML to
    ## self.temp_pieces.

    def	unknown_starttag(self, tag, attrs):
        if self.bodyQ:
            strattrs = "".join([' %s="%s"' % (key, value) for key, value in attrs])
            self.temp_pieces.append("<%(tag)s%(strattrs)s>" % locals())

    def	unknown_endtag(self, tag):
        if self.bodyQ:
            self.temp_pieces.append("</%(tag)s>" % locals())

    def	handle_data(self, data):
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

    def	handle_comment(self, data):
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
        This just sets self.bodyQ to True once we've hit the body tag.

        EXAMPLES:
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
