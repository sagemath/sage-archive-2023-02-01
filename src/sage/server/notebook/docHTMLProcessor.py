"""
Processes SAGE documentation into notebook worksheet format with
evaluatable examples.

This takes in any HTML document, i.e. sage documentation, and returns it in
the editable format (like the notebook edit window). It also returns a
string representing the css link for the document.
The SGML parser is setup to return only the body of the html documentation
page and to re-format sage examples and type-setting.

Note:
This extension of sgmllib.SGMLParser was partly inspired by Mark Pilgrim's 'Dive Into Python' examples.

Author:
    -- Dorian Raymer (2006)


"""

from sgmllib import SGMLParser
from urllib import splittag
from htmlentitydefs import entitydefs

class DocHTMLProcessor(SGMLParser):

    def reset(self):
        """ This function is called by SGMLParser.__init__ so all necessary things
        are initiallized here.
        """
        # flags
        self.bodyQ = False #don't keep anything before the <body> tag
        self.in_verbatim_div = False
        self.in_math_span = False
        self.in_mathdisplay = False

        # lists of what the parser keeps
        self.temp_pieces = []
        # self.all_pieces = []
        self.all_pieces = ''
        self.css_href = None

        # counters
        self.cellcount = 0
        self.allcount = 0

        SGMLParser.reset(self)

    def process_doc_html(self, doc_path, full_path, doc_in):
        """process_doc_html is the only function that needs to be called externally.
        docin should be a properly marked up html file.
        doc_folder tells what part of the documentation (''=main index, ref = reference, tut=tutorial, etc.)
        self.feed() is a SGMLParser method and starts everything off; Most of the functions here
        are extensions to SGMLParser, and may never actually be visibly called here.
        """
        self.doc_path = doc_path
        self.full_path = full_path
        self.feed(doc_in) #SGMLParser call
        self.close()     #SGMLParser call
        self.hand_off_temp_pieces('to_doc_pieces')
        self.all_pieces = self.all_pieces[:-16]  # drop </body></html>
        return self.all_pieces, self.css_href # The goods


    def hand_off_temp_pieces(self, piece_type):
        """ To seperate documentation content from sage examples, everything is split into one of two cell types.
        This function is called to put the current self.temp_pieces into self.all_pieces.
        """
        pieces = "".join(self.temp_pieces)
        pieces = pieces.lstrip()
        if piece_type=='to_doc_pieces':
            # pieces = '%html\n' + pieces
            # self.all_pieces.append(pieces)
            self.all_pieces += pieces
            self.temp_pieces = []
        else:
            pieces = self.process_cell_input_output(pieces)
            # self.all_pieces.append(pieces)
            self.all_pieces += pieces
            self.temp_pieces = []
        self.allcount += 1

    def process_cell_input_output(self, cell_piece):
        """
        All class='verbatim' div's contain code examples.
        Some examples are models of how the function works;
        those begin with INPUT: or something.
        The rest of the examples should have sage:input and
        output. If the example is a model, it is made into a
        div class='usage_model' so it can be stylized.
        If it is actuall input/output, the input is seperated
        from the output according to the Notebook edit format.
        """
        if cell_piece[:5] != 'sage:' and cell_piece[:12] != '&gt;'*3:
            piece = '<div class="verbatim"><pre>'
            piece += cell_piece
            piece = piece.replace('{','{&nbsp;')
            piece = piece.replace('}','}&nbsp;')
            piece += '</pre></div>'
        else:
            # group and format inputs and outputs
            pieces = cell_piece.split('\n')
            output_flag = False
            piece = '{{{\n'
            for p in pieces:
                p = p.lstrip()

                if p[:5] == 'sage:' and not output_flag:
                    piece += p[5:].lstrip() + '\n'
                elif p[:5] == 'sage:' and output_flag:
                    piece += '}}}\n{{{\n' + p[5:].lstrip() + '\n'
                    output_flag = False
                elif p[:12] == '&gt;'*3 and not output_flag:
                    piece += p[12:].lstrip() + '\n'
                elif p[:12] == '&gt;'*3 and output_flag:
                    piece += '}}}\n{{{\n' + p[12:].lstrip() + '\n'
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



    def rewrite_href(self,href_value):
        # hack to make the hrefs work.
        href_value, href_tag = splittag(href_value)
        href_split = href_value.split('/')
        full_path = self.full_path
        if len(href_split) > 1:
            path = '/'.join(href_split[:-1]) + '/'
            full_path += path
            file_name = href_split[-1]
        else:
            file_name = href_value
        # parts = ''
        # for part in href_split:
        #     if part == '..':
        #         poptart = full_path.pop(-1)
        #     else:
        #         parts += '/' + part
        url_path = '/doc_browser?' + full_path + '?'

        if href_tag:
            href_new = url_path + file_name + '#' + href_tag
        else:
            href_new = url_path + file_name

        return href_new

    def rewrite_src(self, src_value):
        # src_split = src_value.split('/')
        # full_path = self.full_path.split('/')
        # for part in src_split:
        return src_value.lstrip('..')


    ##############################################
    ## General tag handlers
    ##

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

    def start_link(self, attrs):
        rel = [value.lower() for key, value in attrs if key=='rel']
        href = [value for key, value in attrs if key=='href']
        if 'stylesheet' in rel:
            self.css_href = href[0]



    def start_body(self, attrs):
        self.bodyQ = True

    def start_a(self, attrs):
        if self.bodyQ:
            count = 0
            for name, value in attrs:
                if name.lower()=='href':
                    href_new = self.rewrite_href(value)
                    attrs[count] = ('href', href_new)
                count += 1
            self.unknown_starttag('a', attrs)


    def start_div(self, attrs):

        for name, value in attrs:
            if name.lower()=='class' and value.lower()=='verbatim':
                self.in_verbatim_div = True
                return
            if name.lower()=='class' and value.lower()=='mathdisplay':
                self.in_mathdisplay = True #left off here
        self.unknown_starttag('div', attrs)


    def end_div(self):
        if self.in_verbatim_div:
            #self.temp_pieces.append(" }}} ")
            self.in_verbatim_div = False
            self.hand_off_temp_pieces('to_cell_pieces')
            return
        self.temp_pieces.append("</div>")

    def start_pre(self, attrs):
        if self.in_verbatim_div:
            self.hand_off_temp_pieces('to_doc_pieces')
            self.cellcount += 1
        else:
            self.unknown_starttag('pre',attrs)

    def end_pre(self):
        if not self.in_verbatim_div:
            self.unknown_endtag('pre')

    def start_span(self, attrs):
        count = 0
        for name, value in attrs:
            if name.lower()=='class' and value.lower()=='math':
                self.in_math_span = True
                attrs[count] = ('class','math')
            count += 1
        self.unknown_starttag('span', attrs)

    def end_span(self):
        if self.in_math_span:
            self.in_math_span = False
        self.unknown_endtag('span')

    def start_img(self, attrs):
        # if in a span with class=math,
        # remove the following img, and just print the alt attribute
        if self.bodyQ:
            if self.in_math_span:
                for name,value in attrs:
                    if name.lower()=='alt':
                        # value = value.replace('$','\\$')
                        tex = value
                        # tex = '\\text{' + value + '}'
                        self.temp_pieces.append(tex)
                        return
            count = 0
            for name, value in attrs:
                if name.lower()=='src':
                    # attrs[count] = ('src',self.doc_path + '/' + value.lstrip('..'))
                    attrs[count] = ('src',self.doc_path + self.full_path + value)
                count += 1
            strattrs = "".join([' %s="%s"' % (key, value) for key, value in attrs])
            self.temp_pieces.append("<img %(strattrs)s>" % locals())



    def end_img(self):
        if self.in_math_span == 1:
            return
        self.unknown_endtag('img')



