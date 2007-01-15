# -*- coding: iso-8859-1 -*-
"""
    MoinMoin - MoinMoin Wiki Markup Parser

    @copyright: 2000, 2001, 2002 by Jürgen Hermann <jh@web.de>
    @license: GNU GPL, see COPYING for details.
"""

import os, re
from MoinMoin import config, wikimacro, wikiutil
from MoinMoin.Page import Page
from MoinMoin.util import web

Dependencies = []

class Parser:
    """
        Object that turns Wiki markup into HTML.

        All formatting commands can be parsed one line at a time, though
        some state is carried over between lines.

        Methods named like _*_repl() are responsible to handle the named regex
        patterns defined in print_html().
    """

    # allow caching
    caching = 1
    Dependencies = []

    # some common strings
    PARENT_PREFIX = wikiutil.PARENT_PREFIX
    attachment_schemas = ["attachment", "inline", "drawing"]
    punct_pattern = re.escape(u'''"\'}]|:,.)?!''')
    url_pattern = (u'http|https|ftp|nntp|news|mailto|telnet|wiki|file|irc|' +
            u'|'.join(attachment_schemas) +
            (config.url_schemas and u'|' + u'|'.join(config.url_schemas) or ''))

    # some common rules
    word_rule = ur'(?:(?<![%(u)s%(l)s])|^)%(parent)s(?:%(subpages)s(?:[%(u)s][%(l)s]+){2,})+(?![%(u)s%(l)s]+)' % {
        'u': config.chars_upper,
        'l': config.chars_lower,
        'subpages': wikiutil.CHILD_PREFIX + '?',
        'parent': ur'(?:%s)?' % re.escape(PARENT_PREFIX),
    }
    url_rule = ur'%(url_guard)s(%(url)s)\:([^\s\<%(punct)s]|([%(punct)s][^\s\<%(punct)s]))+' % {
        'url_guard': u'(^|(?<!\w))',
        'url': url_pattern,
        'punct': punct_pattern,
    }

    ol_rule = ur"^\s+(?:[0-9]+|[aAiI])\.(?:#\d+)?\s"
    dl_rule = ur"^\s+.*?::\s"

    # the big, fat, ugly one ;)
    formatting_rules = ur"""(?P<ent_numeric>&#(\d{1,5}|x[0-9a-fA-F]+);)
(?:(?P<emph_ibb>'''''(?=[^']+'''))
(?P<emph_ibi>'''''(?=[^']+''))
(?P<emph_ib_or_bi>'{5}(?=[^']))
(?P<emph>'{2,3})
(?P<u>__)
(?P<sup>\^.*?\^)
(?P<sub>,,[^,]{1,40},,)
(?P<tt>\{\{\{.*?\}\}\})
(?P<processor>(\{\{\{(#!.*|\s*$)))
(?P<pre>(\{\{\{ ?|\}\}\}))
(?P<small>(\~- ?|-\~))
(?P<big>(\~\+ ?|\+\~))
(?P<strike>(--\(|\)--))
(?P<rule>-{4,})
(?P<comment>^\#\#.*$)
(?P<macro>\[\[(%%(macronames)s)(?:\(.*?\))?\]\]))
(?P<ol>%(ol_rule)s)
(?P<dl>%(dl_rule)s)
(?P<li>^\s+\*\s*)
(?P<li_none>^\s+\.\s*)
(?P<indent>^\s+)
(?P<tableZ>\|\| $)
(?P<table>(?:\|\|)+(?:<[^>]*?>)?(?!\|? $))
(?P<heading>^\s*(?P<hmarker>=+)\s.*\s(?P=hmarker) $)
(?P<interwiki>[A-Z][a-zA-Z]+\:[^\s'\"\:\<\|]([^\s%(punct)s]|([%(punct)s][^\s%(punct)s]))+)
(?P<word>%(word_rule)s)
(?P<url_bracket>\[((%(url)s)\:|#|\:)[^\s\]]+(\s[^\]]+)?\])
(?P<url>%(url_rule)s)
(?P<email>[-\w._+]+\@[\w-]+(\.[\w-]+)+)
(?P<smiley>(?<=\s)(%(smiley)s)(?=\s))
(?P<smileyA>^(%(smiley)s)(?=\s))
(?P<ent_symbolic>&[a-zA-Z]+;)
(?P<ent>[<>&])
(?P<wikiname_bracket>\[".*?"\])
(?P<tt_bt>`.*?`)"""  % {

        'url': url_pattern,
        'punct': punct_pattern,
        'ol_rule': ol_rule,
        'dl_rule': dl_rule,
        'url_rule': url_rule,
        'word_rule': word_rule,
        'smiley': u'|'.join(map(re.escape, config.smileys.keys()))}

    # Don't start p before these
    no_new_p_before = ("heading rule table tableZ tr td "
                       "ul ol dl dt dd li li_none indent "
                       "macro processor pre")
    no_new_p_before = no_new_p_before.split()
    no_new_p_before = dict(zip(no_new_p_before, [1] * len(no_new_p_before)))

    def __init__(self, raw, request, **kw):
        self.raw = raw
        self.request = request
        self.form = request.form
        self._ = request.getText
        self.cfg = request.cfg
        self.line_anchors = kw.get('line_anchors', True)
        self.macro = None
        self.start_line = kw.get('start_line', 0)

        self.is_em = 0
        self.is_b = 0
        self.is_u = 0
        self.is_strike = 0
        self.lineno = 0
        self.in_list = 0 # between <ul/ol/dl> and </ul/ol/dl>
        self.in_li = 0 # between <li> and </li>
        self.in_dd = 0 # between <dd> and </dd>
        self.in_pre = 0
        self.in_table = 0
        self.is_big = False
        self.is_small = False
        self.inhibit_p = 0 # if set, do not auto-create a <p>aragraph
        self.titles = request._page_headings

        # holds the nesting level (in chars) of open lists
        self.list_indents = []
        self.list_types = []

        self.formatting_rules = self.formatting_rules % {'macronames': u'|'.join(wikimacro.getNames(self.cfg))}

    def _close_item(self, result):
        #result.append("<!-- close item begin -->\n")
        if self.in_table:
            result.append(self.formatter.table(0))
            self.in_table = 0
        if self.in_li:
            self.in_li = 0
            if self.formatter.in_p:
                result.append(self.formatter.paragraph(0))
            result.append(self.formatter.listitem(0))
        if self.in_dd:
            self.in_dd = 0
            if self.formatter.in_p:
                result.append(self.formatter.paragraph(0))
            result.append(self.formatter.definition_desc(0))
        #result.append("<!-- close item end -->\n")


    def interwiki(self, url_and_text, **kw):
        # TODO: maybe support [wiki:Page http://wherever/image.png] ?
        if len(url_and_text) == 1:
            url = url_and_text[0]
            text = None
        else:
            url, text = url_and_text

        # keep track of whether this is a self-reference, so links
        # are always shown even the page doesn't exist.
        is_self_reference = 0
        url2 = url.lower()
        if url2.startswith('wiki:self:'):
            url = url[10:] # remove "wiki:self:"
            is_self_reference = 1
        elif url2.startswith('wiki:'):
            url = url[5:] # remove "wiki:"

        tag, tail = wikiutil.split_wiki(url)
        if text is None:
            if tag:
                text = tail
            else:
                text = url
                url = ""
        elif (url.startswith(wikiutil.CHILD_PREFIX) or # fancy link to subpage [wiki:/SubPage text]
              is_self_reference or # [wiki:Self:LocalPage text] or [:LocalPage:text]
              Page(self.request, url).exists()): # fancy link to local page [wiki:LocalPage text]
            return self._word_repl(url, text)

        wikitag, wikiurl, wikitail, wikitag_bad = wikiutil.resolve_wiki(self.request, url)
        href = wikiutil.join_wiki(wikiurl, wikitail)

        # check for image URL, and possibly return IMG tag
        if not kw.get('pretty_url', 0) and wikiutil.isPicture(wikitail):
            return self.formatter.image(src=href)

        # link to self?
        if wikitag is None:
            return self._word_repl(wikitail)

        return (self.formatter.interwikilink(1, tag, tail) +
                self.formatter.text(text) +
                self.formatter.interwikilink(0, tag, tail))

    def attachment(self, url_and_text, **kw):
        """ This gets called on attachment URLs.
        """
        _ = self._
        if len(url_and_text) == 1:
            url = url_and_text[0]
            text = None
        else:
            url, text = url_and_text

        inline = url[0] == 'i'
        drawing = url[0] == 'd'
        url = url.split(":", 1)[1]
        url = wikiutil.url_unquote(url, want_unicode=True)
        text = text or url

        from MoinMoin.action import AttachFile
        if drawing:
            return self.formatter.attachment_drawing(url, text)

        # check for image URL, and possibly return IMG tag
        # (images are always inlined, just like for other URLs)
        if not kw.get('pretty_url', 0) and wikiutil.isPicture(url):
            return self.formatter.attachment_image(url)

        # inline the attachment
        if inline:
            return self.formatter.attachment_inlined(url, text)

        return self.formatter.attachment_link(url, text)

    def _u_repl(self, word):
        """Handle underline."""
        self.is_u = not self.is_u
        return self.formatter.underline(self.is_u)

    def _strike_repl(self, word):
        """Handle strikethrough."""
        # XXX we don't really enforce the correct sequence --( ... )-- here
        self.is_strike = not self.is_strike
        return self.formatter.strike(self.is_strike)

    def _small_repl(self, word):
        """Handle small."""
        if word.strip() == '~-' and self.is_small:
            return self.formatter.text(word)
        if word.strip() == '-~' and not self.is_small:
            return self.formatter.text(word)
        self.is_small = not self.is_small
        return self.formatter.small(self.is_small)

    def _big_repl(self, word):
        """Handle big."""
        if word.strip() == '~+' and self.is_big:
            return self.formatter.text(word)
        if word.strip() == '+~' and not self.is_big:
            return self.formatter.text(word)
        self.is_big = not self.is_big
        return self.formatter.big(self.is_big)

    def _emph_repl(self, word):
        """Handle emphasis, i.e. '' and '''."""
        ##print "#", self.is_b, self.is_em, "#"
        if len(word) == 3:
            self.is_b = not self.is_b
            if self.is_em and self.is_b:
                self.is_b = 2
            return self.formatter.strong(self.is_b)
        else:
            self.is_em = not self.is_em
            if self.is_em and self.is_b:
                self.is_em = 2
            return self.formatter.emphasis(self.is_em)

    def _emph_ibb_repl(self, word):
        """Handle mixed emphasis, i.e. ''''' followed by '''."""
        self.is_b = not self.is_b
        self.is_em = not self.is_em
        if self.is_em and self.is_b:
            self.is_b = 2
        return self.formatter.emphasis(self.is_em) + self.formatter.strong(self.is_b)

    def _emph_ibi_repl(self, word):
        """Handle mixed emphasis, i.e. ''''' followed by ''."""
        self.is_b = not self.is_b
        self.is_em = not self.is_em
        if self.is_em and self.is_b:
            self.is_em = 2
        return self.formatter.strong(self.is_b) + self.formatter.emphasis(self.is_em)

    def _emph_ib_or_bi_repl(self, word):
        """Handle mixed emphasis, exactly five '''''."""
        ##print "*", self.is_b, self.is_em, "*"
        b_before_em = self.is_b > self.is_em > 0
        self.is_b = not self.is_b
        self.is_em = not self.is_em
        if b_before_em:
            return self.formatter.strong(self.is_b) + self.formatter.emphasis(self.is_em)
        else:
            return self.formatter.emphasis(self.is_em) + self.formatter.strong(self.is_b)


    def _sup_repl(self, word):
        """Handle superscript."""
        return self.formatter.sup(1) + \
            self.formatter.text(word[1:-1]) + \
            self.formatter.sup(0)

    def _sub_repl(self, word):
        """Handle subscript."""
        return self.formatter.sub(1) + \
            self.formatter.text(word[2:-2]) + \
            self.formatter.sub(0)


    def _rule_repl(self, word):
        """Handle sequences of dashes."""
        result = self._undent() + self._closeP()
        if len(word) <= 4:
            result = result + self.formatter.rule()
        else:
            # Create variable rule size 1 - 6. Actual size defined in css.
            size = min(len(word), 10) - 4
            result = result + self.formatter.rule(size)
        return result


    def _word_repl(self, word, text=None):
        """Handle WikiNames."""

        # check for parent links
        # !!! should use wikiutil.AbsPageName here, but setting `text`
        # correctly prevents us from doing this for now
        if word.startswith(wikiutil.PARENT_PREFIX):
            if not text:
                text = word
            word = '/'.join(filter(None, self.formatter.page.page_name.split('/')[:-1] + [word[wikiutil.PARENT_PREFIX_LEN:]]))

        if not text:
            # if a simple, self-referencing link, emit it as plain text
            if word == self.formatter.page.page_name:
                return self.formatter.text(word)
            text = word
        if word.startswith(wikiutil.CHILD_PREFIX):
            word = self.formatter.page.page_name + '/' + word[wikiutil.CHILD_PREFIX_LEN:]

        # handle anchors
        parts = word.split("#", 1)
        anchor = ""
        if len(parts)==2:
            word, anchor = parts

        return (self.formatter.pagelink(1, word, anchor=anchor) +
                self.formatter.text(text) +
                self.formatter.pagelink(0, word))

    def _notword_repl(self, word):
        """Handle !NotWikiNames."""
        return self.formatter.nowikiword(word[1:])

    def _interwiki_repl(self, word):
        """Handle InterWiki links."""
        wikitag, wikiurl, wikitail, wikitag_bad = wikiutil.resolve_wiki(self.request, word)
        if wikitag_bad:
            return self.formatter.text(word)
        else:
            return self.interwiki(["wiki:" + word])


    def _url_repl(self, word):
        """Handle literal URLs including inline images."""
        scheme = word.split(":", 1)[0]

        if scheme == "wiki":
            return self.interwiki([word])
        if scheme in self.attachment_schemas:
            return self.attachment([word])

        if wikiutil.isPicture(word):
            word = wikiutil.mapURL(self.request, word)
            # Get image name http://here.com/dir/image.gif -> image
            name = word.split('/')[-1]
            name = ''.join(name.split('.')[:-1])
            return self.formatter.image(src=word, alt=name)
        else:
            return (self.formatter.url(1, word, css=scheme) +
                    self.formatter.text(word) +
                    self.formatter.url(0))


    def _wikiname_bracket_repl(self, word):
        """Handle special-char wikinames."""
        wikiname = word[2:-2]
        if wikiname:
            return self._word_repl(wikiname)
        else:
            return self.formatter.text(word)


    def _url_bracket_repl(self, word):
        """Handle bracketed URLs."""

        # Local extended link?
        if word[1] == ':':
            words = word[2:-1].split(':', 1)
            if len(words) == 1:
                words = words * 2
            words[0] = 'wiki:Self:%s' % words[0]
            return self.interwiki(words, pretty_url=1)
            #return self._word_repl(words[0], words[1])

        # Traditional split on space
        words = word[1:-1].split(None, 1)
        if len(words) == 1:
            words = words * 2

        if words[0][0] == '#':
            # anchor link
            return (self.formatter.url(1, words[0]) +
                    self.formatter.text(words[1]) +
                    self.formatter.url(0))

        scheme = words[0].split(":", 1)[0]
        if scheme == "wiki":
            return self.interwiki(words, pretty_url=1)
        if scheme in self.attachment_schemas:
            return self.attachment(words, pretty_url=1)

        if wikiutil.isPicture(words[1]) and re.match(self.url_rule, words[1]):
            return (self.formatter.url(1, words[0], css='external', do_escape=0) +
                    self.formatter.image(title=words[0], alt=words[0], src=words[1]) +
                    self.formatter.url(0))
        else:
            return (self.formatter.url(1, words[0], css=scheme, do_escape=0) +
                    self.formatter.text(words[1]) +
                    self.formatter.url(0))


    def _email_repl(self, word):
        """Handle email addresses (without a leading mailto:)."""
        return (self.formatter.url(1, "mailto:" + word, css='mailto') +
                self.formatter.text(word) +
                self.formatter.url(0))


    def _ent_repl(self, word):
        """Handle SGML entities."""
        return self.formatter.text(word)
        #return {'&': '&amp;',
        #        '<': '&lt;',
        #        '>': '&gt;'}[word]

    def _ent_numeric_repl(self, word):
        """Handle numeric (decimal and hexadecimal) SGML entities."""
        return self.formatter.rawHTML(word)

    def _ent_symbolic_repl(self, word):
        """Handle symbolic SGML entities."""
        return self.formatter.rawHTML(word)

    def _indent_repl(self, match):
        """Handle pure indentation (no - * 1. markup)."""
        result = []
        if not (self.in_li or self.in_dd):
            self._close_item(result)
            self.in_li = 1
            css_class = None
            if self.line_was_empty and not self.first_list_item:
                css_class = 'gap'
            result.append(self.formatter.listitem(1, css_class=css_class, style="list-style-type:none"))
        return ''.join(result)

    def _li_none_repl(self, match):
        """Handle type=none (" .") lists."""
        result = []
        self._close_item(result)
        self.in_li = 1
        css_class = None
        if self.line_was_empty and not self.first_list_item:
            css_class = 'gap'
        result.append(self.formatter.listitem(1, css_class=css_class, style="list-style-type:none"))
        return ''.join(result)

    def _li_repl(self, match):
        """Handle bullet (" *") lists."""
        result = []
        self._close_item(result)
        self.in_li = 1
        css_class = None
        if self.line_was_empty and not self.first_list_item:
            css_class = 'gap'
        result.append(self.formatter.listitem(1, css_class=css_class))
        return ''.join(result)

    def _ol_repl(self, match):
        """Handle numbered lists."""
        return self._li_repl(match)

    def _dl_repl(self, match):
        """Handle definition lists."""
        result = []
        self._close_item(result)
        self.in_dd = 1
        result.extend([
            self.formatter.definition_term(1),
            self.formatter.text(match[1:-3].lstrip(' ')),
            self.formatter.definition_term(0),
            self.formatter.definition_desc(1),
        ])
        return ''.join(result)


    def _indent_level(self):
        """Return current char-wise indent level."""
        return len(self.list_indents) and self.list_indents[-1]


    def _indent_to(self, new_level, list_type, numtype, numstart):
        """Close and open lists."""
        open = []   # don't make one out of these two statements!
        close = []

        if self._indent_level() != new_level and self.in_table:
            close.append(self.formatter.table(0))
            self.in_table = 0

        while self._indent_level() > new_level:
            self._close_item(close)
            if self.list_types[-1] == 'ol':
                tag = self.formatter.number_list(0)
            elif self.list_types[-1] == 'dl':
                tag = self.formatter.definition_list(0)
            else:
                tag = self.formatter.bullet_list(0)
            close.append(tag)

            del self.list_indents[-1]
            del self.list_types[-1]

            if self.list_types: # we are still in a list
                if self.list_types[-1] == 'dl':
                    self.in_dd = 1
                else:
                    self.in_li = 1

        # Open new list, if necessary
        if self._indent_level() < new_level:
            self.list_indents.append(new_level)
            self.list_types.append(list_type)

            if self.formatter.in_p:
                close.append(self.formatter.paragraph(0))

            if list_type == 'ol':
                tag = self.formatter.number_list(1, numtype, numstart)
            elif list_type == 'dl':
                tag = self.formatter.definition_list(1)
            else:
                tag = self.formatter.bullet_list(1)
            open.append(tag)

            self.first_list_item = 1
            self.in_li = 0
            self.in_dd = 0

        # If list level changes, close an open table
        if self.in_table and (open or close):
            close[0:0] = [self.formatter.table(0)]
            self.in_table = 0

        self.in_list = self.list_types != []
        return ''.join(close) + ''.join(open)


    def _undent(self):
        """Close all open lists."""
        result = []
        #result.append("<!-- _undent start -->\n")
        self._close_item(result)
        for type in self.list_types[::-1]:
            if type == 'ol':
                result.append(self.formatter.number_list(0))
            elif type == 'dl':
                result.append(self.formatter.definition_list(0))
            else:
                result.append(self.formatter.bullet_list(0))
        #result.append("<!-- _undent end -->\n")
        self.list_indents = []
        self.list_types = []
        return ''.join(result)


    def _tt_repl(self, word):
        """Handle inline code."""
        return self.formatter.code(1) + \
            self.formatter.text(word[3:-3]) + \
            self.formatter.code(0)


    def _tt_bt_repl(self, word):
        """Handle backticked inline code."""
        # if len(word) == 2: return "" // removed for FCK editor
        return self.formatter.code(1, css="backtick") + \
            self.formatter.text(word[1:-1]) + \
            self.formatter.code(0)


    def _getTableAttrs(self, attrdef):
        # skip "|" and initial "<"
        while attrdef and attrdef[0] == "|":
            attrdef = attrdef[1:]
        if not attrdef or attrdef[0] != "<":
            return {}, ''
        attrdef = attrdef[1:]

        # extension for special table markup
        def table_extension(key, parser, attrs, wiki_parser=self):
            """ returns: tuple (found_flag, msg)
                found_flag: whether we found something and were able to process it here
                  true for special stuff like 100% or - or #AABBCC
                  false for style xxx="yyy" attributes
                msg: "" or an error msg
            """
            _ = wiki_parser._
            found = False
            msg = ''
            if key[0] in "0123456789":
                token = parser.get_token()
                if token != '%':
                    wanted = '%'
                    msg = _('Expected "%(wanted)s" after "%(key)s", got "%(token)s"') % {
                        'wanted': wanted, 'key': key, 'token': token}
                else:
                    try:
                        dummy = int(key)
                    except ValueError:
                        msg = _('Expected an integer "%(key)s" before "%(token)s"') % {
                            'key': key, 'token': token}
                    else:
                        found = True
                        attrs['width'] = '"%s%%"' % key
            elif key == '-':
                arg = parser.get_token()
                try:
                    dummy = int(arg)
                except ValueError:
                    msg = _('Expected an integer "%(arg)s" after "%(key)s"') % {
                        'arg': arg, 'key': key}
                else:
                    found = True
                    attrs['colspan'] = '"%s"' % arg
            elif key == '|':
                arg = parser.get_token()
                try:
                    dummy = int(arg)
                except ValueError:
                    msg = _('Expected an integer "%(arg)s" after "%(key)s"') % {
                        'arg': arg, 'key': key}
                else:
                    found = True
                    attrs['rowspan'] = '"%s"' % arg
            elif key == '(':
                found = True
                attrs['align'] = '"left"'
            elif key == ':':
                found = True
                attrs['align'] = '"center"'
            elif key == ')':
                found = True
                attrs['align'] = '"right"'
            elif key == '^':
                found = True
                attrs['valign'] = '"top"'
            elif key == 'v':
                found = True
                attrs['valign'] = '"bottom"'
            elif key == '#':
                arg = parser.get_token()
                try:
                    if len(arg) != 6: raise ValueError
                    dummy = int(arg, 16)
                except ValueError:
                    msg = _('Expected a color value "%(arg)s" after "%(key)s"') % {
                        'arg': arg, 'key': key}
                else:
                    found = True
                    attrs['bgcolor'] = '"#%s"' % arg
            return found, self.formatter.rawHTML(msg)

        # scan attributes
        attr, msg = wikiutil.parseAttributes(self.request, attrdef, '>', table_extension)
        if msg:
            msg = '<strong class="highlight">%s</strong>' % msg
        #self.request.log("parseAttributes returned %r" % attr)
        return attr, msg

    def _tableZ_repl(self, word):
        """Handle table row end."""
        if self.in_table:
            result = ''
            # REMOVED: check for self.in_li, p should always close
            if self.formatter.in_p:
                result = self.formatter.paragraph(0)
            result += self.formatter.table_cell(0) + self.formatter.table_row(0)
            return result
        else:
            return self.formatter.text(word)

    def _table_repl(self, word):
        """Handle table cell separator."""
        if self.in_table:
            result = []
            # check for attributes
            attrs, attrerr = self._getTableAttrs(word)

            # start the table row?
            if self.table_rowstart:
                self.table_rowstart = 0
                result.append(self.formatter.table_row(1, attrs))
            else:
                # Close table cell, first closing open p
                # REMOVED check for self.in_li, paragraph should close always!
                if self.formatter.in_p:
                    result.append(self.formatter.paragraph(0))
                result.append(self.formatter.table_cell(0))

            # check for adjacent cell markers
            if word.count("|") > 2:
                if not attrs.has_key('align') and \
                   not (attrs.has_key('style') and 'text-align' in attrs['style'].lower()):
                    # add center alignment if we don't have some alignment already
                    attrs['align'] = '"center"'
                if not attrs.has_key('colspan'):
                    attrs['colspan'] = '"%d"' % (word.count("|")/2)

            # return the complete cell markup
            result.append(self.formatter.table_cell(1, attrs) + attrerr)
            result.append(self._line_anchordef())
            return ''.join(result)
        else:
            return self.formatter.text(word)


    def _heading_repl(self, word):
        """Handle section headings."""
        import sha

        h = word.strip()
        level = 1
        while h[level:level+1] == '=':
            level += 1
        depth = min(5,level)

        # this is needed for Included pages
        # TODO but it might still result in unpredictable results
        # when included the same page multiple times
        title_text = h[level:-level].strip()
        pntt = self.formatter.page.page_name + title_text
        self.titles.setdefault(pntt, 0)
        self.titles[pntt] += 1

        unique_id = ''
        if self.titles[pntt] > 1:
            unique_id = '-%d' % self.titles[pntt]
        result = self._closeP()
        result += self.formatter.heading(1, depth, id="head-"+sha.new(pntt.encode(config.charset)).hexdigest()+unique_id)

        return (result + self.formatter.text(title_text) +
                self.formatter.heading(0, depth))

    def _processor_repl(self, word):
        """Handle processed code displays."""
        if word[:3] == '{{{':
            word = word[3:]

        self.processor = None
        self.processor_name = None
        self.processor_is_parser = 0
        s_word = word.strip()
        if s_word == '#!':
            # empty bang paths lead to a normal code display
            # can be used to escape real, non-empty bang paths
            word = ''
            self.in_pre = 3
            return self._closeP() + self.formatter.preformatted(1)
        elif s_word[:2] == '#!':
            # First try to find a processor for this (will go away in 2.0)
            processor_name = s_word[2:].split()[0]
            self.setProcessor(processor_name)

        if self.processor:
            self.processor_name = processor_name
            self.in_pre = 2
            self.colorize_lines = [word]
            return ''
        elif s_word:
            self.in_pre = 3
            return self._closeP() + self.formatter.preformatted(1) + \
                   self.formatter.text(s_word + ' (-)')
        else:
            self.in_pre = 1
            return ''

    def _pre_repl(self, word):
        """Handle code displays."""
        word = word.strip()
        if word == '{{{' and not self.in_pre:
            self.in_pre = 3
            return self._closeP() + self.formatter.preformatted(self.in_pre)
        elif word == '}}}' and self.in_pre:
            self.in_pre = 0
            self.inhibit_p = 0
            return self.formatter.preformatted(self.in_pre)
        return self.formatter.text(word)


    def _smiley_repl(self, word):
        """Handle smileys."""
        return self.formatter.smiley(word)

    _smileyA_repl = _smiley_repl


    def _comment_repl(self, word):
        # if we are in a paragraph, we must close it so that normal text following
        # in the line below the comment will reopen a new paragraph.
        if self.formatter.in_p:
            self.formatter.paragraph(0)
        self.line_is_empty = 1 # markup following comment lines treats them as if they were empty
        return self.formatter.comment(word)

    def _closeP(self):
        if self.formatter.in_p:
            return self.formatter.paragraph(0)
        return ''

    def _macro_repl(self, word):
        """Handle macros ([[macroname]])."""
        macro_name = word[2:-2]
        self.inhibit_p = 0 # 1 fixes UserPreferences, 0 fixes paragraph formatting for macros

        # check for arguments
        args = None
        if macro_name.count("("):
            macro_name, args = macro_name.split('(', 1)
            args = args[:-1]

        # create macro instance
        if self.macro is None:
            self.macro = wikimacro.Macro(self)
        return self.formatter.macro(self.macro, macro_name, args)

    def scan(self, scan_re, line):
        """ Scans one line

        Append text before match, invoke replace() with match, and add text after match.
        """
        result = []
        lastpos = 0

        ###result.append(u'<span class="info">[scan: <tt>"%s"</tt>]</span>' % line)

        for match in scan_re.finditer(line):
            # Add text before the match
            if lastpos < match.start():

                ###result.append(u'<span class="info">[add text before match: <tt>"%s"</tt>]</span>' % line[lastpos:match.start()])

                if not (self.inhibit_p or self.in_pre or self.formatter.in_p):
                    result.append(self.formatter.paragraph(1, css_class="line862"))
                result.append(self.formatter.text(line[lastpos:match.start()]))

            # Replace match with markup
            if not (self.inhibit_p or self.in_pre or self.formatter.in_p or
                    self.in_table or self.in_list):
                result.append(self.formatter.paragraph(1, css_class="line867"))
            result.append(self.replace(match))
            lastpos = match.end()

        ###result.append('<span class="info">[no match, add rest: <tt>"%s"<tt>]</span>' % line[lastpos:])

        # Add paragraph with the remainder of the line
        if not (self.in_pre or self.in_li or self.in_dd or self.inhibit_p or
                self.formatter.in_p) and lastpos < len(line):
            result.append(self.formatter.paragraph(1, css_class="line874"))
        result.append(self.formatter.text(line[lastpos:]))
        return u''.join(result)

    def replace(self, match):
        """ Replace match using type name """
        result = []
        for type, hit in match.groupdict().items():
            if hit is not None and type != "hmarker":

                ###result.append(u'<span class="info">[replace: %s: "%s"]</span>' % (type, hit))
                if self.in_pre and type not in ['pre', 'ent']:
                    return self.formatter.text(hit)
                else:
                    # Open p for certain types
                    if not (self.inhibit_p or self.formatter.in_p
                            or self.in_pre or (type in self.no_new_p_before)):
                        result.append(self.formatter.paragraph(1, css_class="line891"))

                    # Get replace method and replece hit
                    replace = getattr(self, '_' + type + '_repl')
                    result.append(replace(hit))
                    return ''.join(result)
        else:
            # We should never get here
            import pprint
            raise Exception("Can't handle match " + `match`
                + "\n" + pprint.pformat(match.groupdict())
                + "\n" + pprint.pformat(match.groups()) )

        return ""

    def _line_anchordef(self):
        if self.line_anchors and not self.line_anchor_printed:
            self.line_anchor_printed = 1
            return self.formatter.line_anchordef(self.lineno)
        else:
            return ''

    def format(self, formatter):
        """ For each line, scan through looking for magic
            strings, outputting verbatim any intervening text.
        """
        self.formatter = formatter
        self.hilite_re = self.formatter.page.hilite_re

        # prepare regex patterns
        rules = self.formatting_rules.replace('\n', '|')
        if self.cfg.bang_meta:
            rules = ur'(?P<notword>!%(word_rule)s)|%(rules)s' % {
                'word_rule': self.word_rule,
                'rules': rules,
            }
        self.request.clock.start('compile_huge_and_ugly')
        scan_re = re.compile(rules, re.UNICODE)
        number_re = re.compile(self.ol_rule, re.UNICODE)
        term_re = re.compile(self.dl_rule, re.UNICODE)
        indent_re = re.compile("^\s*", re.UNICODE)
        eol_re = re.compile(r'\r?\n', re.UNICODE)
        self.request.clock.stop('compile_huge_and_ugly')

        # get text and replace TABs
        rawtext = self.raw.expandtabs()

        # go through the lines
        self.lineno = self.start_line
        self.lines = eol_re.split(rawtext)
        self.line_is_empty = 0

        self.in_processing_instructions = 1

        # Main loop
        for line in self.lines:
            self.lineno += 1
            self.line_anchor_printed = 0
            if not self.in_table:
                self.request.write(self._line_anchordef())
            self.table_rowstart = 1
            self.line_was_empty = self.line_is_empty
            self.line_is_empty = 0
            self.first_list_item = 0
            self.inhibit_p = 0

            # ignore processing instructions
            if self.in_processing_instructions:
                found = False
                for pi in ("##", "#format", "#refresh", "#redirect", "#deprecated",
                           "#pragma", "#form", "#acl", "#language"):
                    if line.lower().startswith(pi):
                        self.request.write(self.formatter.comment(line))
                        found = True
                        break
                if not found:
                    self.in_processing_instructions = 0
                else:
                    continue # do not parse this line
            if self.in_pre:
                # TODO: move this into function
                # still looking for processing instructions
                # TODO: use strings for pre state, not numbers
                if self.in_pre == 1:
                    self.processor = None
                    self.processor_is_parser = 0
                    processor_name = ''
                    if (line.strip()[:2] == "#!"):
                        processor_name = line.strip()[2:].split()[0]
                        self.setProcessor(processor_name)

                    if self.processor:
                        self.in_pre = 2
                        self.colorize_lines = [line]
                        self.processor_name = processor_name
                        continue
                    else:
                        self.request.write(self._closeP() +
                                           self.formatter.preformatted(1))
                        self.in_pre = 3
                if self.in_pre == 2:
                    # processing mode
                    endpos = line.find("}}}")
                    if endpos == -1:
                        self.colorize_lines.append(line)
                        continue
                    if line[:endpos]:
                        self.colorize_lines.append(line[:endpos])

                    # Close p before calling processor
                    # TODO: do we really need this?
                    self.request.write(self._closeP())
                    res = self.formatter.processor(self.processor_name,
                                                   self.colorize_lines,
                                                   self.processor_is_parser)
                    self.request.write(res)
                    del self.colorize_lines
                    self.in_pre = 0
                    self.processor = None

                    # send rest of line through regex machinery
                    line = line[endpos+3:]
                    if not line.strip(): # just in the case "}}} " when we only have blanks left...
                        continue
            else:
                # we don't have \n as whitespace any more
                # This is the space between lines we join to one paragraph
                line += ' '

                # Paragraph break on empty lines
                if not line.strip():
                    if self.in_table:
                        self.request.write(self.formatter.table(0))
                        self.request.write(self._line_anchordef())
                        self.in_table = 0
                    # CHANGE: removed check for not self.list_types
                    # p should close on every empty line
                    if self.formatter.in_p:
                        self.request.write(self.formatter.paragraph(0))
                    self.line_is_empty = 1
                    continue

                # Check indent level
                indent = indent_re.match(line)
                indlen = len(indent.group(0))
                indtype = "ul"
                numtype = None
                numstart = None
                if indlen:
                    match = number_re.match(line)
                    if match:
                        numtype, numstart = match.group(0).strip().split('.')
                        numtype = numtype[0]

                        if numstart and numstart[0] == "#":
                            numstart = int(numstart[1:])
                        else:
                            numstart = None

                        indtype = "ol"
                    else:
                        match = term_re.match(line)
                        if match:
                            indtype = "dl"

                # output proper indentation tags
                self.request.write(self._indent_to(indlen, indtype, numtype, numstart))

                # Table mode
                # TODO: move into function?
                if (not self.in_table and line[indlen:indlen + 2] == "||"
                    and line[-3:] == "|| " and len(line) >= 5 + indlen):
                    # Start table
                    if self.list_types and not self.in_li:
                        self.request.write(self.formatter.listitem(1, style="list-style-type:none"))
                        ## CHANGE: no automatic p on li
                        ##self.request.write(self.formatter.paragraph(1))
                        self.in_li = 1

                    # CHANGE: removed check for self.in_li
                    # paragraph should end before table, always!
                    if self.formatter.in_p:
                        self.request.write(self.formatter.paragraph(0))
                    attrs, attrerr = self._getTableAttrs(line[indlen+2:])
                    self.request.write(self.formatter.table(1, attrs) + attrerr)
                    self.in_table = True # self.lineno
                elif (self.in_table and not
                      # intra-table comments should not break a table
                      (line[:2] == "##" or
                       line[indlen:indlen + 2] == "||" and
                       line[-3:] == "|| " and
                       len(line) >= 5 + indlen)):

                    # Close table
                    self.request.write(self.formatter.table(0))
                    self.request.write(self._line_anchordef())
                    self.in_table = 0

            # Scan line, format and write
            formatted_line = self.scan(scan_re, line)
            self.request.write(formatted_line)

            if self.in_pre == 3:
                self.request.write(self.formatter.linebreak())

        # Close code displays, paragraphs, tables and open lists
        self.request.write(self._undent())
        if self.in_pre: self.request.write(self.formatter.preformatted(0))
        if self.formatter.in_p: self.request.write(self.formatter.paragraph(0))
        if self.in_table: self.request.write(self.formatter.table(0))

    # --------------------------------------------------------------------
    # Private helpers

    def setProcessor(self, name):
        """ Set processer to either processor or parser named 'name' """
        cfg = self.request.cfg
        try:
            self.processor = wikiutil.importPlugin(cfg, "processor", name,
                                                   "process")
            self.processor_is_parser = 0
        except wikiutil.PluginMissingError:
            try:
                self.processor = wikiutil.importPlugin(cfg, "parser", name,
                                                   "Parser")
                self.processor_is_parser = 1
            except wikiutil.PluginMissingError:
                self.processor = None


