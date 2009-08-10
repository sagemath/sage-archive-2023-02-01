"""
Sage Notebook Introspection

"""

"""
TODO: - add support for grabbing source code from Pyrex functions
(even if not perfect is better than nothing). - PNG or MathML
output format for docstring
"""

###########################################################################
#       Copyright (C) 2006 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
###########################################################################


def introspect(S, query, format='html'):
    """
    Return introspection from a given query string.

    INPUT:


    -  ``S`` - a Sage0 object, i.e., an interface to a
       running instance of Python with the Sage libraries loaded

    -  ``query`` - a string: - if has no '?' then return
       completion list - if begins or ends in one '?' return docstring -
       if begins or ends in '??' return source code

    -  ``format`` - (string) 'html', 'png', 'none' (only
       html is implemented right now!)
    """
    if format != 'html':
        raise NotImplementedError
    query = query.replace('\n','').strip()
    if query[:9] == '?__last__':
        return get_docstring_last(S, int(query[9:])/15)
    if len(query) > 1:
        if query[:2] == '??':
            return get_source_code(S, query[2:])
        elif query[-2:] == '??':
            return get_source_code(S, query[:-2])
    if len(query) > 0:
        if query[0] == '?':
            return get_docstring(S, query[1:])
        elif query[-1] == '?':
            return get_docstring(S, query[:-1])
    return get_completions(S, query)



def _get_docstring(S, query):
    cmd = '_support_.docstring("%s", globals())'%query
    z = S.eval(cmd)
    z = z.replace('\\n','\n').replace('\\t','        ')[1:-1]
    z = word_wrap(z, ncols=numcols)
    return z

def _get_source_code(S, query):
    cmd = '"".join(_support_.source_code("%s", globals()))'%query
    z = S.eval(cmd)
    z = z.replace('\\n','\n').replace("\\","").replace('\\t','        ')[1:-1]
    return z

def _get_completions(S, query):
    cmd = '"<br>".join(_support_.completions("%s", globals()))'%query
    z = S.eval(cmd)
    _last_ = z
    return z[1:-1]

