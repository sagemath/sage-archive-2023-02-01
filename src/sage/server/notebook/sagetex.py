import os

from notebook_object import notebook

def sagetex(filename, gen=True, **kwds):
    """
    Turn a latex document into an interactive notebook server.

    THIS IS ONLY A PROOF-of-CONCEPT.

    EXAMPLES::

        sage: sagetex('foo.tex')        # not tested
        [pops up web browser with live version of foo.tex.]
    """
    if not os.path.exists(filename):
        raise IOError, "No such file: '%s'"%filename

    if not filename.endswith('.tex'):
        raise ValueError, "File must be a latex file (end in .tex): '%s'"%filename
    if not '\\document' in open(filename).read():
        raise ValueError, "File must be a latex file (contain \\document...): '%s'"%filename

    if gen:
        os.system('latex2html -no_images %s'%filename)

    base = os.path.splitext(os.path.split(filename)[1])[0]
    absp = os.path.abspath(filename)
    path = os.path.splitext(absp)[0]

    for F in os.listdir(base):
        fn = '%s/%s'%(base, F)
        if not fn.endswith('.html'):
            continue
        r = open(fn).read()
        r = r.replace('<PRE>', '<div class="verbatim"><pre>')
        r = r.replace('</PRE>', '</pre></div>')
        open(fn,'w').write(r)

    notebook(sagetex_path=path, start_path="/sagetex/index.html", **kwds)




