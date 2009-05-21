#!/usr/bin/env python
import os, sys, subprocess, shutil, glob, optparse

#We remove the current directory from sys.path right away
#so that we import sage from the proper spot
try:
    sys.path.remove(os.path.realpath(os.getcwd()))
except:
    pass

from sage.misc.cachefunc import cached_method

##########################################
#                Options                 #
##########################################
SAGE_DOC = os.environ['SAGE_DOC']
LANGUAGES = ['en', 'fr']
SPHINXOPTS  = ""
PAPER       = ""

if PAPER == "a4":
    PAPEROPTS = "-D latex_paper_size=a4"
elif PAPER == "letter":
    PAPEROPTS = "-D latex_paper_size=letter"
else:
    PAPEROPTS = ""

#Note that this needs to have the doctrees dir
ALLSPHINXOPTS   = SPHINXOPTS + " " + PAPEROPTS + " . "


##########################################
#          Utility Functions             #
##########################################
def mkdir(path):
    """
    Makes the directory at path if it doesn't exist and returns the
    string path.

    EXAMPLES::

        sage: import os, sys; sys.path.append(os.environ['SAGE_DOC']+'/common/'); import builder
        sage: d = tmp_filename(); d
        '/.../tmp_...'
        sage: os.path.exists(d)
        False
        sage: dd = builder.mkdir(d)
        sage: d == dd
        True
        sage: os.path.exists(d)
        True
    """
    if not os.path.exists(path):
        os.makedirs(path)
    return path

def copytree(src, dst, symlinks=False, ignore_errors=False):
    """
    Recursively copy a directory tree using copy2().

    The destination directory must not already exist.
    If exception(s) occur, an Error is raised with a list of reasons.

    If the optional symlinks flag is true, symbolic links in the
    source tree result in symbolic links in the destination tree; if
    it is false, the contents of the files pointed to by symbolic
    links are copied.

    XXX Consider this example code rather than the ultimate tool.

    """
    names = os.listdir(src)
    mkdir(dst)
    errors = []
    for name in names:
        srcname = os.path.join(src, name)
        dstname = os.path.join(dst, name)
        try:
            if symlinks and os.path.islink(srcname):
                linkto = os.readlink(srcname)
                os.symlink(linkto, dstname)
            elif os.path.isdir(srcname):
                copytree(srcname, dstname, symlinks)
            else:
                shutil.copy2(srcname, dstname)
            # XXX What about devices, sockets etc.?
        except (IOError, os.error), why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except shutil.Error, err:
            errors.extend(err.args[0])
    try:
        shutil.copystat(src, dst)
    except OSError, why:
        errors.extend((src, dst, str(why)))
    if errors and not ignore_errors:
        raise shutil.Error, errors



##########################################
#             Builders                   #
##########################################
def builder_helper(type):
    """
    Returns a function which builds the documentation for
    output type type.
    """
    def f(self):
        output_dir = self._output_dir(type)
        os.chdir(self.dir)

        build_command = 'sphinx-build'
        build_command += ' -b %s -d %s %s %s'%(type, self._doctrees_dir(),
                                               ALLSPHINXOPTS, output_dir)
        print build_command
        subprocess.call(build_command, shell=True)

        print "Build finished.  The built documents can be found in %s"%output_dir

    f.is_output_format = True
    return f


class DocBuilder(object):
    def __init__(self, name, lang='en'):
        """
        INPUT:

        - ``name`` - the name of a subdirectory in SAGE_DOC, such as
          'tutorial' or 'bordeaux_2008'

        - ``lang`` - (default "en") the language of the document.
        """
        if '/' in name:
            lang, name = name.split(os.path.sep)
        self.name = name
        self.lang = lang
        self.dir = os.path.join(SAGE_DOC, lang, name)

        #Make sure the .static and .templates directories are there
        mkdir(os.path.join(self.dir, "static"))
        mkdir(os.path.join(self.dir, "templates"))

    def _output_dir(self, type):
        """
        Returns the directory where the output of type type is stored.
        If the directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: import os, sys; sys.path.append(os.environ['SAGE_DOC']+'/common/'); import builder
            sage: b = builder.DocBuilder('tutorial')
            sage: b._output_dir('html')
            '.../devel/sage/doc/output/html/en/tutorial'
        """
        return mkdir(os.path.join(SAGE_DOC, "output", type, self.lang, self.name))

    def _doctrees_dir(self):
        """
        Returns the directory where the doctrees are stored.  If the
        directory does not exist, then it will automatically be
        created.

        EXAMPLES::

            sage: import os, sys; sys.path.append(os.environ['SAGE_DOC']+'/common/'); import builder
            sage: b = builder.DocBuilder('tutorial')
            sage: b._doctrees_dir()
            '.../devel/sage/doc/output/doctrees/en/tutorial'
        """
        return mkdir(os.path.join(SAGE_DOC, "output", 'doctrees', self.lang, self.name))

    def _output_formats(self):
        """
        Returns a list of the possible output formats.

        EXAMPLES::

            sage: import os, sys; sys.path.append(os.environ['SAGE_DOC']+'/common/'); import builder
            sage: b = builder.DocBuilder('tutorial')
            sage: b._output_formats()
            ['changes', 'html', 'htmlhelp', 'json', 'latex', 'linkcheck', 'pickle', 'web']

        """
        #Go through all the attributes of self and check to
        #see which ones have an 'is_output_format' attribute.  These
        #are the ones created with builder_helper.
        output_formats = []
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                output_formats.append(attr)
        output_formats.sort()
        return output_formats

    def pdf(self):
        """
        Builds the PDF files for this document.  This is done by first
        (re)-building the LaTeX output, going into that LaTeX
        directory, and running 'make all-pdf' there.

        EXAMPLES::

            sage: import os, sys; sys.path.append(os.environ['SAGE_DOC']+'/common/'); import builder
            sage: b = builder.DocBuilder('tutorial')
            sage: b.pdf() #not tested
        """
        self.latex()
        os.chdir(self._output_dir('latex'))
        subprocess.call('make all-pdf', shell=True)

        pdf_dir = self._output_dir('pdf')
        for pdf_file in glob.glob('*.pdf'):
            shutil.move(pdf_file, os.path.join(pdf_dir, pdf_file))

        print "Build finished.  The built documents can be found in %s"%pdf_dir

    def clean(self, *args):
        """
        """
        import shutil
        shutil.rmtree(self._doctrees_dir())
        output_formats = list(args) if args else self._output_formats()
        for format in output_formats:
            shutil.rmtree(self._output_dir(format), ignore_errors=True)

    html = builder_helper('html')
    pickle = builder_helper('pickle')
    web = pickle
    json = builder_helper('json')
    htmlhelp = builder_helper('htmlhelp')
    latex = builder_helper('latex')
    changes = builder_helper('changes')
    linkcheck = builder_helper('linkcheck')

class AllBuilder(object):
    """
    A class used to build all of the documentation.
    """
    def __getattr__(self, attr):
        """
        For any attributes not explicitly defined, we just go through
        all of the documents and call their attr.  For example,
        'AllBuilder().json()' will go through all of the documents
        and call the json() method on their builders.
        """
        from functools import partial
        return partial(self._wrapper, attr)

    def _wrapper(self, name, *args, **kwds):
        """
        This is the function which goes through all of the documents
        and does the actual building.
        """
        for document in self.get_all_documents():
            getattr(get_builder(document), name)(*args, **kwds)

    def get_all_documents(self):
        """
        Returns a list of all of the documents. A document is a directory within one of
        the language subdirectories of SAGE_DOC specified by the global LANGUAGES
        variable.

        EXAMPLES::

            sage: import os, sys; sys.path.append(os.environ['SAGE_DOC']+'/common/'); import builder
            sage: documents = builder.AllBuilder().get_all_documents()
            sage: 'en/tutorial' in documents
            True
        """
        documents = []
        for lang in LANGUAGES:
            for document in os.listdir(os.path.join(SAGE_DOC, lang)):
                documents.append(os.path.join(lang, document))
        return documents

class WebsiteBuilder(DocBuilder):
    def html(self):
        """
        After we've finished building the website index page, we copy
        everything one directory up.
        """
        DocBuilder.html(self)
        html_output_dir = self._output_dir('html')
        copytree(html_output_dir,
                 os.path.realpath(os.path.join(html_output_dir, '..')),
                 ignore_errors=False)

    def clean(self):
        """
        When we clean the output for the website index, we need to
        remove all of the HTML that were placed in the parent
        directory.
        """
        html_output_dir = self._output_dir('html')
        parent_dir = os.path.realpath(os.path.join(html_output_dir, '..'))
        for filename in os.listdir(html_output_dir):
            parent_filename = os.path.join(parent_dir, filename)
            if not os.path.exists(parent_filename):
                continue
            if os.path.isdir(parent_filename):
                shutil.rmtree(parent_filename, ignore_errors=True)
            else:
                os.unlink(parent_filename)

        DocBuilder.clean(self)

class ReferenceBuilder(DocBuilder):
    """
    This the class used to build the reference manual.  It is
    resposible for making sure the auto generated ReST files for the
    Sage library are up to date.

    When building any output, we must first go through and check
    to see if we need to update any of the autogenerated ReST
    files.  There are two cases where this would happen:

    1. A new module gets added to one of the toctrees.

    2. The actual module gets updated and possibly contains a new
       title.
    """
    def __init__(self, *args, **kwds):
        DocBuilder.__init__(self, *args, **kwds)
        self._wrap_builder_helpers()

    def _wrap_builder_helpers(self):
        from functools import partial, update_wrapper
        for attr in dir(self):
            if hasattr(getattr(self, attr), 'is_output_format'):
                f = partial(self._wrapper, attr)
                f.is_output_format = True
                update_wrapper(f, getattr(self, attr))
                setattr(self, attr, f)

    def _wrapper(self, build_type, *args, **kwds):
        """
        This is the wrapper around the builder_helper methods that
        goes through and makes sure things are up to date.
        """
        #Write the .rst files for newly included modules
        for module_name in self.get_newly_included_modules(save=True):
            self.write_auto_rest_file(module_name)

        #Update the .rst files for modified Python modules
        for module_name in self.get_modified_modules():
            self.write_auto_rest_file(module_name.replace(os.path.sep, '.'))

        #Copy over the custom .rst files from _sage
        _sage = os.path.join(self.dir, '_sage')
        if os.path.exists(_sage):
            copytree(_sage, os.path.join(self.dir, 'sage'))

        getattr(DocBuilder, build_type)(self, *args, **kwds)

    def cache_filename(self):
        """
        Returns the filename where the pickle of the dictionary of
        already generated .rst files is stored.
        """
        return os.path.join(self._doctrees_dir(), 'reference.pickle')

    @cached_method
    def get_cache(self):
        """
        Retreive the cache of already generated .rst files.  If it
        doesn't exist, then we just return an empty dictionary.
        """
        filename = self.cache_filename()
        if not os.path.exists(filename):
            return {}

        import cPickle
        file = open(self.cache_filename(), 'rb')
        cache = cPickle.load(file)
        file.close()
        return cache


    def save_cache(self):
        """
        Save the cache of already generated .rst files.
        """
        import cPickle
        file = open(self.cache_filename(), 'wb')
        cPickle.dump(self.get_cache(), file)
        file.close()

    def get_sphinx_environment(self):
        """
        Returns the Sphinx environment for this project.
        """
        from sphinx.environment import BuildEnvironment
        class Foo(object):
            pass
        config = Foo()
        config.values = []

        env_pickle = os.path.join(self._doctrees_dir(), 'environment.pickle')
        try:
            return BuildEnvironment.frompickle(config, env_pickle)
        except IOError:
            pass

    def get_modified_modules(self):
        """
        Returns an iterator for all the modules that have been modified
        since the docuementation was last built.
        """
        env = self.get_sphinx_environment()
        if env is None:
            return
        added, changed, removed = env.get_outdated_files(False)
        for name in changed:
            if name.startswith('sage'):
                yield name

    def print_modified_modules(self):
        """
        Prints a list of all the modules that have been modified since
        the documentation was last built.
        """
        for module_name in self.get_modified_modules():
            print module_name

    def get_all_rst_files(self, exclude_sage=True):
        """
        Returns an iterator for all rst files which are not
        autogenerated.
        """
        for directory, subdirs, files in os.walk(self.dir):
            if exclude_sage and directory.startswith(os.path.join(self.dir, 'sage')):
                continue
            for filename in files:
                if not filename.endswith('.rst'):
                    continue
                yield os.path.join(directory, filename)

    def get_all_included_modules(self):
        """
        Returns an iterator for all modules which are included in the
        reference manual.
        """
        for filename in self.get_all_rst_files():
            for module in self.get_modules(filename):
                yield module

    def get_newly_included_modules(self, save=False):
        """
        Returns an iterator for all modules that appear in the
        toctrees that don't appear in the cache.
        """
        cache = self.get_cache()
        new_modules = []
        for module in self.get_all_included_modules():
            if module not in cache:
                cache[module] = True
                yield module
        if save:
            self.save_cache()

    def print_newly_included_modules(self):
        """
        Prints all of the modules that appear in the toctrees that
        don't appear in the cache.
        """
        for module_name in self.get_newly_included_modules():
            print module_name

    def get_modules(self, filename):
        """
        Given a filename for a .rst file, return an iterator for
        all of the autogenerated rst files that it includes.
        """
        #Create the regular expression used to detect an autogenerated file
        import re
        auto_re = re.compile('^\s*(..\/)*(sage\/[\w\/]*)\s*$')

        #Read the lines
        f = open(filename)
        lines = f.readlines()
        f.close()

        for line in lines:
            match = auto_re.match(line)
            if match:
                yield match.group(2).replace(os.path.sep, '.')

    def get_module_docstring_title(self, module_name):
        """
        Returns the title of the module from its docstring.
        """
        #Try to import the module
        try:
            import sage.all
            __import__(module_name)
        except ImportError, err:
            print "Warning: could not import %s"%module_name
            print err
            return "UNABLE TO IMPORT MODULE"
        module = sys.modules[module_name]

        #Get the docstring
        doc = module.__doc__
        if doc is None:
            doc = module.doc if hasattr(module, 'doc') else ""

        #Extract the title
        i = doc.find('\n')
        if i != -1:
            return doc[i+1:].lstrip().splitlines()[0]
        else:
            return doc

    def write_auto_rest_file(self, module_name):
        """
        Writes the autogenerated .rst file for module_name.
        """
        if not module_name.startswith('sage'):
            return
        filename = self.dir + os.path.sep + module_name.replace('.',os.path.sep) + '.rst'
        mkdir(os.path.dirname(filename))

        outfile = open(filename, 'w')

        title = self.get_module_docstring_title(module_name)

        if title == '':
            print "WARNING: Missing title for", module_name
            title = "MISSING TITLE"

        outfile.write(title + '\n')
        outfile.write('='*len(title) + "\n\n")
        outfile.write('.. This file has been autogenerated.\n\n')
        automodule = '.. automodule:: %s\n   :members:\n   :undoc-members:\n\n'
        outfile.write(automodule%module_name)

        outfile.close()

    def clean_auto(self):
        """
        Remove the cache file for the autogenerated files as well as
        the files themselves.
        """
        if os.path.exists(self.cache_filename()):
            os.unlink(self.cache_filename())

        import shutil
        shutil.rmtree(self.dir + '/sage')

    def get_unincluded_modules(self):
        """
        Returns an iterator for all the modules in the Sage library
        which are not included in the reference manual.
        """
        #Make a dictionary of the included modules
        included_modules = {}
        for module_name in self.get_all_included_modules():
            included_modules[module_name] = True

        base_path = os.path.join(os.environ['SAGE_ROOT'], 'devel', 'sage', 'sage')
        for directory, subdirs, files in os.walk(base_path):
            for filename in files:
                if not (filename.endswith('.py') or
                        filename.endswith('.pyx')):
                    continue

                path = os.path.join(directory, filename)

                #Create the module name
                module_name = path[len(base_path):].replace(os.path.sep, '.')
                module_name = 'sage' + module_name
                module_name = module_name[:-4] if module_name.endswith('pyx') else module_name[:-3]

                #Exclude some ones  -- we don't want init the manual
                if module_name.endswith('__init__') or module_name.endswith('all'):
                    continue

                if module_name not in included_modules:
                    yield module_name

    def print_unincluded_modules(self):
        """
        Prints all of the modules which are not included in the Sage
        reference manual.
        """
        for module_name in self.get_unincluded_modules():
            print module_name

    def print_included_modules(self):
        """
        Prints all of the modules that are included in the Sage reference
        manual.
        """
        for module_name in self.get_all_included_modules():
            print module_name



def get_builder(name):
    """
    Returns a either a AllBuilder or DocBuilder object depending
    on whether ``name`` is 'all' or not.  These are the objects
    which do all the real work in building the documentation.
    """
    if name == 'all':
        return AllBuilder()
    elif name.endswith('reference'):
        return ReferenceBuilder(name)
    elif name.endswith('website'):
        return WebsiteBuilder(name)
    else:
        return DocBuilder(name)


def help_message():
    """
    Returns the help message.
    """
    all_b = AllBuilder()
    docs = all_b.get_all_documents()
    docs = [(d[3:] if d[0:3] == 'en/' else d) for d in docs]
    tut_b = DocBuilder('en/tutorial')
    formats = tut_b._output_formats()
    formats.remove('html')
    help = "Usage: sage -docbuild {document} {format}\n"
    help += "Where {document} is one of:\n    "
    help += "\n    ".join(docs)
    help += "\nor 'all' for all documents, and {format} is one of:\n    "
    help += 'html, pdf, ' + ', '.join(formats)
    help += "\n"
    help += "When building the reference manual, there are several additional\n"
    help += "values for {format}:\n"
    help += "    print_modified_modules: list modules modified since the documentation\n"
    help += "         was last built\n"
    help += "    print_newly_included_modules: list modules added to the\n"
    help += "         documentation since it was last built\n"
    help += "    print_unincluded_modules: list modules not included in the documentation\n"
    help += "    print_included_modules: list modules included in the documentation\n"
    print help


parser = optparse.OptionParser(usage="usage: sage -docbuild [options] name type")
parser.add_option("--jsmath", action="store_true",
                  help="render math using jsMath")
parser.print_help = help_message

if __name__ == '__main__':
    options, args = parser.parse_args()

    if options.jsmath:
        os.environ['SAGE_DOC_JSMATH'] = "True"

    #Get the name of the document we are trying to build
    try:
        name, type = args
    except ValueError:
        print "You must specify the document name and the output format"
        sys.exit(0)

    #Make sure common/static exists
    mkdir(os.path.join(SAGE_DOC, 'common', 'static'))

    #Get the builder and build
    getattr(get_builder(name), type)()
