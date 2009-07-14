#!/usr/bin/env python
import glob, logging, optparse, os, shutil, subprocess, sys, textwrap

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
ALLSPHINXOPTS   = SPHINXOPTS + " " + PAPEROPTS + " "


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
        except (IOError, os.error) as why:
            errors.append((srcname, dstname, str(why)))
        # catch the Error from the recursive copytree so that we can
        # continue with other files
        except shutil.Error as err:
            errors.extend(err.args[0])
    try:
        shutil.copystat(src, dst)
    except OSError as why:
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
        build_command += ' -b %s -d %s %s %s %s'%(type, self._doctrees_dir(),
                                                  ALLSPHINXOPTS, self.dir,
                                                  output_dir)
        logger.warning(build_command)
        subprocess.call(build_command, shell=True)

        logger.warning("Build finished.  The built documents can be found in %s", output_dir)

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

        logger.warning("Build finished.  The built documents can be found in %s", pdf_dir)

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
        # After "sage -clone", refresh the .rst file mtimes in
        # environment.pickle.
        if update_mtimes:
            logger.info("Checking for .rst file mtimes to update...")
            self.update_mtimes()

        #Update the .rst files for modified Python modules
        logger.info("Updating .rst files with modified modules...")
        for module_name in self.get_modified_modules():
            self.write_auto_rest_file(module_name.replace(os.path.sep, '.'))

        #Write the .rst files for newly included modules
        logger.info("Writing .rst files for newly-included modules...")
        for module_name in self.get_newly_included_modules(save=True):
            self.write_auto_rest_file(module_name)

        #Copy over the custom .rst files from _sage
        _sage = os.path.join(self.dir, '_sage')
        if os.path.exists(_sage):
            logger.info("Copying over custom .rst files from %s ...", _sage)
            copytree(_sage, os.path.join(self.dir, 'sage'))

        getattr(DocBuilder, build_type)(self, *args, **kwds)

    def cache_filename(self):
        """
        Returns the filename where the pickle of the dictionary of
        already generated ReST files is stored.
        """
        return os.path.join(self._doctrees_dir(), 'reference.pickle')

    @cached_method
    def get_cache(self):
        """
        Retreive the cache of already generated ReST files.  If it
        doesn't exist, then we just return an empty dictionary.
        """
        filename = self.cache_filename()
        if not os.path.exists(filename):
            return {}
        import cPickle
        file = open(self.cache_filename(), 'rb')
        cache = cPickle.load(file)
        file.close()
        logger.debug("Loaded .rst file cache: %s", filename)
        return cache

    def save_cache(self):
        """
        Save the cache of already generated ReST files.
        """
        import cPickle
        file = open(self.cache_filename(), 'wb')
        cPickle.dump(self.get_cache(), file)
        file.close()
        logger.debug("Saved .rst file cache: %s", self.cache_filename())

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
            env = BuildEnvironment.frompickle(config, env_pickle)
            logger.debug("Opened Sphinx environment: %s", env_pickle)
            return env
        except IOError as err:
            logger.debug("Failed to open Sphinx environment: %s", err)
            pass

    def update_mtimes(self):
        """
        Updates the modification times for ReST files in the Sphinx
        environment for this project.
        """
        env = self.get_sphinx_environment()
        if env is not None:
            import time
            for doc in env.all_docs:
                env.all_docs[doc] = time.time()
            logger.info("Updated %d .rst file mtimes", len(env.all_docs))
            # This is the only place we need to save (as opposed to
            # load) Sphinx's pickle, so we do it right here.
            env_pickle = os.path.join(self._doctrees_dir(),
                                      'environment.pickle')
            env.topickle(env_pickle)
            logger.debug("Saved Sphinx environment: %s", env_pickle)

    def get_modified_modules(self):
        """
        Returns an iterator for all the modules that have been modified
        since the docuementation was last built.
        """
        env = self.get_sphinx_environment()
        if env is None:
            logger.debug("Stopped check for modified modules.")
            return
        try:
            added, changed, removed = env.get_outdated_files(False)
            logger.info("Sphinx found %d modified modules", len(changed))
        except OSError as err:
            logger.debug("Sphinx failed to determine modified modules: %s", err)
            self.clean_auto()
            return
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
        new_modules = 0
        for module in self.get_all_included_modules():
            if module not in cache:
                cache[module] = True
                new_modules += 1
                yield module
        logger.info("Found %d newly included modules", new_modules)
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
        Given a filename for a ReST file, return an iterator for
        all of the autogenerated ReST files that it includes.
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
        except ImportError as err:
            logger.error("Warning: Could not import %s %s", module_name, err)
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
        Writes the autogenerated ReST file for module_name.
        """
        if not module_name.startswith('sage'):
            return
        filename = self.dir + os.path.sep + module_name.replace('.',os.path.sep) + '.rst'
        mkdir(os.path.dirname(filename))

        outfile = open(filename, 'w')

        title = self.get_module_docstring_title(module_name)

        if title == '':
            logger.error("Warning: Missing title for %s", module_name)
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
            logger.debug("Deleted .rst cache file: %s", self.cache_filename())

        import shutil
        try:
            shutil.rmtree(os.path.join(self.dir, 'sage'))
            logger.debug("Deleted auto-generated .rst files in: %s",
                         os.path.join(self.dir, 'sage'))
        except OSError:
            pass

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


def format_columns(lst, align='<', cols=None, indent=4, pad=3, width=80):
    """
    Utility function that formats a list as a simple table and returns
    a Unicode string representation.  The number of columns is
    computed from the other options, unless it's passed as a keyword
    argument.  For help on Python's string formatter, see

    http://docs.python.org/library/string.html#format-string-syntax
    """
    # Can we generalize this (efficiently) to other / multiple inputs
    # and generators?
    size = max(map(len, lst)) + pad
    if cols is None:
        import math
        cols = math.trunc((width - indent) / size)
    s = " " * indent
    for i in xrange(len(lst)):
        if i != 0 and i % cols == 0:
            s += "\n" + " " * indent
        s += "{0:{1}{2}}".format(lst[i], align, size)
    s += "\n"
    return unicode(s)

def help_usage(s=u"", compact=False):
    """
    Appends and returns a brief usage message for the Sage
    documentation builder.  If 'compact' is False, the function adds a
    final newline character.
    """
    s += "sage -docbuild [OPTIONS] DOCUMENT (FORMAT | COMMAND)"
    if not compact:
        s += "\n"
    return s

def help_description(s=u"", compact=False):
    """
    Appends and returns a brief description of the Sage documentation
    builder.  If 'compact' is False, the function adds a final newline
    character.
    """
    s += "Build or return information about Sage documentation.\n"
    s += "    DOCUMENT    name of the document to build\n"
    s += "    FORMAT      document output format\n"
    s += "    COMMAND     document-specific command\n"
    s += "A DOCUMENT and either a FORMAT or a COMMAND are required,\n"
    s += "unless a list of one or more of these is requested."
    if not compact:
        s += "\n"
    return s

def help_examples(s=u""):
    """
    Appends and returns some usage examples for the Sage documentation
    builder.
    """
    s += "Examples:\n"
    s += "    sage -docbuild -FDC all\n"
    s += "    sage -docbuild constructions pdf\n"
    s += "    sage -docbuild reference html -jv3\n"
    s += "    sage -docbuild --jsmath tutorial html\n"
    s += "    sage -docbuild reference print_unincluded_modules\n"
    s += "    sage -docbuild developer -j html --sphinx-opts -q,-aE --verbose 2"
    return s

def get_documents():
    """
    Returns a list of document names the Sage documentation builder
    will accept as command-line arguments.
    """
    all_b = AllBuilder()
    docs = all_b.get_all_documents()
    docs = [(d[3:] if d[0:3] == 'en/' else d) for d in docs]
    return docs

def help_documents(s=u""):
    """
    Appends and returns a tabular list of documents, including a
    shortcut 'all' for all documents, available to the Sage
    documentation builder.
    """
    s += "DOCUMENTs:\n"
    s += format_columns(get_documents() + ['all  (!)'])
    s += "(!) Builds everything.\n"
    return s

def get_formats():
    """
    Returns a list of output formats the Sage documentation builder
    will accept on the command-line.
    """
    tut_b = DocBuilder('en/tutorial')
    formats = tut_b._output_formats()
    formats.remove('html')
    return ['html', 'pdf'] + formats

def help_formats(s=u""):
    """
    Appends and returns a tabular list of output formats available to
    the Sage documentation builder.
    """
    s += "FORMATs:\n"
    s += format_columns(get_formats())
    return s

def help_commands(name='all', s=u""):
    """
    Appends and returns a tabular list of commands, if any, the Sage
    documentation builder can run on the indicated document.  The
    default is to list all commands for all documents.
    """
    # To do: Generate the lists dynamically, using class attributes,
    # as with the Builders above.
    command_dict = { 'reference' : [
        'print_included_modules',   'print_modified_modules       (*)',
        'print_unincluded_modules', 'print_newly_included_modules (*)',
        ] }
    for doc in command_dict:
        if name == 'all' or doc == name:
            s += "COMMANDs for the DOCUMENT '" + doc + "':\n"
            s += format_columns(command_dict[doc])
            s += "(*) Since the last build.\n"
    return s

def help_message_long(option, opt_str, value, parser):
    """
    Prints an extended help message for the Sage documentation builder
    and exits.
    """
    help_funcs = [ help_usage, help_description, help_documents,
                   help_formats, help_commands, parser.format_option_help,
                   help_examples ]
    for f in help_funcs:
        print f()
    sys.exit(0)

def help_message_short(option=None, opt_str=None, value=None, parser=None,
                       error=False):
    """
    Prints a help message for the Sage documentation builder.  The
    message includes command-line usage and a list of options.  The
    message is printed only on the first call.  If error is True
    during this call, the message is printed only if the user hasn't
    requested a list (e.g., documents, formats, commands).
    """
    if not hasattr(parser.values, 'printed_help'):
        if error == True:
            if not hasattr(parser.values, 'printed_list'):
                parser.print_help()
        else:
            parser.print_help()
        setattr(parser.values, 'printed_help', 1)

def help_wrapper(option, opt_str, value, parser):
    """
    A helper wrapper for command-line options to the Sage
    documentation builder that print lists, such as document names,
    formats, and document-specific commands.
    """
    if option.dest == 'commands':
        print help_commands(value),
    if option.dest == 'documents':
        print help_documents(),
    if option.dest == 'formats':
        print help_formats(),
    setattr(parser.values, 'printed_list', 1)


class IndentedHelpFormatter2(optparse.IndentedHelpFormatter, object):
    """
    Custom help formatter class for optparse's OptionParser.
    """
    def format_description(self, description):
        """
        Returns a formatted description, preserving any original
        explicit new line characters.
        """
        if description:
            lines_in = description.split('\n')
            lines_out = [self._format_text(line) for line in lines_in]
            return "\n".join(lines_out) + "\n"
        else:
            return ""

    def format_heading(self, heading):
        """
        Returns a formatted heading using the superclass' formatter.
        If the heading is 'options', up to case, the function converts
        it to ALL CAPS.  This allows us to match the heading 'OPTIONS' with
        the same token in the builder's usage message.
        """
        if heading.lower() == 'options':
            heading = "OPTIONS"
        return super(IndentedHelpFormatter2, self).format_heading(heading)

def setup_parser():
    """
    Sets up and returns a command-line OptionParser instance for the
    Sage documentation builder.
    """
    # Documentation: http://docs.python.org/library/optparse.html
    parser = optparse.OptionParser(add_help_option=False,
                                   usage=help_usage(compact=True),
                                   formatter=IndentedHelpFormatter2(),
                                   description=help_description(compact=True))

    # Standard options. Note: We use explicit option.dest names
    # to avoid ambiguity.
    standard = optparse.OptionGroup(parser, "Standard")
    standard.add_option("-h", "--help",
                        action="callback", callback=help_message_short,
                        help="show a help message and exit")
    standard.add_option("-H", "--help-all",
                        action="callback", callback=help_message_long,
                        help="show an extended help message and exit")
    standard.add_option("-D", "--documents", dest="documents",
                        action="callback", callback=help_wrapper,
                        help="list all available DOCUMENTs")
    standard.add_option("-F", "--formats", dest="formats",
                        action="callback", callback=help_wrapper,
                        help="list all output FORMATs")
    standard.add_option("-C", "--commands", dest="commands",
                        type="string", metavar="DOC",
                        action="callback", callback=help_wrapper,
                        help="list all COMMANDs for DOCUMENT DOC; use 'all' to list all")
    standard.add_option("-j", "--jsmath", dest="jsmath",
                        action="store_true",
                        help="render math using jsMath; FORMATs: html, json, pickle, web")

    standard.add_option("-N", "--no-colors", dest="color", default=True,
                        action="store_false",
                        help="do not color output; does not affect children")
    standard.add_option("-q", "--quiet", dest="verbose",
                        action="store_const", const=0,
                        help="work quietly; same as --verbose=0")
    standard.add_option("-v", "--verbose", dest="verbose",
                        type="int", default=1, metavar="LEVEL",
                        action="store",
                        help="report progress at LEVEL=0 (quiet), 1 (normal), 2 (info), or 3 (debug); does not affect children")
    parser.add_option_group(standard)

    # Advanced options.
    advanced = optparse.OptionGroup(parser, "Advanced",
                                    "Use these options with care.")
    advanced.add_option("-S", "--sphinx-opts", dest="sphinx_opts",
                        type="string", metavar="OPTS",
                        action="store",
                        help="pass comma-separated OPTS to sphinx-build")
    advanced.add_option("-U", "--update-mtimes", dest="update_mtimes",
                        action="store_true",
                        help="before building reference manual, update modification times for auto-generated ReST files")
    parser.add_option_group(advanced)

    return parser

def setup_logger(verbose=1, color=True):
    """
    Sets up and returns a Python Logger instance for the Sage
    documentation builder.  The optional argument sets logger's level
    and message format.
    """
    # Set up colors. Adapted from sphinx.cmdline.
    import sphinx.util.console as c
    if not color or not sys.stdout.isatty() or not c.color_terminal():
        c.nocolor()

    # Available colors: black, darkgray, (dark)red, dark(green),
    # brown, yellow, (dark)blue, purple, fuchsia, turquoise, teal,
    # lightgray, white.  Available styles: reset, bold, faint,
    # standout, underline, blink.

    # Set up log record formats.
    format_std = "%(message)s"
    formatter = logging.Formatter(format_std)

    # format_debug = "%(module)s #%(lineno)s %(funcName)s() %(message)s"
    fields = ['%(module)s', '#%(lineno)s', '%(funcName)s()', '%(message)s']
    colors = ['darkblue', 'darkred', 'brown', 'reset']
    styles = ['reset', 'reset', 'reset', 'reset']
    format_debug = ""
    for i in xrange(len(fields)):
        format_debug += c.colorize(styles[i], c.colorize(colors[i], fields[i]))
        if i != len(fields):
            format_debug += " "

    # Documentation:  http://docs.python.org/library/logging.html
    logger = logging.getLogger('doc.common.builder')

    # Note: There's also Handler.setLevel().  The argument is the
    # lowest severity message that the respective logger or handler
    # will pass on.  The default levels are DEBUG, INFO, WARNING,
    # ERROR, and CRITICAL.  We use "WARNING" for normal verbosity and
    # "ERROR" for quiet operation.  It's possible to define custom
    # levels.  See the documentation for details.
    if verbose == 0:
        logger.setLevel(logging.ERROR)
    if verbose == 1:
        logger.setLevel(logging.WARNING)
    if verbose == 2:
        logger.setLevel(logging.INFO)
    if verbose == 3:
        logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(format_debug)

    handler = logging.StreamHandler()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    return logger


if __name__ == '__main__':
    # Parse the command-line.
    parser = setup_parser()
    options, args = parser.parse_args()

    # Get the name and type (target format) of the document we are
    # trying to build.
    try:
        name, type = args
    except ValueError:
        help_message_short(parser=parser, error=True)
        sys.exit(1)

    # Set up module-wide logging.
    logger = setup_logger(options.verbose, options.color)

    # Process selected options.
    if options.jsmath:
        os.environ['SAGE_DOC_JSMATH'] = "True"

    if options.sphinx_opts:
        ALLSPHINXOPTS += options.sphinx_opts.replace(',', ' ') + " "

    if options.update_mtimes:
        update_mtimes = True
    else:
        update_mtimes = False

    # Make sure common/static exists.
    mkdir(os.path.join(SAGE_DOC, 'common', 'static'))

    # Get the builder and build.
    getattr(get_builder(name), type)()
