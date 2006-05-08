"""
Darcs from SAGE.

These functions make setup and use of darcs with SAGE
easier.
"""

import os, shutil

import package

PAGER='less'   # more doesn't work with darcs!

def darcs_install():
    """
    Download and install a statically binary darcs executable.

    This works on Cygwin, OS X, and intel-based Linux.  If you
    are using something else, obtain darcs yourself and put
    it anywhere in your PATH.
    """
    uname = os.uname()[0]
    if uname[:6] == 'CYGWIN':
        package.install_package('darcs_cygwin')
    elif uname == 'Darwin':
        package.install_package('darcs_darwin')
    elif uname == 'Linux':
        package.install_package('darcs_linux')
    else:
        raise RuntimeError, "No SAGE darcs package available for your platform (this just means you need to get a darcs binary yourself and put it somewhere in your PATH)."

known_installed = False
def darcs_ensure_installed():
    """
    Download and install darcs into your SAGE environment, if
    you do not already have darcs in your PATH.

    The darcs binary is put in <SAGE_ROOT>/bin.
    """
    global known_installed
    if known_installed:
        return
    if os.system('darcs --help 1>/dev/null 2>/dev/null') == 0:
        known_installed = True
        return
    darcs_install()

class Darcs:
    r"""
    This is a darcs repository.

    If you try to use it and don't have darcs installed on your computer,
    darcs will be automatically downloaded and installed (assuming you
    have a net connection).  Also, you do not have to explicitly initialize
    your repository in order to use it.   E.g., if you do
           \code{darcs_src.changes()}
    and you've never used darcs before, darcs will be downloaded and
    installed, then the latest source repository will be downloaded
    and installed.

    The few of the simplest and most useful commands are directly
    provided as member functions.  However, you can use the full
    functionality of darcs by noting that typing, e.g.,
            \code{darcs_src("command line arguments")}
    is \emph{exactly} the same as typing
            \code{cd repo_directory && darcs command line arguments | less}
    """
    def __init__(self, dir, name, url, target=None):
        """
        INPUT:
            dir -- directory that will contain the repository
            name -- a friendly name for the repository (only used for printing)
            url -- a default URL to pull or record sends against (e.g.,
                   this could be a master repository on modular.math.washington.edu)
            target -- if the last part of dir is, e.g., sage-darcs,
                      create a symlink from sage-darcs to target.
                      If target=None, this symlink will not be created.
        """
        self.__dir = os.path.abspath(dir)
        self.__name = name
        self.__url = url
        self.__initialized = False
        self.__target = target

    def __repr__(self):
        return "Darcs repository '%s' in directory %s"%(self.__name,
                                                        self.__dir)

    def __call__(self, cmd, check_initialized=True):
        """
        Run 'darcs cmd' where cmd is an arbitrary string
        in the darcs repository.
        """
	darcs_ensure_installed()
        if check_initialized and not self.__initialized:
            self.initialize()
        darcs_ensure_installed()
        s = 'cd "%s" && darcs %s'%(self.__dir, cmd)
        print s
        return os.system(s)

    def apply(self, patchfile, options=''):
        """
        Apply patches from a darcs patch to the repository.
        """
        patchfile = os.path.abspath(patchfile)
        print "Applying patchfile %s"%patchfile
        self('apply %s "%s"'%(options, patchfile))

    def changes(self, options=''):
        """
        Display the change log for this repository.
        """
        self('changes %s | %s'%(options, PAGER))

    def dir(self):
        """
        Return the directory where this repository is located.
        """
        return self.__dir

    def url(self):
        """
        Return the default 'master url' for this repository.
        """
        return self.__url

    def help(self, cmd=''):
        r"""
        Return help about the given command, or if cmd is omitted
        a list of commands.

        If this darcs object is called darcs_src, then you
        call a command using
             \code{darcs_src('usual darcs command line notation')}
        """
        self('%s --help | %s'%(cmd, PAGER))

    def initialize(self, force=False):
        """
        Create and initialize this darcs repository if you have not
        already done so.
        """
        if force or not os.path.exists('%s/_darcs'%self.__dir):
            if not os.path.exists(self.__dir):
                os.makedirs(self.__dir)
            print "Creating a new darcs repository!  %s"%self.__name
            if self('initialize', check_initialized=False):
                print "WARNING -- problem initializing repository."
                print "Try calling initalize again with the force option?"
                return
            if self.pull('-v -a'):
                print "WARNING -- problem pulling repository."
                print "Try calling initalize again with the force option?"
                return
            n = self.__dir.split('/')[-1]
            if not self.__target is None:
                os.system('cd "%s/.." && ln -snf %s %s'%(self.__dir, n, self.__target))
            if os.path.exists('%s/install'%self.__dir):
                # Darcs pull doesn't preserve permissions.
                os.system('chmod a+x %s/install'%self.__dir)
            self.__initialized = True
            if self.__target == 'sage':
                print ""
                print "Now building the new SAGE libraries"
                os.system('sage -b')
                print "You must restart SAGE in order for the changes to take effect."

    def pull(self, options='', url=None):
        """
        Pull all new patches from the repository at the given url,
        or use the default 'official' repository if no url is
        specified.
        """
        if url is None:
            url = self.__url
        self('pull %s %s'%(options, url))

    def record(self, options=''):
        """
        Interactively record changes as patches between your working
        copy and your repository.

        It's OK to hit control-c and restart if something goes wrong.
        """
        self('record %s'%options)

    def unrecord(self, options=''):
        """
        Remove recorded patches without changing the working copy.
        """
        self('unrecord %s'%options)

    def send(self, filename, options='', url=None):
        """
        Create a darcs patch bundle with the given filename
        against the repository at the given url (which is
        by default the 'official' SAGE repository).

        This is a file that you should probably post to
        sage-devel@lists.sourceforge.net.  It will
        be written to the current directory.

        NOTE: The darcs 'send' command by default tries to email
        patches.  Since email rarely works on users personal machines,
        in SAGE the default is to create a file.
        """
        if url is None:
            url = self.__url
        # We write to a local tmp file, then move, since unders
        # windows darcs has a bug that makes it fail to write
        # to any filename that is at all complicated!
        filename = os.path.abspath(filename) + '.darcs'
        print 'Writing to %s'%filename
        tmpfile = '%s/tmpdarcs'%self.__dir
        if os.path.exists(tmpfile):
            os.unlink(tmpfile)
        self('send %s %s -o tmpdarcs '%(options, url))
        if os.path.exists(tmpfile):
            shutil.move(tmpfile, filename)
            print 'Successfully created darcs patch bundle %s'%filename
        else:
            print 'Problem creating darcs patch bundle %s'%filename

    def what(self, options=''):
        """
        Show all changes between your local darcs repository and
        the working copy of your source code.
        """
        self('what %s | %s'%(options, PAGER))


#############################################################
# Create the default SAGE darcs repositories.
#############################################################

darcs_src = Darcs('%s/devel/sage-darcs'%os.environ['SAGE_ROOT'],
                  'SAGE source code',
        url="http://modular.math.washington.edu/sage/dist/src/sage-darcs",
                  target='sage')

darcs_doc = Darcs('%s/devel/doc-darcs'%os.environ['SAGE_ROOT'],
                  'SAGE documentation',
        url="http://modular.math.washington.edu/sage/dist/src/doc-darcs",
                  target='doc')

darcs_scripts = Darcs('%s/bin/'%os.environ['SAGE_LOCAL'],
                  'SAGE scripts',
                  url='',
                  target='')

