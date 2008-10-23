r"""
Abstract base class for SAGE objects
"""

import cPickle
import os
import sys

# changeto import zlib to use zlib instead; but this
# slows down loading any data stored in the other format
import zlib; comp = zlib
import bz2; comp_other = bz2

base=None

cdef process(s):
    if not base is None and (len(s) == 0 or s[0] != '/'):
        s = base + '/' + s
    if s[-5:] != '.sobj':
        return s + '.sobj'
    else:
        return s


cdef class SageObject:

    #############################################################################
    # Textual representation code
    #############################################################################

    def rename(self, x=None):
        r"""
        Change self so it prints as x, where x is a string.

        \note{This is \emph{only} supported for Python classes that derive
        from SageObject.}

        EXAMPLES:
            sage: x = PolynomialRing(QQ,'x').gen()
            sage: g = x^3 + x - 5
            sage: g
            x^3 + x - 5
            sage: g.rename('a polynomial')
            sage: g
            a polynomial
            sage: g + x
            x^3 + 2*x - 5
            sage: h = g^100
            sage: str(h)[:20]
            'x^300 + 100*x^298 - '
            sage: h.rename('x^300 + ...')
            sage: h
            x^300 + ...

        Real numbers are not Python classes, so rename is not supported:
            sage: a = 3.14
            sage: type(a)
            <type 'sage.rings.real_mpfr.RealLiteral'>
            sage: a.rename('pi')
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support renaming: 3.14000000000000

        \note{The reason C-extension types are not supported by default
        is if they were then every single one would have to carry around
        an extra attribute, which would be slower and waste a lot of
        memory.

        To support them for a specific class, add a \code{cdef public __custom_name}
        attribute.}
        """
        if x is None:
            if hasattr(self, '__custom_name'):
                self.reset_name()
        else:
            try:
                self.__custom_name = str(x)
            except AttributeError:
                raise NotImplementedError, "object does not support renaming: %s"%self

    def reset_name(self):
        if hasattr(self, '__custom_name'):
            del self.__custom_name


    def __repr__(self):
        if hasattr(self, '__custom_name'):
            name = self.__custom_name
            if name is not None:
                return name
        if hasattr(self, '_repr_'):
            return self._repr_()
        return str(type(self))

    def __hash__(self):
        return hash(self.__repr__())


    #############################################################################
    # DATABASE Related code
    #############################################################################

    def version(self):
        r"""
        The version of \sage.

        Call this to save the version of \sage in this object.
        If you then save and load this object it will know in what
        version of \sage it was created.

        This only works on Python classes that derive from SageObject.
        """
        try:
            return self.__version
        except AttributeError:
            import sage.version
            self.__version = sage.version.version
            return self.__version

    def save(self, filename=None, compress=True):
        """
        Save self to the given filename.

        EXAMPLES:
            sage: f = x^3 + 5
            sage: f.save(SAGE_TMP + '/file')
            sage: load(SAGE_TMP + '/file.sobj')
            x^3 + 5
        """
        if filename is None:
            try:
                filename = self._default_filename
            except AttributeError:
                raise RuntimeError, "no default filename, so it must be specified"
        filename = process(filename)
        try:
            self._default_filename = filename
        except AttributeError:
            pass
        open(filename, 'wb').write(self.dumps(compress))

    def dump(self, filename, compress=True):
        """
        Same as self.save(filename, compress)
        """
        return self.save(filename, compress=compress)

    def dumps(self, compress=True):
        """
        Dump self to a string s, which can later be reconstituted
        as self using loads(s).
        """
        # the protocol=2 is very important -- this enables
        # saving extensions classes (with no attributes).
        s = cPickle.dumps(self, protocol=2)
        if compress:
            return comp.compress(s)
        else:
            return s

    def db(self, name, compress=True):
        r"""
        Dumps self into the SAGE database.  Use db(name) by itself to
        reload.

        The database directory is \code{\$HOME/.sage/db}
        """
        #if name is None:
        #    name = self._db_name()
        from sage.misc.all import SAGE_DB
        return self.dump('%s/%s'%(SAGE_DB,name), compress=compress)

##     def _db_name(self):
##         t = str(type(self)).split()[-1][1:-2]
##         try:
##             d = str(self._defining_params_())
##         except AttributeError:
##             d = str(self)
##         d = '_'.join(d.split())
##         from sage.misc.all import SAGE_DB
##         if not os.path.exists('%s/%s'%(SAGE_DB, t)):
##             os.makedirs(t)
##         return '%s/%s'%(t, d)


    #############################################################################
    # Category theory / structure
    #############################################################################

    def category(self):
        from sage.categories.all import Objects
        return Objects()

##     def category(self):
##         try:
##             return self.__category
##         except AttributeError:
##             from sage.categories.all import Objects
##             return Objects()

##     def _set_category(self, C):
##         self.__category = C


    #############################################################################
    # Coercions to interface objects
    #############################################################################

    # SAGE
    def _sage_(self):
        return self

    def _interface_(self, I):
        """
        Return coercion of self to an object of the interface I.

        The result of coercion is cached, unless self is not a C
        extension class or \code{self._interface_is_cached_()} returns
        False.
        """
        c = self._interface_is_cached_()
        if c:
            try:
                X = self.__interface[I]
                X._check_valid()
                return X
            except (AttributeError, TypeError):
                try:
                    self.__interface = {}
                except AttributeError:
                    # do this because C-extension classes won't have
                    # an __interface attribute.
                    pass
            except (KeyError, ValueError):
                pass
        if hasattr(self, '_%s_init_'%I.name()):
            s = self.__getattribute__('_%s_init_'%I.name())()
        elif hasattr(self, '_system_init_'):
            s = self._system_init_(I.name())
        else:
            try:
              s = self._interface_init_()
            except:
                raise NotImplementedError, "coercion of object %s to %s not implemented:\n%s\n%s"%\
                  (repr(self), I)
        X = I(s)
        if c:
            try:
                self.__interface[I] = X
            except AttributeError:
                pass
        return X

    def _interface_init_(self):
        return repr(self)

    def _interface_is_cached_(self):
        """
        Return True if the interface objects are cached.

        If you have an object x and do gp(x), the result is cached if
        this function returns True.
        """
        return True

    def _gap_(self, G=None):
        if G is None:
            import sage.interfaces.gap
            G = sage.interfaces.gap.gap
        return self._interface_(G)

    def _gap_init_(self):
        return self._interface_init_()

    def _gp_(self, G=None):
        if G is None:
            import sage.interfaces.gp
            G = sage.interfaces.gp.gp
        return self._interface_(G)

    def _gp_init_(self):
        return self._pari_init_()

    def _kash_(self, G=None):
        if G is None:
            import sage.interfaces.kash
            G = sage.interfaces.kash.kash
        return self._interface_(G)

    def _kash_init_(self):
        return self._interface_init_()

    def _axiom_(self, G=None):
        if G is None:
            import sage.interfaces.axiom
            G = sage.interfaces.axiom.axiom
        return self._interface_(G)

    def _axiom_init_(self):
        return self._interface_init_()

    def _maxima_(self, G=None):
        if G is None:
            import sage.interfaces.maxima
            G = sage.interfaces.maxima.maxima
        return self._interface_(G)

    def _maxima_init_(self):
        return self._interface_init_()


    def _magma_(self, M=None):
        """
        Given a Magma interpreter M (or None for the default global
        interpreter), return MagmaElement corresponding to self in M.
        This properly keeps tracking of caching, and also checks for
        validity of an element before returning it.  Note that the
        Magma versions of cdef'd elements are not cached.

        This should not be redefined in derived classes.  Instead,
        derived classes should redefine _magma_convert_.

        INPUT:
            M -- a Magma interpreter
        OUTPUT:
            a Magma object

        EXAMPLES:
            sage: n = -3/7
            sage: m2 = Magma()
            sage: n._magma_()                        # optional -- requires magma
            -3/7
            sage: n._magma_().parent()               # optional
            Magma
            sage: n._magma_().parent() is m2         # optional
            False
            sage: n._magma_().parent() is magma      # optional
            True
            sage: n._magma_(m2).parent() is m2       # optional
            True

        This example illustrates caching, which happens automatically
        since K is a Python object:
            sage: K.<a> = NumberField(x^3 + 2)
            sage: K._magma_() is K._magma_()        # optional
            True
            sage: magma2 = Magma()
            sage: K._magma_() is K._magma_(magma2)  # optional
            False
        """
        if M is None:
            import sage.interfaces.magma
            M = sage.interfaces.magma.magma
        c = self._interface_is_cached_()
        if c:
            try:
                X = self.__interface[M]
                X._check_valid()
                return X
            except (AttributeError, TypeError):
                try:
                    self.__interface = {}
                except AttributeError:
                    # do this because C-extension classes won't have
                    # an __interface attribute.
                    c = False
                    pass
            except (KeyError, ValueError):
                pass
        A = self._magma_convert_(M)
        if c:
            self.__interface[M] = A
        return A

    def _magma_init_(self):
        """
        Return an ascii string that evaluates to something equal to
        self in Magma.  Use this for converting very simple things
        (e.g., integers) from Sage to Magma.  For anything much more
        complicated, use _magma_convert_.

        The default coercion for elements from Sage to Magma is to
        call _magma_init_, which just calls the repr method of the
        object.

        OUTPUT:
            string

        EXAMPLES:
            sage: n = -3/7
            sage: n._magma_init_()
            '-3/7'
        """
        return self._interface_init_()

    def _magma_convert_(self, M):
        """
        Given a magma interpreter, this function should return a
        MagmaElement in M.

        One should usually redefine this function in the derived
        class.  There is no need to worry about caching, which is done
        automatically by the infrastructure in sage_object.

        INPUT:
            M -- a Magma interpreter
        OUTPUT:
            a Magma object

        EXAMPLES:
            sage: n = -3/7
            sage: n._magma_convert_(magma)       # optional -- requires magma
            -3/7
        """
        return self._interface_(M)

    def _macaulay2_(self, G=None):
        if G is None:
            import sage.interfaces.macaulay2
            G = sage.interfaces.macaulay2.macaulay2
        return self._interface_(G)

    def _macaulay2_init_(self):
        return self._interface_init_()

    def _maple_(self, G=None):
        if G is None:
            import sage.interfaces.maple
            G = sage.interfaces.maple.maple
        return self._interface_(G)

    def _maple_init_(self):
        return self._interface_init_()

    def _mathematica_(self, G=None):
        if G is None:
            import sage.interfaces.mathematica
            G = sage.interfaces.mathematica.mathematica
        return self._interface_(G)

    def _mathematica_init_(self):
        return self._interface_init_()

    def _octave_(self, G=None):
        if G is None:
            import sage.interfaces.octave
            G = sage.interfaces.octave.octave
        return self._interface_(G)

    def _octave_init_(self):
        return self._interface_init_()

    def _r_init_(self):
        """
        Return default string expression that evaluates in R to this
        object.

        OUTPUT:
            string

        EXAMPLES:
            sage: a = 2/3
            sage: a._r_init_()
            '2/3'
        """
        return self._interface_init_()

    def _singular_(self, G=None, have_ring=False):
        if G is None:
            import sage.interfaces.singular
            G = sage.interfaces.singular.singular
        return self._interface_(G)

    def _singular_init_(self, have_ring=False):
        return self._interface_init_()

    # PARI (slightly different, since is via C library, hence instance is unique)
    def _pari_(self):
        if self._interface_is_cached_():
            try:
                return self.__pari
            except AttributeError:
                pass
        from sage.libs.pari.all import pari
        x = pari(self._pari_init_())
        if self._interface_is_cached_():
            try:
                self.__pari = x
            except AttributeError:
                # do this because C-extension class won't have a __pari attribute.
                pass
        return x

    def _pari_init_(self):
        return self._interface_init_()


##################################################################




def load(filename, compress=True, verbose=True):
    """
    load(filename):

    Load \sage object from the file with name filename, which will
    have an .sobj extension added if it doesn't have one.

    NOTE: There is also a special SAGE command (that is not
    available in Python) called load that you use by typing

                sage: load filename.sage           # not tested

    The documentation below is not for that command.  The documentation
    for load is almost identical to that for attach.  Type attach? for
    help on attach.

    This also loads a ".sobj" file over a network by specifying the full URL.
    (Setting "verbose = False" suppresses the loading progress indicator.)

    EXAMPLE:
        sage: u = 'http://sage.math.washington.edu/home/was/db/test.sobj'  # optional
        sage: s = load(u)                                                  # optional
        Attempting to load remote file: http://sage.math.washington.edu/home/was/db/test.sobj
        Loading: [.]
        sage: s                                                            # optional
        'hello SAGE'
    """

    ## Check if filename starts with "http://" or "https://"
    lower = filename.lower()
    if lower.startswith("http://") or lower.startswith("https://"):
        from sage.misc.remote_file import get_remote_file
        filename = get_remote_file(filename, verbose=verbose)
        tmpfile_flag = True
    elif lower.endswith('.f') or lower.endswith('.f90'):
        globals()['fortran'](filename)
        return
    else:
        tmpfile_flag = False
        filename = process(filename)

    ## Load file by absolute filename
    X = loads(open(filename).read(), compress=compress)
    try:
        X._default_filename = os.path.abspath(filename)
    except AttributeError:
        pass

    ## Delete the tempfile, if it exists
    if tmpfile_flag == True:
        os.unlink(filename)

    return X


def save(obj, filename=None, compress=True, **kwds):
    """
    save(obj, filename=None):

    Save obj to the file with name filename, which will
    have an .sobj extension added if it doesn't have one.
    This will \emph{replace} the contents of filename.

    EXAMPLES:
        sage: a = matrix(2, [1,2,3,-5/2])
        sage: save(a, 'test.sobj')
        sage: load('test')
        [   1    2]
        [   3 -5/2]
        sage: E = EllipticCurve([-1,0])
        sage: P = plot(E)
        sage: save(P, 'test')
        sage: save(P, filename="sage.png", xmin=-2)
        sage: print load('test.sobj')
        Graphics object consisting of 2 graphics primitives
        sage: save("A python string", './test')
        sage: load('./test.sobj')
        'A python string'
        sage: load('./test')
        'A python string'
    """
    # Add '.sobj' if the filename currently has no extension
    if os.path.splitext(filename)[1] == '':
        filename += '.sobj'

    if filename.endswith('.sobj'):
        try:
            obj.save(filename=filename, compress=compress, **kwds)
        except (AttributeError, RuntimeError, TypeError):
            s = cPickle.dumps(obj, protocol=2)
            if compress:
                s = comp.compress(s)
            open(process(filename), 'wb').write(s)
    else:
        # Saving an object to an image file.
        obj.save(filename, **kwds)

def dumps(obj, compress=True):
    """
    dumps(obj):

    Dump obj to a string s.  To recover obj, use loads(s).

    EXAMPLES:
        sage: a = 2/3
        sage: s = dumps(a)
        sage: print len(s)
        49
        sage: loads(s)
        2/3
    """
    if make_pickle_jar:
        picklejar(obj)
    try:
        return obj.dumps(compress)
    except (AttributeError, RuntimeError, TypeError):
        if compress:
            return comp.compress(cPickle.dumps(obj, protocol=2))
        else:
            return cPickle.dumps(obj, protocol=2)

def loads(s, compress=True):
    """
    Recover an object x that has been dumped to a string s
    using s = dumps(x).

    EXAMPLES:
        sage: a = matrix(2, [1,2,3,-4/3])
        sage: s = dumps(a)
        sage: loads(s)
        [   1    2]
        [   3 -4/3]

    One can load uncompressed data even if one messes up
    and doesn't specify compress=False.  This is slightly
    slower though.
        sage: v = [1..10]
        sage: loads(dumps(v, compress=False)) == v
        True
        sage: loads(dumps(v, compress=False), compress=True) == v
        True
        sage: loads(dumps(v, compress=True), compress=False) == v
        True
    """
    if not isinstance(s, str):
        raise TypeError, "s must be a string"
    if compress:
        try:
            return cPickle.loads(comp.decompress(s))
        except Exception, msg1:
            try:
                return cPickle.loads(comp_other.decompress(s))
            except Exception, msg2:
                # Maybe data is uncompressed?
                try:
                    return cPickle.loads(s)
                except Exception, msg3:
                    try: msg1 = str(msg1)
                    except: msg1 = type(msg1)
                    try: msg2 = str(msg2)
                    except: msg2 = type(msg2)
                    try: msg3 = str(msg3)
                    except: msg3 = type(msg3)
                    raise RuntimeError, "%s\n%s\n%s\nUnable to load pickled data."%(msg1,msg2,msg3)
    else:
        try:
            return cPickle.loads(s)
        except:
            # maybe data is compressed anyways??
            return loads(s, compress=True)


cdef bint make_pickle_jar = os.environ.has_key('SAGE_PICKLE_JAR')

def picklejar(obj, dir=None):
    """
    Create pickled sobj of obj in dir, with name the absolute value of
    the hash of the pickle of obj.  This is used in conjection with
    sage.structure.sage_object.unpickle_all.

    To use this to test the whole Sage library right now, set the
    environment variable SAGE_PICKLE_JAR, which will make it so dumps
    will by default call picklejar with the default dir.  Once you do
    that and doctest Sage, you'll find that the SAGE_ROOT/tmp/
    contains a bunch of pickled objects along with corresponding txt
    descriptions of them.  Use the
    sage.structure.sage_object.unpickle_all to see if they unpickle
    later.

    INPUTS:
        obj -- a pickleable object
        dir -- a string or None; if None defaults to
               SAGE_ROOT/tmp/pickle_jar-version

    EXAMPLES:
        sage: dir = tmp_dir()
        sage: sage.structure.sage_object.picklejar(1,dir)
        sage: len(os.listdir(dir))
        2
    """
    if dir is None:
        from sage.version import version
        dir = os.environ['SAGE_ROOT'] + '/tmp/pickle_jar-%s/'%version
    if not os.path.exists(dir):
        os.makedirs(dir)

    s = comp.compress(cPickle.dumps(obj,protocol=2))

    typ = str(type(obj))
    name = ''.join([x if (x.isalnum() or x == '_') else '_' for x in typ])
    base = '%s/%s'%(dir, name)
    if os.path.exists(base):
        i = 0
        while os.path.exists(base + '-%s'%i):
            i += 1
        base += '-%s'%i

    open(base + '.sobj', 'wb').write(s)
    txt = "type(obj) = %s\n"%typ
    import sage.version
    txt += "version = %s\n"%sage.version.version
    txt += "obj =\n'%s'\n"%str(obj)

    open(base + '.txt', 'w').write(txt)

def unpickle_all(dir):
    """
    Unpickle all sobj's in the given directory, reporting failures as
    they occur.  Also printed the number of successes and failure.

    INPUT:
        dir -- string; a directory or name of a .tar.bz2 file that
               decompresses to give a directo pickirectory.

    EXAMPLES:
        sage: dir = tmp_dir()
        sage: sage.structure.sage_object.picklejar('hello', dir)
        sage: sage.structure.sage_object.unpickle_all(dir)
        Successfully unpickled 1 objects.
        Failed to unpickle 0 objects.

    We unpickle the standard pickle jar. This doctest tests that
    all "standard pickles" unpickle.  Every so often the standard pickle jar
    should be updated by running the doctest suite with the environment variable
    SAGE_PICKLE_JAR set, then copying the files from SAGE_ROOT/tmp/pickle_jar*
    into the standard pickle jar.
        sage: std = os.environ['SAGE_DATA'] + '/extcode/pickle_jar/pickle_jar.tar.bz2'
        sage: sage.structure.sage_object.unpickle_all(std)
        Successfully unpickled ... objects.
        Failed to unpickle 0 objects.
    """
    i = 0
    j = 0
    failed = []
    if dir.endswith('.tar.bz2'):
        # create a temporary directory
        from sage.misc.all import tmp_dir
        T = tmp_dir()
        # extract tarball to it
        os.system('cd "%s"; bunzip2 -c "%s" | tar fx - '%(T, os.path.abspath(dir)))
        # Now use the directory in the tarball instead of dir
        dir = T + "/" + os.listdir(T)[0]

    for A in os.listdir(dir):
        if A.endswith('.sobj'):
            try:
                load(dir + '/' + A)
                i += 1
            except Exception, msg:
                j += 1
                print "** failed: ", A
                failed.append(A)

    if len(failed) > 0:
        print "Failed:\n%s"%('\n'.join(failed))
    print "Successfully unpickled %s objects."%i
    print "Failed to unpickle %s objects."%j
