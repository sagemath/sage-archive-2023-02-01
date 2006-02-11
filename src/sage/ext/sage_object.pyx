r"""
Abstract base class for \sage objects
"""

import cPickle
import os

# changeto import zlib to use zlib instead; but this
# slows down loading any data stored in the other format
import zlib; comp = zlib
import bz2; comp_other = bz2

cdef add_ext(s):
    if s[-5:] != '.sobj':
        return s + '.sobj'
    else:
        return s

import sage.interfaces.gap

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
            <type 'mpfr.RealNumber'>
            sage: a.rename('pi')
            Traceback (most recent call last):
            ...
            NotImplementedError: object does not support renaming: 3.1399999999999997


        \note{The reason C-extension types are not supported is if
        they were then every single one would have to carry around an
        extra attribute, which would be slower and waste a lot of
        memory.}
        """
        if x is None and hasattr(self, '__custom_name'):
            self.reset_name()
        else:
            try:
                self.__custom_name = str(x)
            except AttributeError:
                raise NotImplementedError, "object does not support renaming: %s"%self
    def reset_name(self):
        del self.__custom_name

    def __repr__(self):
        try:
            return self.__custom_name
        except AttributeError:
            return self._repr_()

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

    def save(self, filename):
        """
        Save self to the given filename.

        EXAMPLES:
            sage.: f = x^3 + 5
            sage.: f.save('file')
            sage.: load('file')
            x^3 + 5
        """
        open(add_ext(filename), 'w').write(self.dumps())

    def dump(self, filename):
        """
        Same as self.save(filename)
        """
        return self.save(filename)

    def dumps(self):
        """
        Dump self to a string s, which can later be reconstituted
        as self using loads(s).
        """
        # the protocol=2 is very important -- this enables
        # saving extensions classes (with no attributes).
        return comp.compress(cPickle.dumps(self, protocol=2))

    def db(self, name=None):
        r"""
        Dumps self into the SAGE database.  Use db(name) by itself to
        reload.  The default name if none is specified is
        self._db_name(), which incorporates the class name (as
        directory) and self._defining_params() if defined, or
        str(self) if not.

        The database directory is \code{\$HOME/.sage/db}
        """
        if name is None:
            name = self._db_name()
        from sage.misc.all import SAGE_DB
        return self.dump('%s/%s'%(SAGE_DB,name))

    def _db_name(self):
        t = str(type(self)).split()[-1][1:-2]
        try:
            d = str(self._defining_params_())
        except AttributeError:
            d = str(self)
        d = '_'.join(d.split())
        from sage.misc.all import SAGE_DB
        if not os.path.exists('%s/%s'%(SAGE_DB, t)):
            os.makedirs(t)
        return '%s/%s'%(t, d)


    #############################################################################
    # Category theory / structure
    #############################################################################
    def Hom(self, codomain, cat=None):
        r"""
        self.Hom(codomain, cat=None):

        Return the homspace \code{Hom(self, codomain, cat)} of all
        homomorphisms from self to codomain in the category cat.  The
        default category is \code{self.category()}.

        EXAMPLES:
            sage: R, (x,y) = PolynomialRing(Q, 2, 'xy').objgens()
            sage: R.Hom(Q)
            Set of Homomorphisms from Polynomial Ring in x, y over Rational Field to Rational Field

        Homspaces are defined for very general \sage objects, even elements of familiar rings.
            sage: 5.Hom(7)
            Set of Morphisms from 5 to 7 in Category of elements of Integer Ring
            sage: (2/3).Hom(8/1)
            Set of Morphisms from 2/3 to 8 in Category of elements of Rational Field

        This example illustrates the optional third argument:
            sage: QQ.Hom(ZZ, Sets())
            Set of Morphisms from Rational Field to Integer Ring in Category of sets
        """
        from sage.categories.all import Hom
        return Hom(self, codomain, cat)

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
    # Containment testing
    #############################################################################

    def __contains__(self, x):
        r"""
        True if coercion of $x$ into self is possible.  Thus, e.g.,
        $2$ is in the integers and in $\Z/7\Z$ and the element $3$ of
        $\Z/7\Z$ {\em is} in $\Z$.  However, $2/3$ is not in $\Z$.
        Think of this as returning True if $x==y$ for some $y$ in
        self.

        EXAMPLES:
            sage: 2 in Integers(7)
            True
            sage: 2 in IntegerRing()
            True
            sage: Integers(7)(3) in IntegerRing()
            True
        """
        try:
            self(x)
        except TypeError:
            return False
        return True


    #############################################################################
    # Coercions to interface objects
    #############################################################################

    # SAGE
    def _sage_(self):
        return self

    # GAP
    def _gap_(self, G=None):
        if G is None:
            G = sage.interfaces.gap.gap  # default interpreter
        try:
            g = self.__gap
            if g.parent() is G:
                g._check_valid()
                return g
        except (AttributeError, ValueError):
            pass
        g = G(self._gap_init_())
        try:
            self.__gap = g
        except AttributeError:
            # do this because C-extension class won't have a __gap attribute.
            pass
        return g

    def _gap_init_(self):
        return str(self)
        #raise TypeError, "conversion of %s to GAP not yet implemented"%self


    # GP/PARI
    def _gp_(self, G=None):
        if G is None:
            G = sage.interfaces.gp.gp  # default interpreter
        try:
            g = self.__gp
            if g.parent() is G:
                g._check_valid()
                return g
        except (AttributeError, ValueError):
            pass
        g = G(self._gp_init_())
        try:
            self.__gp = g
        except AttributeError:
            # do this because C-extension class won't have a __gp attribute.
            pass
        return g

    def _gp_init_(self):
        return str(self)
        #raise TypeError, "conversion of %s to GP/PARI not yet implemented"%self

    # GP/PARI
    def _pari_(self):
        try:
            return self.__pari
        except AttributeError:
            pass
        from sage.libs.pari.all import pari
        x = pari(self._pari_init_())
        try:
            self.__pari = x
        except AttributeError:
            # do this because C-extension class won't have a __pari attribute.
            pass
        return x

    def _pari_init_(self):
        return self._gp_init_()

    # Singular
    def _singular_(self, G=None):
        if G is None:
            G = sage.interfaces.singular.singular  # default interpreter
        try:
            g = self.__singular
            if g.parent() is G:
                g._check_valid()
                return g
        except (AttributeError, ValueError):
            pass
        g = G(self._singular_init_())
        try:
            self.__singular = g
        except AttributeError:
            # do this because C-extension class won't have a __singular attribute.
            pass
        return g

    def _singular_init_(self):
        return str(self)
        #raise TypeError, "conversion of %s to Singular not yet implemented"%self

    # Maxima
    def _maxima_(self, G=None):
        if G is None:
            G = sage.interfaces.maxima.maxima  # default interpreter
        try:
            g = self.__maxima
            if g.parent() is G:
                g._check_valid()
                return g
        except (AttributeError, ValueError):
            pass
        g = G(self._maxima_init_())
        try:
            self.__maxima = g
        except AttributeError:
            # do this because C-extension class won't have a __maxima attribute.
            pass
        return g

    def _maxima_init_(self):
        return str(self)
        # raise TypeError, "conversion of %s to Maxima not yet implemented"%self

    # Magma
    def _magma_(self, G=None):
        if G is None:
            G = sage.interfaces.magma.magma  # default interpreter
        try:
            g = self.__magma
            if g.parent() is G:
                g._check_valid()
                return g
        except (AttributeError, ValueError):
            pass
        g = G(self._magma_init_())
        try:
            self.__magma = g
        except AttributeError:
            # do this because C-extension class won't have a __magma attribute.
            pass
        return g

    def _magma_init_(self):
        return str(self)
        # raise TypeError, "conversion of %s to Magma not yet implemented"%self




##################################################################


def load(filename):
    """
    load(filename):

    Load \sage object from the file with name filename, which will
    have an .sobj extension added if it doesn't have one.
    """
    return loads(open(add_ext(filename)).read())

def save(obj, filename):
    """
    save(obj, filename):

    Save obj to the file with name filename, which will
    have an .sobj extension added if it doesn't have one.
    This will \emph{replace} the contents of filename.
    """
    try:
        obj.save(add_ext(filename))
    except AttributeError, RuntimeError:
        s = comp.compress(cPickle.dumps(obj, protocol=2))
        open(add_ext(filename), 'w').write(s)

def dumps(obj):
    """
    dumps(obj):

    Dump obj to a string s.  To recover obj, use loads(s).
    """
    try:
        return obj.dumps()
    except (AttributeError, RuntimeError, TypeError):
        return comp.compress(cPickle.dumps(obj, protocol=2))

def loads(s):
    """
    Recover an object x that has been dumped to a string s
    using s = dumps(x).
    """
    if not isinstance(s, str):
        raise TypeError, "s must be a string"
    try:
        return cPickle.loads(comp.decompress(s))
    except:
        return cPickle.loads(comp_other.decompress(s))

