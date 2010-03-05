"""
libSingular: Options.

Singular uses a set of global options to determine verbosity and the
behavior of certain algorithms. We provide an interface to these
options in the most 'natural' python-ic way. Users who do not wish to
deal with Singular functions directly usually do not have to worry
about this interface or Singular options in general since this is
taken care of by higher level functions.

We compute a Groebner basis for Cyclic-5 in two different contexts::

    sage: P.<a,b,c,d,e> = PolynomialRing(GF(127))
    sage: I = sage.rings.ideal.Cyclic(P)
    sage: std = sage.libs.singular.ff.std

By default, tail reductions are performed::

    sage: from sage.libs.singular.option import opt, opt_ctx
    sage: opt['redTail']
    True
    sage: std(I)[-1]
    d^2*e^6 + 28*b*c*d + ...

If we don't want this, we can create an option context, which disables
this::

    sage: with opt_ctx(redTail=False,redSB=False):
    ...      std(I)[-1]
    d^2*e^6 + 8*c^3 + ...

However, this does not affect the global state::

    sage: opt['redTail']
    True

On the other hand, any assignment to an option object will immediately
change the global state::

    sage: opt['redTail'] = False
    sage: opt['redTail']
    False
    sage: opt['redTail'] = True
    sage: opt['redTail']
    True

Assigning values within an option context, only affects this context::

    sage: with opt_ctx:
    ...      opt['redTail'] = False

    sage: opt['redTail']
    True

Option contexts can also be safely stacked::

    sage: with opt_ctx:
    ...       opt['redTail'] = False
    ...       print opt
    ...       with opt_ctx:
    ...           opt['redThrough'] = False
    ...           print opt
    ...
    general options for libSingular (current value 0x00000082)
    general options for libSingular (current value 0x00000002)

    sage: print opt
    general options for libSingular (current value 0x02000082)


The same interface is available for verbosity options::

    sage: from sage.libs.singular.option import opt_verb
    sage: opt_verb['notWarnSB']
    False

AUTHOR:

- Martin Albrecht (2009-08): initial implementation
- Martin Albrecht (2010-01): better interface, verbosity options
"""
#*****************************************************************************
#       Copyright (C) 2010 Martin Albrecht <M.R.Albrecht@rhul.ac.uk>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.libs.singular.decl cimport singular_options, singular_verbose_options

from sage.libs.singular.decl cimport OPT_PROT, OPT_REDSB, OPT_NOT_BUCKETS, OPT_NOT_SUGAR, OPT_SUGARCRIT, OPT_REDTHROUGH
from sage.libs.singular.decl cimport OPT_RETURN_SB, OPT_FASTHC, OPT_OLDSTD, OPT_REDTAIL, OPT_INTSTRATEGY, OPT_NOTREGULARITY
from sage.libs.singular.decl cimport OPT_WEIGHTM, Sy_bit

from sage.libs.singular.decl cimport V_SHOW_MEM, V_YACC, V_REDEFINE, V_READING, V_LOAD_LIB, V_DEBUG_LIB
from sage.libs.singular.decl cimport V_LOAD_PROC, V_DEF_RES, V_SHOW_USE, V_IMAP, V_PROMPT
from sage.libs.singular.decl cimport V_NSB, V_CONTENTSB, V_CANCELUNIT, V_DEG_STOP


cdef class LibSingularOptions_abstract:
    """
    Abstract Base Class for libSingular options.
    """
    cdef unsigned int *global_options
    cdef object name
    cdef object name_map

    def __init__(self, **kwds):
        """
        INPUT:

        - ``**kwds`` - all keyword parameters are immediately applied.

        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: opt = LibSingularOptions(redTail=False)
            sage: opt['redTail']
            False
            sage: opt['redTail'] = True
            sage: opt['redTail']
            True
        """
        for k,v in kwds.iteritems():
            self[k] = v

    def __getitem__(self, name):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt
            sage: opt['redTail']
            True
        """
        try:
            return bool(self.global_options[0] & self.name_map[name])
        except KeyError:
            raise NameError("Option '%s' unknown."%(name,))

    def __setitem__(self, name, value):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx:
            ...      opt['redTail'] = False
            ...      opt['redTail']
            False
        """
        try:
            if value:
                self.global_options[0] = self.global_options[0] | self.name_map[name]
            else:
                self.global_options[0] = self.global_options[0] & ~self.name_map[name]
        except KeyError:
            raise NameError("Option '%s' unknown."%(name,))

    def __int__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt
            sage: hex(int(opt))
            '0x2000082'
        """
        return int(self.global_options[0])

    def load(self, value=None):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt as sopt
            sage: bck = int(sopt); hex(bck)
            '0x2000082'
            sage: sopt['redTail'] = False
            sage: hex(int(sopt))
            '0x82'
            sage: sopt.load(bck)
            sage: sopt['redTail']
            True
        """
        self.global_options[0] = int(value)

    def __repr__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt as sopt
            sage: sopt
            general options for libSingular (current value 0x02000082)
        """
        return "%s options for libSingular (current value 0x%08x)"%(self.name, self.global_options[0])


cdef class LibSingularOptions(LibSingularOptions_abstract):
    """
    Pythonic Interface to libSingular's options.

    Supported options are:

     - ``returnSB`` - the functions syz, intersect, quotient, modulo
       return a standard base instead of a generating set if returnSB
       is set. This option should not be used for lift.

     - ``fastHC`` - tries to find the highest corner of the
       staircase (HC) as fast as possible during a standard basis
       computation (only used for local orderings).

     - ``intStrategy`` - avoids division of coefficients during
       standard basis computations. This option is ring dependent. By
       default, it is set for rings with characteristic 0 and not set
       for all other rings.

     - ``lazy`` - uses a more lazy approach in std computations, which
       was used in SINGULAR version before 2-0 (and which may lead to
       faster or slower computations, depending on the example)

     - ``length`` - select shorter reducers in std computations,

     - ``notRegularity`` - disables the regularity bound for ``res`` and
       ``mres`` (see regularity).

     - ``notSugar`` - disables the sugar strategy during standard
       basis computation.

     - ``notBuckets`` - disables the bucket representation of
       polynomials during standard basis computations. This option
       usually decreases the memory usage but increases the
       computation time. It should only be set for memory-critical
       standard basis computations.

     - ``oldStd`` - uses a more lazy approach in std computations,
       which was used in SINGULAR version before 2-0 (and which may
       lead to faster or slower computations, depending on the
       example)

     - ``prot`` - shows protocol information indicating the progress
       during the following computations: ``facstd``, ``fglm``,
       ``groebner``, ``lres``, ``mres``, ``minres``, ``mstd``,
       ``res``, ``slimgb``,``sres``, ``std``, ``stdfglm``,
       ``stdhilb``, ``syz``.

     - ``redSB`` - computes a reduced standard basis in any standard
       basis computation.

     - ``redTail`` - reduction of the tails of polynomials during
       standard basis computations. This option is ring dependent. By
       default, it is set for rings with global degree orderings and
       not set for all other rings.

     - ``redThrough`` - for inhomogenous input, polynomial reductions
       during standard basis computations are never postponed, but
       always finished through. This option is ring dependent. By
       default, it is set for rings with global degree orderings and
       not set for all other rings.

     - ``sugarCrit`` - uses criteria similar to the homogeneous case
       to keep more useless pairs.

     - ``weightM`` - automatically computes suitable weights for the
       weighted ecart and the weighted sugar method.

    EXAMPLE::

        sage: from sage.libs.singular.option import LibSingularOptions
        sage: libsingular_options = LibSingularOptions()
        sage: libsingular_options
        general options for libSingular (current value 0x02000082)
    """
    def __init__(self, **kwds):
        """
        Create a new option interface.

        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: libsingular_options = LibSingularOptions()
            sage: libsingular_options
            general options for libSingular (current value 0x02000082)
        """
        self.global_options = &singular_options
        self.name = "general"
        self.name_map = {"prot":          Sy_bit(OPT_PROT),
                         "redSB":         Sy_bit(OPT_REDSB),
                         "notBuckets":    Sy_bit(OPT_NOT_BUCKETS),
                         "notSugar":      Sy_bit(OPT_NOT_SUGAR),
                         "sugarCrit":     Sy_bit(OPT_SUGARCRIT),
                         "redThrough":    Sy_bit(OPT_REDTHROUGH),
                         "returnSB":      Sy_bit(OPT_RETURN_SB),
                         "fastHC":        Sy_bit(OPT_FASTHC),
                         "oldStd":        Sy_bit(OPT_OLDSTD),
                         "redTail":       Sy_bit(OPT_REDTAIL),
                         "intStrategy":   Sy_bit(OPT_INTSTRATEGY),
                         "notRegularity": Sy_bit(OPT_NOTREGULARITY),
                         "weightM":       Sy_bit(OPT_WEIGHTM)}
        LibSingularOptions_abstract.__init__(self, **kwds)


#############

cdef class LibSingularVerboseOptions(LibSingularOptions_abstract):
    """
    Pythonic Interface to libSingular's verbosity options.

    Supported options are:

     - ``mem`` -
     - ``yacc`` -
     - ``redefine`` -
     - ``reading`` -
     - ``loadLib`` -
     - ``debugLib`` -
     - ``loadProc`` -
     - ``defRes`` -
     - ``usage`` -
     - ``Imap`` -
     - ``prompt`` -
     - ``notWarnSB`` - do not warn if a basis is not a standard basis
     - ``contentSB`` -
     - ``cancelunit`` -

    EXAMPLE::

        sage: from sage.libs.singular.option import LibSingularVerboseOptions
        sage: libsingular_verbose = LibSingularVerboseOptions()
        sage: libsingular_verbose
        verbosity options for libSingular (current value 0x00002851)
    """
    def __init__(self, **kwds):
        """
        Create a new option interface.

        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularVerboseOptions
            sage: libsingular_verb_options = LibSingularVerboseOptions()
            sage: libsingular_verb_options
            verbosity options for libSingular (current value 0x00002851)
        """
        self.global_options = &singular_verbose_options
        self.name = "verbosity"
        self.name_map = {"mem":      Sy_bit(V_SHOW_MEM),
                         "yacc":     Sy_bit(V_YACC),
                         "redefine": Sy_bit(V_REDEFINE),
                         "reading":  Sy_bit(V_READING),
                         "loadLib":  Sy_bit(V_LOAD_LIB),
                         "debugLib": Sy_bit(V_DEBUG_LIB),
                         "loadProc": Sy_bit(V_LOAD_PROC),
                         "defRes":   Sy_bit(V_DEF_RES),
                         "usage":    Sy_bit(V_SHOW_USE),
                         "Imap":     Sy_bit(V_IMAP),
                         "prompt":   Sy_bit(V_PROMPT),
                         "notWarnSB":Sy_bit(V_NSB),
                         "contentSB":Sy_bit(V_CONTENTSB),
                         "cancelunit":Sy_bit(V_CANCELUNIT),
                         }
        LibSingularOptions_abstract.__init__(self, **kwds)

cdef class LibSingularOptionsContext:
    """
    Option context

    This object localizes changes to options.

    EXAMPLE::

        sage: from sage.libs.singular.option import opt, opt_ctx
        sage: opt
        general options for libSingular (current value 0x02000082)

    ::

        sage: with opt_ctx(redTail=False):
        ...       print opt
        ...       with opt_ctx(redThrough=False):
        ...           print opt
        ...
        general options for libSingular (current value 0x00000082)
        general options for libSingular (current value 0x00000002)

        sage: print opt
        general options for libSingular (current value 0x02000082)
    """
    cdef list bck
    cdef public LibSingularOptions_abstract opt
    cdef object options

    def __init__(self, LibSingularOptions_abstract opt, **kwds):
        """
        Create a new context.

        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptionsContext, opt
            sage: LibSingularOptionsContext(opt)
            general options context for libSingular
        """
        self.bck = []
        self.options = kwds
        self.opt = opt

    def __enter__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx(redTail=False):
            ...     opt['redTail']
            False
        """
        self.bck.append(self.opt.global_options[0])
        opt = self.opt.__class__()
        for k,v in self.options.iteritems():
            opt[k] = v

    def __call__(self, **kwds):
        """
        Return a new option context where ``**kwds`` are applied.

        EXAMPLE::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx(redTail=False):
            ...     opt['redTail']
            False
        """
        new = self.__class__(self.opt, **kwds)
        return new

    def __exit__(self, typ, value, tb):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt, opt_ctx
            sage: opt['redTail']
            True
            sage: with opt_ctx(redTail=False):
            ...     opt['redTail']
            False
        """
        self.opt.global_options[0] = self.bck.pop()

    def __repr__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import opt_ctx
            sage: opt_ctx
            general options context for libSingular
        """
        return "%s options context for libSingular"%(self.opt.name)


opt = LibSingularOptions()
opt_verb = LibSingularVerboseOptions()
opt_ctx = LibSingularOptionsContext(opt)
opt_verb_ctx = LibSingularOptionsContext(opt_verb)
