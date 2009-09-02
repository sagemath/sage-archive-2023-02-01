"""
Pythonic Interface to libSingular's options.

AUTHOR:

- Martin Albrecht (2009-08): initial implementation
"""

from sage.libs.singular.decl cimport singular_options

from sage.libs.singular.decl cimport OPT_PROT, OPT_REDSB, OPT_NOT_BUCKETS, OPT_NOT_SUGAR, OPT_SUGARCRIT, OPT_REDTHROUGH
from sage.libs.singular.decl cimport OPT_RETURN_SB, OPT_FASTHC, OPT_OLDSTD, OPT_REDTAIL, OPT_INTSTRATEGY, OPT_NOTREGULARITY
from sage.libs.singular.decl cimport OPT_WEIGHTM, Sy_bit


name_map = {"prot":          Sy_bit(OPT_PROT),
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

class LibSingularOptions:
    """
    Pythonic Interface to libSingular's options.

    Supported options are:

     - returnSB  the functions syz, intersect, quotient, modulo return a
                 standard base instead of a generating set if returnSB is
                 set. This option should not be used for lift.
     - fastHC    tries to the find the highest corner of the staircase (HC)
                 as fast as possible during a standard basis computation
                 (only used for local orderings).
     - intStrategy  avoids division of coefficients during standard basis
                    computations. This option is ring dependent. By default,
                    it is set for rings with characteristic 0 and not set
                    for all other rings.
     - lazy   uses a more lazy approach in std computations, which was used
              in SINGULAR version before 2-0 (and which may lead to faster
              or slower computations, depending on the example)
     - length   select shorter reduceers in std computations,
     - notRegularity  disables the regularity bound for res and mres (see
                      regularity).
     - notSugar disables the sugar strategy during standard basis
                computation.
     - notBuckets  disables the bucket representation of polynomials during
                   standard basis computations. This option usually
                   decreases the memory usage but increases the computation
                   time. It should only be set for memory-critical standard
                   basis computations.
     - oldStd  uses a more lazy approach in std computations, which was
               used in SINGULAR version before 2-0 (and which may lead to
               faster or slower computations, depending on the example)
     - prot   shows protocol information indicating the progress during the
              following computations: facstd, fglm, groebner, lres, mres,
              minres, mstd, res, slimgb, sres, std, stdfglm, stdhilb,
              syz. See below for more details.
     - redSB  computes a reduced standard basis in any standard basis
              computation.
     - redTail reduction of the tails of polynomials during standard basis
               computations. This option is ring dependent. By default, it
               is set for rings with global degree orderings and not set
               for all other rings.
     - redThrough  for inhomogenous input, polynomial reductions during
                   standard basis computations are never postponed, but
                   always finished through. This option is ring
                   dependent. By default, it is set for rings with global
                   degree orderings and not set for all other rings.
     - sugarCrit  uses criteria similar to the homogeneous case to keep
                  more useless pairs.
    - weightM  automatically computes suitable weights for the weighted
               ecart and the weighted sugar method.

    EXAMPLE::

        sage: from sage.libs.singular.option import LibSingularOptions
        sage: libsingular_options = LibSingularOptions()
        sage: libsingular_options
        libSingular options interface
    """
    def __getitem__(self, name):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: libsingular_options = LibSingularOptions()
            sage: libsingular_options['redTail']
            True
        """
        try:
            return bool(singular_options & name_map[name])
        except KeyError:
            raise NameError("libSingular option '%s' unknown."%(name,))

    def __setitem__(self, name, value):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: libsingular_options = LibSingularOptions()
            sage: libsingular_options['redTail']
            True
            sage: libsingular_options['redTail'] = False
            sage: libsingular_options['redTail']
            False
        """
        global singular_options
        try:
            if value:
                singular_options = singular_options | name_map[name]
            else:
                singular_options = singular_options & ~name_map[name]
        except KeyError:
            raise NameError("libSingular option '%s' unknown."%(name,))

    def __int__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: libsingular_options = LibSingularOptions()
            sage: int(libsingular_options)
            67108994
        """
        return int(singular_options)

    def __call__(self, value):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: libsingular_options = LibSingularOptions()
            sage: bck = int(libsingular_options); bck
            67108994
            sage: libsingular_options['redTail'] = True
            sage: int(libsingular_options)
            100663426
            sage: libsingular_options(bck)
            sage: libsingular_options['redTail']
            False
        """
        global singular_options
        singular_options = int(value)

    def __repr__(self):
        """
        EXAMPLE::

            sage: from sage.libs.singular.option import LibSingularOptions
            sage: libsingular_options = LibSingularOptions()
            sage: libsingular_options
            libSingular options interface
        """
        return "libSingular options interface"
