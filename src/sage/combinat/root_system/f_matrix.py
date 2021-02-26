"""
F-Matrix Factory for FusionRings
"""
# ****************************************************************************
#  Copyright (C) 2019 Daniel Bump <bump at match.stanford.edu>
#                     Guillermo Aboumrad <gh_willieab>
#                     Travis Scrimshaw <tcscrims at gmail.com>
#                     Galit Anikeeva <physicstravels@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from functools import cmp_to_key, partial
from itertools import product
import sage.combinat.root_system.fusion_ring as FusionRing
import sage.graphs
from sage.graphs.generators.basic import EmptyGraph
from sage.matrix.constructor import matrix
from sage.misc.misc import inject_variable
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.ideal import Ideal

from copy import deepcopy
#Import pickle for checkpointing and loading certain variables
try:
    import cPickle as pickle
except:
    import pickle
from multiprocessing import cpu_count, Pool, set_start_method, TimeoutError
import numpy as np
from sage.combinat.root_system.fast_parallel_fmats_methods import *
from sage.combinat.root_system.poly_tup_engine import *
from sage.rings.polynomial.polydict import ETuple
from sage.rings.qqbar import AA, QQbar
from sage.rings.real_double import RDF

from itertools import zip_longest

class FMatrix():
    r"""Return an F-Matrix factory for a FusionRing.

    INPUT:

    - ``FR`` -- a FusionRing.

    We only undertake to compute the F-matrix if the
    FusionRing is *multiplicity free* meaning that
    the Fusion coefficients `N^{ij}_k` are bounded
    by 1. For Cartan Types `X_r` and level `k`,
    the multiplicity-free cases are given by the
    following table.

+------------------------+----------+
| Cartan Type            | `k`      |
+========================+==========+
| `A_1`                  | any      |
+------------------------+----------+
| `A_r, r\geq 2`         | `\leq 2` |
+------------------------+----------+
| `B_r, r\geq 2`         | `\leq 2` |
+------------------------+----------+
| `C_2`                  | `\leq 2` |
+------------------------+----------+
| `C_r, r\geq 3`         | `\leq 1` |
+------------------------+----------+
| `D_r, r\geq 4`         | `\leq 2` |
+------------------------+----------+
| `G_2,F_4,E_r`          | `\leq 2` |
+------------------------+----------+

    Beyond this limitation, computation of the F-matrix
    can involve very large systems of equations. A
    rule of thumb is that this code can compute the
    F-matrix for systems with `\leq 4` primary fields,
    with the exception of `G_2` at level `2`.

    The FusionRing and its methods capture much
    of the structure of the underlying tensor category.
    But an important aspect that is not encoded in the
    fusion ring is the associator, which is a homomorphism
    `(A\otimes B)\otimes C\to A\otimes(B\otimes C)`
    requires an additional tool, the F-matrix or 6j-symbol.
    To specify this, we fix a simple object `D`
    and represent the transformation

    .. MATH::

         \text{Hom}(D,(A\otimes B)\otimes C) \to \text{Hom}(D,A\otimes(B\otimes C))

    by a matrix `F^{ABC}_D`. This depends on a pair of
    additional simple objects `X` and `Y`. Indeed, we can
    get a basis for `\text{Hom}(D,(A\otimes B)\otimes C)`
    indexed by simple objects `X` in which the corresponding
    homomorphism factors through `X\otimes C`, and similarly
    `\text{Hom}(D,A\otimes(B\otimes C))` has a basis indexed
    by `Y`, in which the basis vector factors through `A\otimes Y`.

    See [TTWL2009]_ for an introduction to this topic,
    [EGNO2015]_ Section 4.9 for a precise mathematical
    definition and [Bond2007]_ Section 2.5 for a discussion
    of how to compute the F-matrix. In addition to
    [Bond2007]_ worked out F-matrices may be found in
    [RoStWa2009]_ and [CuWa2015]_.

    The F-matrix is only determined up to a *gauge*. This
    is a family of embeddings `C\to A\otimes B` for
    simple objects `A,B,C` such that `\text{Hom}(C,A\otimes B)`
    is nonzero. Changing the gauge changes the F-matrix though
    not in a very essential way. By varying the gauge it is
    possible to make the F-matrices unitary, or it is possible
    to make them cyclotomic.

    Due to the large number of equations we may fail to find a
    Groebner basis if there are too many variables.


    EXAMPLES::

        sage: I=FusionRing("E8",2,conjugate=True)
        sage: I.fusion_labels(["i0","p","s"],inject_variables=True)
        sage: f = FMatrix(I,inject_variables=True); f
        creating variables fx1..fx14
        Defining fx0, fx1, fx2, fx3, fx4, fx5, fx6, fx7, fx8, fx9, fx10, fx11, fx12, fx13
        F-Matrix factory for The Fusion Ring of Type E8 and level 2 with Integer Ring coefficients

    We've exported two sets of variables to the global namespace.
    We created three variables ``i0, p, s`` to represent the
    primary fields (simple elements) of the FusionRing. Creating
    the FMatrix factory also created variables ``fx1,fx2, ... , fx14``
    in order to solve the hexagon and pentagon equations describing
    the F-matrix. Since  we called ``FMatrix`` with the parameter ``inject_variables``
    set true, these have been exported into the global namespace. This
    is not necessary for the code to work but if you want to
    run the code experimentally you may want access to these
    variables.

    EXAMPLES::

        sage: f.fmatrix(s,s,s,s)
        [fx10 fx11]
        [fx12 fx13]

    The F-matrix has not been computed at this stage, so
    the F-matrix `F^{sss}_s` is filled with variables
    ``fx10``, ``fx11``, ``fx12``, ``fx13``. The task is
    to solve for these.

    As explained above The F-matrix `(F^{ABC}_D)_{X,Y}`
    two other variables `X` and `Y`. We have methods to
    tell us (depending on `A,B,C,D`) what the possibilities
    for these are. In this example with `A=B=C=D=s`
    both `X` and `Y` are allowed to be `i_0` or `s`.

    EXAMPLES::

        sage: f.f_from(s,s,s,s), f.f_to(s,s,s,s)
        ([i0, p], [i0, p])

    The last two statments show that the possible values of
    `X` and `Y` when `A=B=C=D=s` are `i_0` and `p`.

    The F-matrix is computed by solving the so-called
    pentagon and hexagon equations. The *pentagon
    equations* reflect the Mac Lane pentagon axiom in the
    definition of a monoidal category. The hexagon relations
    reflect the axioms of a *braided monoidal category*,
    which are constraints on both the F-matrix and on
    the R-matrix. Finally, orthogonality constraints
    may be imposed to obtain a unitary F-matrix.

    EXAMPLES::

        sage: f.get_pentagons()[1:3]
        [fx1*fx5 - fx7^2, fx5*fx8*fx13 - fx2*fx12]
        sage: f.get_hexagons()[1:3]
        [fx10*fx12 + (-zeta128^32)*fx12*fx13 + (-zeta128^16)*fx12, fx0 - 1]
        sage: f.get_orthogonality_constraints()[1:3]
        [fx1^2 - 1, fx2^2 - 1]

    You may solve these 41+14=55 equations to compute the F-matrix.

    EXAMPLES::

        sage: f.get_solution(output=True)
        Setting up hexagons and pentagons...
        equations: 14
        equations: 37
        Finding a Groebner basis...
        Solving...
        Fixing the gauge...
        adding equation... fx1 - 1
        adding equation... fx11 - 1
        Done!
        {(s, s, s, s, i0, i0): (-1/2*zeta128^48 + 1/2*zeta128^16),
         (s, s, s, s, i0, p): 1,
         (s, s, s, s, p, i0): 1/2,
         (s, s, s, s, p, p): (1/2*zeta128^48 - 1/2*zeta128^16),
         (s, s, p, i0, p, s): (-1/2*zeta128^48 + 1/2*zeta128^16),
         (s, s, p, p, i0, s): (-zeta128^48 + zeta128^16),
         (s, p, s, i0, s, s): 1,
         (s, p, s, p, s, s): -1,
         (s, p, p, s, s, i0): 1,
         (p, s, s, i0, s, p): (-zeta128^48 + zeta128^16),
         (p, s, s, p, s, i0): (-1/2*zeta128^48 + 1/2*zeta128^16),
         (p, s, p, s, s, s): -1,
         (p, p, s, s, i0, s): 1,
         (p, p, p, p, i0, i0): 1}

    We now have access to the values of the F-mstrix using
    the methods :meth:`fmatrix` and :meth:`fmat`.

    EXAMPLES::

        sage: f.fmatrix(s,s,s,s)
        [(-1/2*zeta128^48 + 1/2*zeta128^16)                                  1]
        [                               1/2  (1/2*zeta128^48 - 1/2*zeta128^16)]
        sage: f.fmat(s,s,s,s,p,p)
        (1/2*zeta128^48 - 1/2*zeta128^16)

    """
    def __init__(self, fusion_ring, fusion_label="f", var_prefix='fx', inject_variables=False):
        self._FR = fusion_ring
        if inject_variables:
            self._FR.fusion_labels(fusion_label, inject_variables=True)
        #Set up F-symbols entry by entry
        n_vars = self.findcases()
        self._poly_ring = PolynomialRing(self._FR.field(),n_vars,var_prefix)
        if inject_variables:
            print ("creating variables %s%s..%s%s"%(var_prefix,1,var_prefix,n_vars))
            self._poly_ring.inject_variables()
            # for i in range(self._poly_ring.ngens()):
            #     inject_variable("%s%s"%(var_prefix,i),self._poly_ring.gens()[i])
        self._var_to_sextuple, self._fvars = self.findcases(output=True)
        self._singles = self.singletons()

        #Initialize list of defining equations
        self.ideal_basis = list()

        #Initialize empty set of solved F-symbols
        self.solved = set()

        #New attributes of the FMatrix class
        self._field = self._FR.field()
        r = self._field.defining_polynomial().roots(ring=QQbar,multiplicities=False)[0]
        self._qqbar_embedding = self._field.hom([r],QQbar)
        self._var_to_idx = { var : idx for idx, var in enumerate(self._poly_ring.gens()) }
        self._ks = self.get_known_sq()
        self._nnz = self.get_known_nonz()
        self._known_vals = dict()
        self.symbols_known = False

        #Multiprocessing attributes
        self.mp_thresh = 10000


    #######################
    ### Class utilities ###
    #######################

    def __repr__(self):
        """
        EXAMPLES::

            sage: FMatrix(FusionRing("B2",1))
            F-Matrix factory for The Fusion Ring of Type B2 and level 1 with Integer Ring coefficients
        """
        return "F-Matrix factory for %s"%self._FR

    def remaining_vars(self):
        """
        Return a list of unknown F-symbols (reflects current stage of computation)
        """
        return [var for var in self._poly_ring.gens() if var not in self.solved]

    def clear_equations(self):
        """
        Clear the set of equations to be solved.
        """
        self.ideal_basis = list()

    def clear_vars(self):
        """
        Clear the set of variables.
        """
        self._fvars = { self._var_to_sextuple[key] : key for key in self._var_to_sextuple }
        self.solved = set()

    def fmat(self, a, b, c, d, x, y, data=True):
        """
        Return the F-Matrix coefficient `(F^{a,b,c}_d)_{x,y}`

        EXAMPLES::

            sage: f=FMatrix(FusionRing("G2",1),fusion_label=["i0","t"])
            sage: [f.fmat(t,t,t,t,x,y) for x in f.FR.basis() for y in f.FR.basis()]
            [fx1, fx2, fx3, fx4]
            sage: f.get_solution(output=True)
            Setting up hexagons and pentagons...
            equations: 5
            equations: 10
            Finding a Groebner basis...
            Solving...
            Fixing the gauge...
            adding equation... fx2 - 1
            Done!
            {(t, t, t, i0, t, t): 1,
            (t, t, t, t, i0, i0): (-zeta60^14 + zeta60^6 + zeta60^4 - 1),
            (t, t, t, t, i0, t): 1,
            (t, t, t, t, t, i0): (-zeta60^14 + zeta60^6 + zeta60^4 - 1),
            (t, t, t, t, t, t): (zeta60^14 - zeta60^6 - zeta60^4 + 1)}
            sage: [f.fmat(t,t,t,t,x,y) for x in f.FR.basis() for y in f.FR.basis()]
            [(-zeta60^14 + zeta60^6 + zeta60^4 - 1),
            1,
            (-zeta60^14 + zeta60^6 + zeta60^4 - 1),
            (zeta60^14 - zeta60^6 - zeta60^4 + 1)]

        """
        if self._FR.Nk_ij(a,b,x) == 0 or self._FR.Nk_ij(x,c,d) == 0 or self._FR.Nk_ij(b,c,y) == 0 or self._FR.Nk_ij(a,y,d) == 0:
            return 0

        #Some known zero F-symbols
        if a == self._FR.one():
            if x == b and y == d:
                return 1
            else:
                return 0
        if b == self._FR.one():
            if x == a and y == c:
                return 1
            else:
                return 0
        if c == self._FR.one():
            if x == d and y == b:
                return 1
            else:
                return 0
        if data:
            #Better to use try/except for speed. Somewhat trivial, but worth
            #hours when method is called ~10^11 times
            try:
                return self._fvars[a,b,c,d,x,y]
            except KeyError:
                return 0
        else:
            return (a,b,c,d,x,y)

    def fmatrix(self,a,b,c,d):
        """
        Return the F-Matrix `F^{a,b,c}_d`.

        INPUT:

        - ``a,b,c,d`` -- basis elements of the FusionRing

        EXAMPLES::

            sage: f=FMatrix(FusionRing("A1",2),fusion_label="c")
            sage: f.get_solution(verbose=False);
            equations: 14
            equations: 37
            adding equation... fx4 - 1
            adding equation... fx10 - 1
            sage: f.f_from(c1,c1,c1,c1)
            [c0, c2]
            sage: f.f_to(c1,c1,c1,c1)
            [c0, c2]
            sage: f.fmatrix(c1,c1,c1,c1)
            [ (1/2*zeta32^12 - 1/2*zeta32^4) (-1/2*zeta32^12 + 1/2*zeta32^4)]
            [ (1/2*zeta32^12 - 1/2*zeta32^4)  (1/2*zeta32^12 - 1/2*zeta32^4)]
        """

        X = self.f_from(a,b,c,d)
        Y = self.f_to(a,b,c,d)
        return matrix([[self.fmat(a,b,c,d,x,y) for y in Y] for x in X])

    def field(self):
        return self._field

    def FR(self):
        return self._FR

    def findcases(self,output=False):
        """
        Return unknown F-matrix entries. If run with output=True,
        this returns two dictionaries; otherwise it just returns the
        number of unknown values.

        EXAMPLES::

            sage: f=FMatrix(FusionRing("G2",1),fusion_label=["i0","t"])
            sage: f.findcases()
            5
            sage: f.findcases(output=True)
            ({fx4: (t, t, t, t, t, t),
            fx3: (t, t, t, t, t, i0),
            fx2: (t, t, t, t, i0, t),
            fx1: (t, t, t, t, i0, i0),
            fx0: (t, t, t, i0, t, t)},
            {(t, t, t, i0, t, t): fx0,
            (t, t, t, t, i0, i0): fx1,
            (t, t, t, t, i0, t): fx2,
            (t, t, t, t, t, i0): fx3,
            (t, t, t, t, t, t): fx4})

        """
        i = 0
        if output:
            idx_map = dict()
            ret = dict()
        for (a,b,c,d) in list(product(self._FR.basis(), repeat=4)):
            for x in self.f_from(a, b, c, d):
                for y in self.f_to(a, b, c, d):
                    fm = self.fmat(a, b, c, d, x, y, data=False)
                    if fm is not None and fm not in [0,1]:
                        if output:
                            v = self._poly_ring.gens()[i]
                            ret[(a,b,c,d,x,y)] = v
                            idx_map[v] = (a, b, c, d, x, y)
                        i += 1
        if output:
            return idx_map, ret
        else:
            return i

    def singletons(self):
        """
        Find x_i that are automatically nonzero, because their F-matrix is 1x1
        """
        ret = []
        for (a, b, c, d) in list(product(self._FR.basis(), repeat=4)):
            (ff,ft) = (self.f_from(a,b,c,d),self.f_to(a,b,c,d))
            if len(ff) == 1 and len(ft) == 1:
                v = self._fvars.get((a,b,c,d,ff[0],ft[0]), None)
                if v in self._poly_ring.gens():
                    ret.append(v)
        return ret

    def f_from(self,a,b,c,d):
        r"""
        Return the possible `x` such that there are morphisms
        `d\to x\otimes c\to (a\otimes b)\otimes c`.

        INPUT:

        - ``a,b,c,d`` -- basis elements of the FusionRing.

        EXAMPLES::

            sage: f=FMatrix(FusionRing("A1",3),fusion_label="a")
            sage: f.fmatrix(a1,a1,a2,a2)
            [fx6 fx7]
            [fx8 fx9]
            sage: f.f_from(a1,a1,a2,a2)
            [a0, a2]
            sage: f.f_to(a1,a1,a2,a2)
            [a1, a3]
        """

        return [x for x in self._FR.basis() if self._FR.Nk_ij(a,b,x) != 0 and self._FR.Nk_ij(x,c,d) != 0]

    def f_to(self,a,b,c,d):
        r"""
        Return the possible `y` such that there are morphisms
        `d\to a\otimes y\to a\otimes(b\otimes c)`.

        INPUT:

        - ``a,b,c,d`` -- basis elements of the FusionRing.

        EXAMPLES::

            sage: B=FMatrix(FusionRing("B2",2),fusion_label="b")
            sage: B.fmatrix(b2,b4,b4,b2)
            [fx278 fx279 fx280]
            [fx281 fx282 fx283]
            [fx284 fx285 fx286]
            sage: B.f_from(b2,b4,b4,b2)
            [b1, b3, b5]
            sage: B.f_to(b2,b4,b4,b2)
            [b0, b1, b5]

        """

        return [y for y in self._FR.basis() if self._FR.Nk_ij(b,c,y) != 0 and self._FR.Nk_ij(a,y,d) != 0]

    ####################
    ### Data getters ###
    ####################

    def get_fmats_in_alg_field(self):
        """
        Return F-symbols as elements of the AlgebraicField. This method uses
        self._qqbar_embedding to coerce F-symbols into QQbar.
        """
        return { sextuple : self._qqbar_embedding(fvar) for sextuple, fvar in self._fvars.items() }

    def get_radical_expression(self):
        """
        Return radical expression of F-symbols for easy visualization
        """
        return { sextuple : val.radical_expression() for sextuple, val in get_fmats_in_alg_field().items() }

    def get_known_vals(self):
        """
        Construct a dictionary of idx, known_val pairs for equation substitution
        """
        return { var_idx : self._fvars[self._var_to_sextuple[self._poly_ring.gen(var_idx)]] for var_idx in self.solved }

    def get_known_sq(self,eqns=None):
        """
        Construct a dictionary of known squares. Keys are variable indices and corresponding values are the squares
        """
        if eqns is None:
            eqns = self.ideal_basis
        return { variables(eq_tup)[0] : -eq_tup[-1][1] for eq_tup in eqns if tup_fixes_sq(eq_tup) }

    def get_known_nonz(self):
        """
        Construct an ETuple indicating positions of known nonzero variables.
        MUST be called after self._ks = get_known_sq()
        """
        nonz = { self._var_to_idx[var] : 100 for var in self._singles }
        for idx in self._ks:
            nonz[idx] = 100
        return ETuple(nonz, self._poly_ring.ngens())

    def get_qqbar_embedding(self):
        """
        Return an embedding from the base field containing F-symbols (the
        FusionRing's CyclotomicField, a NumberField, or QQbar) into QQbar
        """
        return self._qqbar_embedding

    def get_coerce_map_from_fr_cyclotomic_field(self):
        """
        Return a coercion map from the FusionRing's cyclotomic field into the
        base field containing all F-symbols (this could be the FusionRing's
        CyclotomicField, a NumberField, or QQbar).
        """
        #If base field is different from associated FusionRing's CyclotomicField,
        #return coercion map
        try:
            return self._coerce_map_from_cyc_field
        #Otherwise, return identity map CyclotomicField <-> CyclotomicField
        except AttributeError:
            F = self._FR.field()
            return F.hom([F.gen()], F)

    #########################
    ### Useful predicates ###
    #########################

    def is_univariate_in_unknown(self,monom_exp):
        """
        Determine if monomial exponent is univariate in a still unknown F-symbol
        """
        return len(monom_exp.nonzero_values()) == 1 and monom_exp.nonzero_positions()[0] not in self.solved

    def is_uni_linear_in_unkwown(self,monom_exp):
        """
        Determine if monomial exponent is univariate and linear in a vstill unknown F-symbol
        """
        return monom_exp.nonzero_values() == [1] and monom_exp.nonzero_positions()[0] not in self.solved

    ##############################
    ### Variables partitioning ###
    ##############################

    def largest_fmat_size(self):
        """
        Get the size of the largest F-matrix F^{abc}_d
        """
        return max(self.fmatrix(*tup).nrows() for tup in product(self._FR.basis(),repeat=4))

    def get_fmats_by_size(self,n):
        """
        Partition the F-symbols according to the size of the F-matrix F^{abc}_d
        they belong to
        """
        fvars_copy = deepcopy(self._fvars)
        solved_copy = deepcopy(self.solved)
        self.clear_vars()
        var_set = set()
        for quadruple in product(self._FR.basis(),repeat=4):
            F = self.fmatrix(*quadruple)
            #Discard trivial 1x1 F-matrix, if applicable
            if F.nrows() == n and F.coefficients() != [1]:
                var_set.update(self._var_to_idx[fx] for fx in F.coefficients())
        self._fvars = fvars_copy
        self.solved = solved_copy
        return var_set

    ############################
    ### Checkpoint utilities ###
    ############################

    def get_fr_str(self):
        """
        Auto-generate filename string for saving results
        """
        return self._FR.cartan_type()[0] + str(self._FR.cartan_type()[1]) + str(self._FR.fusion_level())

    def save_fvars(self,filename):
        """
        Save current variables state
        """
        with open(filename, 'wb') as f:
            pickle.dump([self._fvars, self.solved], f)

    def load_fvars(self,save_dir=""):
        """
        If provided, optional param save_dir should have a trailing forward slash
        """
        with open(save_dir + "saved_fvars_" + self.get_fr_str() + ".pickle", 'rb') as f:
            self._fvars, self.solved = pickle.load(f)

    #################
    ### MapReduce ###
    #################

    def map_triv_reduce(self,mapper,input_iter,worker_pool=None,chunksize=None,mp_thresh=None):
        """
        Apply the given mapper to each element of the given input iterable and
        return the results (with no duplicates) in a list. This method applies the
        mapper in parallel if a worker_pool is provided.

        # INPUT:
        mapper is a string specifying the name of a function defined in the
        fast_parallel_fmats_methods module.

        ##NOTES:
        If worker_pool is not provided, function maps and reduces on a single process.
        If worker_pool is provided, the function attempts to determine whether it should
        use multiprocessing based on the length of the input iterable. If it can't determine
        the length of the input iterable then it uses multiprocessing with the default chunksize of 1
        if chunksize is not explicitly provided.
        """
        if mp_thresh is None:
          mp_thresh = self.mp_thresh
        #Compute multiprocessing parameters
        if worker_pool is not None:
            try:
                n = len(input_iter)
            except:
                n = mp_thresh + 1
            if chunksize is None:
                chunksize = n // (worker_pool._processes**2) + 1
        no_mp = worker_pool is None or n < mp_thresh
        #Map phase. Casting Async Object blocks execution... Each process holds results
        #in its copy of fmats.temp_eqns
        input_iter = zip_longest([],input_iter,fillvalue=(mapper,id(self)))
        if no_mp:
            list(map(executor,input_iter))
        else:
            list(worker_pool.imap_unordered(executor,input_iter,chunksize=chunksize))
        #Reduce phase
        if no_mp:
            results = collect_eqns(0)
        else:
            results = set()
            for child_eqns in worker_pool.imap_unordered(collect_eqns,range(worker_pool._processes)):
                results.update(child_eqns)
            results = list(results)
        return results

    ########################
    ### Equations set up ###
    ########################

    def get_orthogonality_constraints(self,output=True):
        """
        Get equations imposed on the F-matrix by orthogonality.

        If output=True, equations are returned as polynomial objects. Otherwise,
        polynomial generators (stored in the internal tuple representation) are
        appended to self.ideal_basis
        """
        eqns = list()
        for tup in product(self._FR.basis(), repeat=4):
            mat = self.fmatrix(*tup)
            eqns.extend((mat.T * mat - matrix.identity(mat.nrows())).coefficients())
        if output:
            return eqns
        self.ideal_basis.extend([poly_to_tup(eq) for eq in eqns])

    def get_defining_equations(self,option,worker_pool=None,output=True):
        """
        Get the equations defining the ideal generated by the hexagon or
        pentagon relations.

        Use option='hexagons' to get equations imposed on the F-matrix by the hexagon
        relations in the definition of a braided category.

        Use option='pentagons' to get equations imposed on the F-matrix by the pentagon
        relations in the definition of a monoidal category.

        If output=True, equations are returned as polynomial objects. Otherwise,
        polynomial generators (stored in the internal tuple representation) are
        appended to self.ideal_basis
        If a worker_pool is passed, then we use multiprocessing
        """
        n_proc = worker_pool._processes if worker_pool is not None else 1
        params = [(child_id, n_proc) for child_id in range(n_proc)]
        eqns = self.map_triv_reduce('get_reduced_'+option,params,worker_pool=worker_pool,chunksize=1,mp_thresh=0)
        if output:
            return [self.tup_to_fpoly(p) for p in eqns]
        self.ideal_basis.extend(eqns)

    ############################
    ### Equations processing ###
    ############################

    def tup_to_fpoly(self,eq_tup):
        """
        Assemble a polynomial object from its tuple representation
        """
        return tup_to_poly(eq_tup,parent=self._poly_ring)

    def solve_for_linear_terms(self,eqns=None):
        """
        Solve for a linear term occurring in a two-term equation.
        """
        if eqns is None:
            eqns = self.ideal_basis

        linear_terms_exist = False
        for eq_tup in eqns:
            if len(eq_tup) == 1:
                m = eq_tup[0][0]
                if self.is_univariate_in_unknown(m):
                    var = m.nonzero_positions()[0]
                    self._fvars[self._var_to_sextuple[self._poly_ring.gen(var)]] = tuple()
                    self.solved.add(var)
                    linear_terms_exist = True
            if len(eq_tup) == 2:
                monomials = [m for m, c in eq_tup]
                max_var = monomials[0].emax(monomials[1]).nonzero_positions()[0]
                for this, m in enumerate(monomials):
                    other = (this+1)%2
                    if self.is_uni_linear_in_unkwown(m) and m.nonzero_positions()[0] == max_var and monomials[other][m.nonzero_positions()[0]] == 0:
                        var = m.nonzero_positions()[0]
                        rhs_key = monomials[other]
                        rhs_coeff = -eq_tup[other][1] / eq_tup[this][1]
                        self._fvars[self._var_to_sextuple[self._poly_ring.gen(var)]] = ((rhs_key,rhs_coeff),)
                        self.solved.add(var)
                        linear_terms_exist = True
        return linear_terms_exist

    def backward_subs(self):
        """
        Backward substitution step. Traverse variables in reverse lexicographical order.
        """
        for var in reversed(self._poly_ring.gens()):
            sextuple = self._var_to_sextuple[var]
            rhs = self._fvars[sextuple]
            d = { var_idx : self._fvars[self._var_to_sextuple[self._poly_ring.gen(var_idx)]] for var_idx in variables(rhs) if var_idx in self.solved }
            if d:
                kp = compute_known_powers(get_variables_degrees([rhs]), d)
                self._fvars[sextuple] = tuple(subs_squares(subs(rhs,kp),self._ks).items())

    def update_reduction_params(self,eqns=None,worker_pool=None,children_need_update=False):
        """
        Update reduction parameters in all processes
        """
        if eqns is None:
            eqns = self.ideal_basis
        self._ks, self._var_degs = self.get_known_sq(eqns), get_variables_degrees(eqns)
        self._nnz = self.get_known_nonz()
        self._kp = compute_known_powers(self._var_degs,self.get_known_vals())
        if worker_pool is not None and children_need_update:
            #self._nnz and self._kp are computed in child processes to reduce IPC overhead
            n_proc = worker_pool._processes
            new_data = [(self._fvars,self.solved,self._ks,self._var_degs)]*n_proc
            self.map_triv_reduce('update_child_fmats',new_data,worker_pool=worker_pool,chunksize=1,mp_thresh=0)

    def triangular_elim(self,required_vars=None,eqns=None,worker_pool=None,verbose=True):
        """
        Perform triangular elimination of linear terms in two-term equations until no such terms exist
        For optimal usage of TRIANGULAR elimination, pass in a SORTED list of equations
        """
        ret = True
        if eqns is None:
            eqns = self.ideal_basis
            ret = False
        if required_vars is None:
            required_vars = self._poly_ring.gens()
        poly_sortkey = cmp_to_key(poly_tup_cmp)

        #Testing new _fvars representation... Unzip polynomials
        self._fvars = { sextuple : poly_to_tup(rhs) for sextuple, rhs in self._fvars.items() }

        while True:
            linear_terms_exist = self.solve_for_linear_terms(eqns)
            if not linear_terms_exist: break
            self.backward_subs()

            #Support early termination in case only some F-symbols are needed
            req_vars_known = all(self._fvars[self._var_to_sextuple[var]] in self._FR.field() for var in required_vars)
            if req_vars_known: return 1

            #Compute new reduction params, send to child processes if any, and update eqns
            self.update_reduction_params(eqns=eqns,worker_pool=worker_pool,children_need_update=len(eqns)>self.mp_thresh)
            eqns = sorted(self.map_triv_reduce('update_reduce',eqns,worker_pool=worker_pool), key=poly_sortkey)
            if verbose:
                print("Elimination epoch completed... {} eqns remain in ideal basis".format(len(eqns)))

        #Zip up _fvars before exiting
        self._fvars = { sextuple : self.tup_to_fpoly(rhs_tup) for sextuple, rhs_tup in self._fvars.items() }
        if ret:
            return eqns
        self.ideal_basis = eqns

    #####################
    ### Graph methods ###
    #####################

    def equations_graph(self,eqns=None):
        """
        Construct a graph corresponding to equations. The nodes in the graph
        are indices corresponding to variables in the equations and two
        nodes are connected if the corresponding variables appear together in
        a given equation.

        If no list of equations is passed, the graph is built from equations in
        self.ideal_basis
        """
        if eqns is None:
            eqns = self.ideal_basis

        G = sage.graphs.generators.basic.EmptyGraph()
        G.add_vertices([x for eq_tup in eqns for x in variables(eq_tup)])
        for eq_tup in eqns:
            s = [v for v in variables(eq_tup)]
            for x in s:
                for y in s:
                    if y!=x:
                        G.add_edge(x,y)
        return G

    def partition_eqns(self,graph,eqns=None,verbose=True):
        """
        Partition equations corresponding to edges in a disconnected graph
        """
        if eqns is None:
            eqns = self.ideal_basis
        partition = { tuple(c) : [] for c in graph.connected_components() }
        for eq_tup in eqns:
            partition[tuple(graph.connected_component_containing_vertex(variables(eq_tup)[0]))].append(eq_tup)
        if verbose:
            print("Partitioned {} equations into {} components of size:".format(len(eqns),len(graph.connected_components())))
            print(graph.connected_components_sizes())
        return partition

    def par_graph_gb(self,worker_pool=None,eqns=None,term_order="degrevlex",verbose=True):
        """
        Compute a Groebner basis for a set of equations partitioned according to their corresponding graph
        """
        if eqns is None: eqns = self.ideal_basis
        graph = self.equations_graph(eqns)
        small_comps = list()
        temp_eqns = list()

        # #For informative print statement
        # nmax = self.largest_fmat_size()
        # vars_by_size = list()
        # for i in range(nmax+1):
        #     vars_by_size.append(self.get_fmats_by_size(i))

        for comp, comp_eqns in self.partition_eqns(graph,verbose=verbose).items():
            #Check if component is too large to process
            if len(comp) > 60:
                # fmat_size = 0
                # #For informative print statement
                # for i in range(1,nmax+1):
                #     if set(comp).issubset(vars_by_size[i]):
                #         fmat_size = i
                # print("Component of size {} with vars in F-mats of size {} is too large to find GB".format(len(comp),fmat_size))
                temp_eqns.extend(comp_eqns)
            else:
                small_comps.append(comp_eqns)
        input_iter = zip_longest(small_comps,[],fillvalue=term_order)
        small_comp_gb = self.map_triv_reduce('compute_gb',input_iter,worker_pool=worker_pool,chunksize=1,mp_thresh=50)
        ret = small_comp_gb + temp_eqns
        return ret

    def get_component_variety(self,var,eqns):
        """
        Translate equations in each connected component to smaller polynomial rings
        so we can call built-in variety method.
        """
        #Define smaller poly ring in component vars
        R = PolynomialRing(self._FR.field(),len(var),'a',order='lex')

        #Zip tuples into R and compute Groebner basis
        idx_map = { old : new for new, old in enumerate(sorted(var)) }
        nvars = len(var)
        polys = [tup_to_poly(resize(eq_tup,idx_map,nvars),parent=R) for eq_tup in eqns]
        var_in_R = Ideal(sorted(polys)).variety(ring=AA)

        #Change back to fmats poly ring and append to temp_eqns
        inv_idx_map = { v : k for k, v in idx_map.items() }
        return [{ inv_idx_map[i] : value for i, (key, value) in enumerate(sorted(soln.items())) } for soln in var_in_R]

    #######################
    ### Solution method ###
    #######################

    def og_get_numeric_solution(self,eqns=None,verbose=True):
        """
        Get a numeric solution for given equations by using the partitioning the equations
        graph and selecting an arbitrary point on the discrete degrees-of-freedom variety,
        which is computed as a Cartesian product
        """
        if eqns is None:
            eqns = self.ideal_basis
        eqns_partition = self.partition_eqns(self.equations_graph(eqns),verbose=verbose)

        x = self._FR.field()['x'].gen()
        numeric_fvars = dict()
        non_cyclotomic_roots = list()
        must_change_base_field = False
        phi = self._FR.field().hom([self._field.gen()],self._FR.field())
        for comp, part in eqns_partition.items():
            #If component have only one equation in a single variable, get a root
            if len(comp) == 1 and len(part) == 1:
                #Attempt to find cyclotomic root
                univ_poly = tup_to_univ_poly(part[0],x)
                real_roots = univ_poly.roots(ring=AA,multiplicities=False)
                assert real_roots, "No real solution exists... {} has no real roots".format(univ_poly)
                roots = univ_poly.roots(multiplicities=False)
                if roots:
                    numeric_fvars[comp[0]] = roots[0]
                else:
                    non_cyclotomic_roots.append((comp[0],real_roots[0]))
                    must_change_base_field = True
            #Otherwise, compute the component variety and select a point to obtain a numerical solution
            else:
                sols = self.get_component_variety(comp,part)
                assert len(sols) > 1, "No real solution exists... component with variables {} has no real points".format(comp)
                for fx, rhs in sols[0].items():
                    non_cyclotomic_roots.append((fx,rhs))
                must_change_base_field = True

        if must_change_base_field:
            #If needed, find a NumberField containing all the roots
            roots = [self._FR.field().gen()]+[r[1] for r in non_cyclotomic_roots]
            self._field, nf_elts, self._qqbar_embedding = number_field_elements_from_algebraics(roots,minimal=True)
            #Embed cyclotomic field into newly constructed NumberField
            cyc_gen_as_nf_elt = nf_elts.pop(0)
            phi = self._FR.field().hom([cyc_gen_as_nf_elt], self._field)
            numeric_fvars = { k : phi(v) for k, v in numeric_fvars.items() }
            for i, elt in enumerate(nf_elts):
                numeric_fvars[non_cyclotomic_roots[i][0]] = elt

            #Do some appropriate conversions
            new_poly_ring = self._poly_ring.change_ring(self._field)
            nvars = self._poly_ring.ngens()
            self._var_to_idx = { new_poly_ring.gen(i) : i for i in range(nvars) }
            self._var_to_sextuple = { new_poly_ring.gen(i) : self._var_to_sextuple[self._poly_ring.gen(i)] for i in range(nvars) }
            self._poly_ring = new_poly_ring

        #Ensure all F-symbols are known
        self.solved.update(numeric_fvars)
        nvars = self._poly_ring.ngens()
        assert len(self.solved) == nvars, "Some F-symbols are still missing...{}".format([self._poly_ring.gen(fx) for fx in set(range(nvars)).difference(self.solved)])

        #Backward substitution step. Traverse variables in reverse lexicographical order. (System is in triangular form)
        self._fvars = { sextuple : apply_coeff_map(poly_to_tup(rhs),phi) for sextuple, rhs in self._fvars.items() }
        for fx, rhs in numeric_fvars.items():
            self._fvars[self._var_to_sextuple[self._poly_ring.gen(fx)]] = ((ETuple({},nvars),rhs),)
        self.backward_subs()
        self._fvars = { sextuple : constant_coeff(rhs) for sextuple, rhs in self._fvars.items() }
        self.clear_equations()

    def get_numeric_solution(self,eqns=None,verbose=True):
        """
        Get a numeric solution for given equations by using the partitioning the equations
        graph and selecting an arbitrary point on the discrete degrees-of-freedom variety,
        which is computed as a Cartesian product
        """
        if eqns is None:
            eqns = self.ideal_basis
        eqns_partition = self.partition_eqns(self.equations_graph(eqns),verbose=verbose)

        F = self._field
        x = F['x'].gen()
        numeric_fvars = dict()
        non_cyclotomic_roots = list()
        must_change_base_field = False
        phi = F.hom([F.gen()],F)
        for comp, part in eqns_partition.items():
            #If component have only one equation in a single variable, get a root
            if len(comp) == 1 and len(part) == 1:
                #Attempt to find cyclotomic root
                univ_poly = tup_to_univ_poly(part[0],x)
                real_roots = univ_poly.roots(ring=AA,multiplicities=False)
                assert real_roots, "No real solution exists... {} has no real roots".format(univ_poly)
                roots = univ_poly.roots(multiplicities=False)
                if roots:
                    numeric_fvars[comp[0]] = roots[0]
                else:
                    non_cyclotomic_roots.append((comp[0],real_roots[0]))
                    must_change_base_field = True
            #Otherwise, compute the component variety and select a point to obtain a numerical solution
            else:
                sols = self.get_component_variety(comp,part)
                assert len(sols) > 1, "No real solution exists... component with variables {} has no real points".format(comp)
                for fx, rhs in sols[0].items():
                    non_cyclotomic_roots.append((fx,rhs))
                must_change_base_field = True

        if must_change_base_field:
            #Attempt to compute smallest number field containing all the F-symbols
            #If calculation takes too long, we use QQbar as the base field
            proc = Pool(1)
            input_args = ((('get_appropriate_number_field',id(self)),non_cyclotomic_roots),)
            p = proc.apply_async(executor,input_args)
            try:
                self._field, bf_elts, self._qqbar_embedding = p.get(timeout=100)
            except TimeoutError:
                self._field = QQbar
                bf_elts = [self._qqbar_embedding(F.gen())]
                bf_elts += [rhs for fx,rhs in non_cyclotomic_roots]
                self._qqbar_embedding = lambda x : x

            #Embed cyclotomic field into newly constructed base field
            cyc_gen_as_bf_elt = bf_elts.pop(0)
            phi = self._FR.field().hom([cyc_gen_as_bf_elt], self._field)
            self._coerce_map_from_cyc_field = phi
            numeric_fvars = { k : phi(v) for k, v in numeric_fvars.items() }
            for i, elt in enumerate(bf_elts):
                numeric_fvars[non_cyclotomic_roots[i][0]] = elt

            #Do some appropriate conversions
            new_poly_ring = self._poly_ring.change_ring(self._field)
            nvars = self._poly_ring.ngens()
            self._var_to_idx = { new_poly_ring.gen(i) : i for i in range(nvars) }
            self._var_to_sextuple = { new_poly_ring.gen(i) : self._var_to_sextuple[self._poly_ring.gen(i)] for i in range(nvars) }
            self._poly_ring = new_poly_ring

        #Ensure all F-symbols are known
        self.solved.update(numeric_fvars)
        nvars = self._poly_ring.ngens()
        assert len(self.solved) == nvars, "Some F-symbols are still missing...{}".format([self._poly_ring.gen(fx) for fx in set(range(nvars)).difference(self.solved)])

        #Backward substitution step. Traverse variables in reverse lexicographical order. (System is in triangular form)
        self._fvars = { sextuple : apply_coeff_map(poly_to_tup(rhs),phi) for sextuple, rhs in self._fvars.items() }
        for fx, rhs in numeric_fvars.items():
            self._fvars[self._var_to_sextuple[self._poly_ring.gen(fx)]] = ((ETuple({},nvars),rhs),)
        self.backward_subs()
        self._fvars = { sextuple : constant_coeff(rhs) for sextuple, rhs in self._fvars.items() }
        self.clear_equations()

    def find_real_orthogonal_solution(self,use_mp=True,save_dir="",verbose=True):
        """
        Solver
        If provided, optional param save_dir (for saving computed _fvars) should have a trailing forward slash
        Supports "warm" start. Use load_fvars to re-start computation from checkpoint
        """
        #Set multiprocessing parameters. Context can only be set once, so we try to set it
        try:
            set_start_method('fork')
        except RuntimeError:
            pass
        pool = Pool(processes=max(cpu_count()-1,1)) if use_mp else None
        print("Computing F-symbols for {} with {} variables...".format(self._FR, len(self._fvars)))

        #Set up hexagon equations and orthogonality constraints
        poly_sortkey = cmp_to_key(poly_tup_cmp)
        self.get_orthogonality_constraints(output=False)
        self.get_defining_equations('hexagons',worker_pool=pool,output=False)
        if verbose:
            print("Set up {} hex and orthogonality constraints...".format(len(self.ideal_basis)))

        #Set up equations graph. Find GB for each component in parallel. Eliminate variables
        self.ideal_basis = sorted(self.par_graph_gb(worker_pool=pool), key=poly_sortkey)
        print("GB is of length",len(self.ideal_basis))
        self.triangular_elim(worker_pool=pool,verbose=verbose)

        #Report progress and checkpoint!
        if verbose:
            print("Hex elim step solved for {} / {} variables".format(len(self.solved), len(self._poly_ring.gens())))
        filename = save_dir + "saved_fvars_" + self.get_fr_str() + ".pickle"
        self.save_fvars(filename)

        #Update reduction parameters, also in children if any
        self.update_reduction_params(worker_pool=pool,children_need_update=True)

        #Set up pentagon equations in parallel, simplify, and eliminate variables
        self.get_defining_equations('pentagons',worker_pool=pool,output=False)
        if verbose:
            print("Set up {} reduced pentagons...".format(len(self.ideal_basis)))
        self.ideal_basis = sorted(self.ideal_basis, key=poly_sortkey)
        self.triangular_elim(worker_pool=pool,verbose=verbose)

        #Report progress and checkpoint!
        if verbose:
            print("Pent elim step solved for {} / {} variables".format(len(self.solved), len(self._poly_ring.gens())))
        self.save_fvars(filename)

        #Close worker pool to free resources
        if pool is not None: pool.close()

        #Set up new equations graph and compute variety for each component
        self.ideal_basis = sorted(self.par_graph_gb(term_order="lex"), key=poly_sortkey)
        self.triangular_elim(verbose=verbose)
        self.get_numeric_solution(verbose=verbose)
        self.save_fvars(filename)
        self.symbols_known = True
        self._FR._coer = self.get_coerce_map_from_fr_cyclotomic_field()

    #########################
    ### Cyclotomic method ###
    #########################

    def fix_gauge(self, algorithm=''):
        """
        Fix the gauge by forcing F-symbols not already fixed to equal 1.
        This method should be used AFTER adding hex and pentagon eqns to ideal_basis
        """
        while len(self.solved) < len(self._poly_ring.gens()):
            #Get a variable that has not been fixed
            #In ascending index order, for consistent results
            for var in self._poly_ring.gens():
                if var not in self.solved:
                    break

            #Fix var = 1, substitute, and solve equations
            self.ideal_basis.add(var-1)
            print("adding equation...", var-1)
            self.ideal_basis = set(Ideal(list(self.ideal_basis)).groebner_basis(algorithm=algorithm))
            self.substitute_degree_one()
            self.update_equations()

    def substitute_degree_one(self, eqns=None):
        if eqns is None:
            eqns = self.ideal_basis

        new_knowns = set()
        useless = set()
        for eq in eqns:
            #Substitute known value from univariate degree 1 polynomial or,
            #Following Bonderson, p. 37, solve linear equation with two terms
            #for one of the variables
            if eq.degree() == 1 and sum(eq.degrees()) <= 2 and eq.lm() not in self.solved:
                self._fvars[self._var_to_sextuple[eq.lm()]] = -sum(c * m for c, m in zip(eq.coefficients()[1:], eq.monomials()[1:])) / eq.lc()
                #Add variable to set of known values and remove this equation
                new_knowns.add(eq.lm())
                useless.add(eq)

            #Solve equation of the form x_i x_j + k == 0 for x_i
            # print("equation: ", eq, "variables ", eq.variables())
            # if eq.degree() == 2 and max(eq.degrees()) == 1 and len(eq.variables()) == 2 and eq.variable(0) not in self.solved:
            #     self._fvars[self._var_to_sextuple[str(eq.variable(0))]] = - eq.constant_coefficient() / eq.variable(1)
            #     print("Subbed {} for {}".format(- eq.constant_coefficient() / eq.variable(1), eq.variable(0)))
            #     #Add variable to set of known values and remove this equation
            #     new_knowns.add(eq.variable(0))
            #     useless.add(eq)

        #Update fvars depending on other variables
        self.solved.update(new_knowns)
        for sextuple, rhs in self._fvars.items():
            d = { var : self._fvars[self._var_to_sextuple[var]] for var in rhs.variables() if var in self.solved }
            if len(d) == 2: print("THREE TERM LINEAR EQUATION ENCOUNTERED!")
            if d:
                self._fvars[sextuple] = rhs.subs(d)

            # if rhs.variables() and rhs.variable() in self.solved:
            #     assert rhs.is_univariate(), "RHS expression is not univariate"
            #     d = { rhs.variable() :  }
            #     # print("Performing substitution of {} with dictionary {}".format(rhs, d))
            #     self._fvars[sextuple] = rhs.subs(d)

        return new_knowns, useless

    def update_equations(self):
        """
        Update ideal_basis equations by plugging in known values
        """
        special_values = { known : self._fvars[self._var_to_sextuple[known]] for known in self.solved }
        self.ideal_basis = set(eq.subs(special_values) for eq in self.ideal_basis)
        self.ideal_basis.discard(0)

    def find_cyclotomic_solution(self, equations=None, factor=True, verbose=True, prune=True, algorithm='', output=False):
        """
        Solve the the hexagon and pentagon relations to evaluate the F-matrix.
        This method (omitting the orthogonality constraints) produces
        output in the cyclotomic field, but it is very limited in the size
        of examples it can handle: for example, `A_2` at level 2 is
        too large for this method. You may use :meth:`find_real_orthogonal_solution`
        to solve much larger examples.

        INPUT:

        - ``equations`` -- (optional) a set of equations to be
          solved. Defaults to the hexagon and pentagon equations.
        - ``factor`` -- (default: ``True``). Set true to use
          the sreduce method to simplify the hexagon and pentagon
          equations before solving them.
        - ``algorithm`` -- (optional). Algorithm to compute Groebner Basis.
        - ``output`` -- (optional, default False). Output a dictionary of
          F-matrix values. This may be useful to see but may be omitted
          since this information will be available afterwards via the
          :meth:`fmatrix` and :meth:`fmat` methods.

        EXAMPLES::

            sage: f = FMatrix(FusionRing("A2",1),fusion_label="a")
            sage: f.get_solution(verbose=False,output=True)
            equations: 8
            equations: 16
            adding equation... fx4 - 1
            {(a2, a2, a2, a0, a1, a1): 1,
            (a2, a2, a1, a2, a1, a0): 1,
            (a2, a1, a2, a2, a0, a0): 1,
            (a2, a1, a1, a1, a0, a2): 1,
            (a1, a2, a2, a2, a0, a1): 1,
            (a1, a2, a1, a1, a0, a0): 1,
            (a1, a1, a2, a1, a2, a0): 1,
            (a1, a1, a1, a0, a2, a2): 1}

        After you successfully run ``get_solution`` you may check
        the correctness of the F-matrix by running :meth:`hexagon`
        and :meth:`pentagon`. These should return empty lists
        of equations. In this example, we turn off the factor
        and prune optimizations to test all instances.
        
        EXAMPLES::

            sage: f.hexagon(factor=False)
            equations: 0
            []
            sage: f.hexagon(factor=False,side="right")
            equations: 0
            []
            sage: f.pentagon(factor=False,prune=False)
            equations: 0
            []

        """
        if equations is None:
            if verbose:
                print("Setting up hexagons and pentagons...")
            equations = self.get_defining_equations("hexagons")+self.get_defining_equations("pentagons")
            #equations = self.hexagon(verbose=False, factor=factor)+self.pentagon(verbose=False, factor=factor, prune=prune)
        if verbose:
            print("Finding a Groebner basis...")
        self.ideal_basis = set(Ideal(equations).groebner_basis(algorithm=algorithm))
        if verbose:
            print("Solving...")
        self.substitute_degree_one()
        if verbose:
            print("Fixing the gauge...")
        self.fix_gauge(algorithm=algorithm)
        if verbose:
            print("Done!")
        if output:
            return self._fvars


    #####################
    ### Verifications ###
    #####################


    def verify_hexagons(self):
        """
        Ensure the hexagon equations are satisfied
        """
        hex = []
        for a,b,c,d,e,g in product(self._FR.basis(),repeat=6):
            lhs = self._field(self._FR.r_matrix(a,c,e))*self.fmat(a,c,b,d,e,g)*self._field(self._FR.r_matrix(b,c,g))
            rhs = sum(self.fmat(c,a,b,d,e,f)*self._field(self._FR.r_matrix(f,c,d))*self.fmat(a,b,c,d,f,g) for f in self._FR.basis())
            hex.append(lhs-rhs)
        if all(h == self._field.zero() for h in hex):
            print("Success!!! Found valid F-symbols for {}".format(self._FR))
        else:
            print("Ooops... something went wrong... These pentagons remain:")
            print(hex)
            return hex

    def verify_pentagons(self,use_mp=True,prune=False):
        """
        Ensure the pentagon equations are satisfied
        """
        print("Testing F-symbols for {}...".format(self._FR))
        fvars_copy = deepcopy(self._fvars)
        self._fvars = { sextuple : float(RDF(rhs)) for sextuple, rhs in self.get_fmats_in_alg_field().items() }
        if use_mp:
            pool = Pool(processes=cpu_count())
        else:
            pool = None
        n_proc = pool._processes if pool is not None else 1
        params = [(child_id,n_proc) for child_id in range(n_proc)]
        pe = self.map_triv_reduce('pent_verify',params,worker_pool=pool,chunksize=1,mp_thresh=0)
        if np.all(np.isclose(np.array(pe),0,atol=1e-7)):
            print("Success!!! Found valid F-symbols for {}".format(self._FR))
            pe = None
        else:
            print("Ooops... something went wrong... These pentagons remain:")
            print(pe)
        self._fvars = fvars_copy
        return pe

    def fmats_are_orthogonal(self):
        """
        Verify that all F-matrices are real and unitary (orthogonal)
        """
        is_orthog = []
        for a,b,c,d in product(self._FR.basis(),repeat=4):
            mat = self.fmatrix(a,b,c,d)
            is_orthog.append(mat.T * mat == matrix.identity(mat.nrows()))
        return all(is_orthog)
