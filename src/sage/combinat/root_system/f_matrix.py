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

from sage.misc.misc import inject_variable
from sage.matrix.constructor import matrix
from sage.rings.polynomial.all import PolynomialRing
from sage.rings.ideal import Ideal
from sage.combinat.root_system.fusion_ring import FusionRing
import sage.graphs
from sage.graphs.generators.basic import EmptyGraph
from itertools import product
from sage.misc.misc import inject_variable

#Import pickle for checkpointing and loading certain variables
try:
    import cPickle as pickle
except:
    import pickle

from sage.rings.polynomial.polydict import ETuple
# from sage.rings import AlgebraicField as QQbar

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

    The F-matrix is only determined up to a gauge. This
    is a family of embeddings `C\to A\otimes B` for
    simple objects `A,B,C` such that `\text{Hom}(C,A\otimes B)`
    is nonzero. Changing the gauge changes the F-matrix though
    not in a very essential way. By varying the gauge it is
    possible to make the F-matrices unitary, or it is possible
    to make them cyclotomic. We choose the latter.

    Due to the large number of equations we may fail to find a
    Groebner basis if there are too many variables.

    EXAMPLES::

        sage: I=FusionRing("E8",2,conjugate=True)
        sage: I.fusion_labels(["i0","p","s"],inject_variables=True)
        sage: f = FMatrix(I,inject_variables=True); f
        creating variables fx1..fx14
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
    the R-matrix.

    EXAMPLES::

        sage: f.pentagon()[1:3]
        equations: 41
        [-fx0*fx1 + fx1, -fx1*fx2^2 + fx1]
        sage: f.hexagon()[1:3]
        equations: 14
        [fx1*fx5 + fx2, fx2 + 1]

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
        self.FR = fusion_ring
        if self.FR._fusion_labels is None:
            self.FR.fusion_labels(fusion_label, inject_variables=True)
            #Set up F-symbols entry by entry
        n_vars = self.findcases()
        self._poly_ring = PolynomialRing(self.FR.field(),n_vars,var_prefix)
        if inject_variables:
            print ("creating variables %s%s..%s%s"%(var_prefix,1,var_prefix,n_vars))
            for i in range(self._poly_ring.ngens()):
                inject_variable("%s%s"%(var_prefix,i),self._poly_ring.gens()[i])
        self._var_to_sextuple, self._fvars = self.findcases(output=True)
        self._singles = self.singletons()

        #Initialize list of defining equations
        self.ideal_basis = list()

        #Initialize empty set of solved F-symbols
        self.solved = set()

        #New attributes of the FMatrix class
        self._var_to_idx = { var : idx for idx, var in enumerate(self._poly_ring.gens()) }
        self._known_vals = dict()
        self.field = self.FR.field()
        self.temp_eqns = []

    #######################
    ### Class utilities ###
    #######################

    def __repr__(self):
        """
        EXAMPLES::

            sage: FMatrix(FusionRing("B2",1))
            F-Matrix factory for The Fusion Ring of Type B2 and level 1 with Integer Ring coefficients
        """
        return "F-Matrix factory for %s"%self.FR

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
        if self.FR.Nk_ij(a,b,x) == 0 or self.FR.Nk_ij(x,c,d) == 0 or self.FR.Nk_ij(b,c,y) == 0 or self.FR.Nk_ij(a,y,d) == 0:
            return 0

        #Some known zero F-symbols
        if a == self.FR.one():
            if x == b and y == d:
                return 1
            else:
                return 0
        if b == self.FR.one():
            if x == a and y == c:
                return 1
            else:
                return 0
        if c == self.FR.one():
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
        for (a,b,c,d) in list(product(self.FR.basis(), repeat=4)):
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
        for (a, b, c, d) in list(product(self.FR.basis(), repeat=4)):
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

        return [x for x in self.FR.basis() if self.FR.Nk_ij(a,b,x) != 0 and self.FR.Nk_ij(x,c,d) != 0]

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

        return [y for y in self.FR.basis() if self.FR.Nk_ij(b,c,y) != 0 and self.FR.Nk_ij(a,y,d) != 0]

    ####################
    ### Data getters ###
    ####################

    def get_fmats_in_alg_field(self):
        return { sextuple : self._qqbar_embedding(fvar) for sextuple, fvar in self._fvars.items() }

    #Return radical expression for fmats for easy visualization
    def get_radical_expression(self):
        return { sextuple : val.radical_expression() for sextuple, val in get_fmats_in_alg_field().items() }

    #Construct a dictionary of idx, known_val pairs for equation substitution
    def get_known_vals(self):
        return { var_idx : self._fvars[self._var_to_sextuple[self._poly_ring.gen(var_idx)]] for var_idx in self.solved }

    #Construct an ETuple indicating positions of known nonzero variables.
    #MUST be called after fmats._ks = get_known_sq()
    def get_known_nonz(self):
        nonz = { self._var_to_idx[var] : 100 for var in self._singles }
        for idx in self._ks:
            nonz[idx] = 100
        return ETuple(nonz, fmats._poly_ring.ngens())

    ##############################
    ### Variables partitioning ###
    ##############################

    #Get the size of the largest F-matrix F^{abc}_d
    def largest_fmat_size(self):
        return max(self.fmatrix(*tup).nrows() for tup in product(self.FR.basis(),repeat=4))

    #Partition the F-symbols according to the size of the F-matrix F^{abc}_d
    #they belong to
    def get_fmats_by_size(self,n):
        fvars_copy = deepcopy(self._fvars)
        self.clear_vars()
        var_set = set()
        for quadruple in product(self.FR.basis(),repeat=4):
            F = self.fmatrix(*quadruple)
            #Discard trivial 1x1 F-matrix, if applicable
            if F.nrows() == n and F.coefficients() != [1]:
                var_set.update(self._var_to_idx[fx] for fx in F.coefficients())
        self._fvars = fvars_copy
        return var_set

    ############################
    ### Checkpoint utilities ###
    ############################

    #Auto-generate filename string
    def get_fr_str(self):
        return self.FR.cartan_type()[0] + str(self.FR.cartan_type()[1]) + str(self.FR.fusion_level())

    def save_fvars(self,filename):
        with open(filename, 'wb') as f:
            pickle.dump(self._fvars, f)
            # pickle.dump([fmats._fvars, fmats.solved], f)

    #If provided, optional param save_dir should have a trailing forward slash
    def load_fvars(self,save_dir=""):
        with open(save_dir + "saved_fvars_" + self.get_fr_str() + ".pickle", 'rb') as f:
            self._fvars = pickle.load(f)
            # fmats._fvars, fmats.solved = pickle.load(f)

    ########################
    ### Equations set up ###
    ########################

    def poly_to_tup(self,poly):
        return tuple(poly.dict().items())

    def tup_to_poly(self,eq_tup,parent=None):
        if parent is None:
            parent = self._poly_ring
        return parent({ exp : coeff for exp, coeff in eq_tup })

    def get_orthogonality_constraints(self,output=True):
        eqns = list()
        for tup in product(self.FR.basis(), repeat=4):
            mat = self.fmatrix(*tup)
            eqns.extend((mat.T * mat - matrix.identity(mat.nrows())).coefficients())
        if output:
            return eqns
        self.ideal_basis.extend([self.poly_to_tup(eq) for eq in eqns])

    def feq(self, a, b, c, d, e, f, g, k, l, prune=False):
        """
        Return True if the Pentagon axiom ([Bond2007]_ (2.77)) is satisfied.
        """
        lhs = self.fmat(f,c,d,e,g,l)*self.fmat(a,b,l,e,f,k)
        rhs = sum(self.fmat(a,b,c,g,f,h)*self.fmat(a,h,d,e,g,k)*self.fmat(b,c,d,k,h,l) for h in self.FR.basis())
        if lhs != 0 or not prune: # it is believed that if lhs=0, the equation carries no new information
            return lhs - rhs
        else:
            return 0

    def req(self, a, b, c, d, e, g, side="left"):
        """
        Return A hexagon equation (Bond[2007]_ (2.78) or (2.79)).

        INPUT:

        - ``a,b,c,d,e,f`` -- basis elements of the FusionRing
        - ``side`` -- (default left) which hexagon axiom to use

        """
        if side == "left":
            lhs = self.FR.r_matrix(a,c,e)*self.fmat(a,c,b,d,e,g)*self.FR.r_matrix(b,c,g)
            rhs = sum(self.fmat(c,a,b,d,e,f)*self.FR.r_matrix(f,c,d)*self.fmat(a,b,c,d,f,g) for f in self.FR.basis())
        elif side == "right":
            # r(a,b,x) is a root of unity, so its inverse is its complex conjugate
            lhs = self.FR.r_matrix(c,a,e).conjugate()*self.fmat(a,c,b,d,e,g)*self.FR.r_matrix(c,b,g).conjugate()
            rhs = sum(self.fmat(c,a,b,d,e,f)*self.FR.r_matrix(c,f,d).conjugate()*self.fmat(a,b,c,d,f,g) for f in self.FR.basis())
        return lhs-rhs

    def pentagon(self, verbose=False, output=True, factor=False, prune=False):
        """
        Return generators of the ideal of Pentagon equations.

        INPUT:

        - ``verbose`` -- (optional) set True for verbose. Default False
        - ``output`` -- (optional) set True to output a set of equations. Default True
        - ``factor`` -- (optional) set True for sreduce simplified equations.

        In contrast with the hexagon equations, where setting ``factor`` True
        is a big improvement, for the pentagon equations this option produces
        little or no simplification. So the default is False.

        EXAMPLES::

            sage: p = FMatrix(FusionRing("A2",1),fusion_label="c")
            sage: p.pentagon()[-3:]
            equations: 16
            [-fx5*fx6 + fx1, -fx4*fx6*fx7 + fx2, -fx5*fx7^2 + fx3*fx6]
            sage: p.pentagon(factor=True)[-3:]
            equations: 16
            [fx5*fx6 - fx1, fx4*fx6*fx7 - fx2, fx5*fx7^2 - fx3*fx6]

        """
        if output:
            ret = []
        for (a,b,c,d,e,f,g,k,l) in list(product(self.FR.basis(), repeat=9)):
            pd = self.feq(a,b,c,d,e,f,g,k,l,prune=prune)
            if pd != 0:
                if factor:
                    pd = self.sreduce(pd)
                if output:
                    ret.append(pd)
                if verbose:
                    print ("%s,%s,%s,%s,%s,%s,%s,%s,%s : %s"%(a,b,c,d,e,f,g,k,l,pd))
        print ("equations: %s"%len(ret))
        if output:
            return ret

    #####################
    ### Graph methods ###
    #####################

    def equations_graph(self,eqns=None):
        """
        Return a graph whose vertices are variables
        and whose edges correspond to equations
        relating two variables.

        This graph may contain isolated vertices...

        INPUT:

        - ``equations`` -- a list of equations

        EXAMPLES::

            sage: f = FMatrix(FusionRing("A1",3))
            sage: G = f.equation_graph(f.hexagon(factor=True))
            equations: 71
            sage: G.connected_components_number()
            14
            sage: G1=G.connected_components_subgraphs()[0]
            sage: G1.size()
            60
            sage: G1.is_regular()
            True

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

    #######################
    ### Solution method ###
    #######################

    def get_solution(self, equations=None, factor=True, verbose=True, prune=True, algorithm='', output=False):
        """
        Solve the the hexagon and pentagon relations to evaluate the F-matrix.

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
            equations = self.hexagon(verbose=False, factor=factor)+self.pentagon(verbose=False, factor=factor, prune=prune)
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

    #Ensure the hexagon equations are satisfied
    def verify_hexagons(self):
        hex = []
        for a,b,c,d,e,g in product(self.FR.basis(),repeat=6):
            lhs = self.field(self.FR.r_matrix(a,c,e))*self.fmat(a,c,b,d,e,g)*self.field(self.FR.r_matrix(b,c,g))
            rhs = sum(self.fmat(c,a,b,d,e,f)*self.field(self.FR.r_matrix(f,c,d))*self.fmat(a,b,c,d,f,g) for f in self.FR.basis())
            hex.append(lhs-rhs)
        return list(set(hex))

    #Verify that all F-matrices are real and unitary (orthogonal)
    def fmats_are_orthogonal(self):
        is_orthog = []
        for a,b,c,d in product(self.FR.basis(),repeat=4):
            mat = self.fmatrix(a,b,c,d)
            is_orthog.append(mat.T * mat == matrix.identity(mat.nrows()))
        return all(is_orthog)
