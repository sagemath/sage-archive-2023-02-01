r"""
Library interface to Kenzo

Kenzo is a set of lisp functions to compute homology and
homotopy groups of topological spaces.

AUTHORS:

- Miguel Marco, Ana Romero (2019-01): Initial version


For this interface, Kenzo is loaded into ECL which is itself loaded
as a C library in Sage. Kenzo objects in this interface are nothing
but wrappers around ECL objects.
"""
# ****************************************************************************
#       Copyright (C) 2019 Miguel Marco <mmarco@unizar.es>
#                      and Ana Romero <ana.romero@unirioja.es>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************
from __future__ import print_function, absolute_import


from sage.structure.sage_object import SageObject
from sage.homology.homology_group import HomologyGroup
from sage.rings.integer_ring import ZZ
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup
from sage.groups.abelian_gps.abelian_group import AbelianGroup
from sage.categories.commutative_additive_groups import CommutativeAdditiveGroups
from sage.groups.additive_abelian.additive_abelian_group import AdditiveAbelianGroup

from sage.matrix.all import matrix
from sage.homology.chain_complex import ChainComplex
from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet

from sage.libs.ecl import EclObject, ecl_eval, EclListIterator

from sage.env import SAGE_LOCAL
from sys import version

# Redirection of ECL and Maxima stdout to /dev/null
# This is also done in the Maxima library, but we
# also do it here for redundancy.
ecl_eval(r"""(defparameter *dev-null* (make-two-way-stream
              (make-concatenated-stream) (make-broadcast-stream)))""")
ecl_eval("(setf original-standard-output *standard-output*)")
ecl_eval("(setf *standard-output* *dev-null*)")

# Loading and initialization of Kenzo
# Note that it will load kenzo.fas file from $SAGE_LOCAL/lib/ecl/
ecl_eval("(require :kenzo)")
ecl_eval("(in-package :cat)")
ecl_eval("(setf *HOMOLOGY-VERBOSE* nil)")


# defining the auxiliary functions as wrappers over the kenzo ones
chcm_mat = EclObject("chcm-mat")
homologie = EclObject("homologie")
sphere = EclObject("sphere")
crts_prdc = EclObject("crts-prdc")
moore = EclObject("moore")
k_z = EclObject("k-z")
k_z2 = EclObject("k-z2")
k_zp = EclObject("k-zp")
echcm = EclObject("echcm")
loop_space = EclObject("loop-space")
tnsr_prdc = EclObject("tnsr-prdc")
typep = EclObject("typep")
classifying_space = EclObject("classifying-space")
suspension = EclObject("suspension")
homotopy_list = EclObject("homotopy-list")
nth = EclObject("nth")
entry = EclObject("entry")
nlig = EclObject("nlig")
ncol = EclObject("ncol")
array_dimensions = EclObject("array-dimensions")
convertmatrice = EclObject("convertmatrice")
make_array_to_lists = EclObject("make-array-to-lists")
make_array_from_lists = EclObject("make-array-from-lists")
chcm_mat2 = EclObject("chcm-mat2")
build_finite_ss2 = EclObject("build-finite-ss2")
gmsm = EclObject("gmsm")
dgop = EclObject("dgop")
dgop_ext_int = EclObject("dgop-ext-int")
dgop_int_ext = EclObject("dgop-int-ext")
basis_aux1 = EclObject("basis_aux1")
orgn_aux1 = EclObject("orgn_aux1")
dffr_aux1 = EclObject("dffr_aux1")
kabstractsimplex_aux1 = EclObject("kabstractsimplex_aux1")
kchaincomplex_aux1 = EclObject("kchaincomplex_aux1")
sfinitesimplicialset_aux1 = EclObject("sfinitesimplicialset_aux1")
spectral_sequence_group = EclObject("spectral-sequence-group")
spectral_sequence_differential_matrix = EclObject("spectral-sequence-differential-matrix")
eilenberg_moore_spectral_sequence = EclObject("eilenberg-moore-spectral-sequence")


def Sphere(n):
    r"""
    Return the ``n`` dimensional sphere as a Kenzo simplicial set.

    INPUT:

    - ``n`` -- the dimension of the sphere

    OUTPUT:

    - A :class:`KenzoSimplicialSet`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import Sphere # optional - kenzo
        sage: s2 = Sphere(2)                           # optional - kenzo
        sage: s2                                       # optional - kenzo
        [K1 Simplicial-Set]
        sage: [s2.homology(i) for i in range(8)]       # optional - kenzo
        [Z, 0, Z, 0, 0, 0, 0, 0]
    """
    kenzosphere = sphere(n)
    return KenzoSimplicialSet(kenzosphere)


def MooreSpace(m, n):
    r"""
    Return the Moore space ``M(m, n)`` as a Kenzo simplicial set.

    The Moore space ``M(m, n)`` is the space whose n'th homology group
    is isomorphic to the cyclic group of order ``m``, and the rest of the
    homology groups are trivial.

    INPUT:

    - ``m`` - A positive integer. The order of the nontrivial homology group.

    - ``n`` - The dimension in which the homology is not trivial


    OUTPUT:

    - A KenzoSimplicialSet

    EXAMPLES::

        sage: from sage.interfaces.kenzo import MooreSpace   # optional - kenzo
        sage: m24 = MooreSpace(2,4)                          # optional - kenzo
        sage: m24                                            # optional - kenzo
        [K10 Simplicial-Set]
        sage: [m24.homology(i) for i in range(8)]            # optional - kenzo
        [Z, 0, 0, 0, C2, 0, 0, 0]
    """
    kenzomoore = moore(m, n)
    return KenzoSimplicialSet(kenzomoore)


def EilenbergMacLaneSpace(G, n):
    r"""
    Return the Eilenberg-MacLane space ``K(G, n)`` as a Kenzo simplicial group.

    The Eilenberg-MacLane space ``K(G, n)`` is the space whose has n'th homotopy
    group isomorphic to ``G``, and the rest of the homotopy groups are trivial.

    INPUT:

    - ``G`` -- group. Currently only ``ZZ`` and the additive group of two
      elements are supported.

    - ``n`` -- the dimension in which the homotopy is not trivial


    OUTPUT:

    - A :class:`KenzoSimplicialGroup`

    EXAMPLES::

        sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace    # optional - kenzo
        sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
        sage: [e3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, Z, 0, C2, 0, C3]
        sage: f3 = EilenbergMacLaneSpace(AdditiveAbelianGroup([2]), 3)   # optional - kenzo
        sage: [f3.homology(i) for i in range(8)]                         # optional - kenzo
        [Z, 0, 0, C2, 0, C2, C2, C2]
    """
    if G == ZZ:
        kenzospace = k_z(n)
        return KenzoSimplicialGroup(kenzospace)
    elif G == AdditiveAbelianGroup([2]):
        kenzospace = k_z2(n)
        return KenzoSimplicialGroup(kenzospace)
    elif G in CommutativeAdditiveGroups() and G.is_cyclic():
        kenzospace = k_zp(G.cardinality(), n)
        return KenzoSimplicialGroup(kenzospace)
    else:
        raise NotImplementedError("Eilenberg-MacLane spaces are only supported over ZZ and ZZ_n")


class KenzoObject(SageObject):
    r"""
    Wrapper to Kenzo objects

    INPUT:

    - ``kenzo_object`` -- a wrapper around a Kenzo object
      (which is an ecl object).
    """

    def __init__(self, kenzo_object):
        r"""
        Construct the chain complex.

        TESTS::

            sage: from sage.interfaces.kenzo import KenzoObject  # optional -kenzo
            sage: from sage.interfaces.kenzo import sphere       # optional -kenzo
            sage: ks = sphere(2)                                 # optional -kenzo
            sage: ks                                             # optional -kenzo
            <ECL: [K1 Simplicial-Set]>
            sage: s2 = KenzoObject(ks)                           # optional -kenzo
            sage: s2                                             # optional -kenzo
            [K1 Simplicial-Set]
            sage: TestSuite(s2).run(skip='_test_pickling')       # optional -kenzo

        """
        self._kenzo = kenzo_object

    def _repr_(self):
        r"""
        Represent the object.

        It just uses the ecl representation, removing the
        ecl decoration.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace    # optional - kenzo
            sage: m2 = MooreSpace(2,4)                            # optional - kenzo
            sage: m2._repr_()                                     # optional - kenzo
            '[K10 Simplicial-Set]'
        """
        kenzo_string = repr(self._kenzo)
        return kenzo_string[6:-1]

class KenzoSpectralSequence(KenzoObject):
    r"""
    Wrapper around Kenzo spectral sequences
    """

    def group(self, p, i, j):
        r"""
        Return the ``i,j``'th group of the ``p`` page.

        INPUT:

        - ``p`` -- the page to take the group from.

        - ``i`` -- the column where the group is taken from.

        - ``j`` -- the row where the group is taken from.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional -- kenzo
            sage: S2 = Sphere(2)                           # optional -- kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional -- kenzo
            sage: EMS.group(0, -1, 2)                      # optional -- kenzo
            Additive abelian group isomorphic to Z
            sage: EMS.group(0, -1, 3)                      # optional -- kenzo
            Trivial group
        """
        invs = spectral_sequence_group(self._kenzo, p, i ,j).python()
        if not invs:
            invs = []
        return AdditiveAbelianGroup(invs)

    def matrix(self, p, i, j):
        r"""
        Return the matrix that determines the differential from the
        ``i,j``'th group of the ``p``'th page.

        INPUT:

        - ``p`` -- the page.

        - ``i`` -- the column of the differential domain.

        - ``j`` -- the row of the differential domain.
        
        EXAMPLES::
        
            sage: from sage.interfaces.kenzo import Sphere   # optional -- kenzo
            sage: S3 = Sphere(3)                             # optional -- kenzo
            sage: L = S3.loop_space()                        # optional -- kenzo
            sage: EMS = L.em_spectral_sequence()             # optional -- kenzo
            sage: EMS.table(1, -5, -2, 5, 8)                 # optional -- kenzo
              0   Z   Z + Z + Z   Z + Z + Z
              0   0   0           0
              0   0   Z           Z + Z
              0   0   0           0
            sage: EMS.matrix(1, -2 ,8)                       # optional -- kenzo
            [ 3  3  0]
            [-2  0  2]
            [ 0 -3 -3]            
        """
        klist = spectral_sequence_differential_matrix(self._kenzo, p, i, j)
        plist = klist.python()
        if plist is None or plist==[None]:
            i = len(self.group(p, i, j).invariants())
            j = len(self.group(p, i-p, j+p-1).invariants())
            return matrix(i,j)
        return matrix(plist)
    
    def differential(self, p, i, j):
        r"""
        Return the ``(p, i, j)`` differential morphism of the spectral sequence.

        INPUT:

        - ``p`` -- the page.

        - ``i`` -- the column of the differential domain.

        - ``j`` -- the row of the differential domain.
        
        EXAMPLES::
        
            sage: from sage.interfaces.kenzo import Sphere   # optional -- kenzo
            sage: S3 = Sphere(3)                             # optional -- kenzo
            sage: L = S3.loop_space()                        # optional -- kenzo
            sage: EMS = L.em_spectral_sequence()             # optional -- kenzo
            sage: EMS.table(1,-5,-2,5,8)                     # optional -- kenzo
              0   Z   Z + Z + Z   Z + Z + Z
              0   0   0           0
              0   0   Z           Z + Z
              0   0   0           0
            sage: EMS.matrix(1, -3, 8)                       # optional -- kenzo
            [ 2]
            [-2]
            [ 2]
            sage: EMS.differential(1, -3, 8)                 # optional -- kenzo
            Morphism from module over Integer Ring with invariants (0, 0, 0) to module with invariants (0,) that sends the generators to [(2), (-2), (2)]
        """
        domain = self.group(p, i, j)
        codomain = self.group(p, i-p, j+p-1)
        M = self.matrix(p, i, j)
        images = [codomain(r) for r in M.rows()]
        return domain.hom(images, codomain=codomain)

    def table(self, p, i1, i2, j1, j2):
        r"""
        Return a table printing the groups in the ``p`` page.

        INPUT:

        - ``p`` -- the page to print.

        -- ``i1`` -- the first column to print.

        -- ``i2`` -- the last column to print.

        -- ``j1`` -- the first row to print.

        -- ``j2`` -- the last row to print.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional -- kenzo
            sage: S2 = Sphere(2)                           # optional -- kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional -- kenzo
            sage: EMS.table(0, -2, 2, -2, 2)               # optional -- kenzo
              0   Z   0   0   0
              0   0   0   0   0
              0   0   Z   0   0
              0   0   0   0   0
              0   0   0   0   0
        """
        from sage.misc.table import table
        groups = []
        for j in range(j2-j1+1):
            row = []
            for i in range(i1, i2+1):
                group = self.group(p,i,j2-j)
                if group.invariants():
                    row.append(group.short_name())
                else:
                    row.append('0')
            groups.append(row)
        return table(groups)


class KenzoChainComplex(KenzoObject):
    r"""
    Wrapper to Kenzo chain complexes.
    """
    def homology(self, n):
        r"""
        Return the ``n``'th homology group of the kenzo chain complex

        INPUT:

        - ``n`` -- the dimension in which compute the homology

        OUTPUT:

        - An homology group.

        EXAMPLES::

                sage: from sage.interfaces.kenzo import Sphere   # optional - kenzo
                sage: s2 = Sphere(2)                             # optional - kenzo
                sage: s2                                         # optional - kenzo
                [K1 Simplicial-Set]
                sage: s2.homology(2)                             # optional - kenzo
                Z
        """
        echcm1 = echcm(self._kenzo)
        m1 = chcm_mat(echcm1, n)
        m2 = chcm_mat(echcm1, n+1)
        homology = homologie(m1, m2)
        lhomomology = [i for i in EclListIterator(homology)]
        res = []
        for component in lhomomology:
            pair = [i for i in EclListIterator(component)]
            res.append(pair[0].python())
        return HomologyGroup(len(res), ZZ, res)

    def tensor_product(self, other):
        r"""
        Return the tensor product of ``self`` and ``other``.

        INPUT:

        - ``other`` --  The Kenzo object with which to compute the tensor product

        OUTPUT:

        - A :class:`KenzoChainComplex`

        sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
        sage: s2 = Sphere(2)                              # optional - kenzo
        sage: s3 = Sphere(3)                              # optional - kenzo
        sage: p = s2.tensor_product(s3)                   # optional - kenzo
        sage: type(p)                                     # optional - kenzo
        <class 'sage.interfaces.kenzo.KenzoChainComplex'>
        sage: [p.homology(i) for i in range(8)]           # optional - kenzo
        [Z, 0, Z, Z, 0, Z, 0, 0]
        """
        return KenzoChainComplex(tnsr_prdc(self._kenzo, other._kenzo))

    def basis(self, dim):
        r"""
        Return the list of generators of the Kenzo chain complex ``self`` in dimension ``dim``.

        INPUT:

        - ``dim``- An integer number

        OUTPUT:

        - A list of the form ['G"dim"G0', 'G"dim"G1', 'G"dim"G2', ...].

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                      # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)         # optional - kenzo
            sage: kenzo_chcm = KChainComplex(sage_chcm)                                # optional - kenzo
            sage: kenzo_chcm                                                           # optional - kenzo
            [K589 Chain-Complex]
            sage: for i in range(6):                                                   # optional - kenzo
            ....:     print("Basis in dimension %i: %s" % (i, kenzo_chcm.basis(i)))    # optional - kenzo
            Basis in dimension 0: ['G0G0', 'G0G1', 'G0G2']
            Basis in dimension 1: ['G1G0', 'G1G1']
            Basis in dimension 2: None
            Basis in dimension 3: ['G3G0', 'G3G1']
            Basis in dimension 4: ['G4G0', 'G4G1']
            Basis in dimension 5: ['G5G0', 'G5G1', 'G5G2']

        """
        return basis_aux1(self._kenzo, dim).python()

    def differential(self, dim, comb):
        r"""
        Return the differential of a combination.

        INPUT:

        - ``dim``- An integer number

        - ``comb``- A list representing a formal sum of generators in the module
        of dimension ``dim``. For example, to represent G7G12 + 3*G7G0 - 5*G7G3
        we use the list [3, 'G7G0', -5, 'G7G3', 1, 'G7G12']. Note that the
        generators must be in ascending order respect to the number after the
        second G in their representation; the parameter
        ``comb`` = [1, 'G7G12', 3, 'G7G0', -5, 'G7G3'] will produce an error in
        Kenzo.

        OUTPUT:

        - A Kenzo combination representing the differential of the formal
        combination represented by ``comb`` in the chain complex ``self`` in
        dimension ``dim``.

        EXAMPLES::

            sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
            sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
            sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
            sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
            sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)
            sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
            sage: kenzo_chcm                                                      # optional - kenzo
            [K... Chain-Complex]
            sage: kenzo_chcm.basis(4)                                             # optional - kenzo
            ['G4G0', 'G4G1']
            sage: kenzo_chcm.differential(4, [1, 'G4G0'])                                 # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 3}
            <1 * G3G0>
            <3 * G3G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
            sage: kenzo_chcm.basis(5)                                             # optional - kenzo
            ['G5G0', 'G5G1', 'G5G2']
            sage: kenzo_chcm.differential(5, [1, 'G5G0', 2, 'G5G2'])                      # optional - kenzo
            <BLANKLINE>
            ----------------------------------------------------------------------{CMBN 4}
            <6 * G4G0>
            <-3 * G4G1>
            ------------------------------------------------------------------------------
            <BLANKLINE>
        """
        cmbn_list = pairing(comb)
        return KenzoObject(dffr_aux1(self._kenzo, dim, cmbn_list))

    def orgn(self):
        return str(orgn_aux1(self._kenzo))


class KenzoSimplicialSet(KenzoChainComplex):
    r"""
    Wrapper to Kenzo simplicial sets.
    """

    def loop_space(self, n=1):
        r"""
        Return the ``n`` th iterated loop space.

        INPUT:

        - ``n`` -- (default: 1) the number of times to iterate the loop space
          construction

        OUTPUT:

        - A :class:`KenzoSimplicialGroup`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: l2 = s2.loop_space()                        # optional - kenzo
            sage: type(l2)                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: l2 = s2.loop_space()                        # optional - kenzo
            sage: [l2.homology(i) for i in range(8)]          # optional - kenzo
            [Z, Z, Z, Z, Z, Z, Z, Z]
        """
        return KenzoSimplicialGroup(loop_space(self._kenzo, n))

    def cartesian_product(self, other):
        r"""
        Return the cartesian product of ``self`` and ``other``.

        INPUT:

        - ``other`` -- the Kenzo simplicial set with which the product is made

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere    # optional - kenzo
            sage: s2 = Sphere(2)                              # optional - kenzo
            sage: s3 = Sphere(3)                              # optional - kenzo
            sage: p = s2.cartesian_product(s3)                # optional - kenzo
            sage: type(p)                                     # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [p.homology(i) for i in range(6)]           # optional - kenzo
            [Z, 0, Z, Z, 0, Z]
        """
        prod_kenzo = crts_prdc(self._kenzo, other._kenzo)
        return KenzoSimplicialSet(prod_kenzo)

    def suspension(self):
        r"""
        Return the suspension of the simplicial set.

        OUTPUT:

        - A :class:`KenzoSimplicialSet`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import EilenbergMacLaneSpace    # optional - kenzo
            sage: e3 = EilenbergMacLaneSpace(ZZ, 3)                          # optional - kenzo
            sage: s = e3.suspension()                                        # optional - kenzo
            sage: type(s)                                                    # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialSet'>
            sage: [s.homology(i) for i in range(6)]                          # optional - kenzo
            [Z, 0, 0, 0, Z, 0]
        """
        return KenzoSimplicialSet(suspension(self._kenzo))

    def homotopy_group(self, n):
        """
        Return the n'th homotopy group of ``self``
        INPUT:
        - ``n`` - the dimension of the homotopy group to be computed
        EXAMPLES::
            sage: from sage.interfaces.kenzo import Sphere      # optional - kenzo
            sage: s2 = Sphere(2)                                # optional - kenzo
            sage: p = s2.cartesian_product(s2)                  # optional - kenzo
            sage: p.homotopy_group(3)                           # optional - kenzo
            Multiplicative Abelian group isomorphic to Z x Z
        """
        if n not in ZZ or n < 2:
            raise ValueError("homotopy groups can only be computed for dimensions greater than 1")
        lgens = homotopy_list(self._kenzo, n).python()
        if lgens is not None:
            trgens = [0 if i == 1 else i for i in sorted(lgens)]
            return AbelianGroup(trgens)
        else:
            return AbelianGroup([])
        
    def em_spectral_sequence(self):
        r"""
        Return the Eilenberg-Moore spectral sequence of self
        
        EXAMPLES::

            sage: from sage.interfaces.kenzo import Sphere # optional -- kenzo
            sage: S2 = Sphere(2)                           # optional -- kenzo
            sage: EMS = S2.em_spectral_sequence()          # optional -- kenzo
            sage: EMS.table(0, -2, 2, -2, 2)               # optional -- kenzo
              0   Z   0   0   0
              0   0   0   0   0
              0   0   Z   0   0
              0   0   0   0   0
              0   0   0   0   0
        """
        if self.homology(1).invariants():
            raise ValueError("Eilenberg-Moore spectral sequence implemented only for 1-reduced simplicial sets")
        return KenzoSpectralSequence(eilenberg_moore_spectral_sequence(self._kenzo))


class KenzoSimplicialGroup(KenzoSimplicialSet):
    r"""
    Wrapper around Kenzo simplicial groups.
    """

    def classifying_space(self):
        r"""
        Return the classifying space.

        OUTPUT:

        - A :class:`KenzoSimplicialGroup`

        EXAMPLES::

            sage: from sage.interfaces.kenzo import MooreSpace      # optional - kenzo
            sage: m2 = MooreSpace(2,4)                              # optional - kenzo
            sage: l2 = m2.loop_space()                              # optional - kenzo
            sage: c = l2.classifying_space()                        # optional - kenzo
            sage: type(c)                                           # optional - kenzo
            <class 'sage.interfaces.kenzo.KenzoSimplicialGroup'>
            sage: [c.homology(i) for i in range(8)]                 # optional - kenzo
            [Z, 0, 0, 0, C2, 0, 0, 0]
        """
        return KenzoSimplicialGroup(classifying_space(self._kenzo))


def k2s_matrix(kmatrix):
    r"""
    Convert an array of ECL to a matrix of Sage.
    """
    dimensions = array_dimensions(kmatrix).python()
    kmatrix_list = make_array_to_lists(kmatrix).python()
    return matrix(dimensions[0], dimensions[1], kmatrix_list)


def s2k_matrix(smatrix):
    r"""
    Convert a matrix of Sage to an array of ECL.
    """
    initcontents = []
    dimensions = smatrix.dimensions()
    for i in smatrix.rows():
        initcontents.append(i.list())
    return make_array_from_lists(dimensions[0], dimensions[1], initcontents)


def s2k_dictmat(sdictmat):
    r"""
    Convert a dictionary in Sage, whose values are matrices, to an assoc list
    in ECL.
    """
    rslt = EclObject([])
    for k in sdictmat.keys():
        rslt = EclObject(k).cons(s2k_matrix(sdictmat[k])).cons(rslt)
    return rslt


def pairing(slist):
    r"""
    Convert a list of Sage (which has an even length) to an assoc list in ECL.
    """
    rslt = EclObject([])
    for k in range(len(slist) - 1, 0, -2):
        rslt = EclObject(slist[k - 1]).cons(EclObject(slist[k])).cons(rslt)
    return rslt


def KChainComplex(chain_complex):
    r"""
    Construct a KenzoChainComplex from a ChainComplex of degree = -1 in
    Sage.

    INPUT:

    - ``chain_complex`` - A ChainComplex of degree = -1

    OUTPUT:

    - A KenzoChainComplex

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplex                 # optional - kenzo
        sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
        sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
        sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
        sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)    # optional - kenzo
        sage: kenzo_chcm = KChainComplex(sage_chcm)                           # optional - kenzo
        sage: kenzo_chcm                                                      # optional - kenzo
        [K... Chain-Complex]
        sage: kenzo_chcm.homology(5)                                          # optional - kenzo
        Z x Z
    """
    d = chain_complex.differential()
    chcm = s2k_dictmat(d)
    str_orgn = str(d)[1:-1].replace(":", " ").replace(" ", ".").replace("\n", "").replace(",", "")
    return KenzoChainComplex(kchaincomplex_aux1(chcm, str_orgn))


def SChainComplex(kchaincomplex, start=0, end=15):
    r"""
    Convert the KenzoChainComplex ``kchcm`` (between dimensions ``start`` and
    ``end``) to a ChainComplex.

    INPUT:

    - ``kchaincomplex``- A KenzoChainComplex

    - ``start``- An integer number (optional, default 0)

    - ``end``- An integer number greater than or equal to ``start`` (optional, default 15)

    OUTPUT:

    - A ChainComplex

    EXAMPLES::

        sage: from sage.interfaces.kenzo import KChainComplex, SChainComplex    # optional - kenzo
        sage: m1 = matrix(ZZ, 3, 2, [-1, 1, 3, -4, 5, 6])
        sage: m4 = matrix(ZZ, 2, 2, [1, 2, 3, 6])
        sage: m5 = matrix(ZZ, 2, 3, [2, 2, 2, -1, -1, -1])
        sage: sage_chcm = ChainComplex({1: m1, 4: m4, 5: m5}, degree = -1)     # optional - kenzo
        sage: SChainComplex(KChainComplex(sage_chcm)) == sage_chcm             # optional - kenzo
        True
    """
    matrices = {}
    for i in range(start, end):
        dffr_i = chcm_mat2(kchaincomplex._kenzo, i)
        if ((nlig(dffr_i).python() != 0) and (ncol(dffr_i).python() != 0)):
            matrices[i] = k2s_matrix(convertmatrice(dffr_i))
    return ChainComplex(matrices, degree=-1)


def SAbstractSimplex(simplex, dim):
    r"""
    Convert an abstract simplex of Kenzo to an AbstractSimplex.
    INPUT:

    - ``simplex``- An abstract simplex of Kenzo.

    - ``dim``- The dimension of ``simplex``.

    OUTPUT:

    - An AbstractSimplex.

    EXAMPLES::

        sage: from sage.libs.ecl import EclObject, ecl_eval
        sage: from sage.interfaces.kenzo import KenzoObject, SAbstractSimplex    # optional - kenzo
        sage: KAbSm = KenzoObject(ecl_eval("(ABSM 15 'K)"))                      # optional - kenzo
        sage: SAbSm1 = SAbstractSimplex(KAbSm, 2)                                # optional - kenzo
        sage: SAbSm2 = SAbstractSimplex(KAbSm, 7)                                # optional - kenzo
        sage: SAbSm1.degeneracies()                                              # optional - kenzo
        [3, 2, 1, 0]
        sage: SAbSm1.dimension()                                                 # optional - kenzo
        6
        sage: SAbSm2.dimension()                                                 # optional - kenzo
        11
    """
    degeneracies = dgop_int_ext(dgop(simplex._kenzo)).python()
    if degeneracies is None:
        degeneracies = []
    else:
        degeneracies = tuple(degeneracies)
    name = gmsm(simplex._kenzo).python()
    return AbstractSimplex(dim, degeneracies, name=name)


def KAbstractSimplex(simplex):
    r"""
    Convert an AbstractSimplex in Sage to an abstract simplex of Kenzo.
    INPUT:

    - ``simplex``- An AbstractSimplex.

    OUTPUT:

    - An abstract simplex of Kenzo.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex
        sage: from sage.interfaces.kenzo import KAbstractSimplex, SAbstractSimplex    # optional - kenzo
        sage: SAbSm = AbstractSimplex(1, (2,0,3,2,1), name = 'SAbSm')                 # optional - kenzo
        sage: KAbSm = KAbstractSimplex(SAbSm)                                         # optional - kenzo
        sage: SAbSm2 = SAbstractSimplex(KAbSm, 1)                                     # optional - kenzo
        sage: SAbSm.degeneracies() == SAbSm2.degeneracies()                           # optional - kenzo
        True
        sage: SAbSm.dimension() == SAbSm2.dimension()                                 # optional - kenzo
        True
    """
    return KenzoObject(kabstractsimplex_aux1(simplex.degeneracies(), 's'+str(hash(simplex))))


def KFiniteSimplicialSet(sset):
    r"""
    Convert a finite SimplicialSet in Sage to a finite simplicial set of Kenzo.

    INPUT:

    - ``sset``- A finite SimplicialSet.

    OUTPUT:

    - A finite simplicial set of Kenzo.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import AbstractSimplex, SimplicialSet
        sage: from sage.interfaces.kenzo import KFiniteSimplicialSet    # optional - kenzo
        sage: s0 = AbstractSimplex(0, name='s0')
        sage: s1 = AbstractSimplex(0, name='s1')
        sage: s2 = AbstractSimplex(0, name='s2')
        sage: s01 = AbstractSimplex(1, name='s01')
        sage: s02 = AbstractSimplex(1, name='s02')
        sage: s12 = AbstractSimplex(1, name='s12')
        sage: s012 = AbstractSimplex(2, name='s012')
        sage: Triangle = SimplicialSet({s01: (s1, s0), s02: (s2, s0), s12: (s2, s1)}, base_point = s0)
        sage: KTriangle = KFiniteSimplicialSet(Triangle)                # optional - kenzo
        sage: KTriangle.homology(1)                                     # optional - kenzo
        Z
        sage: KTriangle.basis(1)                                        # optional - kenzo
        ['CELL_1_0', 'CELL_1_1', 'CELL_1_2']
        sage: S1 = simplicial_sets.Sphere(1)
        sage: S3 = simplicial_sets.Sphere(3)
        sage: KS1vS3 = KFiniteSimplicialSet(S1.wedge(S3))               # optional - kenzo
        sage: KS1vS3.homology(3)                                        # optional - kenzo
        Z
    """
    from sage.homology.simplicial_set_constructions import ProductOfSimplicialSets
    if isinstance(sset, ProductOfSimplicialSets):
        f0 = KFiniteSimplicialSet(sset.factor(0))
        for f1 in sset.factors()[1:]:
            f0 = f0.cartesian_product(KFiniteSimplicialSet(f1))
        return f0
    else:
        allcells = sset.cells()
        namecells = {c:'cell_{}_{}'.format(d, allcells[d].index(c)) for d in allcells for c in allcells[d]}
        dim = sset.dimension()
        list_rslt = [namecells[i] for i in sset.n_cells(0)]
        if (dim > 0):
            for k in range(1, dim + 1):
                k_cells = sset.n_cells(k)
                if k_cells:
                    list_rslt.append(k)
                    for x in k_cells:
                        list_rslt.append(namecells[x])
                        auxiliar_list = []
                        for z in sset.faces(x):
                            degen_z = z.degeneracies()
                            name = namecells[z.nondegenerate()]
                            degen_z.append(name)
                            auxiliar_list.append(degen_z)
                        list_rslt.append(auxiliar_list)
        return KenzoSimplicialSet(build_finite_ss2(list_rslt))


def SFiniteSimplicialSet(ksimpset, limit):
    r"""Convert the ``limit``-skeleton of a finite simplicial set in Kenzo to a
    finite SimplicialSet in Sage.

    INPUT:

    - ``ksimpset``- A finite simplicial set in Kenzo.

    - ``limit``- A natural number.

    OUTPUT:

    - A finite SimplicialSet.

    EXAMPLES::

        sage: from sage.homology.simplicial_set import SimplicialSet
        sage: from sage.interfaces.kenzo import AbstractSimplex, KFiniteSimplicialSet, SFiniteSimplicialSet, Sphere    # optional - kenzo
        sage: s0 = AbstractSimplex(0, name='s0')                    # optional - kenzo
        sage: s1 = AbstractSimplex(0, name='s1')                    # optional - kenzo
        sage: s2 = AbstractSimplex(0, name='s2')                    # optional - kenzo
        sage: s01 = AbstractSimplex(1, name='s01')                  # optional - kenzo
        sage: s02 = AbstractSimplex(1, name='s02')                  # optional - kenzo
        sage: s12 = AbstractSimplex(1, name='s12')                  # optional - kenzo
        sage: s012 = AbstractSimplex(2, name='s012')                # optional - kenzo
        sage: Triangle = SimplicialSet({s01: (s1, s0), s02: (s2, s0), s12: (s2, s1)}, base_point = s0) # optional - kenzo
        sage: KTriangle = KFiniteSimplicialSet(Triangle)            # optional - kenzo
        sage: STriangle = SFiniteSimplicialSet(KTriangle, 1)        # optional - kenzo
        sage: STriangle.homology()                                  # optional - kenzo
        {0: 0, 1: Z}
        sage: S1 = simplicial_sets.Sphere(1)                        # optional - kenzo
        sage: S3 = simplicial_sets.Sphere(3)                        # optional - kenzo
        sage: KS1vS3 = KFiniteSimplicialSet(S1.wedge(S3))           # optional - kenzo
        sage: SS1vS3 = SFiniteSimplicialSet(KS1vS3, 3)              # optional - kenzo
        sage: SS1vS3.homology()                                     # optional - kenzo
        {0: 0, 1: Z, 2: 0, 3: Z}
    """

    list_orgn = orgn_aux1(ksimpset._kenzo).python()
    if nth(0, list_orgn).python()[0] == 'CRTS-PRDC':
        return SFiniteSimplicialSet(KenzoSimplicialSet(nth(1, list_orgn)),
                                    limit).cartesian_product(SFiniteSimplicialSet(KenzoSimplicialSet(nth(2, list_orgn)), limit))
    rslt = {}
    simplices = []
    faces = []
    bases = []
    names = []
    for k in range(limit + 1):
        basis_k = basis_aux1(ksimpset._kenzo, k)
        names_k = ksimpset.basis(k)
        lbasis_k = [AbstractSimplex(k, name=i) for i in EclListIterator(basis_k)]
        bases.append(lbasis_k)
        names.append(names_k)
    all_simplices = sfinitesimplicialset_aux1(ksimpset._kenzo, limit)
    lall_simplices = [i for i in EclListIterator(all_simplices)]
    dim = 1
    for Kdim in lall_simplices:
        for simp in Kdim:
            index1 = names[dim].index(str(simp.car()))
            lKdim_cdr = []
            for i in EclListIterator(simp.cdr()):
                degenop = dgop_int_ext(dgop(i)).python()
                if degenop is None:
                    degenop = []
                index2 = names[dim - len(degenop) - 1].index(str(gmsm(i)))
                lKdim_cdr.append(bases[dim - len(degenop) - 1][index2].apply_degeneracies(*degenop))
            simplices.append(bases[dim][index1])
            faces.append(tuple(lKdim_cdr))
        dim += 1
    for i in range(len(simplices)):
        rslt[simplices[i]] = faces[i]
    return SimplicialSet(rslt)
