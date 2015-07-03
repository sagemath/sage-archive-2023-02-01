r"""
Bijection between rigged configurations for `B(\infty)` and marginally large tableaux

AUTHORS:

- Travis Scrimshaw (2015-07-01): Initial version

REFERENCES:

.. [RC_MLT] Ben Salisbury and Travis Scrimshaw. *Connecting marginally
   large tableaux and rigged configurations via crystals*.
   Pre-print. :arxiv:`1505.07040`.
"""

#*****************************************************************************
#       Copyright (C) 2015 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations

from sage.combinat.rigged_configurations.bij_type_B import (KRTToRCBijectionTypeB,
                                                            RCToKRTBijectionTypeB)
from sage.combinat.rigged_configurations.bij_type_D import (KRTToRCBijectionTypeD,
                                                            RCToKRTBijectionTypeD)
from sage.combinat.rigged_configurations.bij_type_A import (KRTToRCBijectionTypeA,
                                                            RCToKRTBijectionTypeA)
from sage.combinat.rigged_configurations.bij_type_A2_odd import (KRTToRCBijectionTypeA2Odd,
                                                                 RCToKRTBijectionTypeA2Odd)
from sage.combinat.rigged_configurations.tensor_product_kr_tableaux import TensorProductOfKirillovReshetikhinTableaux
from sage.combinat.crystals.letters import CrystalOfLetters
from sage.combinat.root_system.cartan_type import CartanType
from sage.categories.morphism import Morphism
from sage.categories.homset import Hom
from sage.misc.flatten import flatten

def affine_type_from_classical(ct):
    """
    Return the affine type corresponding to the classical type ``ct``
    used in the bijection.

    EXAMPLES::

        sage: from sage.combinat.rigged_configurations.bij_infinity import affine_type_from_classical
        sage: affine_type_from_classical(CartanType(['A', 6]))
        ['A', 6, 1]
        sage: affine_type_from_classical(CartanType(['B', 6]))
        ['B', 6, 1]
        sage: affine_type_from_classical(CartanType(['C', 6]))
        ['B', 6, 1]^*
        sage: affine_type_from_classical(CartanType(['D', 6]))
        ['D', 6, 1]
    """
    if ct.type() == 'A':
        return CartanType(['A', ct.rank(), 1])
    if ct.type() == 'B':
        return CartanType(['B', ct.rank(), 1])
    if ct.type() == 'C':
        return CartanType(['A', 2*ct.rank()-1, 2])
    if ct.type() == 'D':
        return CartanType(['D', ct.rank(), 1])
    raise NotImplementedError("bijection not implemented for type {}".format(ct))

class FromTableauIsomorphism(Morphism):
    r"""
    Crystal isomorphism of `B(\infty)` in the tableau model to the
    rigged configuration model.
    """
    def _repr_type(self):
        """
        Return the type of morphism of ``self``.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: phi = RC.coerce_map_from(T)
            sage: phi._repr_type()
            'Crystal Isomorphism'
        """
        return "Crystal Isomorphism"

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: phi = RC.coerce_map_from(T)
            sage: ~phi
            Crystal Isomorphism morphism:
              From: The infinity crystal of rigged configurations of type ['A', 3]
              To:   The infinity crystal of tableaux of type ['A', 3]
        """
        return FromRCIsomorphism(Hom(self.codomain(), self.domain()))

    def _call_(self, x):
        r"""
        Return the image of ``x`` in the rigged configuration model
        of `B(\infty)`.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: phi = RC.coerce_map_from(T)
            sage: x = T.an_element().f_string([2,2,1,1,3,2,1,2,1,3])
            sage: y = phi(x); ascii_art(y)
            -4[ ][ ][ ][ ]-2  -3[ ][ ][ ]-1  -1[ ][ ]-1
                              -2[ ]-1
            sage: (~phi)(y) == x
            True
        """
        conj = x.to_tableau().conjugate()
        ct = self.domain().cartan_type()
        act = affine_type_from_classical(ct)
        TP = TensorProductOfKirillovReshetikhinTableaux(act, [[r,1] for r in conj.shape()])
        elt = TP(pathlist=[reversed(row) for row in conj])
        if ct.type() == 'A':
            bij = KRTToRCBijectionTypeA(elt)
        elif ct.type() == 'B':
            bij = MLTToRCBijectionTypeB(elt)
        elif ct.type() == 'C':
            bij = KRTToRCBijectionTypeA2Odd(elt)
        elif ct.type() == 'D':
            bij = MLTToRCBijectionTypeD(elt)
        else:
            raise NotImplementedError("bijection of type {} not yet implemented".format(ct))
        return self.codomain()(bij.run())

class FromRCIsomorphism(Morphism):
    r"""
    Crystal isomorphism of `B(\infty)` in the rigged configuration model
    to the tableau model.
    """
    def _repr_type(self):
        """
        Return the type of morphism of ``self``.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: phi = T.coerce_map_from(RC)
            sage: phi._repr_type()
            'Crystal Isomorphism'
        """
        return "Crystal Isomorphism"

    def __invert__(self):
        """
        Return the inverse of ``self``.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: phi = T.coerce_map_from(RC)
            sage: ~phi
            Crystal Isomorphism morphism:
              From: The infinity crystal of tableaux of type ['A', 3]
              To:   The infinity crystal of rigged configurations of type ['A', 3]
        """
        return FromTableauIsomorphism(Hom(self.codomain(), self.domain()))

    def _call_(self, x):
        r"""
        Return the image of ``x`` in the tableau model of `B(\infty)`.

        EXAMPLES::

            sage: T = crystals.infinity.Tableaux(['A',3])
            sage: RC = crystals.infinity.RiggedConfigurations(['A',3])
            sage: phi = T.coerce_map_from(RC)
            sage: x = RC.an_element().f_string([2,2,1,1,3,2,1,2,1,3])
            sage: y = phi(x); y.pp()
              1  1  1  1  1  2  2  3  4
              2  2  3  4
              3
            sage: (~phi)(y) == x
            True
        """
        lam = [sum(nu)+1 for nu in x]
        ct = self.domain().cartan_type()
        I = ct.index_set()
        if ct.type() == 'B':
            lam[-1] *= 2
        elif ct.type() == 'D':
            lam[-1] = max(lam[-2], lam[-1])
            lam[-2] = lam[-1]
        l = sum([ [[r,1]]*lam[i] for i,r in enumerate(I) ], [])

        RC = RiggedConfigurations(affine_type_from_classical(ct), reversed(l))
        elt = RC(x)
        if ct.type() == 'A':
            bij = RCToKRTBijectionTypeA(elt)
        elif ct.type() == 'B':
            bij = RCToMLTBijectionTypeB(elt)
        elif ct.type() == 'C':
            bij = RCToKRTBijectionTypeA2Odd(elt)
        elif ct.type() == 'D':
            bij = RCToMLTBijectionTypeD(elt)
        else:
            raise NotImplementedError("bijection of type {} not yet implemented".format(ct))
        y = bij.run()

        # Now make the result marginally large
        y = [list(c) for c in y]
        cur = []
        L = CrystalOfLetters(ct)
        for i in I:
            cur.insert(0, L(i))
            c = y.count(cur)
            while c > 1:
                y.remove(cur)
                c -= 1
        return self.codomain()(*flatten(y))

class MLTToRCBijectionTypeB(KRTToRCBijectionTypeB):
    def run(self):
        """
        Run the bijection from a marginally large tableaux to a rigged
        configuration.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['B',4])
            sage: T = crystals.infinity.Tableaux(['B',4])
            sage: Psi = T.crystal_morphism({T.module_generators[0]: RC.module_generators[0]})
            sage: TS = T.subcrystal(max_depth=4)
            sage: all(Psi(b) == RC(b) for b in TS) # long time # indirect doctest
            True
        """
        for cur_crystal in reversed(self.tp_krt):
            # Iterate through the columns
            for col_number, cur_column in enumerate(reversed(cur_crystal.to_array(False))):
                self.cur_path.insert(0, []) # Prepend an empty list

                self.cur_dims.insert(0, [0, 1])

                for letter in reversed(cur_column):
                    self.cur_dims[0][0] += 1

                    val = letter.value # Convert from a CrystalOfLetter to an Integer

                    #print("====================")
                    #print(repr(self.cur_path))
                    #print("--------------------")
                    #print(repr(self.ret_rig_con))
                    #print("--------------------\n")

                    # Build the next state
                    self.cur_path[0].insert(0, [letter]) # Prepend the value
                    if self.cur_dims[0][0] == self.n:
                        # Spinor case, we go from \Lambda_{n-1} -> 2\Lambda_n
                        self.cur_dims.insert(1, [self.n,1])
                        self.cur_path.insert(1, self.cur_path[0])

                    self.next_state(val)

        self.ret_rig_con.set_immutable() # Return it to immutable
        return self.ret_rig_con

class RCToMLTBijectionTypeB(RCToKRTBijectionTypeB):
    def run(self):
        """
        Run the bijection from rigged configurations to a large tableau.

        EXAMPLES::

            sage: RC = crystals.infinity.RiggedConfigurations(['B',4])
            sage: T = crystals.infinity.Tableaux(['B',4])
            sage: Psi = RC.crystal_morphism({RC.module_generators[0]: T.module_generators[0]})
            sage: RCS = RC.subcrystal(max_depth=4)
            sage: all(Psi(nu) == T(nu) for nu in RCS) # long time # indirect doctest
            True
        """
        letters = CrystalOfLetters(self.rigged_con.parent()._cartan_type.classical())

        # This is technically bad, but because the first thing we do is append
        #   an empty list to ret_crystal_path, we correct this. We do it this
        #   way so that we do not have to remove an empty list after the
        #   bijection has been performed.
        ret_crystal_path = []

        while self.cur_dims:
            dim = self.cur_dims[0]
            ret_crystal_path.append([])

            # Assumption: all factors are single columns
            if dim[0] == self.n:
                # Spinor case, since we've done 2\Lambda_n -> \Lambda_{n-1}
                self.cur_dims.pop(1)

            while dim[0] > 0:
                #print("====================")
                #print(repr(self.rigged_con.parent()(*self.cur_partitions, use_vacancy_numbers=True)))
                #print("--------------------")
                #print(ret_crystal_path)
                #print("--------------------\n")

                dim[0] -= 1 # This takes care of the indexing
                b = self.next_state(dim[0])

                # Make sure we have a crystal letter
                ret_crystal_path[-1].append(letters(b)) # Append the rank

            self.cur_dims.pop(0) # Pop off the leading column

        return ret_crystal_path

class MLTToRCBijectionTypeD(KRTToRCBijectionTypeD):
    pass # TODO - placeholder

class RCToMLTBijectionTypeD(RCToKRTBijectionTypeD):
    pass # TODO - placeholder

