r"""
Littelmann paths
"""
#****************************************************************************
#       Copyright (C) 2012 Mark Shimozono
#                          Anne Schilling
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
#****************************************************************************

from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.element_wrapper import ElementWrapper
from sage.structure.parent import Parent
from sage.categories.crystals import Crystals
from sage.categories.finite_crystals import FiniteCrystals
from sage.combinat.root_system.cartan_type import CartanType
from sage.rings.integer import Integer
from sage.combinat.root_system.root_system import RootSystem
from sage.functions.other import floor
from sage.misc.latex import latex


class CrystalOfLSPaths(UniqueRepresentation, Parent):
    r"""
    Crystal graph of LS paths generated from the straight-line path to a given weight.

    INPUT:

    - ``cartan_type`` -- the Cartan type of a finite or affine root system
    - ``starting_weight`` -- a weight given as a list of coefficients of the fundamental weights

    The crystal class of piecewise linear paths in the weight space,
    generated from a straight-line path from the origin to a given
    element of the weight lattice.

    OUTPUT: - a tuple of weights defining the directions of the piecewise linear segments

    EXAMPLES::

        sage: C = CrystalOfLSPaths(['A',2,1],[-1,0,1]); C
        The crystal of LS paths of type ['A', 2, 1] and weight (-1, 0, 1)
        sage: c = C.module_generators[0]; c
        (-Lambda[0] + Lambda[2],)
        sage: [c.f(i) for i in C.index_set()]
        [None, None, (Lambda[1] - Lambda[2],)]

        sage: R = C.R; R
        Root system of type ['A', 2, 1]
        sage: Lambda = R.weight_space().basis(); Lambda
        Finite family {0: Lambda[0], 1: Lambda[1], 2: Lambda[2]}
        sage: b=C(tuple([-Lambda[0]+Lambda[2]]))
        sage: b==c
        True
        sage: b.f(2)
        (Lambda[1] - Lambda[2],)

    For classical highest weight crystals we can also compare the results with the tableaux implementation::

        sage: C = CrystalOfLSPaths(['A',2],[1,1])
        sage: list(set(C.list()))
        [(-Lambda[1] - Lambda[2],), (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2]), (-Lambda[1] + 2*Lambda[2],),
        (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2]), (Lambda[1] - 2*Lambda[2],), (-2*Lambda[1] + Lambda[2],),
        (2*Lambda[1] - Lambda[2],), (Lambda[1] + Lambda[2],)]
        sage: C.cardinality()
        8
        sage: B = CrystalOfTableaux(['A',2],shape=[2,1])
        sage: B.cardinality()
        8
        sage: B.digraph().is_isomorphic(C.digraph())
        True

    TESTS::

        sage: C = CrystalOfLSPaths(['A',2,1],[-1,0,1])
        sage: TestSuite(C).run(skip=['_test_elements', '_test_elements_eq', '_test_enumerated_set_contains', '_test_some_elements'])
        sage: C = CrystalOfLSPaths(['E',6],[1,0,0,0,0,0])
        sage: TestSuite(C).run()

    REFERENCES::

        .. [L] P. Littelmann, Paths and root operators in representation theory. Ann. of Math. (2) 142 (1995), no. 3, 499-525.
    """

    @staticmethod
    def __classcall__(cls, cartan_type, starting_weight):
        """
        cartan_type and starting_weight are lists, which are mutable. The class
        UniqueRepresentation requires immutable inputs. The following code
        fixes this problem.

        TESTS::

            sage: CrystalOfLSPaths.__classcall__(CrystalOfLSPaths,['A',2,1],[-1,0,1])
            The crystal of LS paths of type ['A', 2, 1] and weight (-1, 0, 1)
        """
        cartan_type = CartanType(cartan_type)
        starting_weight = tuple(starting_weight)
        return super(CrystalOfLSPaths, cls).__classcall__(cls, cartan_type, starting_weight)

    def __init__(self, cartan_type, starting_weight):
        """
        EXAMPLES::

            sage: C = CrystalOfLSPaths(['A',2,1],[-1,0,1]); C
            The crystal of LS paths of type ['A', 2, 1] and weight (-1, 0, 1)
            sage: C.R
            Root system of type ['A', 2, 1]
            sage: C.weight
            -Lambda[0] + Lambda[2]
            sage: C.weight.parent()
            Extended weight space over the Rational Field of the Root system of type ['A', 2, 1]
            sage: C.module_generators
            [(-Lambda[0] + Lambda[2],)]
        """
        self._cartan_type = CartanType(cartan_type)
        self.R = RootSystem(cartan_type)

        self._name = "The crystal of LS paths of type %s and weight %s"%(cartan_type,starting_weight)

        if self._cartan_type.is_affine():
            self.extended = True
            Parent.__init__(self, category = Crystals())
        else:
            self.extended = False
            Parent.__init__(self, category = FiniteCrystals())

        Lambda = self.R.weight_space(extended = self.extended).basis()
        offset = self.R.index_set()[Integer(0)]

        zero_weight = self.R.weight_space(extended = self.extended).zero()
        self.weight = sum([zero_weight]+[starting_weight[j-offset]*Lambda[j] for j in self.R.index_set()])

        if self.weight == zero_weight:
            initial_element = self(tuple([]))
        else:
            initial_element = self(tuple([self.weight]))
        self.module_generators = [initial_element]

    def _repr_(self):
        """
        EXAMPLES::

            sage: CrystalOfLSPaths(['B',3],[1,1,0]) # indirect doctest
            The crystal of LS paths of type ['B', 3] and weight (1, 1, 0)
        """
        return self._name

    class Element(ElementWrapper):
        """
        TESTS::

            sage: C = CrystalOfLSPaths(['E',6],[1,0,0,0,0,0])
            sage: c=C.an_element()
            sage: TestSuite(c).run()
        """

        def endpoint(self):
            r"""
            Computes the endpoint of the path.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b.endpoint()
                Lambda[1] + Lambda[2]
                sage: b.f_string([1,2,2,1])
                (-Lambda[1] - Lambda[2],)
                sage: b.f_string([1,2,2,1]).endpoint()
                -Lambda[1] - Lambda[2]
                sage: b.f_string([1,2])
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                sage: b.f_string([1,2]).endpoint()
                0
                sage: b = C([])
                sage: b.endpoint()
                0
            """
            if len(self.value) > 0:
                return sum(self.value)
            return self.parent().R.weight_space(extended = self.parent().extended).zero()

        def compress(self):
            r"""
            Merges consecutive positively parallel steps present in the path.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: Lambda = C.R.weight_space().fundamental_weights(); Lambda
                Finite family {1: Lambda[1], 2: Lambda[2]}
                sage: c = C(tuple([1/2*Lambda[1]+1/2*Lambda[2], 1/2*Lambda[1]+1/2*Lambda[2]]))
                sage: c.compress()
                (Lambda[1] + Lambda[2],)
            """
            def positively_parallel_weights(v, w):
                """
                Checks whether the vectors ``v`` and ``w`` are positive scalar multiples of each other.
                """
                supp = v.support()
                if len(supp) > 0:
                    i = supp[0]
                    if v[i]*w[i] > 0 and v[i]*w == w[i]*v:
                        return True
                return False

            if len(self.value) == 0:
                return self
            q = []
            curr = self.value[0]
            for i in range(1,len(self.value)):
                if positively_parallel_weights(curr,self.value[i]):
                    curr = curr + self.value[i]
                else:
                    q.append(curr)
                    curr = self.value[i]
            q.append(curr)
            return self.parent()(tuple(q))

        def split_step(self, which_step, r):
            r"""
            Splits indicated step into two parallel steps of relative lengths `r` and `1-r`.

            INPUT:
            - ``which_step`` -- a position in the tuple ``self``
            - ``r`` -- a rational number between 0 and 1

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b.split_step(0,1/3)
                (1/3*Lambda[1] + 1/3*Lambda[2], 2/3*Lambda[1] + 2/3*Lambda[2])
            """
            assert which_step in range(len(self.value))
            v = self.value[which_step]
            return self.parent()(self.value[:which_step]+tuple([r*v,(1-r)*v])+self.value[which_step+1:])

        def reflect_step(self, which_step, i):
            r"""
            Apply the `i`-th simple reflection to the indicated step in ``self``.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b.reflect_step(0,1)
                (-Lambda[1] + 2*Lambda[2],)
                sage: b.reflect_step(0,2)
                (2*Lambda[1] - Lambda[2],)
            """
            assert i in self.index_set()
            assert which_step in range(len(self.value))
            return self.parent()(self.value[:which_step]+tuple([self.value[which_step].simple_reflection(i)])+self.value[which_step+1:])

        def _string_data(self, i):
            r"""
            Computes the `i`-string data of ``self``.

            TESTS::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: b = C.module_generators[0]
                sage: b._string_data(1)
                ()
                sage: b._string_data(2)
                ()
                sage: b.f(1)._string_data(1)
                ((0, -1, -1),)
                sage: b.f(1).f(2)._string_data(2)
                ((0, -1, -1),)
            """
            if len(self.value) == 0:
                return ()
            # get the i-th simple coroot
            alv = self.value[0].parent().alphacheck()[i]
            # Compute the i-heights of the steps of vs
            steps = [v.scalar(alv) for v in self.value]
            # Get the wet step data
            minima_pos = []
            ps = 0
            psmin = 0
            for ix in range(len(steps)):
                ps = ps + steps[ix]
                if ps < psmin:
                    minima_pos.append((ix,ps,steps[ix]))
                    psmin = ps
            return tuple(minima_pos)

        def epsilon(self, i):
            r"""
            Returns the distance to the beginning of the `i`-string.

            This method overrides the generic implementation in the category of crystals
            since this computation is more efficient.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: [c.epsilon(1) for c in C]
                [0, 1, 0, 0, 1, 0, 1, 2]
                sage: [c.epsilon(2) for c in C]
                [0, 0, 1, 2, 1, 1, 0, 0]
            """
            return self.e(i,length_only=True)

        def phi(self, i):
            r"""
            Returns the distance to the end of the `i`-string.

            This method overrides the generic implementation in the category of crystals
            since this computation is more efficient.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: [c.phi(1) for c in C]
                [1, 0, 0, 1, 0, 2, 1, 0]
                sage: [c.phi(2) for c in C]
                [1, 2, 1, 0, 0, 0, 0, 1]
            """
            return self.f(i,length_only=True)

        def e(self, i, power=1, to_string_end=False, length_only=False):
            r"""
            Returns the `i`-th crystal raising operator on ``self``.

            INPUT:

            - ``i`` -- element of the index set of the underlying root system
            - ``power`` -- positive integer; specifies the power of the raising operator
               to be applied (default: 1)
            - ``to_string_end`` -- boolean; if set to True, returns the dominant end of the
              `i`-string of ``self``. (default: False)
            - ``length_only`` -- boolean; if set to True, returns the distance to the dominant
               end of the `i`-string of ``self``.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: c = C[2]; c
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                sage: c.e(1)
                sage: c.e(2)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.e(2,to_string_end=True)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.e(1,to_string_end=True)
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                sage: c.e(1,length_only=True)
                0
            """
            assert i in self.index_set()
            data = self._string_data(i)
            # compute the minimum i-height M on the path
            if len(data) == 0:
                M = 0
            else:
                M = data[-1][1]
            max_raisings = floor(-M)
            if length_only:
                return max_raisings
            # set the power of e_i to apply
            if to_string_end:
                p = max_raisings
            else:
                p = power
            if p > max_raisings:
                return None

            # copy the vector sequence into a working vector sequence ws
            #!!! ws only needs to be the actual vector sequence, not some
            #!!! fancy crystal graph element
            ws = self.parent()(self.value)

            ix = len(data)-1
            while ix >= 0 and data[ix][1] < M + p:
            # get the index of the current step to be processed
                j = data[ix][0]
                # find the i-height where the current step might need to be split
                if ix == 0:
                    prev_ht = M + p
                else:
                    prev_ht = min(data[ix-1][1],M+p)
                # if necessary split the step. Then reflect the wet part.
                if data[ix][1] - data[ix][2] > prev_ht:
                    ws = ws.split_step(j,1-(prev_ht-data[ix][1])/(-data[ix][2]))
                    ws = ws.reflect_step(j+1,i)
                else:
                    ws = ws.reflect_step(j,i)
                ix = ix - 1
            #!!! at this point we should return the fancy crystal graph element
            #!!! corresponding to the humble vector sequence ws
            return self.parent()(ws.compress())

        def dualize(self):
            r"""
            Returns dualized path.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: for c in C:
                ...     print c, c.dualize()
                ...
                (Lambda[1] + Lambda[2],) (-Lambda[1] - Lambda[2],)
                (-Lambda[1] + 2*Lambda[2],) (Lambda[1] - 2*Lambda[2],)
                (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2]) (1/2*Lambda[1] - Lambda[2], -1/2*Lambda[1] + Lambda[2])
                (Lambda[1] - 2*Lambda[2],) (-Lambda[1] + 2*Lambda[2],)
                (-Lambda[1] - Lambda[2],) (Lambda[1] + Lambda[2],)
                (2*Lambda[1] - Lambda[2],) (-2*Lambda[1] + Lambda[2],)
                (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2]) (-Lambda[1] + 1/2*Lambda[2], Lambda[1] - 1/2*Lambda[2])
                (-2*Lambda[1] + Lambda[2],) (2*Lambda[1] - Lambda[2],)
            """
            if len(self.value) == 0:
                return self
            dual_path = [-v for v in self.value]
            dual_path.reverse()
            return self.parent()(tuple(dual_path))

        def f(self, i, power=1, to_string_end=False, length_only=False):
            r"""
            Returns the `i`-th crystal lowering operator on ``self``.

            INPUT:

            - ``i`` -- element of the index set of the underlying root system
            - ``power`` -- positive integer; specifies the power of the lowering operator
               to be applied (default: 1)
            - ``to_string_end`` -- boolean; if set to True, returns the anti-dominant end of the
              `i`-string of ``self``. (default: False)
            - ``length_only`` -- boolean; if set to True, returns the distance to the anti-dominant
               end of the `i`-string of ``self``.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: c = C.module_generators[0]
                sage: c.f(1)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.f(1,power=2)
                sage: c.f(2)
                (2*Lambda[1] - Lambda[2],)
                sage: c.f(2,to_string_end=True)
                (2*Lambda[1] - Lambda[2],)
                sage: c.f(2,length_only=True)
                1

                sage: C = CrystalOfLSPaths(['A',2,1],[-1,-1,2])
                sage: c = C.module_generators[0]
                sage: c.f(2,power=2)
                (Lambda[0] + Lambda[1] - 2*Lambda[2],)
            """
            dual_path = self.dualize()
            dual_path = dual_path.e(i, power, to_string_end, length_only)
            if length_only:
                return dual_path
            if dual_path == None:
                return None
            return dual_path.dualize()

        def s(self, i):
            r"""
            Computes the reflection of ``self`` along the `i`-string.

            This method is more efficient than the generic implementation since it uses
            powers of `e` and `f` in the Littelmann model directly.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: c = C.module_generators[0]
                sage: c.s(1)
                (-Lambda[1] + 2*Lambda[2],)
                sage: c.s(2)
                (2*Lambda[1] - Lambda[2],)

                sage: C = CrystalOfLSPaths(['A',2,1],[-1,0,1])
                sage: c = C.module_generators[0]; c
                (-Lambda[0] + Lambda[2],)
                sage: c.s(2)
                (Lambda[1] - Lambda[2],)
                sage: c.s(1)
                (-Lambda[0] + Lambda[2],)
                sage: c.f(2).s(1)
                (Lambda[0] - Lambda[1],)
            """
            ph = self.phi(i)
            ep = self.epsilon(i)
            diff = ph - ep
            if diff >= 0:
                return self.f(i, power=diff)
            else:
                return self.e(i, power=-diff)

        def _latex_(self):
            r"""
            Latex method for ``self``.

            EXAMPLES::

                sage: C = CrystalOfLSPaths(['A',2],[1,1])
                sage: c = C.module_generators[0]
                sage: c._latex_()
                [\Lambda_{1} + \Lambda_{2}]
            """
            return [latex(p) for p in self.value]
