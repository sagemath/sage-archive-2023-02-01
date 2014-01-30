r"""
PQ-Trees

This module implements PQ-Trees and methods to help recognise Interval
Graphs. It is used by :meth:`is_interval
<sage.graphs.generic_graph.GenericGraph.is_interval>`.

Author:

- Nathann Cohen

"""

################################################################################
#      Copyright (C) 2012 Nathann Cohen <nathann.cohen@gail.com>               #
#                                                                              #
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL) #
#                         http://www.gnu.org/licenses/                         #
################################################################################

# Constants, to make the code more readable

FULL           = 2
PARTIAL        = 1
EMPTY          = 0
ALIGNED        = True
UNALIGNED      = False

##########################################################################
# Some Lambda Functions                                                  #
#                                                                        #
# As the elements of a PQ-Tree can be either P-Trees, Q-Trees, or the    #
# sets themselves (the leaves), the following lambda function are        #
# meant to be applied both on PQ-Trees and Sets, and mimic for the       #
# latter the behaviour we expect from the corresponding methods          #
# defined in class PQ                                                    #
##########################################################################

set_contiguous = lambda tree, x : (
    tree.set_contiguous(x) if isinstance(tree, PQ) else
    ((FULL, ALIGNED) if x in tree
     else (EMPTY, ALIGNED)))

new_P = lambda liste : P(liste) if len(liste) > 1 else liste[0]
new_Q = lambda liste : Q(liste) if len(liste) > 1 else liste[0]

flatten = lambda x : x.flatten() if isinstance(x, PQ) else x

impossible_msg = "Impossible"

def reorder_sets(sets):
    r"""
    Reorders a collection of sets such that each element appears on an
    interval.

    Given a collection of sets `C = S_1,...,S_k` on a ground set `X`,
    this function attempts to reorder them in such a way that `\forall
    x \in X` and `i<j` with `x\in S_i, S_j`, then `x\in S_l` for every
    `i<l<j` if it exists.

    INPUT:

    - ``sets`` - a list of instances of ``list, Set`` or ``set``

    ALGORITHM:

    PQ-Trees

    EXAMPLE:

    There is only one way (up to reversal) to represent contiguously
    the sequence ofsets `\{i-1, i, i+1\}`::

        sage: from sage.graphs.pq_trees import reorder_sets
        sage: seq = [Set([i-1,i,i+1]) for i in range(1,15)]

    We apply a random permutation::

        sage: p = Permutations(len(seq)).random_element()
        sage: seq = [ seq[p(i+1)-1] for i in range(len(seq)) ]
        sage: ordered = reorder_sets(seq)
        sage: if not 0 in ordered[0]:
        ...      ordered = ordered.reverse()
        sage: print ordered
        [{0, 1, 2}, {1, 2, 3}, {2, 3, 4}, {3, 4, 5}, {4, 5, 6}, {5, 6, 7}, {8, 6, 7}, {8, 9, 7}, {8, 9, 10}, {9, 10, 11}, {10, 11, 12}, {11, 12, 13}, {12, 13, 14}, {13, 14, 15}]
    """

    if len(sets) == 1:
        return sets

    s = set([])

    for ss in sets:
        for i in ss:
            s.add(i)

    tree = P(sets)


    for i in s:
        tree.set_contiguous(i)
        tree = flatten(tree)

    return tree.ordering()

class PQ:
    r"""
    This class implements the PQ-Tree, used for the recognition of
    Interval Graphs, or equivalently for matrices having the so-caled
    "consecutive ones property".

    Briefly, we are given a collection `C=S_1, ..., S_n` of sets on a
    common ground set `X`, and we would like to reorder the elements
    of `C` in such a way that for every element `x\in X` such that
    `x\in S_i` and `x\in S_j`, `i<j`, we have `x\in S_l` for all
    `i<l<j`. This property could also be rephrased as : the sets
    containing `x` are an interval.

    To achieve it, we will actually compute ALL the orderings
    satisfying such constraints using the structure of PQ-Tree, by
    adding the constraints one at a time.

        * At first, there is no constraint : all the permutations are
          allowed. We will then build a tree composed of one node
          linked to all the sets in our collection (his children). As
          we want to remember that all the permutations of his
          children are allowed, we will label it with "P", making it a
          P-Tree.

        * We are now picking an element `x \in X`, and we want to
          ensure that all the elements `C_x` containing it are
          contiguous. We can remove them from their tree `T_1`, create
          a second tree `T_2` whose only children are the `C_x`, and
          attach this `T_2` to `T_1`. We also make this new tree a
          `P-Tree`, as all the elements of `C_x` can be permuted as
          long as they stay close to each other. Obviously, the whole
          tree `T_2` can be prmuter with the other children of `T_1`
          in any way -- it does not impair the fact that the sequence
          of the children will ensure the sets containing `x` are
          contiguous.

        * We would like to repeat the same procedure for `x' \in X`,
          but we are now encountering a problem : there may be sets
          containing both `x'` and `x`, along with others containing
          only `x` or only `x'`. We can permute the sets containing
          only `x'` together, or the sets containing both `x` and `x'`
          together, but we may NOT permute all the sets containing
          `x'` together as this may break the relationship between the
          sets containing `x`. We need `Q`-Trees. A `Q`-Tree is a tree
          whose children are ordered, even though their order could be
          reversed (if the children of a `Q`-Tree are `c_1c_2
          ... c_k`, we can change it to `c_k ... c_2c_1`). We can now
          express all the orderings satisfying our two constraints the
          following way :

            * We create a tree `T_1` gathering all the elements not
              containing `x` nor `x'`, and make it a `P-Tree`

            * We create 3 `P`-Trees `T_{x, x'}, T_x, T_{x'}`, which
              respectively have for children the elements or our
              collection containing

                * both `x` and `x'`
                * only `x`
                * only `x'`

            * To ensure our constraints on both elements, we create a
              `Q`-tree `T_2` whose children are in order `T_x, T_{x,
              x'}, T_{x'}`

            * We now make this `Q`-Tree `T_2` a children of the
              `P`-Tree `T_1`

    Using these two types of tree, and exploring the different cases
    of intersection, it is possible to represent all the possible
    permutations of our sets satisfying or constraints, or to prove
    that no such ordering exists. This is the whole purpose of this
    class and this algorithm, and is explained with more details in many
    places, for example in the following document from Hajiaghayi [Haj]_.

    REFERENCES:

    .. [Haj] M. Hajiaghayi
      http://www-math.mit.edu/~hajiagha/pp11.ps

    AUTHOR : Nathann Cohen
    """

    def __init__(self, seq):
        r"""
        Construction of a PQ-Tree

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
        """
        from sage.sets.set import Set

        self._children = []
        for e in seq:
            if isinstance(e, list):
                e = Set(e)

            if not e in self._children:
                self._children.append(e)

    def reverse(self):
        r"""
        Recursively reverses ``self`` and its children

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: p.ordering()
            [{1, 2}, {2, 3}, {2, 4}, {8, 2}, {9, 2}]
            sage: p.reverse()
            sage: p.ordering()
            [{9, 2}, {8, 2}, {2, 4}, {2, 3}, {1, 2}]
        """
        for i in self._children:
            if isinstance(i, PQ):
                i.reverse()

        self._children.reverse()

    def __contains__(self, v):
        r"""
        Tests whether there exists an element of ``self`` containing
        an element ``v``

        INPUT:

        - ``v`` -- an element of the ground set

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: 5 in p
            False
            sage: 9 in p
            True

        """
        for i in self:
            if v in i:
                return True
        False

    def split(self, v):
        r"""
        Returns the subsequences of children containing and not
        containing ``v``

        INPUT:

        - ``v`` -- an element of the ground set

        OUTPUT:

        Two lists, the first containing the children of ``self``
        containing ``v``, and the other containing the other children.

        .. NOTE::

           This command is meant to be used on a partial tree, once it
           has be "set continuous" on an element ``v`` and aligned it
           to the right. Hence, none of the list should be empty (an
           exception is raised if that happens, as it would reveal a
           bug in the algorithm) and the sum ``contains +
           does_not_contain`` should be equal to the sequence of
           children of ``self``.

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: p.reverse()
            sage: contains, does_not_contain = p.split(1)
            sage: contains
            [{1, 2}]
            sage: does_not_contain
            [('P', [{9, 2}, {8, 2}, {2, 4}]), {2, 3}]
            sage: does_not_contain + contains == p._children
            True

        """
        contains = []
        does_not_contain = []

        for i in self:
            if v in i:
                contains.append(i)
            else:
                does_not_contain.append(i)

        if not contains or not does_not_contain:
            raise ValueError("None of the sets should be empty !")

        return contains, does_not_contain

    def __iter__(self):
        r"""
        Iterates over the children of ``self``.

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: for i in p:
            ...      print i
            {1, 2}
            {2, 3}
            ('P', [{2, 4}, {8, 2}, {9, 2}])
        """

        for i in self._children:
            yield i

    def cardinality(self):
        r"""
        Returns the number of children of ``self``

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: p.cardinality()
            3
        """
        return len(self._children)

    def ordering(self):
        r"""
        Returns the current ordering given by listing the leaves from
        left to right.

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: p.ordering()
            [{1, 2}, {2, 3}, {2, 4}, {8, 2}, {9, 2}]
        """
        value = []
        for i in self:
            if isinstance(i, PQ):
                value.extend(i.ordering())
            else:
                value.append(i)

        return value

    def __repr__(self):
        r"""
        Succintly represents ``self``.

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([[1,2], [2,3], P([[2,4], [2,8], [2,9]])])
            sage: print p
            ('Q', [{1, 2}, {2, 3}, ('P', [{2, 4}, {8, 2}, {9, 2}])])
        """
        return str((("P" if self.is_P() else "Q"),self._children))

    def simplify(self, v, left = False, right = False):
        r"""
        Returns a simplified copy of self according to the element ``v``

        If ``self`` is a partial P-tree for ``v``, we would like to
        restrict the permutations of its children to permutations
        keeping the children containing ``v`` contiguous. This
        function also "locks" all the elements not containing ``v``
        inside a `P`-tree, which is useful when one want to keep the
        elements containing ``v`` on one side (which is the case when
        this method is called).

        INPUT:

        - ``left, right`` (booleans) -- whether ``v`` is aligned to the
          right or to the left

        - ``v```-- an element of the ground set

        OUTPUT:

        If ``self`` is a `Q`-Tree, the sequence of its children is
        returned. If ``self`` is a `P`-tree, 2 `P`-tree are returned,
        namely the two `P`-tree defined above and restricting the
        permutations, in the order implied by ``left, right`` (if
        ``right =True``, the second `P`-tree will be the one gathering
        the elements containing ``v``, if ``left=True``, the
        opposite).

        .. NOTE::

           This method is assumes that ``self`` is partial for ``v``,
           and aligned to the side indicated by ``left, right``.

        EXAMPLES:

        A `P`-Tree ::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = P([[2,4], [1,2], [0,8], [0,5]])
            sage: p.simplify(0, right = True)
            [('P', [{2, 4}, {1, 2}]), ('P', [{0, 8}, {0, 5}])]

        A `Q`-Tree ::

            sage: q = Q([[2,4], [1,2], [0,8], [0,5]])
            sage: q.simplify(0, right = True)
            [{2, 4}, {1, 2}, {0, 8}, {0, 5}]
        """
        if sum([left, right]) !=1:
            raise ValueError("Exactly one of left or right must be specified")

        if self.is_Q():
            return self._children
        else:

            contains, does_not_contain = self.split(v)

            A = new_P(does_not_contain)
            B = new_P(contains)

            if right:
                return [A, B]
            else:
                return [B, A]

    def flatten(self):
        r"""
        Returns a flattened copy of ``self``

        If self has only one child, we may as well consider its
        child's children, as ``self`` encodes no information. This
        method recursively "flattens" trees having only on PQ-tree
        child, and returns it.

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = Q([P([[2,4], [2,8], [2,9]])])
            sage: p.flatten()
            ('P', [{2, 4}, {8, 2}, {9, 2}])
        """
        if self.cardinality() == 1:
            return flatten(self._children[0])
        else:
            self._children = [flatten(x) for x in self._children]
            return self


    def is_P(self):
        r"""
        Tests whether ``self`` is a `P`-Tree

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: P([[0,1],[2,3]]).is_P()
            True
            sage: Q([[0,1],[2,3]]).is_P()
            False
        """
        return isinstance(self,P)

    def is_Q(self):
        r"""
        Tests whether ``self`` is a `Q`-Tree

        EXAMPLE::

            sage: from sage.graphs.pq_trees import P, Q
            sage: Q([[0,1],[2,3]]).is_Q()
            True
            sage: P([[0,1],[2,3]]).is_Q()
            False
        """
        return isinstance(self,Q)

class P(PQ):
    r"""
    A P-Tree is a PQ-Tree whose children are
    not ordered (they can be permuted in any way)
    """
    def set_contiguous(self, v):
        r"""
        Updates ``self`` so that its sets containing ``v`` are
        contiguous for any admissible permutation of its subtrees.

        This function also ensures, whenever possible,
        that all the sets containing ``v`` are located on an interval
        on the right side of the ordering.

        INPUT:

        - ``v`` -- an element of the ground set

        OUTPUT:

        According to the cases :

            * ``(EMPTY, ALIGNED)`` if no set of the tree contains
              an occurrence of ``v``

            * ``(FULL, ALIGNED)`` if all the sets of the tree contain
              ``v``

            * ``(PARTIAL, ALIGNED)`` if some (but not all) of the sets
              contain ``v``, all of which are aligned
              to the right of the ordering at the end when the function ends

            * ``(PARTIAL, UNALIGNED)`` if some (but not all) of the
              sets contain ``v``, though it is impossible to align them
              all to the right

        In any case, the sets containing ``v`` are contiguous when this
        function ends. If there is no possibility of doing so, the function
        raises a ``ValueError`` exception.

        EXAMPLE:

        Ensuring the sets containing ``0`` are continuous::

            sage: from sage.graphs.pq_trees import P, Q
            sage: p = P([[0,3], [1,2], [2,3], [2,4], [4,0],[2,8], [2,9]])
            sage: p.set_contiguous(0)
            (1, True)
            sage: print p
            ('P', [{1, 2}, {2, 3}, {2, 4}, {8, 2}, {9, 2}, ('P', [{0, 3}, {0, 4}])])

        Impossible situation::

            sage: p = P([[0,1], [1,2], [2,3], [3,0]])
            sage: p.set_contiguous(0)
            (1, True)
            sage: p.set_contiguous(1)
            (1, True)
            sage: p.set_contiguous(2)
            (1, True)
            sage: p.set_contiguous(3)
            Traceback (most recent call last):
            ...
            ValueError: Impossible

        """

        ###############################################################
        # Defining Variables :                                        #
        #                                                             #
        # Collecting the information of which children are FULL of v, #
        # which ones are EMPTY, PARTIAL_ALIGNED and PARTIAL_UNALIGNED #
        #                                                             #
        # Defining variables for their cardinals, just to make the    #
        # code slightly more readable :-)                             #
        ###############################################################

        seq = [set_contiguous(x, v) for x in self]
        self.flatten()
        seq = [set_contiguous(x, v) for x in self]

        f_seq = dict(zip(self, seq))

        set_FULL                = []
        set_EMPTY               = []
        set_PARTIAL_ALIGNED     = []
        set_PARTIAL_UNALIGNED   = []

        sorting = {
            (FULL, ALIGNED)     : set_FULL,
            (EMPTY, ALIGNED)    : set_EMPTY,
            (PARTIAL, ALIGNED)  : set_PARTIAL_ALIGNED,
            (PARTIAL, UNALIGNED) : set_PARTIAL_UNALIGNED
            }

        for i in self:
            sorting[f_seq[i]].append(i)

        n_FULL                  = len(set_FULL)
        n_EMPTY                 = len(set_EMPTY)
        n_PARTIAL_ALIGNED       = len(set_PARTIAL_ALIGNED)
        n_PARTIAL_UNALIGNED     = len(set_PARTIAL_UNALIGNED)

        counts = dict(map(lambda (x,y) : (x,len(y)), sorting.iteritems()))

        # Excludes the situation where there is no solution.
        # read next comment for more explanations

        if (n_PARTIAL_ALIGNED + n_PARTIAL_UNALIGNED > 2 or
            (n_PARTIAL_UNALIGNED >= 1 and n_EMPTY != self.cardinality() -1)):

            raise ValueError(impossible_msg)

        # From now on, there are at most two pq-trees which are partially filled
        # If there is one which is not aligned to the right, all the others are empty

        #########################################################
        # 1/2                                                   #
        #                                                       #
        # Several easy cases where we can decide without paying #
        # attention                                             #
        #########################################################

        # All the children are FULL
        elif n_FULL == self.cardinality():
            return FULL, True

        # All the children are empty
        elif n_EMPTY == self.cardinality():
            return EMPTY, True

        # There is a PARTIAL UNALIGNED element (and all the others are
        # empty as we checked before

        elif n_PARTIAL_UNALIGNED == 1:
            return (PARTIAL, UNALIGNED)

        # If there is just one partial element and all the others are
        # empty, we just reorder the set to put it at the right end

        elif (n_PARTIAL_ALIGNED == 1 and
              n_EMPTY == self.cardinality()-1):

            self._children = set_EMPTY + set_PARTIAL_ALIGNED
            return (PARTIAL, ALIGNED)


        ################################################################
        # 2/2                                                          #
        #                                                              #
        # From now on, there are at most two partial pq-trees and all  #
        # of them have v aligned to their right                        #
        #                                                              #
        # We now want to order them in such a way that all the         #
        # elements containing v are located on the right               #
        ################################################################

        else:

            self._children = []

            # We first move the empty elements to the left, if any

            if n_EMPTY > 0:
                self._children.extend(set_EMPTY)

            # If there is one partial element we but have to add it to
            # the sequence, then add all the full elements

            # We must also make sure these elements will not be
            # reordered in such a way that the elements containing v
            # are not contiguous

            # ==> We create a Q-tree

            if n_PARTIAL_ALIGNED < 2:

                new = []


                # add the partial element, if any
                if n_PARTIAL_ALIGNED == 1:

                    subtree = set_PARTIAL_ALIGNED[0]
                    new.extend(subtree.simplify(v, right = ALIGNED))


                # Then the full elements, if any, in a P-tree (we can
                # permute any two of them while keeping all the
                # elements containing v on an interval

                if n_FULL > 0:

                    new.append(new_P(set_FULL))

                # We lock all of them in a Q-tree

                self._children.append(new_Q(new))

                return PARTIAL, True

            # If there are 2 partial elements, we take care of both
            # ends. We also know it will not be possible to align the
            # interval of sets containing v to the right

            else:

                new = []

                # The second partal element is aligned to the right
                # while, as we want to put it at the end of the
                # interval, it should be aligned to the left
                set_PARTIAL_ALIGNED[1].reverse()

                # 1/3
                # Left partial subtree
                subtree = set_PARTIAL_ALIGNED[0]
                new.extend(subtree.simplify(v, right = ALIGNED))

                # 2/3
                # Center (Full elements, in a P-tree, as they can be
                # permuted)

                if n_FULL > 0:
                    new.append(new_P(set_FULL))

                # 3/3
                # Right partial subtree
                subtree = set_PARTIAL_ALIGNED[1]
                new.extend(subtree.simplify(v, left= ALIGNED))

                # We add all of it, locked in a Q-Tree
                self._children.append(new_Q(new))

                return PARTIAL, False

class Q(PQ):
    r"""
    A Q-Tree is a PQ-Tree whose children are
    ordered up to reversal
    """

    def set_contiguous(self, v):
        r"""
        Updates ``self`` so that its sets containing ``v`` are
        contiguous for any admissible permutation of its subtrees.

        This function also ensures, whenever possible,
        that all the sets containing ``v`` are located on an interval
        on the right side of the ordering.

        INPUT:

        - ``v`` -- an element of the ground set

        OUTPUT:

        According to the cases :

            * ``(EMPTY, ALIGNED)`` if no set of the tree contains
              an occurrence of ``v``

            * ``(FULL, ALIGNED)`` if all the sets of the tree contain
              ``v``

            * ``(PARTIAL, ALIGNED)`` if some (but not all) of the sets
              contain ``v``, all of which are aligned
              to the right of the ordering at the end when the function ends

            * ``(PARTIAL, UNALIGNED)`` if some (but not all) of the
              sets contain ``v``, though it is impossible to align them
              all to the right

        In any case, the sets containing ``v`` are contiguous when this
        function ends. If there is no possibility of doing so, the function
        raises a ``ValueError`` exception.

        EXAMPLE:

        Ensuring the sets containing ``0`` are continuous::

            sage: from sage.graphs.pq_trees import P, Q
            sage: q = Q([[2,3], Q([[3,0],[3,1]]), Q([[4,0],[4,5]])])
            sage: q.set_contiguous(0)
            (1, False)
            sage: print q
            ('Q', [{2, 3}, {1, 3}, {0, 3}, {0, 4}, {4, 5}])

        Impossible situation::

            sage: p = Q([[0,1], [1,2], [2,0]])
            sage: p.set_contiguous(0)
            Traceback (most recent call last):
            ...
            ValueError: Impossible


        """


        #################################################################
        # Guidelines :                                                  #
        #                                                               #
        # As the tree is a Q-Tree, we can but reverse the order in      #
        # which the elements appear. It means that we can but check     #
        # the elements containing v are already contiguous (even        #
        # though we have to take special care of partial elements --    #
        # the endpoints of the interval), and answer accordingly        #
        # (partial, full, empty, aligned..). We also want to align the  #
        # elements containing v to the right if possible.               #
        ################################################################


        ###############################################################
        # Defining Variables :                                        #
        #                                                             #
        # Collecting the information of which children are FULL of v, #
        # which ones are EMPTY, PARTIAL_ALIGNED and PARTIAL_UNALIGNED #
        #                                                             #
        # Defining variables for their cardinals, just to make the    #
        # code slightly more readable :-)                             #
        ###############################################################

        seq = [set_contiguous(x, v) for x in self]
        self.flatten()
        seq = [set_contiguous(x, v) for x in self]

        f_seq = dict(zip(self, seq))

        set_FULL                = []
        set_EMPTY               = []
        set_PARTIAL_ALIGNED     = []
        set_PARTIAL_UNALIGNED   = []

        sorting = {
            (FULL, ALIGNED)     : set_FULL,
            (EMPTY, ALIGNED)    : set_EMPTY,
            (PARTIAL, ALIGNED)  : set_PARTIAL_ALIGNED,
            (PARTIAL, UNALIGNED) : set_PARTIAL_UNALIGNED
            }

        for i in self:
            sorting[f_seq[i]].append(i)

        n_FULL                  = len(set_FULL)
        n_EMPTY                 = len(set_EMPTY)
        n_PARTIAL_ALIGNED       = len(set_PARTIAL_ALIGNED)
        n_PARTIAL_UNALIGNED     = len(set_PARTIAL_UNALIGNED)

        counts = dict(map(lambda (x,y) : (x,len(y)), sorting.iteritems()))

        ###################################################################
        #                                                                 #
        # Picking the good ordering for the children :                    #
        #                                                                 #
        #                                                                 #
        # There is a possibility of aligning to the right iif             #
        # the vector can assume the form (as a regular expression) :      #
        #                                                                 #
        # (EMPTY *) PARTIAL (FULL *) Of course, each of these three       #
        # members could be empty                                          #
        #                                                                 #
        # Hence, in the following case we reverse the vector :            #
        #                                                                 #
        # * if the last element is empty (as we checked the whole         #
        #   vector is not empty                                           #
        #                                                                 #
        # * if the last element is partial, aligned, and all the          #
        #   others are full                                               #
        ###################################################################

        if (f_seq[self._children[-1]] == (EMPTY, ALIGNED) or
            (f_seq[self._children[-1]] == (PARTIAL, ALIGNED) and n_FULL == self.cardinality() - 1)):

            # We reverse the order of the elements in the SET only. Which means that they are still aligned to the right !
            self._children.reverse()


        #########################################################
        # 1/2                                                   #
        #                                                       #
        # Several easy cases where we can decide without paying #
        # attention                                             #
        #########################################################


        # Excludes the situation where there is no solution.
        # read next comment for more explanations

        if (n_PARTIAL_ALIGNED + n_PARTIAL_UNALIGNED > 2 or
            (n_PARTIAL_UNALIGNED >= 1 and n_EMPTY != self.cardinality() -1)):

            raise ValueError(impossible_msg)

        # From now on, there are at most two pq-trees which are partially filled
        # If there is one which is not aligned to the right, all the others are empty

        # First trivial case, no checking neded
        elif n_FULL == self.cardinality():
            return FULL, True

        # Second trivial case, no checking needed
        elif n_EMPTY == self.cardinality():
            return EMPTY, True

        # Third trivial case, no checking needed
        elif n_PARTIAL_UNALIGNED == 1:
            return (PARTIAL, UNALIGNED)

        # If there is just one partial element
        # and all the others are empty, we just reorder
        # the set to put it at the right end

        elif (n_PARTIAL_ALIGNED == 1 and
              n_EMPTY == self.cardinality()-1):

            if set_PARTIAL_ALIGNED[0] == self._children[-1]:
                return (PARTIAL, ALIGNED)

            else:
                return (PARTIAL, UNALIGNED)


        ##############################################################
        # 2/2                                                        #
        #                                                            #
        # We iteratively consider all the children, and check        #
        # that the elements containing v are indeed                  #
        # locate on an interval.                                     #
        #                                                            #
        # We are also interested in knowing whether this interval is #
        # aligned to the right                                       #
        #                                                            #
        # Because of the previous tests, we can assume there are at  #
        # most two partial pq-trees and all of them are aligned to   #
        # their right                                                #
        ##############################################################

        else:

            new_children = []

            # Two variables to remember where we are
            # according to the interval

            seen_nonempty = False
            seen_right_end = False


            for i in self:

                type, aligned = f_seq[i]

                # We met an empty element
                if type == EMPTY:

                    # 2 possibilities :
                    #
                    #  * we have NOT met a non-empty element before
                    #    and it just means we are looking at the
                    #    leading empty elements
                    #
                    #  * we have met a non-empty element before and it
                    #    means we will never met another non-empty
                    #    element again => we have seen the right end
                    #    of the interval

                    new_children.append(i)

                    if seen_nonempty:
                        seen_right_end = True

                # We met a non-empty element
                else:
                    if seen_right_end:
                        raise ValueError(impossible_msg)


                    if type == PARTIAL:

                        # if we see an ALIGNED partial tree after
                        # having seen a nonempty element then the
                        # partial tree must be aligned to the left and
                        # so we have seen the right end

                        if seen_nonempty and aligned:
                            i.reverse()
                            seen_right_end = True

                            # right partial subtree
                            subtree = i
                            new_children.extend(subtree.simplify(v, left = True))

                        # If we see an UNALIGNED partial element after
                        # having met a nonempty element, there is no
                        # solution to the alignment problem

                        elif seen_nonempty and not aligned:
                            raise ValueError(impossible_msg)

                        # If we see an unaligned element but no non-empty
                        # element since the beginning, we are witnessing both the
                        # left and right end

                        elif not seen_nonempty and not aligned:
                            raise ValueError("Bon, ben ca arrive O_o")
                            seen_right_end = True

                        elif not seen_nonempty and aligned:

                            # left partial subtree
                            subtree = i

                            new_children.extend(subtree.simplify(v, right = True))


                    else:
                        new_children.append(i)

                    seen_nonempty = True

            # Setting the updated sequence of children
            self._children = new_children


            # Whether we achieved an alignment to the right is the
            # complement of whether we have seen the right end

            return (PARTIAL, not seen_right_end)

