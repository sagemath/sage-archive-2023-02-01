"""
Braid groups

Braid groups are implemented as a particular case of finitely presented groups,
but with a lot of specific methods for braids.

A braid group can be created by giving the number of strands, and the name of the generators:

EXAMPLES::

    sage: BraidGroup(3)
    Braid group on 3 strands
    sage: BraidGroup(3,'a')
    Braid group on 3 strands
    sage: BraidGroup(3,'a').gens()
    (a0, a1)
    sage: BraidGroup(3,'a,b').gens()
    (a, b)


The elements can be created by operating with the generators, or by passing a list
with the indices of the letters to the group:

EXAMPLES::

    sage: B=BraidGroup(4)
    sage: B.inject_variables()
    Defining s0, s1, s2
    sage: s0*s1*s0
    s0*s1*s0
    sage: B([1,2,1])
    s0*s1*s0

The mapping class action of the braid group over the free group is also implemented.

EXAMPLES::

    sage: B=BraidGroup(4)
    sage: B.inject_variables()
    Defining s0, s1, s2
    sage: F=FreeGroup(4)
    sage: F.inject_variables()
    Defining x0, x1, x2, x3
    sage: A=MappingClassGroupAction(B, F)
    sage: A
    Right action by Braid group on 4 strands on Free Group on generators ('x0', 'x1', 'x2', 'x3')
    sage: A(x0, s1)
    x0
    sage: A(x1, s1)
    x1*x2*x1^-1
    sage: A(x1^-1, s1)
    x1*x2^-1*x1^-1

AUTHOR:

- Miguel Angel Marco Buzunariz
"""

##############################################################################
#       Copyright (C) 2012 Miguel Angel Marco Buzunariz <mmarco@unizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
##############################################################################


from sage.groups import group
from sage.structure.element import Element, MultiplicativeGroupElement
from sage.categories.basic import Groups
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.parent_gens import ParentWithGens
from sage.interfaces.gap import gap
from sage.rings.integer import Integer
from sage.rings.integer_ring import IntegerRing
from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
from sage.matrix.constructor import identity_matrix, matrix
from sage.combinat.permutation import Permutation
from sage.categories.action import Action
from sage.plot.bezier_path import bezier_path
from sage.plot.plot import Graphics, line
from sage.sets.set import Set
from sage.rings.fraction_field import FractionField
from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
from sage.plot.plot3d.shapes2 import bezier3d
from sage.rings.infinity import Infinity
from sage.structure.parent import normalize_names
from sage.groups.perm_gps.permgroup import PermutationGroup
from sage.groups.perm_gps.permgroup_element import PermutationGroupElement
from sage.groups.finitely_presented import FinitelyPresentedGroupElement
from sage.misc.cachefunc import cached_method
import re
import operator


class Braid(FinitelyPresentedGroupElement):
    """
    A class that models elements of the braid group. It is a particular case of element of a finitely
    presented group.

    EXAMPLES::

        sage: B=BraidGroup(4)
        sage: B
        Braid group on 4 strands
        sage: B.inject_variables()
        Defining s0, s1, s2
        sage: s0*s1/s2/s1
        s0*s1*s2^-1*s1^-1
        sage: B((1, 2, -3, -2))
        s0*s1*s2^-1*s1^-1
    """

    def __cmp__(self, other):
        """
        TESTS::

            sage: B=BraidGroup(4)
            sage: b=B([1, 2, 1])
            sage: c=B([2, 1, 2])
            sage: b==c #indirect doctest
            True
            sage: b.__cmp__(c^(-1))
            -1
            sage: B([]).__cmp__(B.one())
            0
        """
        if self.TietzeList()==other.TietzeList():
            return 0
        nfself = map(lambda i: i.TietzeList(), self.left_normal_form())
        nfother = map(lambda i: i.TietzeList(), other.left_normal_form())
        return cmp(nfself, nfother)

    def _latex_(self):
        """
        TESTS::

            sage: B=BraidGroup(4)
            sage: b=B([1, 2, 3, -1, 2, -3])
            sage: b._latex_()
            '\\sigma_{1}\\sigma_{2}\\sigma_{3}\\sigma_{1}^{-1}\\sigma_{2}\\sigma_{3}^{-1}'
        """
        latexrepr=''
        for i in self.TietzeList():
            if i>0:
                latexrepr=latexrepr+"\sigma_{%s}"%i
            if i<0:
                latexrepr=latexrepr+"\sigma_{%s}^{-1}"%(-i)
        return latexrepr

    def strands(self):
        """
        returns the number of strands in the braid.

        EXAMPLES::

            sage: B=BraidGroup(4)
            sage: b=B([1, 2, -1, 3, -2])
            sage: b.strands()
            4
        """
        return self.parent().strands()

    def burau_matrix(self, var='t'):
        """
        Returns the Burau matrix of the braid.

        INPUT:

        - var (default='t'): the name of the variable in the entries of the matrix.

        OUTPUT: The Burau matrix of the braid. It is a matrix whose entries are Laurent polynomials in the variable var.

        EXAMPLES::

            sage: B=BraidGroup(4)
            sage: B.inject_variables()
            Defining s0, s1, s2
            sage: b=s0*s1/s2/s1
            sage: b.burau_matrix()
            [     -t + 1           0    -t^2 + t         t^2]
            [          1           0           0           0]
            [          0           0           1           0]
            [          0        t^-2 t^-1 - t^-2    1 - t^-1]
            sage: s2.burau_matrix('x')
            [     1      0      0      0]
            [     0      1      0      0]
            [     0      0 -x + 1      x]
            [     0      0      1      0]

        REFERENCES:

            http://en.wikipedia.org/wiki/Burau_representation
        """
        R=LaurentPolynomialRing(IntegerRing(), var)
        t=R.gen()
        M=identity_matrix(R, self.strands())
        for i in self.TietzeList():
            A=identity_matrix(R, self.strands())
            if i>0:
                A[i-1, i-1]=1-t
                A[i, i]=0
                A[i, i-1]=1
                A[i-1, i]=t
            if i<0:
                A[-1-i, -1-i]=0
                A[-i, -i]=1-t**(-1)
                A[-1-i, -i]=1
                A[-i, -1-i]=t**(-1)
            M=M*A
        return M

    def permutation(self):
        """
        Returns the permutation induced by the braid in its strands.

        EXAMPLES::

            sage: B=BraidGroup(4)
            sage: B.inject_variables()
            Defining s0, s1, s2
            sage: b=s0*s1/s2/s1
            sage: b.permutation()
            [4, 1, 3, 2]
        """
        per=Permutation((()))
        for i in self.TietzeList():
            j=abs(i)
            per=per*Permutation(((j, j+1)))
        return per

    def plot(self, color='blue', orientation='bottom-top', gap=0.05, aspect_ratio=1, axes=False, **kwds):
        """
        Plots the braid

        The following options are available:

         -``orientation`` - (default: 'bottom-top') determines how the braid is printed. The possible
         values are:

          -'bottom-top', the braid is printed from bottom to top

          -'top-bottom', the braid is printed from top to bottom

          -'left-right', the braid is printed from left to right

         -``gap`` - (default: 0.05) determines the size of the gap left when a strand goes under another.

        EXAMPLES::

            sage: B=BraidGroup(4, 's')
            sage: b=B([1, 2, 3, 1, 2, 1])
            sage: b.plot()

        """
        if orientation=='top-bottom':
            orx=0
            ory=-1
            nx=1
            ny=0
        elif orientation=='left-right':
            orx=1
            ory=0
            nx=0
            ny=-1
        elif orientation=='bottom-top':
            orx=0
            ory=1
            nx=1
            ny=0
        col=color
        br=self.TietzeList()
        n=self.strands()
        a=Graphics()
        op=gap
        for i in range(len(br)):
            m=br[i]
            for j in range(n):
                if m==j+1:
                    a+=bezier_path([[(j*nx+i*orx, i*ory+j*ny), (j*nx+orx*(i+0.25), j*ny+ory*(i+0.25)), (nx*(j+0.5)+orx*(i+0.5), ny*(j+0.5)+ory*(i+0.5))], [(nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)), (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col, **kwds)
                elif m==j:
                    a+=bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)), (nx*(j-0.5+4*op)+orx*(i+0.5-2*op), ny*(j-0.5+4*op)+ory*(i+0.5-2*op)), (nx*(j-0.5+2*op)+orx*(i+0.5-op), ny*(j-0.5+2*op)+ory*(i+0.5-op))]], color=col, **kwds)
                    a+=bezier_path([[(nx*(j-0.5-2*op)+orx*(i+0.5+op), ny*(j-0.5-2*op)+ory*(i+0.5+op)), (nx*(j-0.5-4*op)+orx*(i+0.5+2*op), ny*(j-0.5-4*op)+ory*(i+0.5+2*op)), (nx*(j-1)+orx*(i+0.75), ny*(j-1)+ory*(i+0.75)), (nx*(j-1)+orx*(i+1), ny*(j-1)+ory*(i+1))]], color=col, **kwds)
                elif -m==j+1:
                    a+=bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)), (nx*(j+0.5-4*op)+orx*(i+0.5-2*op), ny*(j+0.5-4*op)+ory*(i+0.5-2*op)), (nx*(j+0.5-2*op)+orx*(i+0.5-op), ny*(j+0.5-2*op)+ory*(i+0.5-op))]], color=col, **kwds)
                    a+=bezier_path([[(nx*(j+0.5+2*op)+orx*(i+0.5+op), ny*(j+0.5+2*op)+ory*(i+0.5+op)), (nx*(j+0.5+4*op)+orx*(i+0.5+2*op), ny*(j+0.5+4*op)+ory*(i+0.5+2*op)), (nx*(j+1)+orx*(i+0.75), ny*(j+1)+ory*(i+0.75)), (nx*(j+1)+orx*(i+1), ny*(j+1)+ory*(i+1))]], color=col, **kwds)
                elif -m==j:
                    a+=bezier_path([[(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+0.25), ny*j+ory*(i+0.25)), (nx*(j-0.5)+orx*(i+0.5), ny*(j-0.5)+ory*(i+0.5))], [(nx*(j-1)+orx*(i+0.75), ny*(j-1)+ory*(i+0.75)), (nx*(j-1)+orx*(i+1), ny*(j-1)+ory*(i+1))]], color=col, **kwds)
                else:
                    a+=line([(nx*j+orx*i, ny*j+ory*i), (nx*j+orx*(i+1), ny*j+ory*(i+1))], color=col, **kwds)
        a.set_aspect_ratio(aspect_ratio)
        a.axes(axes)
        return a

    def plot3d(self):
        """
        Plots the braid in 3d.

        EXAMPLES::

            sage: B=BraidGroup(4, 's')
            sage: b=B([1, 2, 3, 1, 2, 1])
            sage: b.plot3d()

        """
        b=[]
        braid=self.TietzeList()
        n=self.strands()
        for i in range(len(braid)):
            m=braid[i]
            for j in range(n):
                if m==j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j+0.5, i+0.5)], [(0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]]))
                elif -m==j+1:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j+0.5, i+0.5)], [(-0.25, j+1, i+0.75), (0, j+1, i+0.75), (0, j+1, i+1)]]))
                elif m==j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (-0.25, j, i+0.25), (-0.25, j-0.5, i+0.5)], [(-0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]]))
                elif -m==j:
                    b.append(bezier3d([[(0, j, i), (0, j, i+0.25), (0.25, j, i+0.25), (0.25, j-0.5, i+0.5)], [(0.25, j-1, i+0.75), (0, j-1, i+0.75), (0, j-1, i+1)]]))
                else:
                    b.append(bezier3d([[(0, j, i), (0, j, i+1)]]))
        return sum(b)

    @cached_method
    def LKB_matrix(self, variables='x,y'):
        """
        The matrix corresponding to the Lawence-Krammer-Bigelow representation of the braid.
        The matrix is expressed in the basis $\{e_{i, j} \mid 1\\leq i < j \leq n\}$, where
        the indices are ordered lexicographically.
        It is a matrix whose entries are in the ring of Laurent polynomials on the given variables.
        By default, the variables are 'x' and 'y'.

        INPUT:

        -variables (default='x,y'): a string containing the names of the variables, separated by a comma.

        OUTPUT: The matrix corresponding to the Lawence-Krammer-Bigelow representation of the braid.

        EXAMPLES::

            sage: B=BraidGroup(3)
            sage: b=B([1, 2, 1])
            sage: b.LKB_matrix()
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: c=B([2, 1, 2])
            sage: c.LKB_matrix()
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]

        REFERENCES:

        .. [Bigelow] Bigelow, Stephen J. The Lawrence-Krammer representation. arXiv:math/0204057v1

        """
        return self.parent()._LKB_matrix_(tuple(self.TietzeList()), variab=variables)

    @cached_method
    def left_normal_form(self):
        """
        Returns the left normal form of the braid.

        The output is a list of the elements of the left normal form. The first element is a
        power of $\Delta$, and the rest are permutation braids.

        EXAMPLES::


            sage: B=BraidGroup(4)
            sage: b=B([1, 2, 3, -1, 2, -3])
            sage: b.left_normal_form()
            [s0^-1*s1^-1*s2^-1*s0^-1*s1^-1*s0^-1, s0*s1*s2*s1*s0, s0*s2*s1]
            sage: c=B([1])
            sage: c.left_normal_form()
            [<identity ...>, s0]
        """
        lnfp=self._left_normal_form_perm_()
        a=lnfp[0]
        l=lnfp[1:]
        n=self.strands()
        delta=Permutation([n-i for i in range(n)])
        return [self.parent()._permutation_braid(delta).__pow__(a)]+map(lambda i:self.parent()._permutation_braid(i), l)

    @cached_method
    def _left_normal_form_perm_(self):
        """
        Returns the left normal form of the braid, in permutation form.

        The output is a list whose first element is the power of $\Delta$, and the rest are the permutations corresponding to the simple factors.

        EXAMPLES::

            sage: B=BraidGroup(12)
            sage: B([2, 2, 2, 3, 1, 2, 3, 2, 1, -2])._left_normal_form_perm_()
            [-1, [12, 11, 10, 9, 8, 7, 6, 5, 2, 4, 3, 1], [4, 1, 3, 2, 5, 6, 7, 8, 9, 10, 11, 12], [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12], [3, 1, 2, 4, 5, 6, 7, 8, 9, 10, 11, 12], [2, 3, 1, 4, 5, 6, 7, 8, 9, 10, 11, 12]]
            sage: C=BraidGroup(6)
            sage: C([2, 3, -4, 2, 3, -5, 1, -2, 3, 4, 1, -2])._left_normal_form_perm_()
            [-2, [3, 5, 4, 2, 6, 1], [1, 6, 3, 5, 2, 4], [5, 6, 2, 4, 1, 3], [3, 2, 4, 1, 5, 6], [1, 5, 2, 3, 4, 6]]
        """
        n=self.parent().strands()
        delta=0
        Delta=Permutation([n-i for i in range(n)])
        l=list(self.TietzeList())
        if l==[]:
            return [0]
        form=[]
        for i in l:
            if i>0:
                form.append(Permutation((i, i+1)))
            else:
                delta=delta+1
                form=map(lambda a: Delta*a*Delta, form)
                form.append(Delta*Permutation((-i, -i+1)))
        i=0
        j=0
        while j<len(form):
            while i<len(form)-j-1:
                e=form[i].inverse().descents()
                s=form[i+1].descents()
                S=set(s).difference(set(e))
                while S!=set([]):
                    a=list(S)[0]
                    form[i]=form[i]*Permutation((a+1, a+2))
                    form[i+1]=Permutation((a+1, a+2))*form[i+1]
                    e=form[i].inverse().descents()
                    s=form[i+1].descents()
                    S=set(s).difference(set(e))
                if form[i+1].length()==0:
                    form.pop(i+1)
                    i=0
                else:
                    i=i+1
            j=j+1
            i=0
        form=filter(lambda a: a.length()>0, form)
        while form!=[] and form[0]==Delta:
            form.pop(0)
            delta=delta-1
        return [-delta]+form


class BraidGroup(FinitelyPresentedGroup):
    """
    The braid group on n strands.

    EXAMPLES::

        sage: B1=BraidGroup(5)
        sage: B1
        Braid group on 5 strands
        sage: B2=BraidGroup(3)
        sage: B1==B2
        False
        sage: B2 is BraidGroup(3)
        True
    """
    Element=Braid

    def __init__(self, n, names='s'):
        """
        TESTS::


            sage: B1=BraidGroup(5) # indirect doctest
            sage: B1
            Braid group on 5 strands
        """
        if not n in IntegerRing() or n<2:
            raise ValueError, "n must be an integer bigger than one"
        self._freegroup_=FreeGroup(n-1, names)
        self._names_=names
        rels=[]
        for i in range(1, n):
            if i<n-1:
                rels.append(self._freegroup_([i, i+1, i, -i-1, -i, -i-1]))
            for j in range(i+2, n):
                rels.append(self._freegroup_([i, j, -i, -j]))
        self._rels_=tuple(rels)
        lis=libgap([i.gap() for i in rels])
        self._gap_repr_=self._freegroup_.gap()/lis
        #group.Group.__init__(self)
        #self._assign_names(self._freegroup_._gens_str_)
        Parent.__init__(self,gens=self._freegroup_._gens_str_,category=Groups())
        self._nstrands_=n

    def __reduce__(self):
        """
        TESTS::

            sage: B=BraidGroup(3)
            sage: B.__reduce__()
            (<class 'sage.groups.braid.BraidGroup'>, (3,))
            sage: B=BraidGroup(3, 'sigma')
            sage: B.__reduce__()
            (<class 'sage.groups.braid.BraidGroup'>, (3, 'sigma'))

        """
        if self._names_=='s':
            return (BraidGroup, (self._nstrands_, ))
        return (BraidGroup, (self._nstrands_, self._names_))

    def _repr_(self):
        """
        TESTS::

            sage: B1=BraidGroup(5)
            sage: B1 # indirect doctest
            Braid group on 5 strands
        """
        return "Braid group on %s strands"%self._nstrands_

    def size(self):
        """
        TESTS::

            sage: B1=BraidGroup(5)
            sage: B1.size()
            +Infinity
        """
        return Infinity
    cardinality=size
    order=size

    def permutation_group(self):
        """
        Returns an error, since braid groups are infinite.

        TESTS::

            sage: B=BraidGroup(4, 'g')
            sage: try:
            ...     B.permutation_group()
            ... except ValueError as exc:
            ...     print exc
            The group is infinite
        """
        raise ValueError, "The group is infinite"

    def strands(self):
        """
        Returns the number of strands of the braid group

        EXAMPLES::

            sage: B=BraidGroup(4)
            sage: B.strands()
            4
        """
        return self._nstrands_

    def _element_constructor_(self, x):
        """
        TESTS::

            sage: B=BraidGroup(4)
            sage: B([1, 2, 3]) # indirect doctest
            s0*s1*s2
        """
        return Braid(x, parent=self)

    @cached_method
    def _permutation_braid_list(self, p):
        """
        Constructs the braid that corresponds to the given permutation. It is the only braid with
        the following properties:

         - The braid induces the given permutation.

         - The braid is positive (that is, it can be writen without using the inverses of the generators).

         - Every two strands cross each other at most once.

        INPUT:

        -p: a permutation

        Output: the lexicographically smallest word that represents the braid, in Tietze list form.

        EXAMPLES::

            sage: B=BraidGroup(5)
            sage: P=Permutation([5, 3, 1, 2, 4])
            sage: B._permutation_braid_list(P)
            [1, 2, 1, 3, 2, 4]

        """
        if p.length()==0:
            return []
        pl=p
        l=[]
        while pl.length()>0:
            i=1
            while i<max(pl):
                if pl(i)>pl(i+1):
                    l.append(i)
                    pl=Permutation([(i, i+1)])*pl
                    i=1
                else:
                    i=i+1
        return l

    @cached_method
    def _permutation_braid(self, p):
        """
        Constructs the braid that corresponds to the given permutation. It is the only braid with
        the following properties:

         - The braid induces the given permutation.

         - The braid is positive (that is, it can be writen without using the inverses of the generators).

         - Every two strands cross each other at most once.

        INPUT:

        -p: a permutation.

        OUTPUT: The braid that corresponds to the permutation.

        EXAMPLES::

            sage: B=BraidGroup(5)
            sage: P=Permutation([5, 3, 1, 2, 4])
            sage: B._permutation_braid(P)
            s0*s1*s0*s2*s1*s3

        """
        return self(self._permutation_braid_list(p))

    @cached_method
    def _LKB_matrix_(self, braid, variab):
        """
        Returns the matrix corresponding to the Lawrence-Krammer-Bigelow representation of the braid given by the tuple.
        The variables of the matrix must be given.
        This function is included here to improve the time gained by the caching.

        INPUT:

            -braid: a tuple whose entries are the entries of the Tietze list of the braid.

            -variab: the names of the variables that will appear in the matrix. They must be given as a string, separated by a comma

        OUTPUT: The LKB matrix of the braid, with respect to the variables.

        TESTS::

            sage: B=BraidGroup(3)
            sage: B._LKB_matrix_((2, 1, 2), 'x, y')
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: B._LKB_matrix_((1, 2, 1), 'x, y')
            [             0 -x^4*y + x^3*y         -x^4*y]
            [             0         -x^3*y              0]
            [        -x^2*y  x^3*y - x^2*y              0]
            sage: B._LKB_matrix_((-1, -2, -1, 2, 1, 2), 'x, y')
            [1 0 0]
            [0 1 0]
            [0 0 1]

        """
        n=self.strands()
        if len(braid)>1:
            A=self._LKB_matrix_(braid[:1], variab)
            for i in braid[1:]:
                A=A*self._LKB_matrix_((i, ), variab)
            return A
        l=list(Set(range(n)).subsets(2))
        R=LaurentPolynomialRing(IntegerRing(), variab)
        q=R.gens()[0]
        t=R.gens()[1]
        if len(braid)==0:
            return identity_matrix(R, len(l), sparse=True)
        A=matrix(R, len(l), sparse=True)
        if braid[0]>0:
            i=braid[0]-1
            for m in range(len(l)):
                j=min(l[m])
                k=max(l[m])
                if i==j-1:
                    A[l.index(Set([i, k])), m]=q
                    A[l.index(Set([i, j])), m]=q*q-q
                    A[l.index(Set([j, k])), m]=1-q
                elif i==j and not j==k-1:
                    A[l.index(Set([j, k])), m]=0
                    A[l.index(Set([j+1, k])), m]=1
                elif k-1==i and not k-1==j:
                    A[l.index(Set([j, i])), m]=q
                    A[l.index(Set([j, k])), m]=1-q
                    A[l.index(Set([i, k])), m]=(1-q)*q*t
                elif i==k:
                    A[l.index(Set([j, k])), m]=0
                    A[l.index(Set([j, k+1])), m]=1
                elif i==j and j==k-1:
                    A[l.index(Set([j, k])), m]=-t*q*q
                else:
                    A[l.index(Set([j, k])), m]=1
            return A
        else:
            i=-braid[0]-1
            for m in range(len(l)):
                j=min(l[m])
                k=max(l[m])
                if i==j-1:
                    A[l.index(Set([j-1, k])), m]=1
                elif i==j and not j==k-1:
                    A[l.index(Set([j+1, k])), m]=q**(-1)
                    A[l.index(Set([j, k])), m]=1-q**(-1)
                    A[l.index(Set([j, j+1])), m]=t**(-1)*q**(-1)-t**(-1)*q**(-2)
                elif k-1==i and not k-1==j:
                    A[l.index(Set([j, k-1])), m]=1
                elif i==k:
                    A[l.index(Set([j, k+1])), m]=q**(-1)
                    A[l.index(Set([j, k])), m]=1-q**(-1)
                    A[l.index(Set([k, k+1])), m]=-q**(-1)+q**(-2)
                elif i==j and j==k-1:
                    A[l.index(Set([j, k])), m]=-t**(-1)*q**(-2)
                else:
                    A[l.index(Set([j, k])), m]=1
            return A
            #M=self._LKB_matrix_((-braid[0], ), variab).inverse()
            #ML=M.list()
            #MIR=FractionField(PolynomialRing(IntegerRing(), variab))
            #ML=[R(MIR(i).numerator())*R(MIR(i).denominator())**(-1) for i in ML]
            #M=matrix(len(l), len(l), ML)
            #return M


class MappingClassGroupAction(Action):
    r"""
    The action of the braid group the free group as the mapping class group of the punctured disk.
    This action goes as follows:

    .. MATH::

        x_j \cdot \sigma_i=\begin{cases}
        x_{j}\cdot x_{j+1}\cdot {x_j}^{-1} & \text{if $i=j$} \\
        x_{j-1} & \text{if $i=j-1$} \\
        x_{j} & \text{otherwise}
        \end{cases}

    Where $\sigma_i$ are the generators of the braid group on $n$ strands, and $x_j$ the generators
    of the free group of rank $n$.


    EXAMPLES::

        sage: B=BraidGroup(4)
        sage: B.inject_variables()
        Defining s0, s1, s2
        sage: F=FreeGroup(4)
        sage: F.inject_variables()
        Defining x0, x1, x2, x3
        sage: A=MappingClassGroupAction(B, F)
        sage: A
        Right action by Braid group on 4 strands on Free Group on generators ('x0', 'x1', 'x2', 'x3')
        sage: A(x0, s1)
        x0
        sage: A(x1, s1)
        x1*x2*x1^-1
        sage: A(x1^-1, s1)
        x1*x2^-1*x1^-1

    It's also possible to register the action to allow multiplication of free group elements by braids::

        sage: F._unset_coercions_used()
        sage: F.register_action(A)
        sage: x2*s2
        x2*x3*x2^-1
        sage: x2*s1^(-1)
        x2^-1*x1*x2

        """

    def __init__(self, G, M, is_left=0):
        """
        TESTS::

            sage: B=BraidGroup(3)
            sage: G=FreeGroup('a, b, c')
            sage: MappingClassGroupAction(B, G) # indirect doctest
            Right action by Braid group on 3 strands on Free Group on generators ('a', 'b', 'c')
        """
        Action.__init__(self, G, M, is_left, operator.mul)

    def _call_(self, x, b):
        """
        TESTS::

            sage: B=BraidGroup(3)
            sage: G=FreeGroup('a, b, c')
            sage: A=MappingClassGroupAction(B, G)
            sage: A(G.0, B.0) # indirect doctest
            a*b*a^-1
            sage: A(G.1, B.0) # indirect doctest
            a
        """
        t=list(x.TietzeList())
        for j in b.TietzeList():
            s=[]
            for i in t:
                if j==i and i>0:
                    s+=[i, i+1, -i]
                elif j==-i and i<0:
                    s+=[-i, i-1, i]
                elif j==-i and i>0:
                    s+=[i+1]
                elif j==i and i<0:
                    s+=[i-1]
                elif i>0 and j==i-1:
                    s+=[i-1]
                elif i<0 and j==-i-1:
                    s+=[i+1]
                elif i>0 and -j==i-1:
                    s+=[-i, i-1, i]
                elif i<0 and j==i+1:
                    s+=[i, i+1, -i]
                else:
                    s+=[i]
            t=s
        return self.codomain()(list(t))
