#########################################################################
#       Copyright (C) 2011 Cameron Franc and Marc Masdeu
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#########################################################################

from sage.structure.element import ModuleElement
from sage.modules.module import Module
from sage.matrix.constructor import Matrix
from sage.matrix.matrix_space import MatrixSpace
from copy import copy
from sage.rings.finite_rings.integer_mod_ring import Zmod
from sage.rings.all import Integer
from sage.rings.power_series_ring import PowerSeriesRing
from sage.structure.unique_representation import UniqueRepresentation
from sage.rings.rational_field import QQ
from sage.rings.integer_ring import ZZ
from sage.rings.padics.padic_generic import pAdicGeneric
from sage.categories.pushout import pushout

class OCVnElement(ModuleElement):
    r"""
    This class represents elements in an overconvergent coefficient module.

    INPUT:

     - ``parent`` - An overconvergent coefficient module.

     - ``val`` - The value that it needs to store (default: 0). It can be another OCVnElement,
       in which case the values are copied. It can also be a column vector (or something
       coercible to a column vector) which represents the values of the element applied to
       the polynomials `1`, `x`, `x^2`, ... ,`x^n`.

     - ``check`` - boolean (default: True). If set to False, no checks are done and ``val`` is
       assumed to be the a column vector.

    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,parent,val = 0,check = False):
        ModuleElement.__init__(self,parent)
        self._parent=parent
        self._n=self._parent._n
        self._nhalf=Integer(self._n/2)
        self._depth=self._parent._depth
        if check:
            if isinstance(val,self.__class__):
                d=min([val._parent._depth,parent._depth])
                assert(val._parent.weight()==parent.weight())
                self._val=Matrix(self._parent._R,self._depth,1,0)
                for ii in range(d):
                    self._val[ii,0]=val._val[ii,0]
            else:
                try:
                    self._val = Matrix(self._parent._R,self._depth,1,val)
                except:
                    self._val= self._parent._R(val) * MatrixSpace(self._parent._R,self._depth,1)(1)
        else:
            self._val= MatrixSpace(self._parent._R,self._depth,1)(val)
        self.moments = self._val

    def moment(self, i):
        return self.moments[i,0]

    def __getitem__(self,r):
        r"""
        Returns the value of ``self`` on the polynomial `x^r`.

        INPUT:
          - ``r`` - an integer. The power of `x`. 

        EXAMPLES:

        """
        return self._val[r,0]

    def __setitem__(self,r, val):
        r"""
        Sets the value of ``self`` on the polynomial `x^r` to ``val``.

        INPUT:
        - ``r`` - an integer. The power of `x`.
        - ``val`` - a value.

        EXAMPLES:

        """
        self._val[r,0] = val

    def element(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        tmp=self.matrix_rep()
        return [tmp[ii,0] for ii in range(tmp.nrows())]

    def list(self):
        r"""
        EXAMPLES:

        This example illustrates ...

        ::
        """
        return self.element()

    def matrix_rep(self,B=None):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #Express the element in terms of the basis B
        if(B is None):
            B=self._parent.basis()
        A=Matrix(self._parent._R,self._parent.dimension(),self._parent.dimension(),[[b._val[ii,0] for b in B] for ii in range(self._depth)])
        tmp=A.solve_right(self._val)
        return tmp

    def _add_(self,y):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        val=self._val+y._val
        return self.__class__(self._parent,val, check = False)

    def _sub_(self,y):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        val=self._val-y._val
        return self.__class__(self._parent,val, check = False)

    def l_act_by(self,x):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #assert(x.nrows()==2 and x.ncols()==2) #An element of GL2
        return self._l_act_by(x[0,0],x[0,1],x[1,0],x[1,1],extrafactor=x.determinant()**(-self._nhalf))

    def r_act_by(self,x):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #assert(x.nrows()==2 and x.ncols()==2) #An element of GL2
        return self._l_act_by(x[1,1],-x[0,1],-x[1,0],x[0,0],extrafactor=x.determinant()**(-self._nhalf))

    def _l_act_by(self,a,b,c,d,extrafactor=1):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        R=self._parent._R
        if(self._parent.base_ring().is_exact()):
            factor=1
        else:
            t=min([R(x).valuation() for x in [a,b,c,d] if x!=0])
            factor=R.prime()**(-t)
        try:
            x=self._parent._powers[(factor*a,factor*b,factor*c,factor*d)]
            return self.__class__(self._parent,(extrafactor*factor**(-self._n))*(x*self._val), check = False)
        except KeyError:
            tmp = self._parent._get_powers_and_mult(factor*a,factor*b,factor*c,factor*d,extrafactor*factor**(-self._n),self._val)

            return self.__class__(self._parent,tmp)

    def _rmul_(self,a):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #assume that a is a scalar
        return self.__class__(self._parent,a*self._val, check = False)

    def precision_absolute(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #This needs to be thought more carefully...
        if not self._parent.base_ring().is_exact():
            return [self._val[ii,0].precision_absolute() for ii in range(self._depth)]
        else:
            return Infinity

    def precision(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        #This needs to be thought more carefully...
        if not self._parent.base_ring().is_exact():
            return min([self._val[ii,0].precision_absolute() for ii in range(self._depth)])
        else:
            return Infinity

    def precision_relative(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        #This needs to be thought more carefully...
        if not self._parent.base_ring().is_exact():
            return min([self._val[ii,0].precision_relative() for ii in range(self._depth)])
        else:
            return Infinity

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES:

        This example illustrates ...

        ::

        """
        R=PowerSeriesRing(self._parent._R,default_prec=self._depth,name='z')
        z=R.gen()
        s=str(sum([R(self._val[ii,0]*z**ii) for ii in range(self._depth)]))
        return s

    def __cmp__(self,other):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        return cmp(self._val,other._val)

    def __nonzero__(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return self._val!=0

    def evaluate_at_poly(self,P):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        p = self._parent._R.prime()
        try:
            R = pushout(P.parent().base_ring(),self.parent().base_ring())
        except AttributeError:
            R = self.parent().base_ring()

        if hasattr(P,'degree'):
            try:
                r = min([P.degree()+1,self._depth])
                return sum([R(self._val[ii,0])*P[ii] for ii in range(r)])
            except NotImplementedError: pass
        return R(self._val[0,0])*P

    def valuation(self,l=None):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        if not self._parent.base_ring().is_exact():
            if(not l is None and l!=self._parent._R.prime()):
                raise ValueError, "This function can only be called with the base prime"
            return min([self._val[ii,0].valuation() for ii in range(self._depth)])
        else:
            return min([self._val[ii,0].valuation(l) for ii in range(self._depth)])


class OCVn(Module,UniqueRepresentation):
    Element=OCVnElement
    r"""
    This class represents objects in the overconvergent approximation modules used to
    describe overconvergent p-adic automorphic forms. 

    INPUT:

     - ``n`` - integer 

     - ``R`` - ring

     - ``depth`` - integer (Default: None)

     - ``basis`` - (Default: None)


    AUTHORS:

    - Cameron Franc (2012-02-20)
    - Marc Masdeu (2012-02-20)
    """
    def __init__(self,n,R,depth=None,basis=None):
        Module.__init__(self,base=R)
        if basis is not None:
            self._basis=copy(basis)
        self._n=n
        self._R=R
        if R.is_exact():
            self._Rmod=self._R
        else:
            self._Rmod=Zmod(self._R.prime()**(self._R.precision_cap()))

        if depth is None:
            depth=n+1
        if depth != n+1:
            if R.is_exact(): raise ValueError, "Trying to construct an over-convergent module with exact coefficients, how do you store p-adics ??"
        self._depth=depth
        self._PowerSeries=PowerSeriesRing(self._Rmod,default_prec=self._depth,name='z')
        self._powers=dict()
        self._populate_coercion_lists_()

    def is_overconvergent(self):
        return self._depth != self._n+1

    def _an_element_(self):
        r"""
        """
        return OCVnElement(self,Matrix(self._R,self._depth,1,range(1,self._depth+1)), check = False)

    def _coerce_map_from_(self, S):
        r"""

        EXAMPLES:

        ::

        """
        # Nothing coherces here, except OCVnElement
        return False

    def _element_constructor_(self,x,check = True):
        r"""

        EXAMPLES:

        """
        #Code how to coherce x into the space
        #Admissible values of x?
        return OCVnElement(self,x,check)

    def _get_powers_and_mult(self,a,b,c,d,lambd,vect):
        r"""
        Compute the action of a matrix on the basis elements.

        EXAMPLES:

        ::

        """
        R=self._PowerSeries
        r=R([b,a])
        s=R([d,c])
        n=self._n
        if(self._depth==n+1):
            rpows=[R(1)]
            spows=[R(1)]
            for ii in range(n):
                rpows.append(r*rpows[ii])
                spows.append(s*spows[ii])
            x=Matrix(self._Rmod,n+1,n+1,0)
            for ii in range(n+1):
                y=rpows[ii]*spows[n-ii]
                for jj in range(self._depth):
                    x[ii,jj]=y[jj]
        else:
            ratio=r*(s**(-1))
            y=s**n
            x=Matrix(self._Rmod,self._depth,self._depth,0)
            for jj in range(self._depth):
                x[0,jj]=y[jj]
            for ii in range(1,self._depth):
                y*=ratio
                for jj in range(self._depth):
                    x[ii,jj]=y[jj]
        if self._Rmod is self._R:
            xnew=x
        else:
            xnew=x.change_ring(self._R.base_ring())
            xnew=xnew.change_ring(self._R)
        self._powers[(a,b,c,d)]=xnew
        return self._R(lambd) * xnew * vect

    def _repr_(self):
        r"""
        This returns the representation of self as a string.

        EXAMPLES:

        """
        if self.is_overconvergent():
            return "Space of %s-adic distributions with k=%s action and precision cap %s"%(self._R.prime(), self._n, self._depth - 1)
        else:
            if self.base_ring() is QQ:
                V = 'Q^2'
            elif self.base_ring() is ZZ:
                V = 'Z^2'
            elif isinstance(self.base_ring(), pAdicGeneric) and self.base_ring().degree() == 1:
                if self.base_ring().is_field():
                    V = 'Q_%s^2'%(self._R.prime())
                else:
                    V = 'Z_%s^2'%(self._R.prime())
            else:
                V = '(%s)^2'%(self.base_ring())
            return "Sym^%s %s"%(self._n, V)
        # s='Overconvergent coefficient module of weight n = %s over the ring %s and depth %s'%(self._n,self._R,self._depth)
        return s

    def basis(self):
        r"""
        A basis of the module.

        INPUT:

         - ``x`` - integer (default: 1) the description of the
           argument x goes here.  If it contains multiple lines, all
           the lines after the first need to be indented.

         - ``y`` - integer (default: 2) the ...

        OUTPUT:

        integer -- the ...

        EXAMPLES:


        """
        try: return self._basis
        except: pass
        self._basis=[OCVnElement(self,Matrix(self._R,self._depth,1,{(jj,0):1},sparse=False),check = False) for jj in range(self._depth)]
        return self._basis

    def base_ring(self):
        r"""
        This function returns the base ring of the overconvergent element.

        EXAMPLES::

        This example illustrates ...

        ::

        """
        return self._R

    def depth(self):
        r"""
        Returns the depth of the module.
        """
        return self._depth

    def dimension(self):
        r"""
        Returns the dimension (rank) of the module.
        """
        return self._depth

    def precision_cap(self):
        r"""
        Returns the dimension (rank) of the module.
        """
        return self._depth

    def weight(self):
        r"""
        Returns the cohomological weight of the automorphic form.
        """
        return self._n

    def acting_matrix(self,g,d,B=None):
        r"""
        Matrix representation of ``g`` in a given basis.

        """
        if d is None:
            d = self.dimension()
        if B is None:
            B=self.basis()
        A=[(b.l_act_by(g)).matrix_rep(B) for b in B]
        return Matrix(self._R,d,d,[A[jj][ii,0] for ii in range(d) for jj in range(d)]).transpose()


