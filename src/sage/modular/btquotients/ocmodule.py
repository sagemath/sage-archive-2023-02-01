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
from sage.modular.pollack_stevens.sigma0 import Sigma0,Sigma0ActionAdjuster
from sage.categories.action import Action
from sage.modules.free_module_element import free_module_element,vector
from sage.modules.free_module import FreeModule
import operator

# Need this to be pickleable
class _btquot_adjuster(Sigma0ActionAdjuster):
    """
    Callable object that turns matrices into 4-tuples.

    Since the modular symbol and harmonic cocycle code use different
    conventions for group actions, this function is used to make sure
    that actions are correct for harmonic cocycle computations.

    EXAMPLES::

        sage: from sage.modular.btquotients.ocmodule import _btquot_adjuster
        sage: adj = _btquot_adjuster()
        sage: adj(matrix(ZZ,2,2,[1..4]))
        (4, 2, 3, 1)
    """
    def __call__(self, g):
        """
        Turns matrices into 4-tuples.

        INPUT:

        - ``g`` - a 2x2 matrix

        OUTPUT:

        A 4-tuple encoding the entries of ``g``.

        EXAMPLES::

            sage: from sage.modular.btquotients.ocmodule import _btquot_adjuster
            sage: adj = _btquot_adjuster()
            sage: adj(matrix(ZZ,2,2,[1..4]))
            (4, 2, 3, 1)
        """
        a,b,c,d = g.list()
        return tuple([d, b, c, a])

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
                self._val = vector(self._parent._R,self._depth,val._val[:d].list()+[0]*(self._depth-d))
            else:
                try:
                    if hasattr(val,'list'):
                        val = val.list()
                    self._val = vector(self._parent._R,self._depth,val)
                except:
                    self._val= self._parent._R(val) * vector(self._parent._R,self._depth,[1]*self._depth)
        else:
            self._val= FreeModule(self._parent._R,self._depth)(val)
        self._moments = self._val

    def moment(self, i):
        return self._moments[i]

    def __getitem__(self,r):
        r"""
        Returns the value of ``self`` on the polynomial `x^r`.

        INPUT:
          - ``r`` - an integer. The power of `x`. 

        EXAMPLES:

        """
        return self._val[r]

    def __setitem__(self,r, val):
        r"""
        Sets the value of ``self`` on the polynomial `x^r` to ``val``.

        INPUT:
        - ``r`` - an integer. The power of `x`.
        - ``val`` - a value.

        EXAMPLES:

        """
        self._val[r] = val

    def element(self):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return self.matrix_rep().list()

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
        A=Matrix(self._parent._R,self._parent.dimension(),self._parent.dimension(),[[b._val[ii] for b in B] for ii in range(self._depth)])
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
        #return self._l_act_by(x[0,0],x[0,1],x[1,0],x[1,1],extrafactor=x.determinant()**(-self._nhalf))
        return x * self

    def r_act_by(self,x):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::
        """
        return x.adjoint() * self


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
            return [self._val[ii].precision_absolute() for ii in range(self._depth)]
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
            return min([self._val[ii].precision_absolute() for ii in range(self._depth)])
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
            return min([self._val[ii].precision_relative() for ii in range(self._depth)])
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
        s=str(sum([R(self._val[ii]*z**ii) for ii in range(self._depth)]))
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
                return sum([R(self._val[ii])*P[ii] for ii in range(r)])
            except NotImplementedError: pass
        return R(self._val[0])*P

    def valuation(self,l=None):
        r"""

        EXAMPLES:

        This example illustrates ...

        ::

        """
        if not self._parent.base_ring().is_exact():
            if(not l is None and l!=self._parent._R.prime()):
                raise ValueError, "This function can only be called with the base prime"
            return min([self._val[ii].valuation() for ii in range(self._depth)])
        else:
            return min([self._val[ii].valuation(l) for ii in range(self._depth)])


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
        self._act = OCVnWeightKAction(self)
        self._populate_coercion_lists_(action_list = [self._act])

    def is_overconvergent(self):
        return self._depth != self._n+1

    def _an_element_(self):
        r"""
        """
        return OCVnElement(self,vector(self._R,self._depth,range(1,self._depth+1)), check = False)

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
        self._basis=[OCVnElement(self,vector(self._R,self._depth,{jj:1},sparse=False),check = False) for jj in range(self._depth)]
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
        A=[(g * b).matrix_rep(B) for b in B] # b.l_act_by(g)
        return Matrix(self._R,d,d,[A[jj][ii] for ii in range(d) for jj in range(d)]).transpose()


class OCVnWeightKAction(Action):
    r"""

    INPUT:

    - ``Dk`` -- a space of distributions
    - ``character`` -- data specifying a Dirichlet character to apply to the
      top right corner, and a power of the determinant by which to scale.  See
      the documentation of
      :class:`sage.modular.pollack_stevens.distributions.Distributions_factory`
      for more details.
    - ``adjuster`` -- a callable object that turns matrices into 4-tuples.
    - ``on_left`` -- whether this action should be on the left.
    - ``dettwist`` -- a power of the determinant to twist by
    - ``padic`` -- if True, define an action of p-adic matrices (not just integer ones)

    OUTPUT:

    - 

    EXAMPLES::

        sage: from sage.modular.pollack_stevens.distributions import Distributions, Symk
    """
    def __init__(self, Dk):
        r"""
        Initialization.

        """
        self._k = Dk.weight()
        self._dettwist = -ZZ(self._k/2)
        self._Sigma0 = Sigma0(1, base_ring=Dk.base_ring(),adjuster = _btquot_adjuster())
        Action.__init__(self, self._Sigma0, Dk, True, operator.mul)


    def _call_(self, v, g):
        r"""
    
        EXAMPLES:
    
        This example illustrates ...
    
        ::
    
        """
        if self.is_left():
            v,g = g,v

        a,b,c,d = g.matrix().list()
        extrafactor = (a*d - b*c)**self._dettwist
        R=v._parent._R
        if(R.base_ring().is_exact()):
            factor=1
        else:
            t=min([R(x).valuation() for x in [a,b,c,d] if x!=0])
            factor=R.prime()**(-t)
        try:
            x=v._parent._powers[(factor*a,factor*b,factor*c,factor*d)]
            return v.__class__(v._parent,(extrafactor*factor**(-v._n))*(x*v._val), check = False)
        except KeyError:
            tmp = v._parent._get_powers_and_mult(factor*a,factor*b,factor*c,factor*d,extrafactor*factor**(-v._n),v._val)
    
            return v.__class__(v._parent,tmp)
