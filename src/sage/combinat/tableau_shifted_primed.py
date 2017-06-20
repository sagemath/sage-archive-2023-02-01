from six import add_metaclass
import numpy as np

from sage.combinat.partition import Partition,StrictPartitions, OrderedPartitions
from sage.combinat.integer_vector import IntegerVectors

from sage.categories.finite_enumerated_sets import FiniteEnumeratedSets
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass

from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.structure.unique_representation import UniqueRepresentation

# Imports for the crystal

from sage.categories.classical_crystals import ClassicalCrystals
from sage.combinat.root_system.cartan_type import CartanType


@add_metaclass(InheritComparisonClasscallMetaclass)
class ShiftedPrimedTableau(ClonableArray):
    r"""
    A shifted primed tableau with primed elements stored as half-integers 

    A primed tableau is a shifted tableau on the alphabet 
        X' = {1' < 1 < 2' < 2 <...< n' < n} such that 
        1). the entries are weakly increasing along rows and columns
        2). a row can't have two repeated primed elements, and a column
            can't have two repeated non-primed elements
        3). there are only non-primed elements on the main diagonal


    EXAMPLES::

        sage: T = ShiftedPrimedTableaux([4,2])
        sage: T([[1,"2'","3'",3],[2,"3'"]])[1]
        (2.0, 2.5)
        sage: t = ShiftedPrimedTableau([[1,"2p",2.5,3],[0,2,2.5]])
        sage: t[1]
        (2.0, 2.5)
        sage: t = ShiftedPrimedTableau([[1,2,2.5,3],[0,2,2.5]])
        Traceback (most recent call last):
        ...
        ValueError: [(1.0, 2.0, 2.5, 3.0), (2.0, 2.5)] is not an element of Shifted Primed Tableaux
    """
    
    @staticmethod
    def __classcall_private__(cls, T):
        r"""
        This ensures that a shifted tableau is only ever constructed
        as an ``element_class`` call of an appropriate parent.

        EXAMPLES::

            sage: data = [[1,"2'","2",3],[2,"3'"]]
            sage: t = ShiftedPrimedTableau(data)
            sage: T = ShiftedPrimedTableaux(shape=[4,2],weight=(1,3,2))
            sage: t == T(data)
            True
            sage: S = ShiftedPrimedTableaux(shape=[4,2])
            sage: t == S(data)
            True

        """
    
        if isinstance(T, cls):
            return T

        t = preprocessing(T)

        shape = [len(row) for row in t]
        flat = [round(item) for sublist in t for item in sublist]
        if flat == []:
            max_ind = 0
        else:
            max_ind = int(max(flat))
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        
            
        return ShiftedPrimedTableaux(shape = shape, weight = weight)(T)
    

    def __init__(self,parent, T):
        r"""
        Preprocessing of list t and initializing a shifted tableau.

        TESTS::

            sage: s = ShiftedPrimedTableau([[1,"2'","3'",3],[2,"3'"]])
            sage: t = ShiftedPrimedTableaux([4,2])([[1,"2p","3p",3],[0, 2,"3p"]])
            sage: s==t
            True
            sage: t.parent()
            Shifted primed tableaux of shape [4, 2]
            sage: s.parent()
            Shifted Primed Tableaux of weight (1, 2, 3) and shape [4, 2]
            sage: r = ShiftedPrimedTableaux([4, 2])(s); r.parent()
            Shifted primed tableaux of shape [4, 2]
            sage: s is t # identical shifted tableaux are distinct objects
            False

        A shifted primed tableau is shallowly immutable, the rows are represented as tuples.

        ::

            sage: t = ShiftedPrimedTableau([[1,"2p","3p",3],[2,"3p"]])
            sage: t[0][1] = 3
            Traceback (most recent call last):
            ...
            TypeError: 'tuple' object does not support item assignment
        """
        t = preprocessing(T)

        if isinstance(t, ShiftedPrimedTableau):
            # Since we are (supposed to be) immutable, we can share the underlying data
            ClonableArray.__init__(self, parent, t)
            return

        ClonableArray.__init__(self, parent, t)
        # This dispatches the input verification to the :meth:`check`
        # method.

    def __eq__(self, other):
        """
        Check whether ``self`` is equal to ``other``.

        INPUT:

        ``other`` -- the element that ``self`` is compared to

        OUTPUT:

        Boolean.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,"2p"]])
            sage: t == ShiftedPrimedTableaux([2])([1,"2p"])
            True
        """
        if isinstance(other, ShiftedPrimedTableau):
            return list(self) == list(other)
        else:
            return list(self) == preprocessing(other)

    
    def to_matrix(self):
        """
        Return a 2-dimensional numpy.array representation of a shifted tableau

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p'],[3]])
            sage: mat = t.to_matrix()
            sage: mat
            array([[ 1. ,  1.5,  2. ,  2. ],
                   [ 0. ,  2. ,  2.5,  0. ],
                   [ 0. ,  0. ,  3. ,  0. ]])
            sage: t == ShiftedPrimedTableau(mat)
            True
        """
        array = []
        m = len(self[0])
        for i in range(len(self)):
            array.append([0]*i + list(self[i]) + [0]*(m-i-len(self[i])))
        array = np.array(array, dtype = 'float')
        return array
    
    def check(self):
        """
        Check that ``self`` is a valid primed tableaux.

        EXAMPLES::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: t = T([[1,'2p',2,2],[2,'3p']])
            sage: t.check()
        """
        if [len(_) for _ in self] not in StrictPartitions():
            raise ValueError('shape must be a strict partition')
        for i,row in enumerate(self):
            if i > 0:
                if not all(val > self[i-1][j+1] for j,val in enumerate(row) if int(val)==val):
                    raise ValueError('column is not strictly increasing in non-primes')
                if not all(val >= self[i-1][j+1] for j,val in enumerate(row) if int(val)!=val):
                    raise ValueError('column is not weakly increasing in primes')
            if not all(row[j] <= row[j+1] for j in range(len(row)-1) if int(row[j])==row[j]):
                raise ValueError('row is not weakly increasing in non-primes')
            if not all(row[j] < row[j+1] for j in range(len(row)-1) if int(row[j])!=row[j]):
                raise ValueError('row is not strictly increasing in primes')
        if not all(int(row[0])==row[0] for row in self):
            raise ValueError('diagonal elements must be non-primes')

    def _repr_(self):
        """
        Represent Shifted Primed Tableau as a list of rows, rows are represented as tuples of half-integers.

        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t
            [(1.0, 1.5, 2.0, 2.0), (2.0, 2.5)]
        """
        return repr([tuple(_) for _ in self])
    
    def pp(self):
        """
        Print out a nice version of ``self``.
        
        EXAMPLES::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.pp()
            1  2' 2  2  
               2  3'
            sage: t = ShiftedPrimedTableau([[10,'11p',11,11],[11,'12']])
            sage: t.pp()
            10  11' 11  11  
                11  12
        """
        if [item for sublist in self for item in sublist] == []:
            max_ind = 0
        else:
            max_ind = max([item for sublist in self for item in sublist])
        max_len = len(str(round(max_ind)))
        string_list = ''
        for i,row in enumerate(self):
            string_list += ' '
            string_list += ' ' * i * (max_len)
            for val in row:
                if int(val)==val:
                    string_list += str(int(val)) + ' ' * (max_len-len(str(int(val))))
                else:
                    string_list += str(int(val+.5)) + "'" + ' '*(max_len-1- len(str(int(val+.5))))
            string_list += '\n'
        string_list = string_list[:-2]
        print(string_list)
        return

    def _latex_(self):
        r"""
        Return LaTex code for ``self``.

        EXAMPLES::

            sage: T = ShiftedPrimedTableaux([4,2])
            sage: latex(T([[1,"2p",2,"3p"],[2,3]]))
            {\def\lr#1{\multicolumn{1}{|@{\hspace{.6ex}}c@{\hspace{.6ex}}|}{\raisebox{-.3ex}{$#1$}}}
            \raisebox{-.6ex}{$\begin{array}[b]{*{4}c}\cline{1-4}
            \lr{1}&\lr{2'}&\lr{2}&\lr{3'}\\\cline{1-4}
            &\lr{2}&\lr{3}\\\cline{2-3}
            \end{array}$}
            }
        """
        from sage.combinat.output import tex_from_array
        L = list()
        for i,row in enumerate(self):
            num_list = [None]*i
            for let in row:
                if int(let)==let:
                    num_list.append(int(let))
                else:
                    num_list.append(str(int(let+0.5))+"'")
            L.append(num_list)
                
        return tex_from_array(L)

    
    def shape(self):
        """
        Return the shape of the underlying partition of ``self`` in list format.
        
        EXAMPLES::
           
            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.shape()
            [4, 2]
        """
        return ([len(row) for row in self])
    
    def __call__(self,*cell):
        """
        Function call of ``self``.

        INPUT:

        - ``cell`` -- a pair of integers, tuple, or list specifying a cell in
          the tableau.

        OUTPUT:

        - the element in the corresponding cell. if the element is primed, 
          returnes half-integer value.

        EXAMPLES::

            sage: t  = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t(1,0)
            2.0
            sage: t((1,2))
            Traceback (most recent call last):
            ...
            IndexError: the cell (1,2) is not contained in the shifted tableau 
            [(1.0, 1.5, 2.0, 2.0), (2.0, 2.5)]
            sage: t((1,1))
            2.5
        """
        
        try:
            i,j = cell
        except ValueError:
            i,j = cell[0]

        try:
            return self[i][j]
        except IndexError:
            raise IndexError("the cell (%d,%d) is not contained in the shifted tableau \n%s"%(i, j, repr(self)))

        
    def weight(self):
        """
        Returnes the weight of ``self`` in tuple format.
        The weight of a shifted primed tableau is defined to be the vector with i-th component
        equal to the number of letters i and i' in the tableau.
        
        EXAMPLE::
        
           sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
           sage: t.weight()
           (1, 4, 1)
        """
        flat = [round(item) for sublist in self for item in sublist]
        if flat == []:
            max_ind = 0
        else:
            max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return tuple(weight)

    def reading_word_with_positions(self):
        """
        Returnes the reading word of ``self`` together with positions of the corresponding 
        letters in ``self``.
        The reading word of a shifted primed tableau is constructed as follows:
        (1) List all primed letters in the tableau, column by column, in decreasing order 
        within each column, moving from the rightmost column to the left, and with all the 
        primes removed (i.e. all letters are increased by half a unit).
        (2) Then list all unprimed elements, row by row, in increasing order within each row, 
        moving from the bottommost row to the top. 

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.reading_word_with_positions()
            [((1, 2), 3), ((0, 1), 2), ((1, 1), 2), ((0, 0), 1), ((0, 2), 2), ((0, 3), 2)]
        """
        mat = self.to_matrix()
        list_with_positions = []
        for (i,j), x in np.ndenumerate(mat[:,::-1].T):
            if int(x) != x:
                list_with_positions.append(((j,mat.shape[1]-i-1),int(x+0.5)))
        for (i,j), x in np.ndenumerate(mat[::-1,:]):
            if int(x) == x  and int(x) != 0:
                list_with_positions.append(((mat.shape[0]-i-1,j),int(x)))
        return list_with_positions

    def reading_word(self):
        """
        Returnes the reading word of ``self``.
        The reading word of a shifted primed tableau is constructed as follows:
        (1) List all primed letters in the tableau, column by column, in decreasing order 
        within each column, moving from the rightmost column to the left, and with all the 
        primes removed (i.e. all letters are increased by half a unit).
        (2) Then list all unprimed elements, row by row, in increasing order within each row, 
        moving from the bottommost row to the top. 

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([[1,'2p',2,2],[2,'3p']])
            sage: t.reading_word()
            [3, 2, 2, 1, 2, 2]
        """
        return [tup[1] for tup in self.reading_word_with_positions()]


    def f(self, ind):
        """
        A function to compute the action of the crystal operator $f_i$ on a Shifted Primed Tableau 
        using cases from the paper [GPS.17].

        INPUT:

        self -- shifted primed tableau
        ind -- index of the crystal operator $f_i$.
    
        OUTPUT:
    
        Primed tableau or 'None'.

        EXAMPLES::
    
            sage: t = ShiftedPrimedTableau([[1,1,1,1,2.5],[2,2,2,2.5],[3,3]])
            sage: t.pp()
            1  1  1  1  3' 
               2  2  2  3' 
                  3  3  
            sage: s = t.f(2)
            sage: print(s)
            None
            sage: t = ShiftedPrimedTableau([[1,1,1,1.5,2.5],[0,2,2,3,3],[0,0,3,4]])
            sage: t.pp()
            1  1  1  2' 3' 
               2  2  3  3  
                  3  4  
            sage: s = t.f(2)
            sage: s.pp()
            1  1  1  2' 3'    
               2  3' 3  3     
                  3  4     

        """
        if self is None:
            return None
        
        T = self.to_matrix()
       
        read_word = self.reading_word_with_positions()
        read_word = [num for num in read_word if num[1]==ind or num[1]==ind+1]

        element_to_change = None
        count=0
        
        for element in read_word:
            if element[1] == ind+1:
                count += 1
            elif count == 0:
                element_to_change = element
            else:
                count -= 1

        if element_to_change == None:
            return None

        (r,c), elt =  element_to_change

        if  T[r,c] == ind- .5:
            T = T.T + .5
            r,c = c,r
        h,l = T.shape

        if (c+1==l or T[r,c+1]>=ind+1 or T[r,c+1]<1):

            (tp_r,tp_c)=(r,c)
            while True:

                if tp_r+1==h or T[tp_r+1,tp_c]>ind+1 or T[tp_r+1,tp_c]<1:
                    break
                if tp_r<=tp_c and T[tp_r+1,tp_r+1]==ind+1:
                    tp_r += 1
                    tp_c = tp_r
                    break
                if  (ind+.5 not in T[tp_r+1]):
                    break
                tp_r += 1
                tp_c = np.where(T[tp_r]==ind+.5)[0]


            if tp_r == r:
                T[r,c] += 1

            elif tp_r == tp_c:
                T[r,c] += .5

            else:
                T[r,c] += .5
                T[tp_r,tp_c] += .5

        elif T[r,c+1]==ind+.5:
            T[r,c+1] += .5
            T[r,c] += .5

        if r>c:
            T = T.T - .5

        return(ShiftedPrimedTableau(T))

    def e(self,ind):
        """
        A function to compute an action of the crystal operator $e_i$ on a Primed Tableau 
        using cases from the paper [GPS.17].

        INPUT:

        self -- shifted primed tableau
        ind -- index of the crystal operator $e_i$.
    
        OUTPUT:
    
        Primed tableau or 'None'.
    
        EXAMPLES::
    
            sage: t = ShiftedPrimedTableau([[1,1,1,1.5,2.5],[0,2,"3p",3,3],[0,0,3,4]])
            sage: t.pp()
            1  1  1  2' 3' 
               2  3' 3  3  
                  3  4  
            sage: s = t.e(2)
            sage: s.pp()
            1  1  1  2' 3' 
               2  2  3  3  
                  3  4 
            sage: t == s.f(2)
            True
    
        """
        if self is None:
            return None
        T = self.to_matrix()
       
        read_word = self.reading_word_with_positions()
        read_word = [num for num in read_word if num[1]==ind or num[1]==ind+1]

        element_to_change = None
        count=0
        
        for element in read_word[::-1]:
            if element[1] == ind:
                count += 1
            elif count == 0:
                element_to_change = element
            else:
                count -= 1

        if element_to_change == None:
            return None
        (r,c), elt =  element_to_change

        if  T[r,c] == ind + .5:
            T = T.T + .5
            r,c = c,r
        h,l = T.shape

        if (c==0 or T[r,c-1]<=ind or T[r,c-1]<1):

            (tp_r,tp_c)=(r,c)
            while True:

                if tp_r==0 or T[tp_r-1,tp_c]<ind or T[tp_r-1,tp_c]<1:
                    break
                if  (ind+.5 not in T[tp_r-1]):
                    break
                tp_r -= 1
                tp_c = np.where(T[tp_r]==ind+.5)[0]


            if tp_r == r:
                T[r,c] -= 1

            elif tp_r == tp_c:
                T[r,c] -= .5

            else:
                T[r,c] -= .5
                T[tp_r,tp_c] -= .5

        elif T[r,c-1]==ind+.5:
            T[r,c-1] -= .5
            T[r,c] -= .5

        if r>c:
            T = T.T - .5

        return(ShiftedPrimedTableau(T))

    def epsilon(self,i):
        r"""
        Compute value of the crystal function $\varepsilon_i$ applied to a shifted primed tableau ``self``.
        The function $\varepsilon_i$ is defined to be the  maximal number of operators $e_i$ one can apply to
        ``self`` before vanishing into ``None``.

        INPUT:
        
        self -- shifted primed tableau
        ind -- index of the function $\varepsilon_i$ associated with a crystal operator $e_i$.
    
        OUTPUT:
    
        Value of the function $\varepsilon_i$.
    
        EXAMPLES::
    
            sage: t = ShiftedPrimedTableau([(1.0, 1.0, 1.0, 1.5, 2.5), (2.0, 2.5, 3.0, 3.0), (3.0, 4.0)])
            sage: t.epsilon(2)
            3
            sage: s = t.e(2).e(2).e(2)
            sage: s.epsilon(2)
            0
        """
        b = self
        count = -1
        while b is not None:
            b = b.e(i)
            count +=1
        return count

    def phi(self,i):
        r"""
        Compute value of the crystal function $\phi_i$ applied to a shifted primed tableau ``self``.
        The function $\phi_i$ is defined to be the  maximal number of operators $f_i$ one can apply to
        ``self`` before vanishing into ``None``.

        INPUT:
        
        self -- shifted primed tableau
        ind -- index of the function $\varepsilon_i$ associated with a crystal operator $f_i$.
    
        OUTPUT:
    
        Value of the function $\varepsilon_i$.
    
        EXAMPLES::
    
            sage: t = ShiftedPrimedTableau([(1.0, 1.0, 1.0, 1.5, 2.0), (2.0, 2.0, 2.0, 2.5), (3.0, 4.0)])
            sage: t.phi(2)
            3
            sage: s = t.f(2).f(2).f(2)
            sage: s.phi(2)
            0
        """
        b = self
        count = -1
        while b is not None:
            b = b.f(i)
            count +=1
        return count

    def is_highest_weight(self):
        """
        Check wether the shifted primed tableau ``self`` is a highest weight element of the crystal.
        Highest weight element defined to be vanishing under any crystal operator $e_i$.

        EXAMPLE::

            sage: t = ShiftedPrimedTableau([(1.0, 1.0, 1.0, 1.0, 1.0), (2.0, 2.0, 2.0, 2.5), (3.0, 3.0)])
            sage: print(t.e(1))
            None
            sage: print(t.e(2))
            None
            sage: t.is_highest_weight()
            True

        """
        read_w = self.reading_word()
        count = {}
        for l in read_w[::-1]:
            try:
                count[l] +=1
            except KeyError:
                count[l] = 1
            if l>1:
                try:
                    if count[l]>count[l-1]:
                        return False
                except KeyError:
                    return False
        return True
    
        
                    

class ShiftedPrimedTableaux(UniqueRepresentation, Parent):
    """
    A factory for the various classes of shifted standard tableaux.

    INPUT:

    - a weight and/or a partition

    OUTPUT:

    - with no argument, the class of all primed tableaux

    - with a list or partition argument, the class of all primed tableaux of that
        shape (infinite set if we don't specify the weight or maximum element)
        - with an additional integer argument ``max_element``, the class of all primed
            tableaux of a given shape and a given maximum element
        
    - with a tuple argument, the class of all primed tableaux of that
        weight (finite set)

    A primed tableau is a shifted tableau on the alphabet 
        X' = {1' < 1 < 2' < 2 <...< n' < n} such that 
        1). the entries are weakly increasing along rows and columns
        2). a row can't have two repeated primed elements, and a column
            can't have two repeated non-primed elements
        3). there are only non-primed elements along the main diagonal
        
    The weight of a tableau is defined to be the vector with i-th component
        equal to the number of letters i and i' in the tableau.
    The sum of the entries in the weight vector must be equal to the number 
        of boxes in the partition.

    All classes of Shifted Primed Tableaux are not iterateble.

    EXAMPLES::

        sage: SPT = ShiftedPrimedTableaux(weight=(1,2,2), shape=[3,2]); SPT
        Shifted Primed Tableaux of weight (1, 2, 2) and shape [3, 2]
        sage: SPT.list()
        [[(1.0, 2.0, 2.0), (3.0, 3.0)],
         [(1.0, 1.5, 2.5), (2.0, 3.0)],
         [(1.0, 1.5, 2.5), (2.0, 2.5)],
         [(1.0, 1.5, 2.0), (3.0, 3.0)]]

        sage: SPT = ShiftedPrimedTableaux((1,2)); SPT
        Shifted Primed Tableaux of weight (1, 2)
        sage: list(SPT)
        [[(1.0, 2.0, 2.0)], [(1.0, 1.5, 2.0)], [(1.0, 1.5), (2.0,)]]

        sage: SPT = ShiftedPrimedTableaux([3,2], max_element = 2); SPT
        Shifted Primed Tableaux of shape [3, 2] and maximum element 2
        sage: list(SPT)
        [[(1.0, 1.0, 1.0), (2.0, 2.0)], [(1.0, 1.0, 1.5), (2.0, 2.0)]]
        sage: SPT.list_highest_weight()
        [[(1.0, 1.0, 1.0), (2.0, 2.0)]]

    .. SEEALSO::

        - :class:`ShiftedPrimedTableau`
    """

    @staticmethod
    def __classcall_private__(cls, *args, **kwargs):
        r"""
        This is a factory class which returns the appropriate parent based on
        arguments.  See the documentation for :class:`ShiftedPrimedTableaux` for
        more information.


        TESTS::

            sage: ShiftedPrimedTableaux([])
            Shifted Primed Tableaux of shape []

            sage: ShiftedPrimedTableaux(3)
            Traceback (most recent call last):
            ...
            ValueError: weight argument must be a tuple and shape argument must be a strictly increasing partition

            sage: ShiftedPrimedTableaux(weight=(2,2,2), shape=[3,2])
            Traceback (most recent call last):
            ...
            ValueError: the sum of weights should match the sum of parts of the shape

            sage: ShiftedPrimedTableaux([[1]])
            Traceback (most recent call last):
            ...
            ValueError: [[1]] is not a partition

            sage: ShiftedPrimedTableaux(weight=(2,2,2), max_element=2)
            Traceback (most recent call last):
            ...
            ValueError: maximum element can not be smaller then the length of the weight vector

        """
        weight = None
        shape = None
        max_element = None

        if 'max_elt' in kwargs:
            max_element = int(kwargs['max_elt'])

        if 'max_element' in kwargs:
            max_element = int(kwargs['max_element'])

        if 'max' in kwargs:
            max_element = int(kwargs['max'])
        
        if 'size' in kwargs and isinstance(kwargs['size'],(list,tuple,Partition)):
            shape = Partition(kwargs['size'])
        
        if 'shape' in kwargs:
            shape = Partition(kwargs['shape'])

        if 'weight' in kwargs:
            weight = tuple(kwargs['weight'])

        if args:
            if isinstance(args[0],tuple) and weight==None:
                weight = args[0]
                if len(args)>1:
                    if (isinstance(args[1],(list,Partition)) or args[1]==None) and shape==None:
                        shape = args[1]
                    else:
                        raise ValueError('weight argument must be a tuple and shape argument must be a strictly increasing partition')
            
            elif isinstance(args[0],(list,Partition)) and shape==None:
                shape = args[0]
                if len(args)>1:
                    if (isinstance(args[1],tuple) or args[1]==None) and weight==None:
                        weight = args[1]
                    else:
                        raise ValueError('weight argument must be a tuple and shape argument must be a strictly increasing partition')
            else:
                raise ValueError('weight argument must be a tuple and shape argument must be a strictly increasing partition')

        if shape is not None:
            try:
                shape = Partition(shape)
            except ValueError:
                raise ValueError('{} is not a partition'.format(shape))
            
        if weight is not None:
            while weight[-1]==0:
                weight = weight[:-1]

        if max_element is not None and weight is not None:
            if len(weight)!=max_element:
                raise ValueError("maximum element can not be smaller then the length of the weight vector")

        if max_element is not None and shape is not None:
            if max_element<len(shape):
                raise ValueError("maximum element can not be smaller then the number of rows")

    
        if shape is None and weight is None:
            if max_element is not None:
                raise ValueError("maximum element can not be the only argument, specify shape or weight")
            return ShiftedPrimedTableaux_all()        

        elif weight is None and shape in StrictPartitions():
            return ShiftedPrimedTableaux_shape(Partition(shape), max_element)
        
        elif shape is None:
            return ShiftedPrimedTableaux_weight(weight)

        if shape not in StrictPartitions():
            raise ValueError("{} is not a strict partition".format(shape))
            
        if sum(shape) != sum(weight):
            raise ValueError("the sum of weights should match the sum of parts of the shape")

        return ShiftedPrimedTableaux_weight_shape(weight,shape)

    def __contains__(self, T):
        """
        EXAMPLES::

            sage: [(1,'2p',2,2),(2,'3p')] in ShiftedPrimedTableaux()
            True
            sage: [(1,1),(2,2)]           in ShiftedPrimedTableaux()
            False
            sage: [(1,1),('2p',)]         in ShiftedPrimedTableaux()
            False
            sage: [(1,1,'3p'),(2,'3p')]   in ShiftedPrimedTableaux()
            True
            sage: [(1,1,2),(2,2)]         in ShiftedPrimedTableaux()
            False
            sage: [1,1,1] in ShiftedPrimedTableaux()
            True

        """
        if isinstance(T, ShiftedPrimedTableau) or T == []:
            return True

        t=preprocessing(T)

        if [len(_) for _ in t] not in StrictPartitions():
            return False # must have strict partition shape

        for i,row in enumerate(t):
            if i > 0:
                if not all(val > t[i-1][j+1] for j,val in enumerate(row) if int(val) == val):
                    return False
                if not all(val >= t[i-1][j+1] for j,val in enumerate(row) if int(val) != val):
                    return False
            if not all(row[j] <= row[j+1] for j in range(len(row)-1) if int(row[j]) == row[j]):
                return False
            if not all(row[j] < row[j+1] for j in range(len(row)-1) if int(row[j]) != row[j]):
                return False
        if not all(int(row[0]) == row[0] for row in t):
            return False

        return True

    _is_a = __contains__



class ShiftedPrimedTableaux_all(ShiftedPrimedTableaux):

    Element = ShiftedPrimedTableau

    def __init__(self):
        Parent.__init__(self, category=FiniteEnumeratedSets())

    def _repr_(self):
        return "Shifted Primed Tableaux"

    def _element_constructor_(self, T):
        t = preprocessing(T)
        if not t in self:
            raise ValueError("{} is not an element of {}".format(t, self))
        
        return self.element_class(self, t)

    def __contains__(self, x):
        
        if isinstance(x, ShiftedPrimedTableau):
            return True
        return (ShiftedPrimedTableaux.__contains__(self, x))
    

    
class ShiftedPrimedTableaux_shape(ShiftedPrimedTableaux):
    """
    Shifted Primed tableaux of a fixed shape.
    """
    Element = ShiftedPrimedTableau

    def __init__(self, shape, max_elt):
        r"""
        Initializes the class of semistandard tableaux of a fixed shape.

        """
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._max_elt = max_elt
        self._shape = shape

    def _repr_(self):
        """
        Return a string representation of ``self``.

        TESTS::

            sage: ShiftedPrimedTableaux([3,2,1])    # indirect doctest
            Shifted Primed tableaux of shape [3, 2, 1]
        """
        if self._max_elt is None:
            return "Shifted Primed Tableaux of shape {}".format(self._shape)
        if self._max_elt is not None:
            return "Shifted Primed Tableaux of shape {} and maximum element {}".format(self._shape, self._max_elt)

    def __contains__(self, T):
        """
        TESTS::

           sage: [[1,'2p',2,2],[2,'3p']] in ShiftedPrimedTableaux([4,2])
           True
        """
        
        if isinstance(T, ShiftedPrimedTableau):
            return [len(row) for row in T] == self._shape

        t = preprocessing(T)

        if self._max_elt is not None:
            return (ShiftedPrimedTableaux.__contains__(self, t)
                and [len(row) for row in t] == self._shape
                and max(flatten(t)) <= self._max_elt)
        
        return (ShiftedPrimedTableaux.__contains__(self, t)
                and [len(row) for row in t] == self._shape)

    def _element_constructor_(self, T):
        r"""
        Constructs an object from ``t`` as an element of ``self``, if
        possible.

        INPUT:

        - ``t`` -- data which can be interpreted as a primed tableau

        OUTPUT:

        - the corresponding primed tableau object

        TESTS::

            sage: ShiftedPrimedTableaux([3])([[1,2,3]]).parent() is ShiftedPrimedTableaux([3])
            True
            sage: ShiftedPrimedTableaux([3])([[1,2]])
            Traceback (most recent call last):
            ...
            ValueError: [[1, 2]] is not an element of Shifted Primed tableaux of shape [3]
        """
        t = preprocessing(T)
        if not t in self:
            raise ValueError("{} is not an element of {}".format(t, self))
        
        return self.element_class(self, t)

    def shape(self):
        """
        Return the shape of the shifted tableaux ``self``.

        EXAMPLES::

            sage: ShiftedPrimedTableaux([6,4,3,1]).shape()
            [6, 4, 3, 1]
        """
        return self._shape

    def __iter__(self):
        if self._max_elt is None:
            raise ValueError("set is infinite")
        for weight in OrderedPartitions(sum(self._shape)+self._max_elt,k=self._max_elt):
            weight_n = tuple([w-1 for w in weight])
            for tab in ShiftedPrimedTableaux(shape = self._shape, weight = weight_n):
                yield (tab)

    def list_decreasing_weight(self):
        list_dw = []
        if self._max_elt is None:
            max_element = sum(self._shape)
        else:
            max_element = self._max_elt
        for weight in Partition(self._shape).dominated_partitions(rows=max_element):
            list_dw.extend(list(ShiftedPrimedTableaux(weight=tuple(weight),
                                                     shape = self._shape)))
        return list_dw

    def list_highest_weight(self):
        return [tab for tab in self.list_decreasing_weight() if tab.is_highest_weight()]



class ShiftedPrimedTableaux_weight(ShiftedPrimedTableaux):
    """
    Shifted Primed Tableaux of fixed weight
    """
    Element = ShiftedPrimedTableau

    def __init__(self, weight):
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._weight = weight

    def _repr_(self):
        return "Shifted Primed Tableaux of weight {}".format(self._weight)

    def __contains__(self, T):
        if isinstance(x, ShiftedPrimedTableau):
            return T.weight() == self._weight

        t = preprocessing(T)

        flat = [round(item) for sublist in t for item in sublist]
        if flat == []:
            max_ind = 0
        else:
            max_ind = max(flat)
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return (ShiftedPrimedTableaux.__contains__(self, t)
                and weight == self._weight)
    
    def _element_constructor_(self, T):
        t = preprocessing(T)
        if not t in self:
            raise ValueError("{} is not an element of {}".format(t, self))
        return self.element_class(self, t)

    def __iter__(self):
        #for size in Partition(sorted(list(self._weight), key=int, reverse=True)).dominated_partitions():
        return (tab for shape_ in Partitions(sum(self._weight))
                if all(shape_[i] > shape_[i+1] for i in range(len(shape_)-1))
                for tab in ShiftedPrimedTableaux(shape = shape_, weight = self._weight))
    


class ShiftedPrimedTableaux_weight_shape(ShiftedPrimedTableaux):
    Element = ShiftedPrimedTableau

    def __init__(self, weight, shape):
        Parent.__init__(self, category=FiniteEnumeratedSets())
        self._weight = weight
        self._shape = shape

    def _repr_(self):
        return ("Shifted Primed Tableaux of weight {} and shape {}" .format(self._weight,self._shape))

    def __contains__(self, T):
        if isinstance(T, ShiftedPrimedTableau):
            return (T.weight() == self._weight and T.shape() == self._shape)

        t = preprocessing(T)

        flat = [round(item) for sublist in t for item in sublist]
        if flat == []:
            max_ind = 0
        else:
            max_ind = int(max(flat))
        weight = tuple([flat.count(i+1) for i in range(max_ind)])
        return(ShiftedPrimedTableaux.__contains__(self, t) and weight == self._weight
               and [len(row) for row in t] == self._shape)

    def _element_constructor_(self, T):
        t = preprocessing(T)
        if not t in ShiftedPrimedTableaux():
            raise ValueError("{!r} is not an element of Shifted Primed Tableaux" .format(t))
        if not t in self:
            raise ValueError("{} is not an element of {}".format(t, self))
        
        return self.element_class(self, t)

    def __iter__(self):
        if not self._shape.dominates(Partition(sorted(list(self._weight), key=int, reverse=True))):
            return
            yield
        #TODO: More efficient algorithm with generators
        full_shape = self._shape
        sub_tab = []
        tab_list_new = [[]]
        for i,w in enumerate(self._weight):
            tab_list_old = tab_list_new
            tab_list_new = []
            for sub_tab in tab_list_old:        
                sub_shape = [len(sub_tab[r]) for r in range(len(sub_tab))] 
                for strip in add_strip(sub_shape, full_shape, w):
                    l = int(len(strip)/2)
                    if len(sub_shape)<len(full_shape):
                        new_tab = [sub_tab[r] + [float(i+.5)]*int(strip[r])
                                   + [float(i+1)]*int(strip[-r-1]) for r in range(l-1)]
                        if strip[l]!=0:
                            new_tab.append([float(i+1)]*int(strip[l]))
                    else:
                        new_tab = [sub_tab[r] + [float(i+.5)]*int(strip[r])
                                   +  [float(i+1)]*int(strip[-r-1]) for r in range(l)]
                    tab_list_new.append(new_tab)
        for tab in tab_list_new:
            yield(ShiftedPrimedTableau(tab))



###########
# Crystal #
###########

class SPTCrystal(UniqueRepresentation, Parent):
    @staticmethod
    def __classcall_private__(cls,*args,**kwargs):
        shape = None
        n= None
        if 'shape' in kwargs:
            shape = tuple(kwargs['shape'])
        if 'n' in kwargs:
            n = int(kwargs['n'])
        if args:
            if isinstance(args[0],(list,tuple,Partition)):
                shape = tuple(args[0])
                if len(args)>1 and isinstance(args[1],int):
                    n = args[1]
            else:
                if isinstance(args[0],int):
                    n = args[0]
                    if len(args)>1 and isinstance(args[1],(list,tuple,Partition)):
                        shape = tuple(args[1])
        if shape is None:
            raise ValueError('input size of tableaux as a list or a tuple')

        if n is None:
            n = sum(shape)-1
            
        if n+1 < len(shape):
            raise ValueError('crystal index should be greater or equal to the number of rows minus one')
        return ShiftedPrimedTableauxCrystal(shape =shape, n=n)
   

class ShiftedPrimedTableauxCrystal(SPTCrystal):  
     def __init__(self, shape=[4,2],n=2):
        Parent.__init__(self, category = ClassicalCrystals())
        self.n = n
        self._shape = shape
        self._cartan_type = CartanType(['A',n])
        self.module_generators = ShiftedPrimedTableaux(shape=shape,max_element=n+1).list_decreasing_weight()

     def _repr_(self):
        return ("Crystal of Shifted Primed Tableaux of type A_%s of shape "%(self.n) + str(self._shape))

####################
# Helper functions #
####################

def add_strip(sub_tab, full_tab, length):
    if sum(sub_tab)+length > sum(full_tab):
        raise ValueError("strip does not fit")

    if len(sub_tab)==0:
        cliff_list = []
    else:
        cliff_list = [int(sub_tab[0]!=full_tab[0])]

    for row in range(1,len(sub_tab)):
        if sub_tab[row] == full_tab[row]:
            cliff_list.append(0)
        elif sub_tab[row-1]-1 == sub_tab[row]:
            cliff_list[-1] += 1
        else:
            cliff_list.append(1)
    
    if len(sub_tab)<len(full_tab):
        cliff_list.append(0)
        
    for primes_num in range(min(sum(cliff_list),length)+1):
        for primed_list in IntegerVectors(n=primes_num, k=len(cliff_list), outer=cliff_list):
            row=0
            primed_strip = list()
            for i,cliff in enumerate(cliff_list):
                if cliff == 0:
                    row += 1
                    primed_strip.append(0)
                    pass
                primed_strip.extend([int(primed_list[i] > j) for j in range(cliff)])
                row += cliff
            plat_list = list()

            if len(sub_tab)<len(full_tab) and len(sub_tab)!=0:
                plat_list.append(min(sub_tab[-1] + primed_strip[-2] - 1, full_tab[len(sub_tab)]))
            for row in range(1,len(sub_tab))[::-1]:
                plat_list.append(min(sub_tab[row-1]+primed_strip[row-1]-1, full_tab[row])
                                 -sub_tab[row]-primed_strip[row])
            if len(sub_tab)>0:
                plat_list.append(full_tab[0] - sub_tab[0]- primed_strip[0])
            else:
                plat_list.append(full_tab[0])
           
            if sum(plat_list)<length - primes_num:
                pass
            for non_primed_strip in IntegerVectors(n=length-primes_num, k=len(plat_list), outer = plat_list):
                yield (list(primed_strip) + list(non_primed_strip))


# Helper function, preprocessing of the input array to fit Shifted Primed Tableau format. 
# Right now it requires 4 preprocessing calls to initialize ShiftedPrimedTableau element.
# Donno how to fix.

def preprocessing(T):
    
    if isinstance(T,ShiftedPrimedTableau):
        return T
    elif isinstance(T, np.ndarray):
        t = T.tolist()
    elif not isinstance(T,(list,tuple)):
        raise ValueError ("input must be list, tuple or tableau")
    elif all(isinstance(row, (list,tuple)) for row in T):
        t = []
        # Preprocessing list T for primes and other symbols
        for row in T:
            row_new = []
            for element in row:
                try:
                    if int(element*2) == element*2:
                        row_new.append(float(element))
                        continue
                    else:
                        raise ValueError("all numbers must be half-integers")
                except (TypeError,ValueError):
                    pass
                if isinstance(element,str):
                    if element[-1] == "'" and element[:-1].isdigit() == True:
                        row_new.append(float(float(element[:-1]) - .5))
                        continue
                    if element[-1] == "p" and element[:-1].isdigit() == True:
                        row_new.append(float(float(element[:-1]) - .5))
                        continue
                    try:
                        row_new.append(float(element))
                    except ValueError:
                        raise ValueError("primed elements should be half-integers or have p or ' at the end")
            t.append(row_new)
    else:
        row_new=[]
        for element in T:
            try:
                if int(element*2) == element*2:
                    row_new.append(float(element))
                    continue
                else:
                    raise ValueError("all numbers must be half-integers")
            except (ValueError,TypeError):
                pass
            if isinstance(element,str):
                if element[-1] == "'" and element[:-1].isdigit() == True:
                    row_new.append(float(float(element[:-1]) - .5))
                    continue
                if element[-1] == "p" and element[:-1].isdigit() == True:
                    row_new.append(float(float(element[:-1]) - .5))
                    continue
                try:
                    row_new.append(float(element))
                except ValueError:
                    raise ValueError("primed elements should be half-integers or have p or ' at the end")
        t=[row_new]


    #Accounting for zeros at the beginning and at the end of a row
    i=0
    while i < len(t):
        row = t[i]
        try:
            while row[0]==0:
                row.pop(0)
            while row[-1]==0:
                row.pop(-1)
        except IndexError:
            t.remove(t[i])
            continue
        t[i] = row
        i += 1

    # Normalize t to be a list of tuples.
    t = [tuple(_) for _ in t]
    return t
