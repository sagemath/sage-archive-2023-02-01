from sage.structure.element_wrapper import ElementWrapper
from sage.structure.sage_object import SageObject
from sage.structure.parent import Parent
from sage.structure.parent_gens import normalize_names
from sage.rings.integer_ring import ZZ
 
class ClusterAlgebraElement(ElementWrapper):

    # This function needs to be removed once AdditiveMagmas.Subobjects
    # implements _add_
    def _add_(self, other):
        return self.parent().retract(self.lift() + other.lift())

class ClusterAlgebraSeed(SageObject):
    
    def __init__(self, B, C, G, parent):
        self._B = copy(B)
        self._C = copy(C)
        self._G = copy(G)
        self._parent = parent

    def __copy__(self):
        other = type(self).__new__(type(self))
        other._B = copy(self._B)
        other._C = copy(self._C)
        other._G = copy(self._G)
        other._parent = self._parent
        return other

    def parent(self):
        return self._parent

    def F_polynomial(self, j):
        gvect = tuple(self._G.column(j))
        return self.parent().F_polynomial(gvect)
    
    #maybe this alias sould be removed
    F = F_polynomial
    
    def cluster_variable(self, j):
        gvect = tuple(self._G.column(j))
        return self.parent().cluster_variable(gvect)

    def g_vector(self, j):
        return tuple(self._G.column(j))

    def g_matrix(self):
        return copy(self._G)

    def c_vector(self, j):
        return tuple(self._C.column(j))

    def c_matrix(self):
        return copy(self._C)
        
    def mutate(self, k, inplace=True, mutating_F=True):
        if inplace:
            seed = self
        else:
            seed = copy(self)
        
        n = self.parent().rk()        
    
        if k not in xrange(n):
            raise ValueError('Cannot mutate in direction' + str(k) + '.')
        
        # find sign of k-th c-vector
        if any(x > 0 for x in seed._C.column(k)):
            eps = +1
        else:
            eps = -1
        
        # store the g-vector to be mutated in case we are mutating also F-polynomials
        old_gvect = tuple(seed._G.column(k))
 
        # mutate G-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[j,k] += max(0, -eps*seed._B[j,k])
        J[k,k] = -1
        seed._G = seed._G*J
 
        # F-polynomials
        if mutating_F:
            gvect = tuple(seed._G.column(k))
            if not self.parent()._F_dict.has_key(gvect):
                self.parent()._F_dict.setdefault(gvect, self._mutated_F(k, old_gvect))
        
        # C-matrix
        J = identity_matrix(n)
        for j in xrange(n):
            J[k,j] += max(0, eps*seed._B[k,j])
        J[k,k] = -1
        seed._C = seed._C*J
        
        # B-matrix
        seed._B.mutate(k)
            
        if not inplace:
            return seed

    def _mutated_F(self, k, old_gvect):
        alg = self.parent()
        pos = alg._U(1)
        neg = alg._U(1)
        for j in xrange(alg.rk()):
            if self._C[j,k] > 0:
                pos *= alg._U.gen(j)**self._C[j,k]
            else:
                neg *= alg._U.gen(j)**(-self._C[j,k])
            if self._B[j,k] > 0:
                pos *= self.F_polynomial(j)**self._B[j,k]
            elif self._B[j,k] <0:
                neg *= self.F_polynomial(j)**(-self._B[j,k])
        return (pos+neg)//alg.F_polynomial(old_gvect)

    def mutation_sequence(self, sequence, inplace=True, mutating_F=True):
        seq = iter(sequence)
        
        if inplace:  
            seed = self
        else:                                                                                                                                                                                                                                                             
            seed = self.mutate(seq.next(), inplace=False, mutating_F=mutating_F)
       
        for k in seq:
            seed.mutate(k, inplace=True, mutating_F=mutating_F)

        if not inplace:
            return seed
        

class ClusterAlgebra(Parent):
 
    Element = ClusterAlgebraElement
 
    def __init__(self, B0, scalars=ZZ):
        self._B0 = copy(B0)
        self._n = self._B0.ncols()
        self._m = self._B0.nrows()-self._n
        
        # maybe here scalars can be replaced with just ZZ
        self._U = PolynomialRing(scalars,['u%s'%i for i in xrange(self._n)])
        self._F_dict = dict([ (tuple(v) ,self._U(1)) for v in identity_matrix(self._n).columns() ])
        
        self.Seed = ClusterAlgebraSeed(self._B0[:self._n,:self._n], identity_matrix(self._n), identity_matrix(self._n), self)

        names = normalize_names(self._n, 'x')
        base = LaurentPolynomialRing(scalars,names)
        names += normalize_names(self._m, 'y')
        self._ambient = LaurentPolynomialRing(scalars,names)
        Parent.__init__(self, base=base, category=CommutativeAlgebras(scalars).Subobjects())
 
    def _repr_(self):
        return "Cluster Algebra of rank %s"%self.rk()

    def rk(self):
        return self._n

    def F_polynomial(self, gvect):
        gvect = tuple(gvect)
        try:
            return self._F_dict[gvect]
        except:
            # TODO: improve this error message
            raise ValueError("This F-polynomial has not been computed yet. Did you explore the tree with compute_F=False ?")
    
    # maybe this alias could be removed
    F = F_polynomial

    def cluster_variable(self, gvect):
        raise NotImplementedError("I will do this, I swear it on the Beatles.")

    def ambient(self):
        return self._ambient
 
    def lift(self, x):
        return x.value
 
    def retract(self, x):
        return self(x)
 

#from sage.structure.sage_object import SageObject
#from sage.matrix.constructor import identity_matrix
#
#class ClusterAlgebra(SageObject):
#
#    def __init__(self, data):
#        if isinstance(data, Matrix):
#            self._B0 = copy(data) 
#            self._m = self._B0.nrows() 
#            self._n = self._B0.ncols()
#            B = copy(self._B0[:self._n,:self._n])
#            if not B.is_skew_symmetrizable(positive=True):
#                raise ValueError('data must have skew-symmetrizable principal part.')
#            self.current_seed = Seed(B)
#            self._U = PolynomialRing(QQ,['u%s'%i for i in xrange(self._n)])
#            self._F = dict([ (self.current_seed.g_vector(i),self._U(1)) for i in xrange(self._n) ])
#            self._R = LaurentPolynomialRing(QQ,['x%s'%i for i in xrange(self._m)])
#
#            def col_to_y(j):
#                return  prod([self._R.gen(i)**self._B0[i,j] for i in xrange(self._n,self._m)])
#            self._y = dict([ (self._U.gen(j),col_to_y(j)) for j in xrange(self._n)]) 
#
#            def col_to_yhat(j):
#                return prod([self._R.gen(i)**self._B0[i,j] for i  in xrange(self._m)])
#            self._yhat = dict([ (self._U.gen(j),col_to_yhat(j)) for j in xrange(self._n)])
#            
#            # should keep track of the seed too ?
#
#        # at the moment we only deal with b_matrices
#        else:
#            raise NotImplementedError('At the moment we only deal with matrices.')
#
#    def __copy__(self):
#        other = ClusterAlgebra(self._B0)
#        other.current_seed = copy(self.current_seed)
#        other._F = copy(self._F)
#
#    def mutate(self, sequence, inplace=True, mutate_F=True):
#        if not isinstance(mutate_F, bool):
#            raise ValueError('mutate_F must be boolean.')
#
#        if not isinstance(inplace, bool):
#            raise ValueError('inplace must be boolean.')
#        if inplace:
#            seed = self.current_seed
#        else:
#            seed = copy(self.current_seed)
#            
#        if sequence in xrange(seed._n):
#            seq = iter([sequence])
#        else:
#            seq = iter(sequence)
#            
#        for k in seq:
#            if k not in xrange(seed._n):
#                raise ValueError('Cannot mutate in direction' + str(k) + '.')
#
#            # G-matrix
#            J = identity_matrix(seed._n)
#            if any(x > 0 for x in seed._C.column(k)):
#                eps = +1
#            else:
#                eps = -1
#            for j in xrange(seed._n):
#                J[j,k] += max(0, -eps*seed._B[j,k])
#            J[k,k] = -1
#            seed._G = seed._G*J
#
#            # F-polynomials
#            if mutating_F:
#                gvect = tuple(seed._G.column(k))
#                if not self._F_dict.has_key(gvect):
#                    self._F_dict.setdefault(gvect,self._mutated_F(k))
#                seed._F[k] = self._F_dict[gvect]
#            
#            # C-matrix
#            J = identity_matrix(seed._n)
#            if any(x > 0 for x in seed._C.column(k)):
#                eps = +1
#            else:
#                eps = -1
#            for j in xrange(seed._n):
#                J[k,j] += max(0, eps*seed._B[k,j])
#            J[k,k] = -1
#            seed._C = seed._C*J
#            
#            # B-matrix
#            seed._B.mutate(k)
#            
#        if not inplace:
#            return seed
#
#class Seed(SageObject):
#
#    def __init__(self, B):
#        n = B.ncols()
#        self._B = copy(B)
#        self._C = identity_matrix(n)
#        self._G = identity_matrix(n)
#        self._path = []
#
#    def __copy__(self):
#        other = Seed(self._B)
#        other._C = copy(self._C)
#        other._G = copy(self._G)
#        other._path = copy(self._path)
#        return other
#
#    def g_vectors(self):
#        return tuple(self._G.columns())
#
#    def g_vector(self, i):
#        return self._G.columns()[i]
#
#
#
#
######################## OLD ###################
#
#import time
#from sage.matrix.matrix import Matrix
#
#
#
#    def _mutated_F(self, k):
#        pos = self._U(1)
#        neg = self._U(1)
#        for j in xrange(self._n):
#            if self._C[j,k] > 0:
#                pos *= self._U.gen(j)**self._C[j,k]
#            else:
#                neg *= self._U.gen(j)**(-self._C[j,k])
#            if self._B[j,k] > 0:
#                pos *= self._F[j]**self._B[j,k]
#            else:
#                neg *= self._F[j]**(-self._B[j,k])
#        return (pos+neg)//self._F[k]
#
#
#    def find_cluster_variables(self, glist_tofind=[], num_mutations=infinity):
#        MCI = self.mutation_class_iter()
#        mutation_counter = 0
#        ## the last check should be done more efficiently
#        while mutation_counter < num_mutations and (glist_tofind == [] or not Set(glist_tofind).issubset(Set(self._F_dict.keys()))):
#            try:
#                MCI.next()
#            except:
#                break
#            mutation_counter += 1
#        print "Found after "+str(mutation_counter)+" mutations."
#
#
#    def mutation_class_iter(self, depth=infinity, show_depth=False, return_paths=False, up_to_equivalence=True, only_sink_source=False):
#        depth_counter = 0
#        n = self._n
#        timer = time.time()
#        if return_paths:
#            yield (self,[])
#        else:
#            yield self
#        if up_to_equivalence:
#            cl = Set(self._G.columns())
#        else:
#            cl = tuple(self._G.columns())
#        clusters = {}
#        clusters[ cl ] = [ self, range(n), [] ]
#        gets_bigger = True
#        if show_depth:
#            timer2 = time.time()
#            dc = str(depth_counter)
#            dc += ' ' * (5-len(dc))
#            nr = str(len(clusters))
#            nr += ' ' * (10-len(nr))
#            print "Depth: %s found: %s Time: %.2f s"%(dc,nr,timer2-timer)
#        while gets_bigger and depth_counter < depth:
#            gets_bigger = False
#            keys = clusters.keys()
#            for key in keys:
#                sd = clusters[key]
#                while sd[1]:
#                    i = sd[1].pop()
#                    if not only_sink_source or all( entry >= 0 for entry in sd[0]._B.row( i ) ) or all( entry <= 0 for entry in sd[0]._B.row( i ) ):
#                        sd2  = sd[0].mutate( i, inplace=False )
#                        if up_to_equivalence:
#                            cl2 = Set(sd2._G.columns())
#                        else:
#                            cl2 = tuple(sd2._G.columns())
#                        if cl2 in clusters:
#                            if not up_to_equivalence and i in clusters[cl2][1]:
#                                clusters[cl2][1].remove(i)
#                        else:
#                            gets_bigger = True
#                            if only_sink_source:
#                                orbits = range(n)
#                            else:
#                                orbits = [ index for index in xrange(n) if index > i or sd2._B[index,i] != 0 ]
#
#                            clusters[ cl2 ] = [ sd2, orbits, clusters[key][2]+[i] ]
#                            if return_paths:
#                                yield (sd2,clusters[cl2][2])
#                            else:
#                                yield sd2
#            depth_counter += 1
#            if show_depth and gets_bigger:
#                timer2 = time.time()
#                dc = str(depth_counter)
#                dc += ' ' * (5-len(dc))
#                nr = str(len(clusters))
#                nr += ' ' * (10-len(nr))
#                print "Depth: %s found: %s Time: %.2f s"%(dc,nr,timer2-timer)
#
#
#    def cluster_variable(self, k):
#        if k in range(self._n):
#            g_mon = prod([self._R.gen(i)**self._G[i,k] for i in xrange(self._n)])
#            F_std = self._F[k].subs(self._yhat)
#            F_trop = self._F[k].subs(self._y).denominator()
#            return g_mon*F_std*F_trop
#        if isinstance(k,tuple):
#            if not k in self._F_dict.keys():
#                self.find_cluster_variables([k])
#            g_mon = prod([self._R.gen(i)**k[i] for i in xrange(self._n)])
#            F_std = self._F_dict[k].subs(self._yhat)
#            F_trop = self._F_dict[k].subs(self._y).denominator()
#            return g_mon*F_std*F_trop
#
#        
#
#         
