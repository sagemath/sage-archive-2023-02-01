r"""
Index notation for tensors.

AUTHORS:

- Eric Gourgoulhon, Michal Bejger (2014): initial version

"""
#******************************************************************************
#       Copyright (C) 2014 Eric Gourgoulhon <eric.gourgoulhon@obspm.fr>
#       Copyright (C) 2014 Michal Bejger <bejger@camk.edu.pl>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#******************************************************************************

from sage.structure.sage_object import SageObject

class TensorWithIndices(SageObject):
    r"""
    Index notation for tensors.
    
    This is a technical class to allow one to write some tensor operations 
    (contractions and symmetrizations) in index notation.
    
    INPUT:
    
    - ``tensor`` -- a tensor (or a tensor field)
    - ``indices`` -- string containing the indices, as single letters; the
      contravariant indices must be stated first and separated from the 
      covariant indices by the character '_'. 
    
    EXAMPLES:
    
    Index representation of tensors on a rank-3 free module::
    
        sage: M = FiniteRankFreeModule(QQ, 3, name='M')
        sage: e = M.basis('e')
        sage: a = M.tensor((2,0), name='a')
        sage: a[:] = [[1,2,3], [4,5,6], [7,8,9]]
        sage: b = M.tensor((0,2), name='b')
        sage: b[:] = [[-1,2,-3], [-4,5,6], [7,-8,9]]
        sage: t = a*b ; t.set_name('t') ; t
        type-(2,2) tensor t on the rank-3 free module M over the Rational Field
        sage: from sage.tensor.modules.tensor_with_indices import TensorWithIndices
        sage: T = TensorWithIndices(t, '^ij_kl') ; T
        t^ij_kl

    The :class:`TensorWithIndices` object is returned by the square 
    bracket operator acting on the tensor and fed with the string specifying
    the indices::
    
        sage: a['^ij']
        a^ij
        sage: type(a['^ij'])
        <class 'sage.tensor.modules.tensor_with_indices.TensorWithIndices'>
        sage: b['_ef']
        b_ef
        sage: t['^ij_kl']
        t^ij_kl
        
    The symbol '^' may be omitted, since the distinction between covariant 
    and contravariant indices is performed by the index position relative to 
    the symbol '_'::
    
        sage: t['ij_kl']
        t^ij_kl

    Also, LaTeX notation may be used::
     
        sage: t['^{ij}_{kl}']
        t^ij_kl

    If some operation is asked in the index notation, the resulting tensor
    is returned, not a :class:`TensorWithIndices` object; for instance, for 
    a symmetrization::
    
        sage: s = t['^(ij)_kl'] ; s  # the symmetrization on i,j is indicated by parentheses
        type-(2,2) tensor on the rank-3 free module M over the Rational Field
        sage: s.symmetries()
        symmetry: (0, 1);  no antisymmetry
        sage: s == t.symmetrize((0,1))
        True
    
    The letters denoting the indices can be chosen freely; since they carry no
    information, they can even be replaced by dots::
    
        sage: t['^(..)_..'] == t.symmetrize((0,1))
        True

    Similarly, for an antisymmetrization::
    
        sage: s = t['^ij_[kl]'] ; s # the symmetrization on k,l is indicated by square brackets
        type-(2,2) tensor on the rank-3 free module M over the Rational Field
        sage: s.symmetries()
        no symmetry;  antisymmetry: (2, 3)
        sage: s == t.antisymmetrize((2,3))
        True

    Another example of an operation indicated by indices is a contraction::
    
        sage: s = t['^ki_kj'] ; s  # contraction on the repeated index k
        endomorphism on the rank-3 free module M over the Rational Field
        sage: s == t.self_contract(0,2)
        True

    Indices not involved in the contraction may be replaced by dots::
    
        sage: s == t['^k._k.']
        True

    The contraction of two tensors is indicated by repeated indices and 
    the * operator::
    
        sage: s = a['^ik']*b['_kj'] ; s
        endomorphism on the rank-3 free module M over the Rational Field
        sage: s == a.contract(1, b, 0)
        True
        sage: s = t['^.k_..']*b['_.k'] ; s
        type-(1,3) tensor on the rank-3 free module M over the Rational Field
        sage: s == t.contract(1, b, 1)
        True
        sage: t['^{ik}_{jl}']*b['_{mk}'] == s # LaTeX notation
        True
    
    """
    def __init__(self, tensor, indices):
        self._tensor = tensor # may be changed below
        self._changed = False # indicates whether self contains an altered 
                              # version of the original tensor (True if 
                              # symmetries or contractions are indicated in the 
                              # indices)
        # Suppress all '{' and '}' comming from LaTeX notations:
        indices = indices.replace('{','').replace('}','')
        # Suppress the first '^':
        if indices[0] == '^':
            indices = indices[1:]
        if '^' in indices:
            raise IndexError("The contravariant indices must be placed first.")
        con_cov = indices.split('_')
        con = con_cov[0]
        # Contravariant indices
        # ---------------------
        #  Search for (anti)symmetries:
        if '(' in con:
            sym1 = con.index('(')
            sym2 = con.index(')')-2
            if con.find('(', sym1+1) != -1 or '[' in con:
                raise NotImplementedError("Multiple symmetries are not " + 
                                          "treated yet.")
            self._tensor = self._tensor.symmetrize(range(sym1, sym2+1))
            self._changed = True # self does no longer contain the original tensor
            con = con.replace('(','').replace(')','')
        if '[' in con:
            sym1 = con.index('[')
            sym2 = con.index(']')-2
            if con.find('[', sym1+1) != -1 or '(' in con:
                raise NotImplementedError("Multiple symmetries are not " + 
                                          "treated yet.")
            self._tensor = self._tensor.antisymmetrize(range(sym1, sym2+1))
            self._changed = True # self does no longer contain the original tensor
            con = con.replace('[','').replace(']','')
        if len(con) != self._tensor._tensor_type[0]:
            raise IndexError("Number of contravariant indices not compatible " + 
                            "with the tensor type.")
        self._con = con
        # Covariant indices
        # -----------------
        if len(con_cov) == 1:
            if tensor._tensor_type[1] != 0:
                raise IndexError("Number of covariant indices not compatible " + 
                                "with the tensor type.")
            self._cov = ''
        elif len(con_cov) == 2:
            cov = con_cov[1]
            #  Search for (anti)symmetries:
            if '(' in cov:
                sym1 = cov.index('(')
                sym2 = cov.index(')')-2
                if cov.find('(', sym1+1) != -1 or '[' in cov:
                    raise NotImplementedError("Multiple symmetries are not " + 
                                              "treated yet.")
                csym1 = sym1 + self._tensor._tensor_type[0]
                csym2 = sym2 + self._tensor._tensor_type[0]
                self._tensor = self._tensor.symmetrize(range(csym1, csym2+1))
                self._changed = True # self does no longer contain the original 
                                     # tensor
                cov = cov.replace('(','').replace(')','')
            if '[' in cov:
                sym1 = cov.index('[')
                sym2 = cov.index(']')-2
                if cov.find('[', sym1+1) != -1 or '(' in cov:
                    raise NotImplementedError("Multiple symmetries are not " + 
                                              "treated yet.")
                csym1 = sym1 + self._tensor._tensor_type[0]
                csym2 = sym2 + self._tensor._tensor_type[0]
                self._tensor = self._tensor.antisymmetrize(range(csym1, csym2+1))
                self._changed = True # self does no longer contain the original 
                                     # tensor
                cov = cov.replace('[','').replace(']','')
            if len(cov) != tensor._tensor_type[1]:
                raise IndexError("Number of covariant indices not compatible " + 
                                "with the tensor type.")
            self._cov = cov
        else:
            raise IndexError("Two many '_' in the list of indices.")
        # Treatment of possible self-contractions:
        # ---------------------------------------
        contraction_pairs = []
        for ind in self._con:
            if ind != '.' and ind in self._cov:
                pos1 = self._con.index(ind)
                pos2 = self._tensor._tensor_type[0] + self._cov.index(ind)
                contraction_pairs.append((pos1, pos2))
        if len(contraction_pairs) > 1:
            raise NotImplementedError("Multiple contractions are not " + 
                                      "implemented yet.")
        if len(contraction_pairs) == 1:
            pos1 = contraction_pairs[0][0]
            pos2 = contraction_pairs[0][1]
            self._tensor = self._tensor.self_contract(pos1, pos2)
            self._changed = True # self does no longer contain the original 
                                 # tensor
            ind = self._con[pos1]
            self._con = self._con.replace(ind, '')
            self._cov = self._cov.replace(ind, '')


    def _repr_(self):
        r"""
        String representation of the object.
        """
        if self._tensor._name is not None:
            name = self._tensor._name
        else:
            name = 'X'
        if self._con == '':
            return name + '_' + self._cov
        elif self._cov == '':
            return name + '^' + self._con
        else:
            return name + '^' + self._con + '_' + self._cov

    def update(self):
        r"""
        Return the tensor contains in ``self`` if it differs from that used
        for creating ``self``, otherwise return ``self``. 
        """
        if self._changed:
            return self._tensor
        else:
            return self

    def __mul__(self, other):
        r"""
        Tensor contraction on specified indices
        """
        if not isinstance(other, TensorWithIndices):
            raise TypeError("The second item of * must be a tensor with " + 
                            "specified indices.")
        contraction_pairs = []
        for ind in self._con:
            if ind != '.':
                if  ind in other._cov:
                    pos1 = self._con.index(ind)
                    pos2 = other._tensor._tensor_type[0] + other._cov.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if  ind in other._con:
                    raise IndexError("The index " + str(ind) + " appears twice "
                                    + "in a contravariant position.")
        for ind in self._cov:
            if ind != '.':
                if ind in other._con:
                    pos1 = self._tensor._tensor_type[0] + self._cov.index(ind)
                    pos2 = other._con.index(ind)
                    contraction_pairs.append((pos1, pos2))
                if ind in other._cov:
                    raise IndexError("The index " + str(ind) + " appears twice "
                                    + "in a covariant position.")
        if contraction_pairs == []:
            # No contraction is performed: the tensor product is returned
            return self._tensor * other._tensor
        if len(contraction_pairs) > 1:
            raise NotImplementedError("Multiple contractions are not " + 
                                      "implemented yet.")
        pos1 = contraction_pairs[0][0]
        pos2 = contraction_pairs[0][1]
        return self._tensor.contract(pos1, other._tensor, pos2)

        
