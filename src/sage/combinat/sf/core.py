r"""
A place for things that may be useful in core sage (not specific to k-combinatorics)
"""
# ^*^ sphinx insert ^*^


def summands(poly):
    r""" Iterate through the summands of a symmetric function.

    For example, ``(s[2, 1] + s[3])`` has summands ``s[2, 1]`` and ``s[3]``.
    """
    parent_basis = poly.parent()
    return (coeff * parent_basis(index) for index, coeff in poly)

def prod(lis):
    return reduce(operator.mul, lis, 1)

def is_k_schur(obj):
    # checks if obj is a k-schur function (coming from the 'kSchur_with_category' class)
    try:
        classname = obj.parent().__class__.__name__
        return classname == 'kSchur_with_category'
    except:
        return False


# class InfiniteDimensionalFreeRing (CommutativeRing, InfiniteDimensionalFreeAlgebra):
#   pass

# base_ring=IntegerRing()
# algebras = Algebras(base_ring.category()).WithBasis()
# commutative_rings = CommutativeRings()
# F = ForgetfulFunctor(algebras, commutative_rings)
# InfiniteDimensionalFreeRing = F(InfiniteDimensionalFreeAlgebra)

# idea: coerce

# idea: manually change the category
# found '_init_category_' '_initial_coerce_list' '_initial_convert_list' '_unset_category' 'category' 'categories' 'coerce' 'hom' 'is_ring'

# Let us declare a coercion from `\ZZ[x]` to `\ZZ[z]`::
#  |
#  |                  sage: Z.<z> = ZZ[]
#  |                  sage: phi = Hom(X, Z)(z)
#  |                  sage: phi(x^2+1)
#  |                  z^2 + 1
#  |                  sage: phi.register_as_coercion()
#  |
#  |              Now we can add elements from `\ZZ[x]` and `\ZZ[z]`, because
#  |              the elements of the former are allowed to be implicitly
#  |              coerced into the later::
#  |
#  |                  sage: x^2 + z
#  |                  z^2 + z

# idea: patch SymmetricFunctions to accept Algebras, not just 'commutative rings'
