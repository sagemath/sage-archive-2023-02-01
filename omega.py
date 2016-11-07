
from sage.groups.indexed_free_group import IndexedFreeAbelianGroup
from six import iteritems

class OmegaGroupElement(IndexedFreeAbelianGroup.Element):

    def __init__(self, parent, x, normalize=True):
        if normalize:
            from sage.misc.misc_c import prod
            L = parent.factor_ring()
            constant = prod(
                (f.constant_coefficient()**e for f, e in iteritems(x)), L.one())
            x = {f/f.constant_coefficient(): e for f, e in iteritems(x)}
        super(OmegaGroup.Element, self).__init__(parent, x)


class OmegaGroup(IndexedFreeAbelianGroup):

    Element = OmegaGroupElement

    @staticmethod
    def __classcall__(cls, base_ring, names):
        from sage.structure.unique_representation import UniqueRepresentation
        names = LaurentPolynomialRing(QQ, names).variable_names()
        return UniqueRepresentation.__classcall__(cls, base_ring, names)


    def __init__(self, base_ring, names):
        r"""
        EXAMPLES::

            sage: G = OmegaGroup(QQ, 'x, y'); G
            OmegaGroup in x, y over Rational Field
        """
        from sage.rings.polynomial.laurent_polynomial_ring import LaurentPolynomialRing
        L = LaurentPolynomialRing(base_ring, names)
        self._factor_ring_ = L
        super(OmegaGroup, self).__init__(indices=L, prefix='',
                                         names=names,
                                         bracket='(', latex_bracket='(')


    def factor_ring(self):
        return self._factor_ring_


    def _repr_(self):
        return 'OmegaGroup in {} over {}'.format(
            ', '.join(self.variable_names()),
            self.factor_ring().base_ring())


    def _element_constructor_(self, x=None):
        r"""
        TESTS::

            sage: G = OmegaGroup(QQ, 'mu, x, y')
            sage: var('mu, x, y')
            (mu, x, y)
            sage: G(1 / (1-mu*x) / (1-y/mu))
            (1 - mu^-1*y)^-1*(-mu*x + 1)^-1
       """
        if isinstance(x, Expression):
            import operator
            from sage.symbolic.operators import mul_vararg
            factors = []
            if x.operator() == mul_vararg:
                x = list(
                    (f.operands() if f.operator() == operator.pow else (f, 1))
                    for f in x.operands())
            elif x.operator() == operator.pow:
                x = [x.operands()]
            else:
                x = [(x, 1)]
            L = self.factor_ring()
            x = list((L(f), e) for f, e in x)
        return super(OmegaGroup, self)._element_constructor_(x)
