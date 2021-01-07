import sage.all
from sage.manifolds.differentiable.symplectic_form import SymplecticForm
from sage.manifolds.differentiable.examples.symplectic_vector_space import SymplecticVectorSpace
import pytest


class TestR2VectorSpace:

    @pytest.fixture
    def M(self):
        return SymplecticVectorSpace(2, 'R2', symplectic_name='omega')

    @pytest.fixture
    def omega(self, M: SymplecticVectorSpace):
        return M.symplectic_form()

    def test_repr(self, M: SymplecticVectorSpace):
        assert M._repr_() == '2-dimensional symplectic vector space R2'

    def test_display(self, omega: SymplecticForm):
        assert str(omega.display()) == r'omega = -dq/\dp'


class TestR4VectorSpace:

    @pytest.fixture
    def M(self):
        return SymplecticVectorSpace(4, 'R4', symplectic_name='omega')

    @pytest.fixture
    def omega(self, M: SymplecticVectorSpace):
        return M.symplectic_form()

    def test_repr(self, M: SymplecticVectorSpace):
        assert M._repr_() == '4-dimensional symplectic vector space R4'

    def test_display(self, omega: SymplecticForm):
        assert str(omega.display()) == r'omega = -dq1/\dp1 - dq2/\dp2'
