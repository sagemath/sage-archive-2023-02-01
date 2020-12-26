import pytest
from .generic_backend import GenericBackend

class GenericBackendTests:

    @pytest.fixture
    def backend(self) -> GenericBackend:
        raise NotImplementedError

    def test_ncols_nonnegative(self, backend: GenericBackend):
        assert backend.ncols() >= 0

    def test_sage_unittest_testsuite(self, backend: GenericBackend):
        # TODO: Remove this test as soon as all old test methods are migrated
        from sage.misc.sage_unittest import TestSuite
        TestSuite(backend).run(verbose=True, raise_on_failure=True, skip="_test_pickling")
