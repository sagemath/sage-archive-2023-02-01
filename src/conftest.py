from __future__ import annotations
from typing import Any
import pytest

# Ignore a few test files that are (not yet) using pytest
collect_ignore = [
    "sage/libs/gap/test_long.py",
    "sage/structure/test_factory.py",
    "sage/misc/nested_class_test.py",
    "sage/repl/rich_output/backend_test.py"
]


@pytest.fixture(autouse=True)
def add_imports(doctest_namespace: dict[str, Any]):
    """
    Add global imports for doctests.

    See `pytest documentation <https://docs.pytest.org/en/stable/doctest.html#doctest-namespace-fixture>`.
    """
    import sage.all  # type: ignore # implicitly used below by calling locals()
    doctest_namespace.update(**locals())
