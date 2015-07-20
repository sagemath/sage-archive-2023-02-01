r"""
Index of encoders

The ``codes.encoders`` object may be used to access the encoders that Sage can build.

**Generic encoders**

- :func:`linear_code.LinearCodeGeneratorMatrixEncoder <sage.coding.linear_code.LinearCodeGeneratorMatrixEncoder>`
- :func:`linear_code.LinearCodeParityCheckEncoder <sage.coding.linear_code.LinearCodeParityCheckEncoder>`

**Generalized Reed-Solomon code encoders**

- :func:`grs.GRSEvaluationVectorEncoder <sage.coding.grs.GRSEvaluationVectorEncoder>`
- :func:`grs.GRSEvaluationPolynomialEncoder <sage.coding.grs.GRSEvaluationPolynomialEncoder>`

.. NOTE::

    To import these names into the global namespace, use:

        sage: from sage.coding.encoders_catalog import *
"""

from linear_code import (LinearCodeGeneratorMatrixEncoder, LinearCodeParityCheckEncoder)
from grs import (GRSEvaluationVectorEncoder, GRSEvaluationPolynomialEncoder)
