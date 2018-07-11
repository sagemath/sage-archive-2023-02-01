# -*- coding: utf-8 -*-
r"""
Testing for databases at runtime
"""

import os

from . import StaticFile
from sage.env import SAGE_SHARE

class DatabaseCremona(StaticFile):
    r"""
    A :class:`Feature` which describes the presence of John Cremona's full
    database of elliptic curves.

    EXAMPLES::

        sage: from sage.features.databases import DatabaseCremona
        sage: DatabaseCremona().is_present()  # optional: database_cremona_ellcurve
        FeatureTestResult("Cremona's database of elliptic curves", True)
    """
    def __init__(self, name="cremona", spkg="database_cremona_ellcurve"):
        r"""
        TESTS::

            sage: from sage.features.databases import DatabaseCremona
            sage: isinstance(DatabaseCremona(), DatabaseCremona)
            True
        """
        StaticFile.__init__(self, "Cremona's database of elliptic curves",
            filename='{}.db'.format(name.replace(' ', '_')),
            search_path=[ os.path.join(SAGE_SHARE, 'cremona') ],
            search_path_override='CREMONA_DATA_DIR',
            spkg=spkg,
            url="https://github.com/JohnCremona/ecdata")
