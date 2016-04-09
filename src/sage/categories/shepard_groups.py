r"""
Shepard Groups
"""
#*****************************************************************************
#  Copyright (C) 2016 Travis Scrimshaw <tscrim at ucdavis.edu>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************
from sage.categories.category_singleton import Category_singleton
from sage.categories.generalized_coxeter_groups import GeneralizedCoxeterGroups
from sage.misc.cachefunc import cached_method


class ShepardGroups(Category_singleton):
    r"""
    The category of Shepard groups.

    EXAMPLES::

        sage: from sage.categories.shepard_groups import ShepardGroups
        sage: C = ShepardGroups(); C
        Category of shepard groups

    TESTS::

        sage: TestSuite(C).run()
    """
    @cached_method
    def super_categories(self):
        """
        EXAMPLES::

            sage: from sage.categories.shepard_groups import ShepardGroups
            sage: ShepardGroups().super_categories()
            [Category of finite well generated generalized coxeter groups]
        """
        return [GeneralizedCoxeterGroups().Finite()]
