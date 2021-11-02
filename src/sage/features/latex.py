# -*- coding: utf-8 -*-
r"""
Check for pdflatex and equivalent programs
"""
# ****************************************************************************
#       Copyright (C) 2021 Sebastien Labbe <slabqc@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ****************************************************************************

from . import Executable, FeatureTestResult

class latex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``latex``

    EXAMPLES::

        sage: from sage.features.latex import latex
        sage: latex().is_present()             # optional: latex
        FeatureTestResult('latex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import latex
            sage: isinstance(latex(), latex)
            True
        """
        Executable.__init__(self, "latex", executable="latex",
                            url="https://www.latex-project.org/")

    def is_functional(self):
        r"""
        Return whether `latex` in the path is functional.

        EXAMPLES:

            sage: from sage.features.latex import latex
            sage: latex().is_functional()             # optional: latex
            FeatureTestResult('latex', True)
        """
        from sage.misc.latex import _run_latex_, _latex_file_
        from sage.misc.temporary_file import tmp_filename
        try:
            f = tmp_filename(ext='.tex')
            O = open(f, 'w')
            O.write(_latex_file_('2+3'))
            O.close()
            _run_latex_(f)
            return FeatureTestResult(self, True)
        except Exception:
            return FeatureTestResult(self, False, reason="Running latex on a sample file raised an exception")

class pdflatex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``pdflatex``

    EXAMPLES::

        sage: from sage.features.latex import pdflatex
        sage: pdflatex().is_present()             # optional: pdflatex
        FeatureTestResult('pdflatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import pdflatex
            sage: isinstance(pdflatex(), pdflatex)
            True
        """
        Executable.__init__(self, "pdflatex", executable="pdflatex",
                            url="https://www.latex-project.org/")

class xelatex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``xelatex``

    EXAMPLES::

        sage: from sage.features.latex import xelatex
        sage: xelatex().is_present()             # optional: xelatex
        FeatureTestResult('xelatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import xelatex
            sage: isinstance(xelatex(), xelatex)
            True
        """
        Executable.__init__(self, "xelatex", executable="xelatex",
                            url="https://www.latex-project.org/")

class lualatex(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``lualatex``

    EXAMPLES::

        sage: from sage.features.latex import lualatex
        sage: lualatex().is_present()             # optional: lualatex
        FeatureTestResult('lualatex', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.latex import lualatex
            sage: isinstance(lualatex(), lualatex)
            True
        """
        Executable.__init__(self, "lualatex", executable="lualatex",
                            url="https://www.latex-project.org/")
