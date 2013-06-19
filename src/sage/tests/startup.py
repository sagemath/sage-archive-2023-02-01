r"""
Ensure that certain modules are not loaded on startup.

EXAMPLES::

    sage: 'numpy' in sys.modules
    False
    sage: 'sage.libs.gap.libgap' in sys.modules
    False
"""
