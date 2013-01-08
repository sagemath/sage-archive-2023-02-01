r"""
Ensure that certain modules are not loaded on startup::

    sage: sage0("'numpy' in sys.modules")  # long time
    False
"""
