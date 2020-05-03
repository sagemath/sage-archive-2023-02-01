numpy
=====

Description
-----------

This package adds numerical linear algebra and other numerical computing
capabilities to python.

.. _upstream_contact:

Upstream Contact
----------------

-  Travis Oliphant
-  Fernando Perez
-  Brian Granger

Dependencies
------------

-  GNU patch
-  Python
-  Lapack
-  Blas
-  Atlas
-  Fortran

.. _special_updatebuild_instructions:

Special Update/Build Instructions
---------------------------------

-  Scipy uses numpy's distutils to control its compilation of fortran
   code.

   Whenever numpy is updated it is necessary to make sure that scipy
   still builds ok.
