database_knotinfo: Tables of Knots and Links from the KnotInfo Databases
========================================================================

Description
-----------

Database for named knots and links provided at

https://knotinfo.math.indiana.edu/

and

https://linkinfo.sitehost.iu.edu'

Dependencies
------------

- Sage library


Upstream Contact
----------------

- Charles Livingston <livingst@indiana.edu>
- Allison H. Moore <moorea14@vcu.edu>

Update Instructions
-------------------

- See the Python script ``create_knotinfo_tarball.py`` in the current directory.

Changelog
---------

- 20200713 (Sebastian Oehms, 13 Juli 2020, :trac:`30352`, initial version)

  The tarball has been created from the both download files at the
  given date:

  ``knotinfo_data_complete.xls``
  ``linkinfo_data_complete.xlsx``

  exporting them to CSV via LibreOffice.

  The second file has been changed manually deleting one character:
  a trailing "}" occuring in the homfly_polynomial column of the last
  link ``L11n459{1,1,1}``.

- 20210201 (Sebastian Oehms, 1 February 2021, :trac:`30352`, upgrade)

  Three new columns have been added to ``knotinfo_data_complete.xls``
  (``montesinos_notation``, ``boundary_slopes`` and ``pretzel_notation``).
  ``linkinfo_data_complete.xlsx`` remains unchanged.

  The tarball has been created using ``create_knotinfo_tarball.py``.
  The fix concerning the misplaced character for ``L11n459{1,1,1}``
  is performed in :meth:`KnotInfoBase.homfly_polynomial`, now.
