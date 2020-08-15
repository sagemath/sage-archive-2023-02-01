=  KnotInfo Database =

== Description ==

Database for named knots and links provided at

https://knotinfo.math.indiana.edu/

and

https://linkinfo.sitehost.iu.edu'

== Dependencies ==

 * Sage library

== Changelog ==

=== knotinfo-20200713.tar.bz2 (Sebastian Oehms, 13 Juli 2020) ===

 * #30352: Initial version

 The tarball has been created from the both download files at the
 given date:

 `knotinfo_data_complete.xls`
 `linkinfo_data_complete.xlsx` 

 exporting them to CSV via LibreOffice.

 The second file has been changed manually deleting one character:
 a trailing "}" occuring in the homfly_polynomial column of the last
 link `L11n459{1,1,1}`.
