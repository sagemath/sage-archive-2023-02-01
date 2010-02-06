Inclusion Procedure for New Packages
====================================

For a package to become part of Sage's standard distribution, it
must meet the following requirements:

- **License**. The license must be a GPL version 2+ compatible
  license.

- **Build Support**. The code must build on the following supported
  architectures and compilers (and intended port targets):

  - Linux: x86, x86_64, Itanium, ppc, ppc64, Sparc (gcc 3.4--4.3)

  - Apple Mac OS X: ppc, ppc64, x86, x86_64 (Xcode 2.5+)

  - Microsoft Windows: x86, x86_64 MSVC 2005/Intel Fortran (MinGW or
    Cygwin support is insufficient!)

  - Solaris 10: Sparc, x86, x86_64 (Sun Forte 12)

  Remarks:

  - Some Sage developers are willing to help you port to OS X, Solaris
    and Windows. But this is no guarantee and you or your project are
    expected to do the heavy lifting and also support those ports
    upstream if there is no Sage developer who is willing to share the
    burden.

  Potential future ports include FreeBSD (x86, x86_64), OpenBSD (x86,
  x86_64), HPUX (Itanium), AIX (PPC64), and ARM (OS X).

- **Quality**. The code should be "better" than any other available
  code (that passes the two above criteria), and the authors need to
  justify this. The comparison should be made to both Python and other
  software. Criteria in passing the quality test include:

  - Speed

  - Documentation

  - Usability

  - Memory leaks

  - Maintainable

  - Reasonable build time, size, dependencies

-  **Refereeing**. The code must be refereed, as discussed in
   :ref:`chapter-trac`.
