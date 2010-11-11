Inclusion Procedure for New Packages
====================================

For a package to become part of Sage's standard distribution, it
must meet the following requirements:

- **License**. The license must be a GPL version 2+ compatible
  license.


- **Build Support**. The code **must** build on all the `fully supported platforms <http://wiki.sagemath.org/SupportedPlatforms#Fullysupported-SageisALWAYScheckedonALLtheseplatformsBEFOREareleaseismade/>`_.

  A standard package should also work on all the platforms where Sage is
  `expected to work <http://wiki.sagemath.org/SupportedPlatforms#Expectedtowork-Sagewillprobablywork.2Cbutitisnotalwaystested.>`_, but since we don't fully
  support these platforms and often lack the resources to test on them, you
  are not expected to confirm your packages works on those platforms.
  However, if you can, it is better to do so. As noted
  `here <http://wiki.sagemath.org/SupportedPlatforms#Expectedtowork-Sagewillprobablywork.2Cbutitisnotalwaystested.>`_, a failure of Sage to work on a
  platform where it is expected to work, will be considered a bug.

  There is no need to worry too much about platforms where Sage will
  `probably not work <http://wiki.sagemath.org/SupportedPlatforms#Probablywillnotwork-Portingworkmaybeongoing>`_ though if it's clear that there is
  significant effort taking place to port Sage to a platform, then you should
  aim to ensure your package does not cause unnecessary headaches to those
  working on the port.

  If it's clear that a port is stagnent, with nobody working on
  it, then you can safely ignore it.

  Remarks:

  - Some Sage developers are willing to help you port to OS X, Solaris
    and Windows. But this is no guarantee and you or your project are
    expected to do the heavy lifting and also support those ports
    upstream if there is no Sage developer who is willing to share the
    burden.
  - One of the best ways to ensure your code works on multiple platforms
    is to only use commands which are defined by `POSIX.1-2008 <http://www.opengroup.org/onlinepubs/9699919799/>`_ and only use options which are defined
    in the POSIX standard. For example, do not use the -p option to `uname <http://www.opengroup.org/onlinepubs/9699919799/utilities/uname.html>`_ as
    the '-p' option is not defined by the POSIX standard, so is not portable.
    If you must use a non-POSIX command, or a option which is not defined
    by POSIX, then ensure the code only gets executed on the platform(s)
    where that command and/or option will be acceptable.


- **Quality**. The code should be "better" than any other available
  code (that passes the two above criteria), and the authors need to
  justify this. The comparison should be made to both Python and other
  software. Criteria in passing the quality test include:

  - Speed

  - Documentation

  - Usability

  - Memory leaks

  - Maintainable

  - Portability

  - Reasonable build time, size, dependencies


- **Previously an optional package**. Usually a new standard package must have spent some time as an optional package. However, sometimes this is not possible, if for example a new library is needed to permit an updated version of a standard package to function.

-  **Refereeing**. The code must be refereed, as discussed in
   :ref:`chapter-trac`.
