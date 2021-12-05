# -*- conf-unix -*-
[metadata]
name = sagemath-environment
version = file: VERSION.txt
description = Sage: Open Source Mathematics Software: System and software environment
long_description = file: README.rst
long_description_content_type = text/x-rst
license = GNU General Public License (GPL) v2 or later
author = The Sage Developers
author_email = sage-support@googlegroups.com
url = https://www.sagemath.org

classifiers =
    Development Status :: 6 - Mature
    Intended Audience :: Education
    Intended Audience :: Science/Research
    License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
    Programming Language :: Python :: 3 :: Only
    Programming Language :: Python :: 3.7
    Programming Language :: Python :: 3.8
    Programming Language :: Python :: 3.9
    Programming Language :: Python :: Implementation :: CPython
    Topic :: Scientific/Engineering :: Mathematics

[options]
python_requires = >=3.7, <3.10
install_requires =
    esyscmd(`sage-get-system-packages install-requires \
        | sed "2,\$s/^/    /;"')dnl

py_modules =
    sage.all__sagemath_environment
    sage.env
    sage.version
    sage.misc.all__sagemath_environment
    sage.misc.package
    sage.misc.temporary_file
    sage.misc.viewer

packages =
    sage.features

[options.extras_require]
# sage.env can optionally use sage_conf
sage_conf = esyscmd(`sage-get-system-packages install-requires sage_conf')
