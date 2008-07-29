#! /usr/bin/env python

# To use:
#       python setup.py install

import os, sys
from distutils.core import setup, Extension

def src(*args):
  path = args[0]
  for d in args[1:]:
    path = os.path.join(path, d)
  names = filter(lambda name: name[-4:] == ".cpp", os.listdir(path))
  return map(lambda name: os.path.join(path, name), names)

def dep(*args):
  path = args[0]
  for d in args[1:]:
    path = os.path.join(path, d)
  names = filter(lambda name: name[-2:] == ".h", os.listdir(path))
  return map(lambda name: os.path.join(path, name), names)

util_src = src("ginv", "util")
util_dep = dep("ginv", "util")

monom_src = src("ginv", "monom")
monom_dep = dep("ginv", "monom")

coeff_src = src("ginv", "coeff")
coeff_dep = dep("ginv", "coeff")

poly_src = src("ginv", "poly")
poly_dep = dep("ginv", "poly")

criteria_src = src("ginv", "criteria")
criteria_dep = dep("ginv", "criteria")

division_src = src("ginv", "division")
division_dep = dep("ginv", "division")

algorithm_src = src("ginv", "algorithm")
algorithm_dep = dep("ginv", "algorithm")

gauss_src = src("ginv", "gauss")
gauss_dep = dep("ginv", "gauss")

modules_src = src("modules")
modules_dep = dep("modules")

sources = util_src + monom_src + coeff_src + poly_src + criteria_src + division_src + algorithm_src + gauss_src + modules_src
depends = util_dep + monom_dep + coeff_dep + poly_dep + criteria_dep + division_dep + algorithm_dep + gauss_dep + modules_dep

if sys.byteorder == 'little':
  define_macros = [('WORDS_LITTLEENDIAN', None)]
elif sys.byteorder == 'big':
  define_macros = [('WORDS_BIGENDIAN', None)]

undef_macros = None
include_dirs = None
library_dirs = None

if sys.platform in ['win32']:
  gmpcxx_src = src(os.path.join("gmp", "gmpcxx"))
  gmpcxx_dep = dep(os.path.join("gmp", "gmpcxx"))

  sources += gmpcxx_src
  depends += gmpcxx_dep

  include_dirs = ['gmp']
  library_dirs = ['gmp']
  libraries = ['gcc', 'gmp']
  define_macros.append(('WIN32', None))
  extra_compile_args = ['-O2']
elif sys.platform in ['darwin']:
  include_dirs = ['/sw/include']
  library_dirs = ['/sw/lib']
  libraries = ['gmp', 'gmpxx']
  extra_compile_args = None
else:
  libraries = ['gmp', 'gmpxx']
  #extra_compile_args = None
#  extra_compile_args = ['-O0', '-g']
#  undef_macros = ['NDEBUG']
  #extra_compile_args = ['-O3', '-march=k8', '-pipe', '-mmmx', '-msse', '-msse2', '-m3dnow', '-fomit-frame-pointer']
  extra_compile_args = ['-O3', '-pipe', '-mmmx', '-msse', '-m3dnow', '-fomit-frame-pointer']

# SAGE specific
  extra_compile_args.append("-I%s/include"%os.environ["SAGE_LOCAL"])


if sys.version_info[0] == 2 and sys.version_info[1] < 5:
  # PY_SSIZE_T_CLEAN
  define_macros.append(('Py_ssize_t', 'int'))
  define_macros.append(('lenfunc', 'inquiry'))
  define_macros.append(('ssizeargfunc', 'intargfunc'))

ext_modules = [
    Extension('ginv',
               sources = sources,
               depends = depends,
               extra_compile_args = extra_compile_args,
               define_macros = define_macros,
               #define_macros=[('NDEBUG', '1')],
               undef_macros = undef_macros,
               #undef_macros=['NDEBUG'],
               include_dirs = include_dirs,
               library_dirs = library_dirs,
               libraries = libraries)
    ]

setup(name="ginv",
      version = "1.9",
      description = "ginv",
      author = "Blinkov Yu.A.",
      author_email = "BlinkovUA@info.sgu.ru",
      url = "http://invo.jinr.ru",
      ext_modules = ext_modules,
     )
