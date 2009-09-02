from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

setup(
   cmdclass = {'build_ext': build_ext},
   ext_modules = [Extension("sage.numerical.mipGlpk",
           ["patch/mipGlpk.pyx"],
         include_dirs=["../../../local/include/","../../../devel/sage/c_lib/include/"],
         language='c++',
         libraries=["csage","stdc++","glpk"])]
   )
