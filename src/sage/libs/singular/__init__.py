from sage.libs.singular.function import singular_function, lib

from sage.libs.singular.option import LibSingularOptions

libsingular_options = LibSingularOptions()

## We predefine a couple of often used functions here to avoid the fetch overhead ##
groebner = singular_function('groebner')

