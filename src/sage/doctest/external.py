"""
Available External Softwares

This module makes up a list of external softwares that Sage interfaces. Availability 
of each software is tested only when necessary.

AUTHORS:

- Kwankyu Lee (2016-03-09) -- initial version, based on the codes 
                              by Robert Bradshaw and Nathann Cohen 
"""

#*****************************************************************************
#       Copyright (C) 2016 KWANKYU LEE <ekwankyu@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from util import LazySet

# Functions in this module whose name is of the form 'has_xxx' tests if the 
# software xxx is available to Sage.
prefix = 'has_'   

def has_internet():
    from six.moves import urllib
    try:
        urllib.request.urlopen("http://sagemath.org")
        return True
    except Exception:
        return False

def has_latex():
    from sage.misc.latex import view
    try:
        view('$e^{\pi}+1=0$')
        return True
    except Exception:
        return False

def has_magma():
    from sage.interfaces.magma import magma
    try:
        magma('2+3')
        return True
    except Exception:
        return False 

def has_matlab():
    from sage.interfaces.matlab import matlab
    try:
        matlab('2+3')
        return True
    except Exception:
        return False 

def has_mathematica():
    from sage.interfaces.mathematica import mathematica
    try:
        mathematica('2+3')
        return True
    except Exception:
        return False  

def has_maple():
    from sage.interfaces.maple import maple
    try:
        maple('2+3')
        return True
    except Exception:
        return False   

def has_macaulay2():
    from sage.interfaces.macaulay2 import macaulay2
    try:
        macaulay2('2+3')
        return True
    except Exception:
        return False   

def has_octave():
    from sage.interfaces.octave import octave
    try:
        octave('2+3')
        return True
    except Exception:
        return False  

def has_scilab():
    from sage.interfaces.scilab import scilab
    try:
        scilab('2+3')
        return True
    except Exception:
        return False

def has_cplex():
    from sage.numerical.mip import MixedIntegerLinearProgram
    try:
        MixedIntegerLinearProgram(solver='cplex')
        return True
    except Exception:
        return False

def has_gurobi():
    from sage.numerical.mip import MixedIntegerLinearProgram
    try:
        MixedIntegerLinearProgram(solver='gurobi')
        return True
    except Exception:
        return False   

def external_softwares():
    """
    Return the list of external softwares supported by this module.
    """
    supported = list()
    for func in globals():
        if func.startswith(prefix):
            supported.append(func[len(prefix):])
    return supported

external_softwares = external_softwares()

def lookup(software):
    """
    Return true if the software is available on the system.
    """
    if software in external_softwares:
        func = globals().get(prefix + software)
        return func()
    else:
        return False

available_softwares = LazySet(lookup)
