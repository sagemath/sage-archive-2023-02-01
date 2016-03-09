"""
Tests for Available Nonpackages

This module contains tests to make up a list of available nonpackages (external softwares) that Sage interfaces.

AUTHORS:

- Kwankyu Lee (2016-03-09) -- initial version.
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

# The test for a nonpackage xxx has the name `test_xxx` and checks 
# the system has the software xxx installed and available, and 
# returns `True` if affirmative. If you write a new test, then insert 
# also an appropriate line in the function `test_nonpackages` at the end.

def test_latex():
    from sage.misc.latex import have_latex
    if have_latex():
        return True
    return False

def test_magma():
    from sage.interfaces.magma import magma
    try:
        magma(1)
        return True
    except Exception:
        return False 

def test_matlab():
    from sage.interfaces.matlab import matlab
    try:
        matlab(1)
        return True
    except Exception:
        return False 

def test_mathematica():
    from sage.interfaces.mathematica import mathematica
    try:
        mathematica(1)
        return True
    except Exception:
        return False  

def test_macaulay2():
    from sage.interfaces.macaulay2 import macaulay2
    try:
        macaulay2(1)
        return True
    except Exception:
        return False   

def test_octave():
    from sage.interfaces.octave import octave
    try:
        octave(1)
        return True
    except Exception:
        return False  

def test_scilab():
    from sage.interfaces.scilab import scilab
    try:
        scilab(1)
        return True
    except Exception:
        return False      

def test_maple():
    from sage.interfaces.maple import maple
    try:
        maple(1)
        return True
    except Exception:
        return False 

def test_cplex():
    from sage.numerical.mip import MixedIntegerLinearProgram
    try:
        MixedIntegerLinearProgram(solver='cplex')
        return True
    except Exception:
        return False

def test_gurobi():
    from sage.numerical.mip import MixedIntegerLinearProgram
    try:
        MixedIntegerLinearProgram(solver='gurobi')
        return True
    except Exception:
        return False   

def test_nonpackages():
    pkgs = list()
    if test_latex(): pkgs.append('latex')
    if test_magma(): pkgs.append('magma')
    if test_matlab(): pkgs.append('matlab')     
    if test_mathematica(): pkgs.append('mathematica')  
    if test_macaulay2(): pkgs.append('macaulay2')  
    if test_octave(): pkgs.append('octave')  
    if test_scilab(): pkgs.append('scilab')  
    if test_maple(): pkgs.append('maple')
    if test_cplex(): pkgs.append('cplex')  
    if test_gurobi(): pkgs.append('gurobi')  
    return pkgs

available_nonpackages = test_nonpackages()



