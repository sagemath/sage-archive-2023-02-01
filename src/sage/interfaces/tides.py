r"""
This module contains tools to write the .c files needed for TIDES [TI]_ .

Tides is an integration engine based on the Taylor method. It is implemented
as a c library. The user must translate its IVP into a pair of .c files that
will then be compiled and linked against the TIDES library. The reulting binary
will produce the desored output. The tools in this module can be used to
automate the generation of these files from the symbolic expression of the
differential equation.

##########################################################################
#  Copyright (C) 2014 Miguel Marco <mmarco@unizar.es>, Marcos Rodriguez
#   <marcos@uunizar.es>
#
#  Distributed under the terms of the GNU General Public License (GPL):
#
#                  http://www.gnu.org/licenses/
##########################################################################

AUTHORS:

- Miguel Marco (06-2014) - Implementation of tides solver

- Marcos Rodriguez (06-2014) - Implementation of tides solver

- Alberto Abad (06-2014) - tides solver

- Roberto Barrio (06-2014) - tides solver

REFERENCES:

.. [ALG924] A. Abad, R. Barrio, F. Blesa, M. Rodriguez. Algorithm 924. *ACM
Transactions on Mathematical Software*, *39*(1), 1–28.

.. [TI](http://www.unizar.es/acz/05Publicaciones/Monografias/MonografiasPublicadas/Monografia36/IndMonogr36.htm)
A. Abad, R. Barrio, F. Blesa, M. Rodriguez.
TIDES tutorial: Integrating ODEs by using the Taylor Series Method.
"""



from  sage.rings.real_mpfr import RealField
import shutil
import os
from sage.calculus.all import symbolic_expression
from sage.misc.flatten import flatten
from sage.ext.fast_callable import fast_callable
from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.semirings.non_negative_integer_semiring', 'NN')
from sage.misc.functional import N
from sage.functions.log import log
from sage.functions.other import floor, sqrt
from sage.functions.trig import sin, cos, arcsin, arctan, arccos


def subexpressions_list(f, parameters=[]):
    """
    Construct the lists with the intermediate steps on the evaluation of the
    function.

    INPUT:

    - ``f`` -- a symbollic function of several components.

    - ``paramters`` -- a list of the parameters that appear in the function
    this should be the symbollic constants that appear in f but are not
    arguments.

    OTUPUT:

    - a list of the intermediate subexpressions that appear in the evaluation
    of f.

    - a list with the operations used to construct each of the subexpressions.
    each element of this list is a tuple, formed by a string describing the
    operation made, and the operands.

    For the trigonometric functions, some extra expressions will be added.
    These extra expressions will be used later to compute their derivatives.


    EXAMPLES::

        sage: from sage.interfaces.tides import subexpressions_list
        sage: var('x,y')
        (x, y)
        sage: f(x,y) = [x^2+y, cos(x)/log(y)]
        sage: subexpressions_list(f)
        ([x^2, x^2 + y, sin(x), cos(x), log(y), cos(x)/log(y)],
        [('mul', x, x),
        ('add', y, x^2),
        ('sin', x),
        ('cos', x),
        ('log', y),
        ('div', log(y), cos(x))])

    ::

        sage: f(a)=[cos(a), arctan(a)]
        sage: from sage.interfaces.tides import subexpressions_list
        sage: subexpressions_list(f)
        ([sin(a), cos(a), a^2, a^2 + 1, arctan(a)],
        [('sin', a), ('cos', a), ('mul', a, a), ('add', 1, a^2), ('atan', a)])


    """
    variables = f[0].arguments()
    varpar = list(parameters) + list(variables)
    F = symbolic_expression([i(*variables) for i in f]).function(*varpar)
    lis = flatten([fast_callable(i,vars=varpar).op_list() for i in F], max_level=1)
    deflist = []
    stack = []
    const =[]
    stackcomp=[]
    detail=[]
    for i in lis:
        if i[0] == 'load_arg':
            stack.append(varpar[i[1]])
        elif i[0] == 'ipow':
            if i[1] in NN:
                basis = stack[-1]
                for j in range(i[1]-1):
                    a=stack.pop(-1)
                    detail.append(('mul', a, basis))
                    stack.append(a*basis)
                    stackcomp.append(stack[-1])
            else:
                detail.append(('pow',stack[-1],i[1]))
                stack[-1]=stack[-1]**i[1]
                stackcomp.append(stack[-1])

        elif i[0] == 'load_const':
            const.append(i[1])
            stack.append(i[1])
        elif i == 'mul':
            a=stack.pop(-1)
            b=stack.pop(-1)
            detail.append(('mul', a, b))
            stack.append(a*b)
            stackcomp.append(stack[-1])

        elif i == 'div':
            a=stack.pop(-1)
            b=stack.pop(-1)
            detail.append(('div', a, b))
            stack.append(b/a)
            stackcomp.append(stack[-1])

        elif i == 'add':
            a=stack.pop(-1)
            b=stack.pop(-1)
            detail.append(('add',a,b))
            stack.append(a+b)
            stackcomp.append(stack[-1])

        elif i == 'pow':
            a=stack.pop(-1)
            b=stack.pop(-1)
            detail.append(('pow', b, a))
            stack.append(b**a)
            stackcomp.append(stack[-1])

        elif i[0] == 'py_call' and str(i[1])=='log':
            a=stack.pop(-1)
            detail.append(('log', a))
            stack.append(log(a))
            stackcomp.append(stack[-1])

        elif i[0] == 'py_call' and str(i[1])=='exp':
            a=stack.pop(-1)
            detail.append(('exp', a))
            stack.append(exp(a))
            stackcomp.append(stack[-1])

        elif i[0] == 'py_call' and str(i[1])=='sin':
            a=stack.pop(-1)
            detail.append(('sin', a))
            detail.append(('cos', a))
            stackcomp.append(sin(a))
            stackcomp.append(cos(a))
            stack.append(sin(a))

        elif i[0] == 'py_call' and str(i[1])=='cos':
            a=stack.pop(-1)
            detail.append(('sin', a))
            detail.append(('cos', a))
            stackcomp.append(sin(a))
            stackcomp.append(cos(a))
            stack.append(cos(a))

        elif i[0] == 'py_call' and str(i[1])=='tan':
            a=stack.pop(-1)
            b = sin(a)
            c = cos(a)
            detail.append(('sin', a))
            detail.append(('cos', a))
            detail.append(('div', b, c))
            stackcomp.append(b)
            stackcomp.append(c)
            stackcomp.append(b/c)
            stack.append(b/c)

        elif i[0] == 'py_call' and str(i[1])=='arctan':
            a=stack.pop(-1)
            detail.append(('mul', a, a))
            detail.append(('add', 1, a*a))
            detail.append(('atan', a))
            stackcomp.append(a*a)
            stackcomp.append(1+a*a)
            stackcomp.append(arctan(a))
            stack.append(arctan(a))

        elif i[0] == 'py_call' and str(i[1])=='arcsin':
            a=stack.pop(-1)
            detail.append(('mul', a, a))
            detail.append(('mul', -1, a*a))
            detail.append(('add', 1, -a*a))
            detail.append(('pow', 1- a*a, 0.5))
            detail.append(('asin', a))
            stackcomp.append(a*a)
            stackcomp.append(-a*a)
            stackcomp.append(1-a*a)
            stackcomp.append(sqrt(1-a*a))
            stackcomp.append(arcsin(a))
            stack.append(arcsin(a))

        elif i[0] == 'py_call' and str(i[1])=='arccos':
            a=stack.pop(-1)
            detail.append(('mul', a, a))
            detail.append(('mul', -1, a*a))
            detail.append(('add', 1, -a*a))
            detail.append(('pow', 1- a*a, 0.5))
            detail.append(('mul', -1, sqrt(1-a*a)))
            detail.append(('acos', a))
            stackcomp.append(a*a)
            stackcomp.append(-a*a)
            stackcomp.append(1-a*a)
            stackcomp.append(sqrt(1-a*a))
            stackcomp.append(-sqrt(1-a*a))
            stackcomp.append(arccos(a))
            stack.append(arccos(a))

        elif i[0] == 'py_call' and 'sqrt' in str(i[1]):
            a=stack.pop(-1)
            detail.append(('pow', a, 0.5))
            stackcomp.append(sqrt(a))
            stack.append(sqrt(a))


        elif i == 'neg':
            a = stack.pop(-1)
            detail.append(('mul', -1, a))
            stack.append(-a)
            stackcomp.append(-a)

    return stackcomp,detail



def remove_repeated(l1, l2):
    """
    Given two lists, remove the repeated elements in l1, and the elements
    in l2 that are on the same position.
    positions.

    EXAMPLES::

        sage: from sage.interfaces.tides import (subexpressions_list, remove_repeated)
        sage: f(a)=[1 + a^2, arcsin(a)]
        sage: l1, l2 = subexpressions_list(f)
        sage: l1, l2
        ([a^2, a^2 + 1, a^2, -a^2, -a^2 + 1, sqrt(-a^2 + 1), arcsin(a)],
        [('mul', a, a),
        ('add', 1, a^2),
        ('mul', a, a),
        ('mul', -1, a^2),
        ('add', 1, -a^2),
        ('pow', -a^2 + 1, 0.5),
        ('asin', a)])
        sage: remove_repeated(l1, l2)
        sage: l1, l2
        ([a^2, a^2 + 1, -a^2, -a^2 + 1, sqrt(-a^2 + 1), arcsin(a)],
        [('mul', a, a),
        ('add', 1, a^2),
        ('mul', -1, a^2),
        ('add', 1, -a^2),
        ('pow', -a^2 + 1, 0.5),
        ('asin', a)])


    """
    for i in range(len(l1)-1):
        j=i+1
        while j<len(l1):
            if l1[j] == l1[i]:
                l1.pop(j)
                l2.pop(j)
            else:
                j+=1



def remove_constants(l1,l2):
    """
    Given two lists, remove the entries in the first that are real constants,
    and also the corresponding elements in the second one.

        sage: from sage.interfaces.tides import subexpressions_list, remove_constants
        sage: f(a)=[1+cos(7)*a]
        sage: l1, l2 = subexpressions_list(f)
        sage: l1, l2
        ([sin(7), cos(7), a*cos(7), a*cos(7) + 1],
        [('sin', 7), ('cos', 7), ('mul', cos(7), a), ('add', 1, a*cos(7))])
        sage: remove_constants(l1,l2)
        sage: l1, l2
        ([a*cos(7), a*cos(7) + 1], [('mul', cos(7), a), ('add', 1, a*cos(7))])

    """
    i=0
    while i < len(l1):
        if l1[i] in RealField():
            l1.pop(i)
            l2.pop(i)
        else:
            i+=1



def genfiles_mintides(integrator, driver, f, ics, initial, final, delta,
                      tolrel=1e-16, tolabs=1e-16, output = ''):
    r"""
    Generate the needed files for the min_tides library.

    INPUT:

    - ``integrator`` -- the name of the integrator file.

    - ``driver`` -- the name of the driver file.

    - ``f`` -- the function that determines the differential equation.

    - ``ics`` -- a list or tuple with the initial conditions.

    - ``initial`` -- the initial time for the integration.

    - ``final`` -- the final time for the integration.

    - ``delta`` -- the step of the output.

    - ``tolrel`` -- the relative tolerance.

    - ``tolabs`` -- the absolute tolerance.

    -  ``output`` -- the name of the file that the compiled integrator will write to

    This function creates two files, integrator and driver, that can be used
    later with the min_tides library [TI]_.


    TESTS::

        sage: from tempfile import mkdtemp
        sage: from sage.interfaces.tides import genfiles_mintides
        sage: import os
        sage: import shutil
        sage: var('t,x,y,X,Y')
        (t, x, y, X, Y)
        sage: f(t,x,y,X,Y)=[X, Y, -x/(x^2+y^2)^(3/2), -y/(x^2+y^2)^(3/2)]
        sage: tempdir = mkdtemp()
        sage: intfile = tempdir + '/integrator.c'
        sage: drfile = tempdir + '/driver.c'
        sage: genfiles_mintides(intfile, drfile, f, [1,0, 0, 0.2], 0, 10, 0.1, output = 'out')
        sage: fileint = open(intfile)
        sage: l = fileint.readlines()
        sage: fileint.close()
        sage: filter(lambda a: len(a)>2, l[28:])
        ['#include "minc_tides.h"\n',
        'void    mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO)\n',
        '\tint VAR,PAR,TT,i,j, inext;\n',
        '\tVAR = 5;\n',
        '\tPAR = 0;\n',
        '\tTT = 12;\n',
        '\tdouble XX[TT+1][MO+1];\n',
        '\tfor(j=0; j<=TT; j++)\n',
        '\t\tfor(i=0; i<=ORDER; i++)\n',
        '\t\t\tXX[j][i] = 0.e0;\n',
        '\tXX[0][0] = t;\n',
        '\tXX[0][1] = 1.e0;\n',
        '\tfor(i=1;i<=VAR;i++) {\n',
        '\t\tXX[i][0] = v[i-1];\n',
        '\t}\n',
        '\tfor(i=0;i<ORDER;i++) {\n',
        '\t\tXX[5][i] = mul_mc(XX[1],XX[1],i);\n',
        '\t\tXX[6][i] = mul_mc(XX[2],XX[2],i);\n',
        '\t\tXX[7][i] = XX[6][i] + XX[5][i];\n',
        '\t\tXX[8][i] = pow_mc_c(XX[7],-1.50000000000000,XX[8], i);\n',
        '\t\tXX[9][i] = mul_mc(XX[1],XX[8],i);\n',
        '\t\tXX[10][i] = -1.00000000000000*XX[9][i];\n',
        '\t\tXX[11][i] = mul_mc(XX[2],XX[8],i);\n',
        '\t\tXX[12][i] = -1.00000000000000*XX[11][i];\n',
        '\t\tXX[1][i+1] = XX[3][i] / (i+1.0);\n',
        '\t\tXX[2][i+1] = XX[4][i] / (i+1.0);\n',
        '\t\tXX[3][i+1] = XX[10][i] / (i+1.0);\n',
        '\t\tXX[4][i+1] = XX[12][i] / (i+1.0);\n',
        '\t}\n',
        '\tfor(j=0; j<=VAR; j++)\n',
        '\t\tfor(i=0; i<=ORDER; i++)\n',
        '\t\t\tXVAR[i][j] = XX[j][i];\n']
        sage: filedr = open(drfile)
        sage: l = filedr.readlines()
        sage: filedr.close()
        sage: filter(lambda a: len(a)>2, l[28:])
        ['#include "minc_tides.h"\n',
        'int main() {\n',
        '    int  i, VARS, PARS; \n',
        '\tVARS = 4 ;\n',
        '\tPARS = 1;\n',
        '\tdouble tolrel, tolabs, tini, tend, dt; \n',
        '\tdouble v[VARS], p[PARS]; \n',
        '\tv[0] = 1 ; \n',
        '\tv[1] = 0 ; \n',
        '\tv[2] = 0 ; \n',
        '\tv[3] = 0.200000000000000 ; \n',
        '\ttini = 0 ;\n',
        '\ttend = 10 ;\n',
        '\tdt   = 0.100000000000000 ;\n',
        '\ttolrel = 1e-16 ;\n',
        '\ttolabs = 1e-16 ;\n',
        '\textern char ofname[500];\tstrcpy(ofname, "out");\n',
        '\tminc_tides(v,VARS,p,PARS,tini,tend,dt,tolrel,tolabs);\n',
        '\treturn 0; \n']
        sage: shutil.rmtree(tempdir)



    """
    from sage.misc.misc import SAGE_ROOT
    RR = RealField()

    l1, l2 = subexpressions_list(f)

    remove_repeated(l1, l2)
    remove_constants(l1, l2)
    #generate the corresponding c lines

    l3=[]
    var = f[0].arguments()
    for i in l2:
        oper = i[0]
        if oper in ["log", "exp", "sin", "cos"]:
            a = i[1]
            if a in var:
                l3.append((oper, 'XX[{}]'.format(var.index(a))))
            elif a in l1:
                l3.append((oper, 'XX[{}]'.format(l1.index(a)+len(var))))

        else:
            a=i[1]
            b=i[2]
            consta=False
            constb=False

            if a in var:
                aa = 'XX[{}]'.format(var.index(a))
            elif a in l1:
                aa = 'XX[{}]'.format(l1.index(a)+len(var))
            else:
                consta=True
                aa = str(a)
            if b in var:
                bb = 'XX[{}]'.format(var.index(b))
            elif b in l1:
                bb = 'XX[{}]'.format(l1.index(b)+len(var))
            else:
                constb=True
                bb = str(b)
            if consta:
                oper += '_c'
                if not oper=='div':
                    bb, aa = aa, bb
            elif constb:
                oper += '_c'
            l3.append((oper, aa, bb))



    n = len(var)
    res = []
    for i in range(len(l3)):
        el = l3[i]
        string = "XX[{}][i] = ".format(i + n)
        if el[0] == 'add':
            string += el[1] + "[i] + " + el[2] +"[i];"
        elif el[0] == 'add_c':
            string += "(i==0)? {}+".format(N(el[2])) + el[1] + "[0] : "+ el[1]+ "[i];"
        elif el[0] == 'mul':
            string += "mul_mc("+el[1]+","+el[2]+",i);"
        elif el[0] == 'mul_c':
            string += str(N(el[2])) + "*"+ el[1] + "[i];"
        elif el[0] == 'pow_c':
            string += "pow_mc_c("+el[1]+","+str(N(el[2]))+",XX[{}], i);".format(i+n)
        elif el[0] == 'div':
            string += "div_mc("+el[2]+","+el[1]+",XX[{}], i);".format(i+n)
        elif el[0] == 'div_c':
            string += "inv_mc("+str(N(el[2]))+","+el[1]+",XX[{}], i);".format(i+n)
        elif el[0] == 'log':
            string += "log_mc("+el[1]+",XX[{}], i);".format(i+n)
        elif el[0] == 'exp':
            string += "exp_mc("+el[1]+",XX[{}], i);".format(i+n)
        elif el[0] == 'sin':
            string += "sin_mc("+el[1]+",XX[{}], i);".format(i+n+1)
        elif el[0] == 'cos':
            string += "cos_mc("+el[1]+",XX[{}], i);".format(i+n-1)


        res.append(string)

    l1 = list(var)+l1
    indices = [l1.index(i(*var))+n for i in f]
    for i in range (1, n):
        res.append("XX[{}][i+1] = XX[{}][i] / (i+1.0);".format(i,indices[i-1]-n))


    code = res


    outfile = open(integrator, 'a')
    outfile.write('/****************************************************************************\n')
    outfile.write('\tThis file has been created by SageTIDES (1.00)\n')
    outfile.write('\n')
    outfile.write('\tCopyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Marco, M. Rodriguez\n')
    outfile.write('\tGrupo de Mecanica Espacial\n')
    outfile.write('\tUniversity of Zaragoza\n')
    outfile.write('\tSPAIN\n')
    outfile.write('\n')
    outfile.write('\thttp://gme.unizar.es/software/tides\n')
    outfile.write('\tContact: <tides@unizar.es>\n')
    outfile.write('\n')
    outfile.write('\tThis file is part of TIDES.\n')
    outfile.write('\n')
    outfile.write('\tTIDES is free software: you can redistribute it and/or modify\n')
    outfile.write('\tit under the terms of the GNU General Public License as published by\n')
    outfile.write('\tthe Free Software Foundation, either version 3 of the License, or\n')
    outfile.write('\t(at your option) any later version.\n')
    outfile.write('\n')
    outfile.write('\tTIDES is distributed in the hope that it will be useful,\n')
    outfile.write('\tbut WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    outfile.write('\tMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    outfile.write('\tGNU General Public License for more details.\n')
    outfile.write('\n')
    outfile.write('\tYou should have received a copy of the GNU General Public License\n')
    outfile.write('\talong with TIDES.  If not, see <http://www.gnu.org/licenses/>.\n')
    outfile.write('\n')
    outfile.write('*****************************************************************************/\n')
    outfile.write('\n')
    outfile.write('\n')
    outfile.write('\n')
    outfile.write('#include "minc_tides.h"\n')
    outfile.write('\n')
    outfile.write('void    mincseries(double t,double *v, double *p, double **XVAR,int ORDER, int MO)\n')
    outfile.write('{\n')
    outfile.write('\tint VAR,PAR,TT,i,j, inext;\n')
    outfile.write('\n')
    outfile.write("\tVAR = {};\n".format(n))
    outfile.write("\tPAR = {};\n".format(0))
    outfile.write("\tTT = {};\n".format(len(res)))

    outfile.write('\tdouble XX[TT+1][MO+1];\n')
    outfile.write('\n')
    outfile.write('\tfor(j=0; j<=TT; j++)\n')
    outfile.write('\t\tfor(i=0; i<=ORDER; i++)\n')
    outfile.write('\t\t\tXX[j][i] = 0.e0;\n')
    outfile.write('\tXX[0][0] = t;\n')
    outfile.write('\tXX[0][1] = 1.e0;\n')
    outfile.write('\tfor(i=1;i<=VAR;i++) {\n')
    outfile.write('\t\tXX[i][0] = v[i-1];\n')
    outfile.write('\t}\n')
    outfile.write('\n')
    outfile.write('\tfor(i=0;i<ORDER;i++) {\n')
    outfile.write('\n')

    outfile.writelines(["\t\t"+i+"\n" for i in code])

    outfile.write('\t}\n')
    outfile.write('\n')
    outfile.write('\tfor(j=0; j<=VAR; j++)\n')
    outfile.write('\t\tfor(i=0; i<=ORDER; i++)\n')
    outfile.write('\t\t\tXVAR[i][j] = XX[j][i];\n')
    outfile.write('}\n')
    outfile.write('\n')


    outfile = open(driver, 'a')

    outfile.write('/****************************************************************************\n')
    outfile.write('    Driver file of the minc_tides program\n')
    outfile.write('    This file has been created by SageTIDES\n')
    outfile.write('\n')
    outfile.write('    Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez\n')
    outfile.write('    Grupo de Mecanica Espacial\n')
    outfile.write('    University of Zaragoza\n')
    outfile.write('    SPAIN\n')
    outfile.write('\n')
    outfile.write('    http://gme.unizar.es/software/tides\n')
    outfile.write('    Contact: <tides@unizar.es>\n')
    outfile.write('\n')
    outfile.write('    This file is part of TIDES.\n')
    outfile.write('\n')
    outfile.write('    TIDES is free software: you can redistribute it and/or modify\n')
    outfile.write('    it under the terms of the GNU General Public License as published by\n')
    outfile.write('    the Free Software Foundation, either version 3 of the License, or\n')
    outfile.write('    (at your option) any later version.\n')
    outfile.write('\n')
    outfile.write('    TIDES is distributed in the hope that it will be useful,\n')
    outfile.write('    but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    outfile.write('    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    outfile.write('    GNU General Public License for more details.\n')
    outfile.write('\n')
    outfile.write('    You should have received a copy of the GNU General Public License\n')
    outfile.write('    along with TIDES.  If not, see <http://www.gnu.org/licenses/>.\n')
    outfile.write('\n')
    outfile.write('*****************************************************************************/\n')
    outfile.write('\n')
    outfile.write('#include "minc_tides.h"\n')
    outfile.write('\n')
    outfile.write('int main() {\n')
    outfile.write('\n')
    outfile.write('    int  i, VARS, PARS; ')


    outfile.write('\n\tVARS = {} ;\n'.format(n-1))
    outfile.write('\tPARS = 1;\n')
    outfile.write('\tdouble tolrel, tolabs, tini, tend, dt; \n')
    outfile.write('\tdouble v[VARS], p[PARS]; \n')
    for i in range(len(ics)):
        outfile.write('\tv[{}] = {} ; \n'.format(i, ics[i]))
    outfile.write('\ttini = {} ;\n'.format(initial))
    outfile.write('\ttend = {} ;\n'.format(final))
    outfile.write('\tdt   = {} ;\n'.format(delta))
    outfile.write('\ttolrel = {} ;\n'.format(tolrel))
    outfile.write('\ttolabs = {} ;\n'.format(tolabs))
    outfile.write('\textern char ofname[500];')
    outfile.write('\tstrcpy(ofname, "'+ output +'");\n')
    outfile.write('\tminc_tides(v,VARS,p,PARS,tini,tend,dt,tolrel,tolabs);\n')
    outfile.write('\treturn 0; \n }')
    outfile.close()

def genfiles_mpfr(integrator, driver, f, ics, initial, final, delta,
                  parameters=[], parameter_values =[], dig = 20, tolrel=1e-16,
                  tolabs=1e-16, output = ''):
    r"""
        Generate the needed files for the mpfr module of the tides library.

    INPUT:

    - ``integrator`` -- the name of the integrator file.

    - ``driver`` -- the name of the driver file.

    - ``f`` -- the function that determines the differential equation.

    - ``ics`` -- a list or tuple with the initial conditions.

    - ``initial`` -- the initial time for the integration.

    - ``final`` -- the final time for the integration.

    - ``delta`` -- the step of the output.

    - ``parameters`` -- the variables inside the function that should be treated
    as parameters.

    - ``parameter_values`` -- the values of the parameters for the particular IVP

    - ``dig`` -- the number of digits of precission that will be used in the integration

    - ``tolrel`` -- the relative tolerance.

    - ``tolabs`` -- the absolute tolerance.

    -  ``output`` -- the name of the file that the compiled integrator will write to

    This function creates two files, integrator and driver, that can be used
    later with the tides library ([TI]_).


    TESTS::

        sage: from tempfile import mkdtemp
        sage: from sage.interfaces.tides import genfiles_mpfr
        sage: import os
        sage: import shutil
        sage: var('t,x,y,X,Y')
        (t, x, y, X, Y)
        sage: f(t,x,y,X,Y)=[X, Y, -x/(x^2+y^2)^(3/2), -y/(x^2+y^2)^(3/2)]
        sage: tempdir = mkdtemp()
        sage: drfile = tempdir + '/driver.c'
        sage: intfile = tempdir + '/integrator.c'
        sage: genfiles_mpfr(intfile, drfile, f, [1,0, 0, 0.2], 0, 10, 0.1, output = 'out', dig = 50)
        sage: fileint = open(intfile)
        sage: l = fileint.readlines()
        sage: fileint.close()
        sage: l[56]
        '\t\tmpfrts_add_t(itd, link[1], link[0], link[2], i);\n'
        sage: filedr = open(drfile)
        sage: l = filedr.readlines()
        sage: filedr.close()
        sage: l[45]
        '\tfor(i=0; i<nvar; i++) mpfr_init2(v[i], TIDES_PREC);\n'
        sage: l[56]
        '\tmpfr_init2(tini, TIDES_PREC); \n'
        sage: shutil.rmtree(tempdir)


    """
    from sage.misc.misc import SAGE_ROOT
    RR = RealField()
    l1, l2 = subexpressions_list(f, parameters)
    remove_repeated(l1, l2)
    remove_constants(l1, l2)
    l3=[]
    var = f[0].arguments()
    for i in l2:
        oper = i[0]
        if oper in ["log", "exp", "sin", "cos", "atan", "asin", "acos"]:
            a = i[1]
            if a in var:
                l3.append((oper, 'var[{}]'.format(var.index(a))))
            elif a in parameters:
                l3.append((oper, 'par[{}]'.format(parameters.index(a))))
            else:
                l3.append((oper, 'link[{}]'.format(l1.index(a))))

        else:
            a=i[1]
            b=i[2]
            consta=False
            constb=False

            if a in var:
                aa = 'var[{}]'.format(var.index(a))
            elif a in l1:
                aa = 'link[{}]'.format(l1.index(a))
            elif a in parameters:
                aa = 'par[{}]'.format(parameters.index(a))
            else:
                consta=True
                aa = str(a)
            if b in var:
                bb = 'var[{}]'.format(var.index(b))
            elif b in l1:
                bb = 'link[{}]'.format(l1.index(b))
            elif b in parameters:
                bb = 'par[{}]'.format(parameters.index(b))
            else:
                constb=True
                bb = str(b)
            if consta:
                oper += '_c'
                if not oper=='div':
                    bb, aa = aa,bb
            elif constb:
                oper += '_c'
            l3.append((oper, aa, bb))


    n = len(var)
    code = []


    l1 = list(var)+l1
    indices = [l1.index(i(*var))+n for i in f]
    for i in range (1, n):
        aux = indices[i-1]-n
        if aux < n:
            code.append('mpfrts_var_t(itd, var[{}], var[{}], i);'.format(aux, i))
        else:
            code.append('mpfrts_var_t(itd, link[{}], var[{}], i);'.format(aux-n, i))

    for i in range(len(l3)):
        el = l3[i]
        string = "mpfrts_"
        if el[0] == 'add':
            string += 'add_t(itd, ' + el[1] + ', ' + el[2] + ', link[{}], i);'.format(i)
        elif el[0] == 'add_c':
            string += 'add_t_c(itd, "' + str(N(el[2], digits=dig)) + '", ' + el[1] + ', link[{}], i);'.format(i)
        elif el[0] == 'mul':
            string += 'mul_t(itd, ' + el[1] + ', ' + el[2] + ', link[{}], i);'.format(i)
        elif el[0] == 'mul_c':
            string += 'mul_t_c(itd, "' + str(N(el[2], digits=dig)) + '", ' + el[1] + ', link[{}], i);'.format(i)
        elif el[0] == 'pow_c':
            string += 'pow_t_c(itd, ' + el[1] + ', "' + str(N(el[2],digits=dig)) + '", link[{}], i);'.format(i)
        elif el[0] == 'div':
            string += 'div_t(itd, ' + el[2] + ', ' + el[1] + ', link[{}], i);'.format(i)
        elif el[0] == 'div_c':
            string += 'div_t_cv(itd, "' + str(N(el[2], digits=dig)) + '", ' + el[1] + ', link[{}], i);'.format(i)
        elif el[0] == 'log':
            string += 'log_t(itd, ' + el[1]  + ', link[{}], i);'.format(i)
        elif el[0] == 'exp':
            string += 'exp_t(itd, ' + el[1]  + ', link[{}], i);'.format(i)
        elif el[0] == 'sin':
            string += 'sin_t(itd, ' + el[1]  + ', link[{}], link[{}], i);'.format(i+1, i)
        elif el[0] == 'cos':
            string += 'cos_t(itd, ' + el[1]  + ', link[{}], link[{}], i);'.format(i-1, i)
        elif el[0] == 'atan':
            indarg = l1.index(1+l2[i][1]**2)-n
            string += 'atan_t(itd, ' + el[1] + ', link[{}], link[{}], i);'.format(indarg, i)
        elif el[0] == 'asin':
            indarg = l1.index(sqrt(1-l2[i][1]**2))-n
            string += 'asin_t(itd, ' + el[1] + ', link[{}], link[{}], i);'.format(indarg, i)
        elif el[0] == 'acos':
            indarg = l1.index(-sqrt(1-l2[i][1]**2))-n
            string += 'acos_t(itd, ' + el[1] + ', link[{}], link[{}], i);'.format(indarg, i)
        code.append(string)

    VAR = n-1
    PAR = len(parameters)
    TT =  len(code)+1-VAR

    outfile = open(integrator, 'a')


    outfile.write('/****************************************************************************\n')
    outfile.write('\tThis file has been created by SageTIDES (1.0)\n')
    outfile.write('\n')
    outfile.write('\tCopyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez\n')
    outfile.write('\tGrupo de Mecanica Espacial\n')
    outfile.write('\tUniversity of Zaragoza\n')
    outfile.write('\tSPAIN\n')
    outfile.write('\n')
    outfile.write('\thttp://gme.unizar.es/software/tides\n')
    outfile.write('\tContact: <tides@unizar.es>\n')
    outfile.write('\n')
    outfile.write('\tThis file is part of TIDES.\n')
    outfile.write('\n')
    outfile.write('\tTIDES is free software: you can redistribute it and/or modify\n')
    outfile.write('\tit under the terms of the GNU General Public License as published by\n')
    outfile.write('\tthe Free Software Foundation, either version 3 of the License, or\n')
    outfile.write('\t(at your option) any later version.\n')
    outfile.write('\n')
    outfile.write('\tTIDES is distributed in the hope that it will be useful,\n')
    outfile.write('\tbut WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    outfile.write('\tMERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    outfile.write('\tGNU General Public License for more details.\n')
    outfile.write('\n')
    outfile.write('\tYou should have received a copy of the GNU General Public License\n')
    outfile.write('\talong with TIDES.  If not, see <http://www.gnu.org/licenses/>.\n')
    outfile.write('\n')
    outfile.write('*****************************************************************************/\n')
    outfile.write('\n')
    outfile.write('#include "mp_tides.h"\n')
    outfile.write('\n')
    outfile.write('\n')
    outfile.write('\n')
    outfile.write('long  function_iteration(iteration_data *itd, mpfr_t t, mpfr_t v[], mpfr_t p[], int ORDER, mpfr_t *cvfd)\n')
    outfile.write('{\n')
    outfile.write('\n')
    outfile.write('\tint i;\n')
    outfile.write('    int NCONST = 0;\n')
    outfile.write('    mpfr_t ct[0];\n')
    outfile.write('\n')
    outfile.write('\n')


    outfile.write("\n\tstatic int VARIABLES = {};\n".format(VAR))
    outfile.write("\tstatic int PARAMETERS = {};\n".format(PAR))
    outfile.write("\tstatic int LINKS = {};\n".format(TT))
    outfile.write('\tstatic int   FUNCTIONS        = 0;\n')
    outfile.write('\tstatic int   POS_FUNCTIONS[1] = {0};\n')
    outfile.write('\n\tinitialize_mp_case();\n')
    outfile.write('\n\tfor(i=0;  i<=ORDER; i++) {\n')
    for i in code:
        outfile.write('\t\t'+i+'\n')

    outfile.write('\t}\n\twrite_mp_solution();\n\n')
    outfile.write('\tclear_vpl();\n\tclear_cts();\n')
    outfile.write('\treturn NUM_COLUMNS;\n}')
    outfile.close()


    npar = len(parameter_values)
    outfile = open(driver, 'a')

    outfile.write('/****************************************************************************\n')
    outfile.write('Driver file of the mp_tides program\n')
    outfile.write('This file has been created by SageTIDES (1.0)\n')
    outfile.write('\n')
    outfile.write('    Copyright (C) 2010 A. Abad, R. Barrio, F. Blesa, M. Rodriguez\n')
    outfile.write('    Grupo de Mecanica Espacial\n')
    outfile.write('    University of Zaragoza\n')
    outfile.write('    SPAIN\n')
    outfile.write('\n')
    outfile.write('    http://gme.unizar.es/software/tides\n')
    outfile.write('    Contact: <tides@unizar.es>\n')
    outfile.write('\n')
    outfile.write('    This file is part of TIDES.\n')
    outfile.write('\n')
    outfile.write('    TIDES is free software: you can redistribute it and/or modify\n')
    outfile.write('    it under the terms of the GNU General Public License as published by\n')
    outfile.write('    the Free Software Foundation, either version 3 of the License, or\n')
    outfile.write('    (at your option) any later version.\n')
    outfile.write('\n')
    outfile.write('    TIDES is distributed in the hope that it will be useful,\n')
    outfile.write('    but WITHOUT ANY WARRANTY; without even the implied warranty of\n')
    outfile.write('    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the\n')
    outfile.write('    GNU General Public License for more details.\n')
    outfile.write('\n')
    outfile.write('    You should have received a copy of the GNU General Public License\n')
    outfile.write('    along with TIDES.  If not, see <http://www.gnu.org/licenses/>.\n')
    outfile.write('\n')
    outfile.write('    *****************************************************************************/\n')
    outfile.write('\n')
    outfile.write('    #include "mpfr.h"\n')
    outfile.write('    #include "mp_tides.h"\n')
    outfile.write('    long  function_iteration(iteration_data *itd, mpfr_t t, mpfr_t v[], mpfr_t p[], int ORDER, mpfr_t *cvfd);\n')
    outfile.write('\n')
    outfile.write('    int main() {\n')
    outfile.write('\n')
    outfile.write('        int i;\n')
    outfile.write('\n')
    outfile.write('\n')

    outfile.write('\tint nfun = 0;\n')
    outfile.write('\tset_precision_digits({});'.format(dig))
    outfile.write('\n\tint npar = {};\n'.format(npar))
    outfile.write('\tmpfr_t p[npar];\n')
    outfile.write('\tfor(i=0; i<npar; i++) mpfr_init2(p[i], TIDES_PREC);\n')

    for i in range(npar):
        outfile.write('\tmpfr_set_str(p[{}], "{}", 10, TIDES_RND);\n'.format(i,N(parameter_values[i],digits=dig)))
    outfile.write('\tint nvar = {};\n\tmpfr_t v[nvar];\n'.format(VAR))
    outfile.write('\tfor(i=0; i<nvar; i++) mpfr_init2(v[i], TIDES_PREC);\n')
    for i in range(len(ics)):
        outfile.write('\tmpfr_set_str(v[{}], "{}", 10, TIDES_RND);\n'.format(i,N(ics[i],digits=dig)))
    outfile.write('\tmpfr_t tolrel, tolabs;\n')
    outfile.write('\tmpfr_init2(tolrel, TIDES_PREC); \n')
    outfile.write('\tmpfr_init2(tolabs, TIDES_PREC); \n')
    outfile.write('\tmpfr_set_str(tolrel, "{}", 10, TIDES_RND);\n'.format(N(tolrel,digits=dig)))
    outfile.write('\tmpfr_set_str(tolabs, "{}", 10, TIDES_RND);\n'.format(N(tolabs,digits=dig)))

    outfile.write('\tmpfr_t tini, dt; \n')
    outfile.write('\tmpfr_init2(tini, TIDES_PREC); \n')
    outfile.write('\tmpfr_init2(dt, TIDES_PREC); \n')


    outfile.write('\tmpfr_set_str(tini, "{}", 10, TIDES_RND);;\n'.format(N(initial,digits=dig)))
    outfile.write('\tmpfr_set_str(dt, "{}", 10, TIDES_RND);\n'.format(N(delta,digits=dig)))
    outfile.write('\tint nipt = {};\n'.format(floor((final-initial)/delta)))
    outfile.write('\tFILE* fd = fopen("' + output + '", "w");\n')
    outfile.write('\tmp_tides_delta(function_iteration, NULL, nvar, npar, nfun, v, p, tini, dt, nipt, tolrel, tolabs, NULL, fd);\n')
    outfile.write('\tfclose(fd);\n\treturn 0;\n}')
    outfile.close()
