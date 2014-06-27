from  sage.rings.real_mpfr import RealField
import shutil
import os
from sage.calculus.all import symbolic_expression
from sage.misc.flatten import flatten
from sage.ext.fast_callable import fast_callable
from sage.misc.lazy_import import lazy_import
lazy_import('sage.rings.semirings.non_negative_integer_semiring', 'NN')
from sage.misc.functional import N


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
    operation made, and the operands. For the cosines (resp. sines), the
    corresponding sine (resp. cosine) is also added.


    EXAMPLES::

        sage: from sage.calculus.tides.file_generator import subexpressions_list
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
    later with the min_tides library.


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


    shutil.copy(SAGE_ROOT+'/src/sage/calculus/tides/seriesFile00.txt', integrator)
    outfile = open(integrator, 'a')
    outfile.write("\tVAR = {};\n".format(n))
    outfile.write("\tPAR = {};\n".format(0))
    outfile.write("\tTT = {};\n".format(len(res)))
    infile = open(SAGE_ROOT+'/src/sage/calculus/tides/seriesFile01.txt')
    for i in infile:
        outfile.write(i)
    infile.close()
    outfile.writelines(["\t\t"+i+"\n" for i in code])

    infile = open(SAGE_ROOT+'/src/sage/calculus/tides/seriesFile02.txt')
    for i in infile:
        outfile.write(i)
    outfile.close()


    shutil.copy(SAGE_ROOT+'/src/sage/calculus/tides/driverFile00.txt', driver)
    outfile = open(driver, 'a')
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


