r"""
Generic interface to LattE integrale programs
"""
#*****************************************************************************
#       Copyright (C) 2017 Vincent Delecroix <vincent.delecroix@gmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def count(arg, ehrhart_polynomial=False, multivariate_generating_function=False, raw_output=False, verbose=False, **kwds):
    r"""
    Call to the program count from LattE integrale

    INPUT:

    - ``arg`` -- a cdd or LattE description string

    - ``ehrhart_polynomial``, ``multivariate_generating_function``  -- to
      compute Ehrhart polynomial or multivariate generating function instead of
      just counting points

    - ``raw_output`` -- if ``True`` then return directly the output string from LattE

    - For all other options of the count program, consult the LattE manual

    OUTPUT:

    Either a string (if ``raw_output`` if set to ``True``) or an integer (when
    counting points), or a polynomial (if ``ehrhart_polynomial`` is set to
    ``True``) or a multivariate THING (if ``multivariate_generating_function``
    is set to ``True``)

    EXAMPLES::

        sage: from sage.interfaces.latte import count
        sage: P = 2 * polytopes.cube()

    Counting integer points from either the H or V representation::

        sage: count(P.cdd_Hrepresentation(), cdd=True)   # optional - latte_int
        125
        sage: count(P.cdd_Vrepresentation(), cdd=True)   # optional - latte_int
        125

    Ehrhart polynomial::

        sage: count(P.cdd_Hrepresentation(), cdd=True, ehrhart_polynomial=True)  # optional - latte_int
        64*t^3 + 48*t^2 + 12*t + 1

    Multivariate generating function currently only work with ``raw_output=True``::

        sage: opts = {'cdd': True,
        ....:         'multivariate_generating_function': True,
        ....:         'raw_output': True}
        sage: cddin = P.cdd_Hrepresentation()
        sage: print(count(cddin, **opts))  # optional - latte_int
        x[0]^2*x[1]^(-2)*x[2]^(-2)/((1-x[1])*(1-x[2])*(1-x[0]^(-1)))
         + x[0]^(-2)*x[1]^(-2)*x[2]^(-2)/((1-x[1])*(1-x[2])*(1-x[0]))
         + x[0]^2*x[1]^(-2)*x[2]^2/((1-x[1])*(1-x[0]^(-1))*(1-x[2]^(-1)))
         + x[0]^(-2)*x[1]^(-2)*x[2]^2/((1-x[1])*(1-x[0])*(1-x[2]^(-1)))
         + x[0]^2*x[1]^2*x[2]^(-2)/((1-x[2])*(1-x[0]^(-1))*(1-x[1]^(-1)))
         + x[0]^(-2)*x[1]^2*x[2]^(-2)/((1-x[2])*(1-x[0])*(1-x[1]^(-1)))
         + x[0]^2*x[1]^2*x[2]^2/((1-x[0]^(-1))*(1-x[1]^(-1))*(1-x[2]^(-1)))
         + x[0]^(-2)*x[1]^2*x[2]^2/((1-x[0])*(1-x[1]^(-1))*(1-x[2]^(-1)))

    TESTS:

    Testing raw output::

        sage: from sage.interfaces.latte import count
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: count(cddin, cdd=True, raw_output=True)  # optional - latte_int
        '19'
        sage: count(cddin, cdd=True, raw_output=True, ehrhart_polynomial=True) # optional - latte_int
        ' + 1 * t^0 + 10/3 * t^1 + 8 * t^2 + 20/3 * t^3'
        sage: count(cddin, cdd=True, raw_output=True, multivariate_generating_function=True) # optional - latte_int
        'x[0]^(-1)*x[1]^(-1)/((1-x[0]*x[2])*(1-x[0]^(-1)*x[1])*...x[0]^(-1)*x[2]^(-1)))\n'

    Testing the ``verbose`` option::

        sage: n = count(cddin, cdd=True, verbose=True, raw_output=True)  # optional - latte_int
        This is LattE integrale ...
        ...
        Invocation: count '--redundancy-check=none' --cdd /dev/stdin
        ...
        Total Unimodular Cones: ...
        Maximum number of simplicial cones in memory at once: ...
        <BLANKLINE>
        ****  The number of lattice points is:   ****
        Total time: ... sec
    """
    from subprocess import Popen, PIPE
    from sage.misc.misc import SAGE_TMP
    from sage.rings.integer import Integer

    args = ['count']
    if ehrhart_polynomial and multivariate_generating_function:
        raise ValueError
    if ehrhart_polynomial:
        args.append('--ehrhart-polynomial')
    elif multivariate_generating_function:
        args.append('--multivariate-generating-function')

    if 'redundancy_check' not in kwds:
        args.append('--redundancy-check=none')

    for key,value in kwds.items():
        if value is None or value is False:
            continue

        key = key.replace('_','-')
        if value is True:
            args.append('--{}'.format(key))
        else:
            args.append('--{}={}'.format(key, value))

    if multivariate_generating_function:
        from sage.misc.temporary_file import tmp_filename
        filename = tmp_filename()
        with open(filename, 'w') as f:
            f.write(arg)
        args += [filename]
    else:
        args += ['/dev/stdin']

    try:
        # The cwd argument is needed because latte
        # always produces diagnostic output files.
        latte_proc = Popen(args,
                           stdin=PIPE, stdout=PIPE,
                           stderr=(None if verbose else PIPE),
                           cwd=str(SAGE_TMP))
    except OSError:
        from sage.misc.package import PackageNotFoundError
        raise PackageNotFoundError('latte_int')

    ans, err = latte_proc.communicate(arg)
    ret_code = latte_proc.poll()
    if ret_code:
        if err is None:
            err = ", see error message above"
        else:
            err = ":\n" + err
        raise RuntimeError("LattE integrale program failed (exit code {})".format(ret_code) + err.strip())


    if ehrhart_polynomial:
        ans = ans.splitlines()[-2]
        if raw_output:
            return ans
        else:
            from sage.rings.polynomial.polynomial_ring_constructor import PolynomialRing
            from sage.rings.rational_field import QQ
            R = PolynomialRing(QQ, 't')
            return R(ans)
    elif multivariate_generating_function:
        with open(filename + '.rat') as f:
            ans = f.read()
        if raw_output:
            return ans
        else:
            raise NotImplementedError("there is no Sage object to handle multivariate series from LattE, use raw_output=True")
    else:
        ans = ans.splitlines()[-1]
        if not ans:
            # opening a file is slow (30e-6s), so we read the file
            # numOfLatticePoints only in case of a IndexError above
            with open(SAGE_TMP+'/numOfLatticePoints', 'r') as f:
                ans = f.read()

        if raw_output:
            return ans
        else:
            from sage.rings.integer import Integer
            return Integer(ans)

def integrate(arg, polynomial, raw_output=False, verbose=False, **kwds):
    r"""
    Call to the function integrate from LattE integrale, to compute the integral of a polynomial over a polytope.

    INPUT:

    - ``arg`` -- a cdd or LattE description string

    - ``polynomial`` -- the polynomial to be integrated, as

    - ``raw_output`` -- if ``True`` then return directly the output string from LattE

    - ``verbose`` -- if ``True`` then return directly verbose output from LattE

    - For all other options of the count program, consult the LattE manual

    OUTPUT:

    Either a string (if ``raw_output`` if set to ``True``) or a rational.

    EXAMPLES::

        sage: from sage.interfaces.latte import integrate
        sage: P = 2 * polytopes.cube()
        sage: x, y, z = polygen(QQ, 'x, y, z')

    Integrating over a polynomial over a polytope in either the H or V representation::

        sage: integrate(P.cdd_Hrepresentation(), x^2*y^2*z^2, cdd=True)   # optional - latte_int
        4096/27
        sage: integrate(P.cdd_Vrepresentation(), x^2*y^2*z^2, cdd=True)   # optional - latte_int
        4096/27

    TESTS:

    Testing raw output::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: x, y, z = polygen(QQ, 'x, y, z')
        sage: f = 3*x^2*y^4*z^6 + 7*y^3*z^5
        sage: integrate(cddin, f, cdd=True, raw_output=True)  # optional - latte_int
        '629/47775'

    Testing the ``verbose`` option::

        sage: ans = integrate(cddin, f, cdd=True, verbose=True, raw_output=True)  # optional - latte_int
        This is LattE integrale ...
        ...
        Invocation: integrate --valuation=integrate --triangulate --redundancy-check=none --cdd --monomials=... /dev/stdin
        ...
    """
    from subprocess import Popen, PIPE
    from sage.misc.misc import SAGE_TMP
    from sage.rings.integer import Integer

    args = ['integrate']
    args.append('--valuation=integrate')
    args.append('--triangulate')

    if 'redundancy_check' not in kwds:
        args.append('--redundancy-check=none')

    for key,value in kwds.items():
        if value is None or value is False:
            continue

        key = key.replace('_','-')
        if value is True:
            args.append('--{}'.format(key))
        else:
            args.append('--{}={}'.format(key, value))

    # transform polynomial to LattE description
    coefficients_list = polynomial.coefficients()
    exponents_list = [list(exponent_vector_i) for exponent_vector_i in polynomial.exponents()]

    # assuming that the order in coefficients() and exponents() methods match
    monomials_list = zip(coefficients_list, exponents_list)
    monomials_list = [list(monomial_i) for monomial_i in monomials_list]
    monomials_list = str(monomials_list)

    # write the polynomial description to a temp file
    from sage.misc.temporary_file import tmp_filename
    filename_polynomial = tmp_filename()

    with open(filename_polynomial, 'w') as f:
        f.write(monomials_list)
    args += ['--monomials=' + filename_polynomial]

    args += ['/dev/stdin']

    try:
        # The cwd argument is needed because latte
        # always produces diagnostic output files.
        latte_proc = Popen(args,
                           stdin=PIPE, stdout=PIPE,
                           stderr=(None if verbose else PIPE),
                           cwd=str(SAGE_TMP))
    except OSError:
        from sage.misc.package import PackageNotFoundError
        raise PackageNotFoundError('latte_int')

    ans, err = latte_proc.communicate(arg)
    ret_code = latte_proc.poll()
    if ret_code:
        if err is None:
            err = ", see error message above"
        else:
            err = ":\n" + err
        raise RuntimeError("LattE integrale program failed (exit code {})".format(ret_code) + err.strip())

    ans = ans.splitlines()
    ans = ans[-5].split()
    assert(ans[0]=='Answer:')
    ans = ans[1]

    if raw_output:
        return ans
    else:
        from sage.rings.rational import Rational
        return Rational(ans)

def volume(arg, algorithm='triangulate', raw_output=False, verbose=False, **kwds):
    r"""
    Call to the function integrate from LattE integrale, to compute the volume of a polytope.

    INPUT:

    - ``arg`` -- a cdd or LattE description string

    - ``algorithm`` -- (default: 'triangulate') the integration method. Use 'triangulate' for polytope triangulation
    or 'cone-decompose' for tangent cone decomposition method.

    - ``raw_output`` -- if ``True`` then return directly the output string from LattE

    - ``verbose`` -- if ``True`` then return directly verbose output from LattE

    - For all other options of the count program, consult the LattE manual

    OUTPUT:

    Either a string (if ``raw_output`` if set to ``True``) or a rational.

    EXAMPLES::

        sage: from sage.interfaces.latte import volume
        sage: P = 2 * polytopes.cube()

    Computing the volume of a polytope in either the H or V representation::

        sage: volume(P.cdd_Hrepresentation(), cdd=True)   # optional - latte_int
        64
        sage: volume(P.cdd_Vrepresentation(), cdd=True)   # optional - latte_int
        64

    TESTS:

    Testing triangulate algorithm::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: volume(cddin, algorithm='triangulate', cdd=True)  # optional - latte_int
        20/3

    Testing convex decomposition algorithm::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: volume(cddin, algorithm='cone-decompose', cdd=True)  # optional - latte_int
        20/3

    Testing raw output::

        sage: from sage.interfaces.latte import integrate
        sage: P = polytopes.cuboctahedron()
        sage: cddin = P.cdd_Vrepresentation()
        sage: volume(cddin, cdd=True, raw_output=True)  # optional - latte_int
        '20/3'

    Testing the ``verbose`` option::

        sage: ans = volume(cddin, cdd=True, verbose=True, raw_output=True)  # optional - latte_int
        This is LattE integrale ...
        ...
        Invocation: integrate --valuation=volume --triangulate --redundancy-check=none --cdd /dev/stdin
        ...
    """
    from subprocess import Popen, PIPE
    from sage.misc.misc import SAGE_TMP
    from sage.rings.integer import Integer

    args = ['integrate']
    args.append('--valuation=volume')

    if algorithm=='triangulate':
        args.append('--triangulate')
    elif algorithm=='cone-decompose':
        args.append('--cone-decompose')

    if 'redundancy_check' not in kwds:
        args.append('--redundancy-check=none')

    for key,value in kwds.items():
        if value is None or value is False:
            continue

        key = key.replace('_','-')
        if value is True:
            args.append('--{}'.format(key))
        else:
            args.append('--{}={}'.format(key, value))

    args += ['/dev/stdin']

    try:
        # The cwd argument is needed because latte
        # always produces diagnostic output files.
        latte_proc = Popen(args,
                           stdin=PIPE, stdout=PIPE,
                           stderr=(None if verbose else PIPE),
                           cwd=str(SAGE_TMP))
    except OSError:
        from sage.misc.package import PackageNotFoundError
        raise PackageNotFoundError('latte_int')

    ans, err = latte_proc.communicate(arg)
    ret_code = latte_proc.poll()
    if ret_code:
        if err is None:
            err = ", see error message above"
        else:
            err = ":\n" + err
        raise RuntimeError("LattE integrale program failed (exit code {})".format(ret_code) + err.strip())

    ans = ans.splitlines()

    ans = ans[-5].split()
    assert(ans[0]=='Answer:')
    ans = ans[1]

    if raw_output:
        return ans
    else:
        from sage.rings.rational import Rational
        return Rational(ans)
