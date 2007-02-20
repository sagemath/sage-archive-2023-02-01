from sage.all import *

verbose = False

timeout = 5

def report(F, title):
    systems = ['sage', 'magma']
    print '\n\n'
    print '='*70
    print ' '*10 + title
    print '='*70
    os.system('uname -a')
    print '\n'
    for f in F:
        print "-"*70
        print f.__doc__.strip()
        print ('%15s'*len(systems))%tuple(systems)
        w = []
        for s in systems:
            alarm(timeout)
            try:
                t = f(system=s)
            except KeyboardInterrupt:
                t = -timeout
            alarm(0)
            w.append(t)
        w = tuple(w)
        print ('%15.3f'*len(w))%w


#######################################################################
# Dense Benchmarks over ZZ
#######################################################################

def report_ZZ():
    F = [rank_ZZ, rank2_ZZ, nullspace_ZZ, charpoly_ZZ, smithform_ZZ,
         matrix_multiply_ZZ, det_ZZ, det_QQ]
    title = 'Dense benchmarks over ZZ'
    report(F, title)

# Integer Nullspace

def nullspace_ZZ(n=300, min=0, max=10, system='sage'):
    """
    Nullspace over ZZ:
    Given a n+1 x n (with n=300) matrix over ZZ with random entries
    between min=0 and max=10, compute the nullspace.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n+1, n, x=min, y=max+1).change_ring(QQ)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RationalField(), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system


def charpoly_ZZ(n=100, min=0, max=9, system='sage'):
    """
    Characteristic polynomial over ZZ:
    Given a n x n (with n=100) matrix over ZZ with random entries
    between min=0 and max=9, compute the charpoly.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        v = A.charpoly()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
K := CharacteristicPolynomial(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system




def rank_ZZ(n=500, min=0, max=9, system='sage'):
    """
    Rank over ZZ:
    Given a n x (n+10) (with n=500) matrix over ZZ with random entries
    between min=0 and max=9, compute the rank.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n+10, x=min, y=max+1)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(IntegerRing(), n, n+10)![Random(%s,%s) : i in [1..n*(n+10)]];
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system

def rank2_ZZ(n=500, min=0, max=9, system='sage'):
    """
    Rank over ZZ:
    Given a (n + 10) x n (with n=500) matrix over ZZ with random entries
    between min=0 and max=9, compute the rank.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n+10, n, x=min, y=max+1)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(IntegerRing(), n+10, n)![Random(%s,%s) : i in [1..n*(n+10)]];
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system

# Smith Form

def smithform_ZZ(n=128, min=0, max=9, system='sage'):
    """
    Smith Form over ZZ:
    Given a n x n (with n=128) matrix over ZZ with random entries
    between min=0 and max=9, compute the Smith normal form.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        v = A.elementary_divisors()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
K := ElementaryDivisors(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system


def matrix_multiply_ZZ(n=200, min=-9, max=9, system='sage', times=1):
    """
    Matrix multiplication over ZZ
    Given an n x n (with n=200) matrix A over ZZ with random entries
    between min=-9 and max=9, inclusive, compute A * (A+1).
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        B = A + 1
        t = cputime()
        for z in range(times):
            v = A * B
        return cputime(t)/times
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
"""%(n,min,max,times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError, 'unknown system "%s"'%system

def matrix_add_ZZ(n=200, min=-9, max=9, system='sage', times=10):
    """
    Matrix addition over ZZ
    Given an n x n (with n=200) matrix A and B over ZZ with random entries
    between min=-9 and max=9, inclusive, compute A + B.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        B = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        for z in range(times):
            v = A + B
        return cputime(t)/times
    elif system == 'magma':
        code = """
n := %s;
min := %s;
max := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(min,max) : i in [1..n^2]];
B := MatrixAlgebra(IntegerRing(), n)![Random(min,max) : i in [1..n^2]];
t := Cputime();
for z in [1..%s] do
    K := A + B;
end for;
s := Cputime(t);
"""%(n,min,max,times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError, 'unknown system "%s"'%system


def det_ZZ(n=400, min=1, max=100, system='sage'):
    """
    Dense integer determinant over ZZ.
    Given an n x n (with n=400) matrix A over ZZ with random entries
    between min=1 and max=100, inclusive, compute det(A).
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
d := Determinant(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system


def det_QQ(n=300, num_bound=10, den_bound=10, system='sage'):
    """
    Dense rational determinant over QQ.
    Given an n x n (with n=300) matrix A over QQ with random entries
    with numerator and denominator between min=-10 and 10,
    inclusive, compute det(A).
    """
    if system == 'sage':
        A = random_matrix(QQ, n, n, num_bound=num_bound, den_bound=den_bound)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(RationalField(), n)![Random(%s,%s)/Random(1,%s) : i in [1..n^2]];
t := Cputime();
d := Determinant(A);
s := Cputime(t);
"""%(n,-num_bound, num_bound, den_bound)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system

#######################################################################
# Dense Benchmarks over GF(p), for small p.
#######################################################################

def report_GF(p=16411):
    F = [rank_GF, rank2_GF, nullspace_GF, charpoly_GF,
         matrix_multiply_GF, det_GF]
    title = 'Dense benchmarks over GF with prime %i' % p
    report(F, title)

# Nullspace over GF

def nullspace_GF(n=300, p=16411, system='sage'):
    """
    Given a n+1 x n (with n=300) matrix over GF(p) p=16411 with random
    entries, compute the nullspace.
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n+1)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(RMatrixSpace(GF(%s), n, n+1));
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(n,p)
        if verbose: print code
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError, 'unknown system "%s"'%system


# Characteristic Polynomial over GF

def charpoly_GF(n=100, p=16411, system='sage'):
    """
    Given a n x n (with n=100) matrix over GF with random entries,
    compute the charpoly.
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n)
        t = cputime()
        v = A.charpoly()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := CharacteristicPolynomial(A);
s := Cputime(t);
"""%(n,p)
        if verbose: print code
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError, 'unknown system "%s"'%system

def matrix_add_GF(n=1000, p=16411, system='sage',times=100):
    """
    Given two n x n (with n=1000) matrix over GF with random entries,
    add them.
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n)
        B = random_matrix(GF(p), n, n)
        t = cputime()
        for n in range(times):
            v = A + B
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
B := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
for z in [1..%s] do
    K := A + B;
end for;
s := Cputime(t);
"""%(n,p,p,times)
        if verbose: print code
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError, 'unknown system "%s"'%system



# Matrix multiplication over GF(p)

def matrix_multiply_GF(n=100, p=16411, system='sage', times=3):
    """
    Given an n x n (with n=100) matrix A over GF(p) with random
    entries, compute A * (A+1).
    """
    if system == 'sage':
        A = random_matrix(GF(p), n)
        B = A + 1
        t = cputime()
        for n in range(times):
            v = A * B
        return cputime(t) / times
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
"""%(n,p,times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError, 'unknown system "%s"'%system


def rank_GF(n=500, p=16411, min=0, max=9, system='sage'):
    """
    Rank over GF:
    Given a n x (n+10) (with n=500) matrix over ZZ with random entries
    between min=0 and max=9, compute the rank.
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n+10)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(n,p)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system

def rank2_GF(n=500, p=16411, min=0, max=9, system='sage'):
    """
    Rank over GF(p):
    Given a (n + 10) x n (with n=500) matrix over GF(p) with random entries
    between min=0 and max=9, compute the rank.
    """
    if system == 'sage':
        A = random_matrix(GF(p), n+10, n)
        t = cputime()
        v = A.rank()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
K := Rank(A);
s := Cputime(t);
"""%(n,p)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system

def det_GF(n=400, p=16411 ,min=1, max=100, system='sage'):
    """
    Dense integer determinant over GF.
    Given an n x n (with n=400) matrix A over GF with random entries
    between min=1 and max=100, inclusive, compute det(A).
    """
    if system == 'sage':
        A = random_matrix(GF(p), n, n, x=min, y=max+1)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := Random(MatrixAlgebra(GF(%s), n));
t := Cputime();
d := Determinant(A);
s := Cputime(t);
"""%(n,p)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system





#######################################################################
# Dense Benchmarks over QQ
#######################################################################

def hilbert_matrix(n):
  A = Matrix(QQ,n,n)
  for i in range(A.nrows()):
      for j in range(A.ncols()):
          A[i,j] =  QQ(1)/((i+1)+(j+1)-1)
  return A

# Reduced row echelon form over QQ

def echelon_QQ(n=100, min=0, max=9, system='sage'):
    """
    Given a n x (2*n) (with n=100) matrix over QQ with random integer entries
    between min=0 and max=9, compute the reduced row echelon form.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, 2*n, x=min, y=max+1).change_ring(QQ)
        t = cputime()
        v = A.echelon_form()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RationalField(), n, 2*n)![Random(%s,%s) : i in [1..n*2*n]];
t := Cputime();
K := EchelonForm(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system

# Invert a matrix over QQ.

def inverse_QQ(n=100, min=0, max=9, system='sage'):
    """
    Given a n x n (with n=100) matrix over QQ with random integer entries
    between min=0 and max=9, compute the reduced row echelon form.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1).change_ring(QQ)
        t = cputime()
        v = A**(-1)
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(RationalField(), n)![Random(%s,%s) : i in [1..n*n]];
t := Cputime();
K := A^(-1);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system


# Matrix multiplication over QQ
def matrix_multiply_QQ(n=100, bnd=2, system='sage', times=1):
    """
    Given an n x n (with n=100) matrix A over QQ with random entries
    whose numerators and denominators are bounded by b, compute A *
    (A+1).
    """
    if system == 'sage':
        A = random_matrix(QQ, n, n, num_bound=bnd, den_bound=bnd)
        B = A + 1
        t = cputime()
        for z in range(times):
            v = A * B
        return cputime(t)/times
    elif system == 'magma':
        A = magma(random_matrix(QQ, n, n, num_bound=bnd, den_bound=bnd))
        code = """
n := %s;
A := %s;
B := A + 1;
t := Cputime();
for z in [1..%s] do
    K := A * B;
end for;
s := Cputime(t);
"""%(n, A.name(), times)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError, 'unknown system "%s"'%system


# Determinant of Hilbert matrix
def det_hilbert_QQ(n=80, system='sage'):
    if system == 'sage':
        A = hilbert_matrix(n)
        t = cputime()
        d = A.determinant()
        return cputime(t)
    elif system == 'magma':
        code = """
h := HilbertMatrix(%s);
tinit := Cputime();
d := Determinant(h);
s := Cputime(tinit);
delete h;
"""%n
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))

# inverse of Hilbert matrix
def invert_hilbert_QQ(n=40, system='sage'):
    if system == 'sage':
        A = hilbert_matrix(n)
        t = cputime()
        d = A**(-1)
        return cputime(t)
    elif system == 'magma':
        code = """
h := HilbertMatrix(%s);
tinit := Cputime();
d := h^(-1);
s := Cputime(tinit);
delete h;
"""%n
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))


#######################################################################
# Dense Benchmarks over machine reals
# Note that the precision in reals for MAGMA is base 10, while in
# sage it is in base 2
#######################################################################

def report_RR():
    F = [rank_RR, rank2_RR, nullspace_RR, charpoly_RR, smithform_RR,
         matrix_multiply_RR, det_RR]
    title = 'Dense benchmarks over RR'
    report(F, title)

# Real Nullspace

def nullspace_RR(n=300, min=0, max=10, system='sage'):
    """
    Nullspace over RR:
    Given a n+1 x n (with n=300) matrix over RR with random entries
    between min=0 and max=10, compute the nullspace.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n+1, n, x=min, y=max+1).change_ring(RR)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RealField(16), n+1,n)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(n,min,max)
        if verbose: print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system


