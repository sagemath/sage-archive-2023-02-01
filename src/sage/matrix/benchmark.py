from sage.all import *

#######################################################################
# Dense Benchmarks over ZZ
#######################################################################



# Integer Nullspace

def nullspace_ZZ(n=300, min=0, max=10, system='sage'):
    """
    Given a n+1 x n (with n=300) matrix over ZZ with random entries
    between min=0 and max=10, compute the nullspace.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n+1, x=min, y=max+1).change_ring(QQ)
        t = cputime()
        v = A.kernel()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := RMatrixSpace(RationalField(), n, n+1)![Random(%s,%s) : i in [1..n*(n+1)]];
t := Cputime();
K := Kernel(A);
s := Cputime(t);
"""%(n,min,max)
        print code
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError, 'unknown system "%s"'%system


# Characteristic Polynomial over ZZ

def charpoly_ZZ(n=100, min=0, max=9, system='sage'):
    """
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
        print code
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError, 'unknown system "%s"'%system




# Smith Form

def smithform_ZZ(n=100, min=0, max=9, system='sage'):
    """
    Given a n x n (with n=100) matrix over ZZ with random entries
    between min=0 and max=9, compute the charpoly.
    """
    if system == 'sage':
        A = random_matrix(ZZ, n, n, x=min, y=max+1)
        t = cputime()
        v = A.smith_form()
        return cputime(t)
    elif system == 'magma':
        code = """
n := %s;
A := MatrixAlgebra(IntegerRing(), n)![Random(%s,%s) : i in [1..n^2]];
t := Cputime();
K := SmithForm(A);
s := Cputime(t);
"""%(n,min,max)
        print code
        magma.eval(code)
        return float(magma.eval('s'))
    else:
        raise ValueError, 'unknown system "%s"'%system


# Matrix multiplication over ZZ
def matrix_multiply_ZZ(n=100, min=-9, max=9, system='sage', times=1):
    """
    Given an n x n (with n=100) matrix A over ZZ with random entries
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
        print code
        magma.eval(code)
        return float(magma.eval('s'))/times
    else:
        raise ValueError, 'unknown system "%s"'%system



#######################################################################
# Dense Benchmarks over GF(p), for small p.
#######################################################################

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
        print code
        magma.eval(code)
        return magma.eval('s')
    else:
        raise ValueError, 'unknown system "%s"'%system


# Characteristic Polynomial over GF

def charpoly_GF(n=100, p=16411, system='sage'):
    """
    Given a n x n (with n=100) matrix over ZZ with random entries,
    compute the charpoly.
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
        print code
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
        print code
        magma.eval(code)
        return float(magma.eval('s'))/times
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
        print code
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
        print code
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
        print code
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
        print code
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
        print code
        magma.eval(code)
        return float(magma.eval('s'))

