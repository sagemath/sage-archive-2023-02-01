"""
Generic Strassen-Winograd matrix multiplication method

CURRENTLY -- half way done and broken -- ignore...

"""

import constructor

# NOTATION:
#   We store a submatrix by giving the ul (=upper left) and lr (=lower right)
#   position of the matrix in a bigger matrix.

def decompose(submatrix):
    """
    Given positions that define a submatrix, return the positions that
    define the upper left, upper right, lower left, lower right.
    """
    xmin, ymin = submatrix[0]
    xmax, ymax = submatrix[1]
    xmid = (xmax + xmin)//2
    ymid = (ymax + ymin)//2

    A00 = ((xmin, ymin), (xmid,ymid))
    A01 = ((xmid, ymin), (xmax,ymid))
    A10 = ((xmin, ymid), (xmid,ymax))
    A11 = ((xmid, ymid), (xmax,ymax))
    return A00, A01, A10, A11

def mat_range(M):
    return ((0,0),(M.nrows()-1,M.ncols()-1))

def internal_add(A,B, AR, BR):
    """
    internal_add is adding submatrices together, and
    constructing a matrix to return.
    """
    A_xmin, A_ymin = AR[0]
    A_xmax, A_ymax = AR[1]
    B_xmin, B_ymin = BR[0]
    B_xmax, B_ymax = BR[1]

    C = constructor.Matrix(A.base_ring(), A_xmax - A_xmin, A_ymax - A_ymin)
    xRange = range(AR[3]-AR[1]);
    yRange = range(AR[2]-AR[0]);
    for i in xRange:
            for j in yRange:
                    C[j,i] = A[AR[0]+j,AR[1]+i]+B[BR[0]+j,BR[1]+i]
    return C

def internal_sub(A,B,AR,BR):
    """
    For now, assume that A and B are the same size.
    """
    C = constructor.Matrix(A.base_ring(),AR[2]-AR[0],BR[3]-BR[1])
    xRange = range(AR[3]-AR[1]);
    yRange = range(AR[2]-AR[0]);
    for i in xRange:
            for j in yRange:
                    C[j,i] = A[AR[0]+j,AR[1]+i]-B[BR[0]+j,BR[1]+i]
    return C

def strassen_internal_add(C,A,B,CR):
    xRange = range(CR[3]-CR[1]);
    yRange = range(CR[2]-CR[0]);
    for i in xRange:
            for j in yRange:
                    C[CR[0]+j,CR[1]+i] = A[j,i]+B[j,i]


def strassen_internal_mult(A,B,AR,BR):
    """
    This is the crux of Strassun matrix multiply.
    Assume for now that the sizes are 2^n.
    AR = [AMinY,AMinX,AMaxY,AMaxX]. This is a bit counter intuitive,
    but I've made it to look like the notes that I have.
    """
    if AR[2]-AR[0]  == 1:
        C = constructor.Matrix(A.base_ring(),1,1,[A[AR[0],AR[1]]*B[BR[0],BR[1]]])
        return C

    A11R,A12R,A21R,A22R = decompose(AR)
    B11R,B12R,B21R,B22R = decompose(BR)

    """
    We need to create a matrix of the right type.
    I will figure out how to do this soon. Note that
    I am taking the ring as defined for A.
    """
    C = constructor.Matrix(A.base_ring(),BR[2]-BR[0],AR[3]-AR[1])
    [C11R,C12R,C21R,C22R] = decompose(mat_range(C));
    """
    this is what I'm doing:
    S1 = A21 + A22,  T1 = B12 - B11
    S2 = S1 - A11,     T2 = B22 - T1
    S3 = A11 - A21,  T3 = B22 - B12
    S4 = A12 - S2,    T4 = B21 - T2

    """
    S1 = internal_add(A,A,A21R,A22R)
    S2 = internal_sub(S1,A,mat_range(S1),A11R)
    S3 = internal_sub(A,A,A11R,A21R)
    S4 = internal_sub(A,S2,A12R,mat_range(S2))

    T1 = internal_sub(B,B,B12R,B11R)
    T2 = internal_sub(B,T1,B22R,mat_range(T1))
    T3 = internal_sub(B,B,B22R,B12R)
    T4 = internal_sub(B,T2,B21R,mat_range(T2))

    """
    Now we do the matrix multiplications:
    P1 = A11*B11
    P2 = A12*B21
    P3 = S1*T1
    P4 = S2*T2
    P5 = S3*T3
    P6 =  S4*B22
    P7 = A22*T4
    """
    P1 = strassen_internal_mult(A,B,A11R,B11R)
    P2 = strassen_internal_mult(A,B,A12R,B21R)
    P3 = strassen_internal_mult(S1,T1,mat_range(S1),mat_range(T1))
    P4 = strassen_internal_mult(S2,T2,mat_range(S2),mat_range(T2))
    P5 = strassen_internal_mult(S3,T3,mat_range(S3),mat_range(T3))
    P6 = strassen_internal_mult(S4,B,mat_range(S4),B22R)
    P7 = strassen_internal_mult(A,T4,A22R,mat_range(T4))

    """
    Now bunch more additions...
    U1 = P1 + P2
    U2 = P1 + P4
    U3 = U2 + P5
    U4 = U3 + P7
    U5 = U3 + P3
    U6 = U2 + P3
    U7 = U6 + P6
    C11 = U1, C12 = U7, C21 = U5, C22 = U6
    """
    strassen_internal_add(C,P1,P2,C11R)
    U1 = P1+P2
    U2 = P1+P4
    U3 = U2+P5
    U4 = U3+P7
    U5 = U3+P3
    U6 = U2+P3
    U7 = U6+P6
    strassen_internal_add(C,U6,P6,C12R)
    strassen_internal_add(C,U3,P7,C21R)
    strassen_internal_add(C,U3,P3,C22R)

    print S1,S2,S3,S4
    print T1,T2,T3,T4
    print P1,P2,P3,P4,P5,P6,P7
    print U1,U2,U3,U4,U5,U6,U7

    return C


def strassen_mult(A, B):
    return strassen_internal_mult(A, B, mat_range(A), mat_range(B))

