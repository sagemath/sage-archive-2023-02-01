cdef extern from "eclib/moddata.h":
    pass
cdef extern from "eclib/symb.h":
    pass
cdef extern from "eclib/cusp.h":
    pass

cdef extern from "eclib/homspace.h":

    # From mat.h
    ctypedef int scalar   # TODO: int or long??

    ctypedef struct mat "mat":
        scalar* get_entries()   # TODO: possibly not int --
        scalar sub(long,long)

    long nrows(mat M)
    long ncols(mat M)
    mat addscalar(mat M, scalar)
    long rank(mat M)

    # Constructors
    mat *new_mat "new mat" (mat m)

    # General C++ stuff
    void delete_mat "delete "(mat* m)

cdef class Matrix:
    cdef mat* M

    cdef set(self, mat*  M)

cdef class MatrixFactory:
    cdef new_matrix(self, mat M)




