#*****************************************************************************
#      Copyright (C) 2007 Robert Bradshaw <robertwb@math.washington.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************


# Utility functions for operating on the point_c struct.
# These would be macros, but inline functions should do the same job.


cdef extern from "math.h":
    double stdmath_sqrt "sqrt" (double)

cdef inline bint point_c_set(point_c* res, P) except -2:
    res.x, res.y, res.z = P[0], P[1], P[2]

cdef inline bint point_c_eq(point_c P, point_c Q):
    return P.x == Q.x and P.y == Q.y and P.z == Q.z

cdef inline int point_c_cmp(point_c P, point_c Q):
    """
    Lexographic order
    """
    if P.x == Q.x:
        if P.y == Q.y:
            if P.z == Q.z:
                return 0
            elif P.z < Q.z:
                return -1
            else:
                return 1
        elif P.y < Q.y:
            return -1
        else:
            return 1
    elif P.x < Q.x:
        return -1
    else:
        return 1

cdef inline void point_c_lower_bound(point_c* res, point_c P, point_c Q):
    res.x = P.x if P.x < Q.x else Q.x
    res.y = P.y if P.y < Q.y else Q.y
    res.z = P.z if P.z < Q.z else Q.z

cdef inline void point_c_upper_bound(point_c* res, point_c P, point_c Q):
    res.x = P.x if P.x > Q.x else Q.x
    res.y = P.y if P.y > Q.y else Q.y
    res.z = P.z if P.z > Q.z else Q.z

cdef inline void point_c_add(point_c* res, point_c P, point_c Q):
    res.x = P.x + Q.x
    res.y = P.y + Q.y
    res.z = P.z + Q.z

cdef inline void point_c_sub(point_c* res, point_c P, point_c Q):
    res.x = P.x - Q.x
    res.y = P.y - Q.y
    res.z = P.z - Q.z

cdef inline void point_c_mul(point_c* res, point_c P, double a):
    res.x = a * P.x
    res.y = a * P.y
    res.z = a * P.z

cdef inline double point_c_dot(point_c P, point_c Q):
    return P.x*Q.x + P.y*Q.y + P.z*Q.z

cdef inline double point_c_len(point_c P):
    return stdmath_sqrt(point_c_dot(P, P))

cdef inline void point_c_transform(point_c* res, double* M, point_c P):
    """
    M is a flattened 4x4 matrix, row major, representing a Euclidean Transformation.
    Operate on P as a point.
    """
    res.x = M[0]*P.x + M[1]*P.y + M[2]*P.z + M[3]
    res.y = M[4]*P.x + M[5]*P.y + M[6]*P.z + M[7]
    res.z = M[8]*P.x + M[9]*P.y + M[10]*P.z + M[11]

cdef inline void point_c_stretch(point_c* res, double* M, point_c P):
    """
    M is a flattened 4x4 matrix, row major, representing a Euclidean Transformation.
    Operate on P as a vector (i.e. ignore the translation component)
    """
    res.x = M[0]*P.x + M[1]*P.y + M[2]*P.z
    res.y = M[4]*P.x + M[5]*P.y + M[6]*P.z
    res.z = M[8]*P.x + M[9]*P.y + M[10]*P.z

