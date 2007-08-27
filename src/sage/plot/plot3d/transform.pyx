from math import atan2, sin, cos

from sage.rings.real_double import RDF
from sage.misc.functional import sqrt, atan, acos

from sage.matrix.constructor import matrix
from sage.modules.free_module_element import vector

pi = RDF.pi()


class Transformation:
    def __init__(self, scale=(1,1,1),
                       rot=None,
                       trans=[0,0,0],
                       m=None):

        if scale is None:
            scale = (1,1,1)
        if trans is None:
            trans = [0,0,0]

        # TODO: determine for sure if x3d does scale or rotation first
        if m is None:
            m = matrix(RDF, 3, 3,
                      [scale[0], 0, 0, 0, scale[1], 0, 0, 0, scale[2]])

        if rot is not None:
            # rotate about v by theta
            vx, vy, vz, theta = rot
            t = atan2(vy,vz) + pi/2
            m = self.rotX(t) * m
            new_y = vy*cos(t) - vz*sin(t)
            # v = [vx, new_y, 0]
            s = atan2(vx,new_y) + pi/2
            m = self.rotZ(s) * m
            # v = [new_x, 0, 0]
            m = self.rotX(theta) * m
            # now put back to our former reference frame
            m = self.rotZ(-s) * m
            m = self.rotX(-t) * m

        self.matrix = m.augment(matrix(RDF, 3, 1, list(trans))) \
                       .stack(matrix(RDF, 1, 4, [0,0,0,1]))

    def is_skew(self, eps=1e-5):
        dx, dy, dz = self.matrix.submatrix(0,0,3,3).columns()
        return abs(dx.dot_product(dy)) + abs(dx.dot_product(dz)) + abs(dy.dot_product(dz)) > eps

    def is_uniform(self, eps=1e-5):
        cols = self.matrix.submatrix(0,0,3,3).columns()
        lens = [col.dot_product(col) for col in cols]
        return abs(lens[0] - lens[1]) + abs(lens[0] - lens[2]) < eps

    def is_uniform_on(self, basis, eps=1e-5):
        basis = [vector(RDF, self.transform_vector(b)) for b in basis]
        a = basis.pop()
        len_a = a.dot_product(a)
        return max([len_a - b.dot_product(b) for b in basis]) < eps

    def rotX(self, theta):
        return matrix(RDF, 3, 3, [1, 0, 0,
                                  0, cos(theta), -sin(theta),
                                  0, sin(theta), cos(theta)])

    def rotZ(self, theta):
        return matrix(RDF, 3, 3, [cos(theta), -sin(theta), 0,
                                  sin(theta), cos(theta), 0,
                                  0, 0, 1])

    def transform_point(self, x):
        Tx = self.matrix * vector(RDF, [x[0], x[1], x[2], 1])
        return (Tx[0], Tx[1], Tx[2])

    def transform_vector(self, x):
        Tx = self.matrix * vector(RDF, [x[0], x[1], x[2], 0])
        return (Tx[0], Tx[1], Tx[2])

    def __mul__(self, other):
        T = Transformation()
        T.matrix = self.matrix * other.matrix
        return T

    def __call__(self, p):
        res = self.matrix * vector(RDF, [p[0], p[1], p[2], 1])
        return tuple(res[:3])
