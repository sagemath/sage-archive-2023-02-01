r"""
Methods for creating the five platonic solids.

EXAMPLES:
    sage: from sage.plot.plot3d.platonic import *
    sage: all = Tetrahedron(color='blue').translate(0,-3.5,0) + Cube(color=(.25,0,.5)).translate(0,-2,0) + Octahedron(color='red') + Dodecahedron(color='orange').translate(0,2,0) + Icosahedron(color='yellow').translate(0,4,0)
    sage: all.scale(.25).show()

    sage: Icosahedron().stickers(['red','blue'], .075, .1).show()

AUTHOR:
    -- Robert Bradshaw 2007-08: initial version

"""


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



from sage.rings.all import RDF
from sage.matrix.constructor import matrix

from shapes import Box
from index_face_set import IndexFaceSet


def Tetrahedron(**kwds):
    """
    sage: S = Tetrahedron(color='yellow') + Sphere(.7, color='red')
    sage: S.show()
    """
    RR = RDF
    one = RR(1)
    sqrt2 = RR(2).sqrt()
    sqrt6 = RR(6).sqrt()
    point_list = [(0,0,1),
                  (2*sqrt2/3,        0, -one/3),
                  ( -sqrt2/3,  sqrt6/3, -one/3),
                  ( -sqrt2/3, -sqrt6/3, -one/3)]
    face_list = [[0,1,2],[1,3,2],[0,2,3],[0,3,1]]
    return IndexFaceSet(face_list, point_list, enclosed=True, **kwds)


def Cube(**kwds):
    RR = RDF
    sqrt3 = RR(3).sqrt()
    return Box(1/sqrt3, 1/sqrt3, 1/sqrt3, **kwds)


def Octahedron(**kwds):
    return Box(1,1,1).dual(**kwds)


def Dodecahedron(**kwds):
    """
    CONSTRUCTION:
        We let one point be $Q = (0,1,1)$.

        Now there are three points spaced equally on a circle
        around the north pole. The other requirement is that
        the angle between them be the angle of a pentagon, namely
        $3\pi/5$. This is enough to determine them. Placing one
        on the $xz$-plane we have.

        $P_1 = \left(t, 0, \sqrt{1-t^2}\right)$

        $P_2 = \left(-\frac{1}{2}t, \frac{\sqrt{3}}{2}t, \sqrt{1-t^2}\right)$

        $P_3 = \left(-\frac{1}{2}t, \frac{\sqrt{3}}{2}t, \sqrt{1-t^2}\right)$

        $Solving $\frac{(P_1-Q) \cdot (P_2-Q)}{|P_1-Q||P_2-Q|} = \cos(3\pi/5)$ we get $t = 2/3$.

        Now we have 6 points $R_1, ..., R_6$ to close the three top pentagons.
        These can be found by mirroring $P_2$ and $P_3$ by the $yz$-plane and
        rotating around the $y$-axis by the angle $\theta$ from $Q$ to $P_1$.
        Note that $\cos(\theta) = t = 2/3$ and so $\sin(\theta) = \sqrt{5}/3$.
        Rotation gives us the other four.

        Now we reflect through the origin for the bottom half.


    AUTHOR:
        -- Robert Bradshaw
    """
    RR = RDF
    one = RR(1)
    sqrt3 = RR(3).sqrt();
    sqrt5 = RR(5).sqrt()
    R3 = RR**3
    rot = matrix(RR, [[  -one/2,-sqrt3/2, 0],
                      [ sqrt3/2,  -one/2, 0],
                      [       0,       0, 1]])
    rot2 = rot*rot

    # The top
    Q = R3([0,0,1])
    # The first ring
    P1 = R3([2*one/3, 0, sqrt5/3])
    # The second ring
    R1 = R3([sqrt5/3,  1/sqrt3, one/3])
    R2 = R3([sqrt5/3, -1/sqrt3, one/3])

    top = [Q, P1, rot*P1, rot2*P1, R1, rot*R2, rot*R1, rot2*R2, rot2*R1, R2]
    point_list = top + [-p for p in reversed(top)]

    top_faces = [[0,1,4,5,2],
                 [0,2,6,7,3],
                 [0,3,8,9,1],
                 [1,9,13,12,4],
                 [2,5,11,10,6],
                 [3,7,15,14,8]]
    face_list = top_faces + [reversed([19-p for p in f]) for f in top_faces]

    return IndexFaceSet(face_list, point_list, enclosed=True, **kwds)

#    if style == 'vertices' or style == 'edges':
#        from sage.plot.plot import rainbow
#        colors = rainbow(len(vs), 'rgbtuple')
#        #vertex_spheres = [Box(.05, .05, .05, color=color).translate(p) for p in vs]
#        vertex_spheres = [Box(.05, .05, .05, color=c).translate(p) for p,c in zip(vs,colors)]
#        faces = IndexFaceSet([[tuple(vs[i]) for i in f] for f in face_list])
#        vertex_spheres += [faces.stickers(['red','yellow','blue','purple','black','orange'], .1, .1)] # [faces]
#        return Graphics3dGroup(vertex_spheres)


def Icosahedron(**kwds):
    return Dodecahedron().dual(**kwds)