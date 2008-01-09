r"""
Methods for creating the five platonic solids.

EXAMPLES:
The five platonic solids in a row;
    sage: G = tetrahedron((0,-3.5,0), color='blue') + cube((0,-2,0),color=(.25,0,.5)) + octahedron(color='red') + dodecahedron((0,2,0), color='orange') + icosahedron(center=(0,4,0), color='yellow')
    sage: G.show(aspect_ratio=[1,1,1])

All the platonic solids in the same place:
    sage: G = tetrahedron(color='blue',opacity=0.7) + cube(color=(.25,0,.5), opacity=0.7) + octahedron(color='red', opacity=0.7) + dodecahedron(color='orange', opacity=0.7) + icosahedron(opacity=0.7)
    sage: G.show(aspect_ratio=[1,1,1])

Display nice faces only:
    sage: icosahedron().stickers(['red','blue'], .075, .1)

AUTHOR:
    -- Robert Bradshaw 2007-08: initial version
    -- William Stein

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

from shapes import Box, ColorCube
from shapes2 import frame3d
from index_face_set import IndexFaceSet

def index_face_set(face_list, point_list, enclosed, **kwds):
    if kwds.has_key('center'):
        center = kwds['center']
        del kwds['center']
    else:
        center = (0,0,0)
    if kwds.has_key('size'):
        size = kwds['size']
        del kwds['size']
    else:
        size = 1
    I = IndexFaceSet(face_list, point_list, enclosed=enclosed, **kwds)
    return prep(I, center, size, kwds)

def prep(G, center, size, kwds):
    if center != (0,0,0):
        G = G.translate(center)
    if size != 1:
        G = G.scale(size)
    G._set_extra_kwds(kwds)
    return G

def tetrahedron(center=(0,0,0), size=1, **kwds):
    """
    A 3d tetrahedron.

    INPUTS:
        center -- (default: (0,0,0))
        size -- (default: 1)
        color -- a word that describes a color
        rgbcolor -- (r,g,b) with r, g, b between 0 and 1 that describes a color
        opacity -- (default: 1) if less than 1 then is transparent

    EXAMPLES:
    A default colored tetrahedron at the origin:
        sage: tetrahedron()

    A transparent green tetrahedron in front of a solid red one:
        sage: tetrahedron(opacity=0.8, color='green') + tetrahedron((-2,1,0),color='red')

    A translucent tetrahedron sharing space with a sphere:
        sage: tetrahedron(color='yellow',opacity=0.7) + sphere(r=.5, color='red')


    A big tetrahedron:
        sage: tetrahedron(size=10)

    A wide tetrahedron:
        sage: tetrahedron(aspect_ratio=[1,1,1]).scale((4,4,1))

    A red and blue tetrahedron touching noses:
        sage: tetrahedron(color='red') + tetrahedron((0,0,-2)).scale([1,1,-1])

    AUTHOR:
        -- Robert Bradshaw and William Stein
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
    return index_face_set(face_list, point_list, enclosed=True, center=center, size=size, **kwds)

def cube(center=(0,0,0), size=1, color=None, frame_thickness=0, frame_color=None, **kwds):
    """
    A 3d cube centered at the origin with default side lengths 1.

    INPUTS:
        center -- (default: (0,0,0))
        size -- (default: 1) the side lengths of the cube
        color -- a string that describes a color; this can also be a list of 3-tuples or
                 strings length 6 or 3, in which case the faces (and oppositive faces)
                 are colored.
        frame_thickness -- (default: 0) if positive, then thickness of the frame
        frame_color -- (default: None) if given, gives the color of the frame
        opacity -- (default: 1) if less than 1 then is transparent

    EXAMPLES:
    A simple cube:
        sage: cube()

    A red cube:
        sage: cube(color="red")

    A transparent grey cube that contains a red cube:
        sage: cube(opacity=0.8, color='grey') + cube(size=3/4)

    A transparent colored cube:
        sage: cube(color=['red', 'green', 'blue'], opacity=0.5)

    A bunch of random cubes:
        sage: sum([cube((10*random(), 10*random(), 10*random()), size=random()/3, color=(random(),random(), random())) for _ in [1..30]])

    Non-square cubes (boxes):
        sage: cube(aspect_ratio=[1,1,1]).scale([1,2,3])
        sage: cube(color=['red', 'blue', 'green'],aspect_ratio=[1,1,1]).scale([1,2,3])

    And one that is colored:
        sage: cube(color=['red', 'blue', 'green', 'black', 'white', 'orange'], aspect_ratio=[1,1,1]).scale([1,2,3])

    A nice translucent color cube with a frame:
        sage: c = cube(color=['red', 'blue', 'green'], frame=False, frame_thickness=2, frame_color='brown', opacity=0.8)
        sage: c

    A raytraced color cube with frame and transparency:
        sage: c.show(viewer='tachyon')

    AUTHOR:
        -- William Stein
    """
    if isinstance(color, (list, tuple)) and len(color) > 0 and isinstance(color[0], (list,tuple,str)):
        B = ColorCube(size=[0.5,0.5,0.5], colors=color, **kwds)
    else:
        if color is not None:
            kwds['color'] = color
        B = Box(0.5,0.5,0.5, **kwds)
    if frame_thickness > 0:
        if frame_color is None:
            B += frame3d((-0.5,-0.5,-0.5),(0.5,0.5,0.5), thickness=frame_thickness)
        else:
            B += frame3d((-0.5,-0.5,-0.5),(0.5,0.5,0.5), thickness=frame_thickness, color=frame_color)
    return prep(B, center, size, kwds)

def octahedron(center=(0,0,0), size=1, **kwds):
    """
    Return an octahedron.

    INPUTS:
        center -- (default: (0,0,0))
        size -- (default: 1)
        color -- a string that describes a color; this can also be a list of 3-tuples or
                 strings length 6 or 3, in which case the faces (and oppositive faces)
                 are colored.
        opacity -- (default: 1) if less than 1 then is transparent

    EXAMPLES:
        sage: octahedron((1,4,3), color='orange') + octahedron((0,2,1), size=2, opacity=0.6)
    """
    return prep(Box(1,1,1).dual(**kwds), center, size, kwds)

def dodecahedron(center=(0,0,0), size=1, **kwds):
    """
    A dodecahedron.

    INPUTS:
        center -- (default: (0,0,0))
        size -- (default: 1)
        color -- a string that describes a color; this can also be a list of 3-tuples or
                 strings length 6 or 3, in which case the faces (and oppositive faces)
                 are colored.
        opacity -- (default: 1) if less than 1 then is transparent

    EXAMPLES:
    A plain Dodecahedron:
        sage: dodecahedron()

    A translucent dodecahedron that contains a black sphere:
        sage: dodecahedron(color='orange', opacity=0.8) + sphere(size=0.5, color='black')

    CONSTRUCTION:
        This is how we construct a dodecahedron.   We let one point be $Q = (0,1,0)$.

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
        -- Robert Bradshaw, William Stein
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

    return index_face_set(face_list, point_list, enclosed=True, center=center, size=size, **kwds)

#    if style == 'vertices' or style == 'edges':
#        from sage.plot.plot import rainbow
#        colors = rainbow(len(vs), 'rgbtuple')
#        #vertex_spheres = [Box(.05, .05, .05, color=color).translate(p) for p in vs]
#        vertex_spheres = [Box(.05, .05, .05, color=c).translate(p) for p,c in zip(vs,colors)]
#        faces = IndexFaceSet([[tuple(vs[i]) for i in f] for f in face_list])
#        vertex_spheres += [faces.stickers(['red','yellow','blue','purple','black','orange'], .1, .1)] # [faces]
#        return Graphics3dGroup(vertex_spheres)


def icosahedron(center=(0,0,0), size=1, **kwds):
    """
    An icosahedron.

    INPUTS:
        center -- (default: (0,0,0))
        size -- (default: 1)
        color -- a string that describes a color; this can also be a list of 3-tuples or
                 strings length 6 or 3, in which case the faces (and oppositive faces)
                 are colored.
        opacity -- (default: 1) if less than 1 then is transparent

    EXAMPLES:
        sage: icosahedron()

    Two icosahedrons at different positions of different sizes.
        sage: icosahedron((-1/2,0,1), color='orange') + icosahedron((2,0,1), size=1/2, aspect_ratio=[1,1,1])
    """
    return prep(dodecahedron().dual(**kwds), center, size, kwds)
