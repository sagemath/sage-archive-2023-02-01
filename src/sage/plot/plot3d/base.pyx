r"""
Base classes for 3D Graphics objects and plotting.

EXAMPLES:
    sage: from sage.plot.graphics3d import *
    sage: S = ColorCube(.35, ['green', 'yellow', 'blue']) + Sphere(.2, color='red').translate(.4,.4,.4)
    sage: S.show()

    sage: from sage.plot.plot3d.plot3d import plot3d
    sage: from sage.plot.plot3d.shapes import *
    sage: S = Sphere(.5, color='yellow')
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)
    sage: S += Sphere(.1, color='white').translate(.45,-.1,.15) + Sphere(.05, color='black').translate(.51,-.1,.17)
    sage: S += Sphere(.1, color='white').translate(.45, .1,.15) + Sphere(.05, color='black').translate(.51, .1,.17)
    sage: S += Sphere(.1, color='yellow').translate(.5, 0, -.2)
    sage: def f(x,y): return math.exp(x/5)*math.cos(y)
    sage: P = plot3d(f,(-5,5),(-5,5), ['red','yellow'], max_depth=10)
    sage: cape_man = P.scale(.2)+S.translate(1,0,0)
    sage.: cape_man.show()


AUTHOR:
    -- Robert Bradshaw 2007-02: inital version
    -- Robert Bradshaw 2007-08: cythonization, much optimization

TODO:
    -- finish integrating tachyon
    -- good default lights, camera
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


include "../../ext/python_list.pxi"


from math import atan2

from sage.modules.free_module_element import vector

from sage.rings.real_double import RDF
from sage.misc.functional import sqrt, atan, acos
#from sage.functions.all import *
from texture import Texture, is_Texture
from transform import Transformation
pi = RDF.pi()

from sage.interfaces.tachyon import tachyon_rt


default_texture = Texture()


cdef class Graphics3d(SageObject):


    def __add__(self, other):
        return Graphics3dGroup([self, other])

    def transform(self, **kwds):
        return TransformGroup([self], **kwds)

    def translate(self, *x):
        if isinstance(x[0], (tuple, list)):
            x = x[0]
        return self.transform(trans=x)

    def scale(self, *x):
        if isinstance(x[0], (tuple, list)):
            x = x[0]
        return self.transform(scale=x)

    def rotate(self, v, theta):
        vx, vy, vz = v
        return self.transform(rot=[vx, vy, vz, theta])

    def rotateX(self, theta):
        return self.rotate((1,0,0), theta)

    def rotateY(self, theta):
        return self.rotate((0,1,0), theta)

    def rotateZ(self, theta):
        return self.rotate((0,0,1), theta)


    def viewpoint(self):
        return Viewpoint(0,0,6)

    def default_render_params(self):
        return RenderParams(ds=.075)

    def x3d(self):
        return """
<X3D version='3.0' profile='Immersive' xmlns:xsd='http://www.w3.org/2001/XMLSchema-instance' xsd:noNamespaceSchemaLocation=' http://www.web3d.org/specifications/x3d-3.0.xsd '>
<head>
<meta name='title' content='sage3d'/>
</head>
<Scene>
%s
%s
</Scene>
</X3D>
"""%(self.viewpoint().x3d_str(), self.x3d_str())

    def tachyon(self):
        render_params = self.default_render_params()
        # switch from LH to RH coords to be consistant with java rendition
        render_params.push_transform(Transformation(scale=[1,-1,1]))
        return """
begin_scene
resolution 400 400


         camera
            zoom 1.0
            aspectratio 1.0
            antialiasing 1
            raydepth 8
            center  2.0 1.0 0.75
            viewdir  -2.0 -1.0 -0.5
            updir  0.0 0.0 1.0
         end_camera


      light center  4.0 3.0 2.0
            rad 0.2
            color  1.0 1.0 1.0

      plane
        center -2000 -1000 -500
        normal -2.0 -1.0 -0.5
        TEXTURE
            AMBIENT 1.0 DIFFUSE 1.0 SPECULAR 1.0 OPACITY 1.0
            COLOR 1.0 1.0 1.0
            TEXFUNC 0

    %s

    %s

end_scene""" % (
               "\n".join([t.tachyon_str() for t in self.texture_set()]),
               self.tachyon_str(render_params))

    def tachyon_str(self, render_params):
        """
        DO NOT override this method, override tachyon_repr instead.
        """
        return "\n".join(flatten_list(self.tachyon_repr(render_params)))

    def obj(self):
        return self.obj_str(self.default_render_params())

    def obj_str(self, render_params):
        """
        DO NOT override this method, override obj_repr instead.
        """
        return "\n".join(flatten_list([self.obj_repr(render_params), ""]))

    def texture_set(self):
        return set()

    def mtl_str(self):
        return "\n\n".join([t.mtl_str() for t in self.texture_set()]) + "\n"

    def flatten(self, T=None):
        if T is None:
            return self
        else:
            return self.transform(T=T)

    def show(self, interactive=True, filename="shape", verbosity=0):
        tachyon_rt(self.tachyon(), filename+".png", verbosity, True, '')
        if interactive:
            f = open(filename+".obj", "w")
            f.write("mtllib %s.mtl\n" % filename)
            f.write(self.obj())
            f.close()
            f = open(filename+".mtl", "w")
            f.write(self.mtl_str())
            f.close()


class Graphics3dGroup(Graphics3d):
    def __init__(self, all=[]):
        self.all = all

    def transform(self, **kwds):
        return TransformGroup(self.all, **kwds)

    def tachyon_repr(self, render_params):
        return [g.tachyon_repr(render_params) for g in self.all]

    def x3d_str(self):
        return "\n".join([g.x3d_str() for g in self.all])

    def obj_repr(self, render_params):
        return [g.obj_repr(render_params) for g in self.all]

    def texture_set(self):
        return reduce(set.union, [g.texture_set() for g in self.all])

    def flatten(self, T=None):
        if len(self.all) == 1:
            return self.all[0].flatten(T)
        all = []
        for g in self.all:
            g = g.flatten(T)
            if type(g) is Graphics3dGroup:
                all += g.all
            else:
                all.append(g)
        return Graphics3dGroup(all)



class TransformGroup(Graphics3dGroup):

    def __init__(self, all=[], rot=None, trans=None, scale=None, T=None):
        Graphics3dGroup.__init__(self, all)
        self._rot = rot
        self._trans = trans
        if scale is not None and len(scale) == 1:
            if isinstance(scale, (tuple, list)):
                scale = scale[0]
            scale = (scale, scale, scale)
        self._scale = scale
        if T is not None:
            self.T = T

    def transform(self, **kwds):
        # TODO: flatten right here
        return TransformGroup([self], **kwds)

    def x3d_str(self):
        s = "<Transform"
        if self._rot is not None:
            s += " rotation='%s %s %s %s'"%tuple(self._rot)
        if self._trans is not None:
            s += " translation='%s %s %s'"%tuple(self._trans)
        if self._scale is not None:
            s += " scale='%s %s %s'"%tuple(self._scale)
        s += ">\n"
        s += Graphics3dGroup.x3d_str(self)
        s += "\n</Transform>"
        return s

    def tachyon_repr(self, render_params):
        render_params.push_transform(self.get_transformation())
        rep = [g.tachyon_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def obj_repr(self, render_params):
        render_params.push_transform(self.get_transformation())
        rep = [g.obj_repr(render_params) for g in self.all]
        render_params.pop_transform()
        return rep

    def get_transformation(self):
        try:
            return self.T
        except AttributeError:
            self.T = Transformation(self._scale, self._rot, self._trans)
            return self.T

    def flatten(self, T=None):
        assert False, "broken"
        all = []
        for g in self.all:
            g = g.flatten().transform(T=self.get_transformation())
            if type(g) is Graphics3dGroup:
                all += g.all
            else:
                all.append(g)
        self.all = all




class Viewpoint(Graphics3d):

    def __init__(self, *x):
        if isinstance(x[0], (tuple, list)):
            x = tuple(x[0])
        self.pos = x

    def x3d_str(self):
        return "<Viewpoint position='%s %s %s'/>"%self.pos



cdef class PrimativeObject(Graphics3d):
    def __init__(self, **kwds):
        try:
            self.texture = kwds['texture']
            if not is_Texture(self.texture):
                self.texture = Texture(self.texture)
        except KeyError:
            try:
                self.texture = kwds['color']
                if not is_Texture(self.texture):
                    self.texture = Texture(self.texture)
            except KeyError:
                self.texture = default_texture

    def set_texture(self, texture):
        if not is_Texture(texture):
            texture = Texture(texture)
        self.texture = texture

    def x3d_str(self):
        return "<Shape>" + self.x3d_geometry() + self.texture.x3d_str() + "</Shape>\n"

    def tachyon_repr(self, render_params):
        return self.triangulation().tachyon_repr(render_params)

    def obj_repr(self, render_params):
        return self.triangulation().obj_repr(render_params)

    def texture_set(self):
        return set([self.texture])



class BoundingSphere(SageObject):
    def __init__(self, cen, r):
        self.cen = vector(RDF, cen)
        self.r = r
    def __repr__(self):
        return "Center %s radius %s" % (self.cen, self.r)
    def __add__(self, other):
        if self.cen == other.cen:
            return self if self.r > other.r else other
        diff = other.cen - self.cen
        dist = (diff[0]*diff[0] + diff[1]*diff[1] + diff[2]*diff[2]).sqrt()
        diam = dist + self.r + other.r
        off  = diam/2 - self.r
        return BoundingSphere(self.cen + (off/dist)*diff, diam/2)
    def transform(self, T):
        return BoundingSphere(T.transform(self.cen), self.r * T.max_scale())


class RenderParams(SageObject):
    """
    This class is a container for all parameters that may be
    needed to render triangulate/render an object to a certain
    format. It can contain both cumulative and global parameters.
    """
    def __init__(self, **kwds):
        self._uniq_counter = 0
        self.obj_vertex_offset = 1
        self.transform_list = []
        self.transform = None
        self.ds = 1
        self.__dict__.update(kwds)

    def push_transform(self, T):
        self.transform_list.append(self.transform)
        if self.transform is None:
            self.transform = T
        else:
            self.transform = self.transform * T

    def pop_transform(self):
        self.transform = self.transform_list.pop()

    def unique_name(self, desc="name"):
        self._uniq_counter += 1
        return "%s_%s" % (desc, self._uniq_counter)

def flatten_list(L):
    """
    This is an optimized routine to turn a list of lists (of lists ...) into a single
    list. We generate data in a non-flat format to avoid multiple data copying, and
    then concatinate it all at the end.

    This is NOT recursive, otherwise there would be a lot of redundant copying (which
    we are trying to avoid in the first place, though at least it would be just the
    pointers).
    """
    if not PyList_CheckExact(L):
        return [L]
    flat = []
    L_stack = []; L_pop = L_stack.pop
    i_stack = []; i_pop = i_stack.pop
    cdef Py_ssize_t i = 0
    while i < PyList_GET_SIZE(L) or PyList_GET_SIZE(L_stack) > 0:
        while i < PyList_GET_SIZE(L):
            tmp = <object>PyList_GET_ITEM(L, i)
            if PyList_CheckExact(tmp):
                PyList_Append(L_stack, L)
                L = tmp
                PyList_Append(i_stack, i)
                i = 0
            else:
                PyList_Append(flat, tmp)
                i += 1
        if PyList_GET_SIZE(L_stack) > 0:
            L = L_pop()
            i = i_pop()
            i += 1
    return flat
