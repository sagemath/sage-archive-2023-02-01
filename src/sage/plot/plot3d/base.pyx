r"""
Base classes for 3D Graphics objects and plotting.

EXAMPLES:
    sage: from sage.plot.graphics3d import *
    sage: S = ColorCube(.35, ['green', 'yellow', 'blue']) + Sphere(.2, color='red').translate(.4,.4,.4)
    sage: S.show()

AUTHOR:
    -- Robert Bradshaw

TODO:
   -- finish integrating tachyon
   -- good default lights, camera
"""
from random import randint
from math import atan2


from sage.rings.real_double import RDF
from sage.misc.functional import sqrt, atan, acos
#from sage.functions.all import *
from texture import Texture, is_Texture
from transform import Transformation
pi = RDF.pi()

from sage.interfaces.tachyon import tachyon_rt


default_texture = Texture()


cdef class Graphics3d_new(SageObject):


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
   self.tachyon_str(Transformation(scale=[1,-1,1])))  # switch from LH to RH coords to be consistant with java rendition

    def viewpoint(self):
        return Viewpoint(0,0,6)

#    def show(self, **kwds):
#        f = open("test.x3d", "w")
#        f.write(self.x3d())
#        f.close()

    def mtl_str(self):
        return "\n\n".join([t.mtl_str() for t in self.texture_set()]) + "\n"

    def texture_set(self):
        return set()

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
            f.write(self.obj_str(None))
            f.close()
            f = open(filename+".mtl", "w")
            f.write(self.mtl_str())
            f.close()


class Graphics3dGroup(Graphics3d_new):
    def __init__(self, all=[]):
        self.all = all

    def transform(self, **kwds):
        return TransformGroup(self.all, **kwds)

    def tachyon_str(self, transform):
        return "\n".join([g.tachyon_str(transform) for g in self.all])

    def x3d_str(self):
        return "\n".join([g.x3d_str() for g in self.all])

    def obj_str(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        return "\n\n".join([g.obj_str(transform, point_list) for g in self.all]) + "\n"

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

    def tachyon_str(self, transform):
        if transform is None:
            composite_transform = self.get_transformation()
        else:
            composite_transform = transform * self.get_transformation()
        return "\n".join([g.tachyon_str(composite_transform) for g in self.all])

    def obj_str(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        if transform is None:
            composite_transform = self.get_transformation()
        else:
            composite_transform = transform * self.get_transformation()
        return "\n\n".join([g.obj_str(composite_transform, point_list) for g in self.all])

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




class Viewpoint(Graphics3d_new):

    def __init__(self, *x):
        if isinstance(x[0], (tuple, list)):
            x = tuple(x[0])
        self.pos = x

    def x3d_str(self):
        return "<Viewpoint position='%s %s %s'/>"%self.pos



cdef class PrimativeObject(Graphics3d_new):
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

    def x3d_str(self):
        return "<Shape>" + self.x3d_geometry() + self.texture.x3d_str() + "</Shape>\n"

    def tachyon_str(self, transform):
        try:
            return self.tachyon_geometry(transform) + "\n" + self.texture.tachyon_str()
        except AttributeError:
            return self.triangulation().tachyon_str(transform)

    def tachyon_geometry(self, transform):
        return self.triangulation().tachyon_geometry(transform)

    def obj_str(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        return "g obj%s\n\nusemtl "%randint(0,10000000) + self.texture.id + "\n" + self.obj_geometry(transform, point_list)

    def obj_geometry(self, transform, point_list=None):
        if point_list is None:
            point_list = []
        return self.triangulation().obj_geometry(transform, point_list)

    def texture_set(self):
        return set([self.texture])
