r"""
Base classes for 3D Graphics objects and plotting.

EXAMPLES:
    sage: from sage.plot.plot3d.shapes import *
    sage: from sage.plot.plot3d.plot3d import plot3d
    sage: S = Sphere(.5, color='yellow')
    sage: S += Cone(.5, .5, color='red').translate(0,0,.3)
    sage: S += Sphere(.1, color='white').translate(.45,-.1,.15) + Sphere(.05, color='black').translate(.51,-.1,.17)
    sage: S += Sphere(.1, color='white').translate(.45, .1,.15) + Sphere(.05, color='black').translate(.51, .1,.17)
    sage: S += Sphere(.1, color='yellow').translate(.5, 0, -.2)
    sage: def f(x,y): return math.exp(x/5)*math.cos(y)
    sage: P = plot3d(f,(-5,5),(-5,5), ['red','yellow'], max_depth=10)
    sage: cape_man = P.scale(.2)+S.translate(1,0,0)
    sage: cape_man.show()


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


import os
from math import atan2
from random import randint
import zipfile
from cStringIO import StringIO

import sage.misc.misc

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
        # Use == not "other is 0" here, since e.g., Sage integer zero is not 0.
        if other == 0 or other is None:
            return self
        elif self == 0 or self is None:
            return other
        return Graphics3dGroup([self, other])

    def bounding_box(self):
        raise NotImplementedError

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
            center  2.3 2.4 2.0
            viewdir  -2.3 -2.4 -2.0
            updir  0.0 0.0 1.0
         end_camera


      light center  4.0 3.0 2.0
            rad 0.2
            color  1.0 1.0 1.0

      plane
        center -2000 -1000 -500
        normal 2.3 2.4 2.0
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

    def export_jmol(self, filename='jmol_shape.jmol', force_reload=False,
                    zoom=100, spin=False, background=(1,1,1), stereo=False,
                    perspective_depth = True,
                    orientation = (-764,-346,-545,76.39)):
                    # orientation chosen to look same as tachyon
        render_params = self.default_render_params()
        render_params.output_file = filename
        render_params.force_reload = render_params.randomize_counter = force_reload
        render_params.output_archive = zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED, True)

        f = StringIO()

        # Set the scene background color
        f.write('background [%s,%s,%s]\n'%tuple([int(a*255) for a in background]))
        if spin:
            f.write('spin ON\n')
        else:
            f.write('spin OFF\n')
        if stereo:
            if stereo is True: stereo = "redblue"
            f.write('stereo %s\n' % stereo)
        if orientation:
            f.write('moveto 0 %s %s %s %s\n'%tuple(orientation))

        f.write('zoom %s\n'%zoom)

        if perspective_depth:
            f.write('set perspectivedepth ON\n')
        else:
            f.write('set perspectivedepth OFF\n')

        # Put the rest of the object in
        f.write("\n".join(flatten_list([self.jmol_repr(render_params), ""])))

        render_params.output_archive.writestr('SCRIPT', f.getvalue())
        render_params.output_archive.close()

    def jmol_repr(self, render_params):
        raise NotImplementedError

    def texture_set(self):
        return set()

    def mtl_str(self):
        return "\n\n".join([t.mtl_str() for t in self.texture_set()]) + "\n"

    def flatten(self, T=None):
        if T is None:
            return self
        else:
            return self.transform(T=T)

    def _rescale_for_aspect_ratio_and_zoom(self, b, aspect_ratio, zoom):
        if aspect_ratio is None:
            return (b*zoom,b*zoom,b*zoom), (-b*zoom,-b*zoom,-b*zoom)
        box = [b*w for w in aspect_ratio]
        # Now take the maximum length in box and rescale to b.
        s = b / max(box)
        box_max = tuple([s*w*zoom for w in box])
        box_min = tuple([-w*zoom for w in box_max])
        return box_min, box_max

    def _prepare_for_jmol(self, frame, axes, aspect_ratio, zoom):
        box_min, box_max = self._rescale_for_aspect_ratio_and_zoom(6.0, aspect_ratio, zoom)
        return self._transform_to_bounding_box(box_min, box_max, frame=frame,
                                               axes=axes, thickness=1)

    def _prepare_for_tachyon(self, frame, axes, aspect_ratio, zoom):
        box_min, box_max = self._rescale_for_aspect_ratio_and_zoom(1.0, aspect_ratio, zoom)
        A = self._transform_to_bounding_box(box_min, box_max,
                                            frame=frame, axes=axes, thickness=0.5)
        return A

    def _transform_to_bounding_box(self, xyz_min, xyz_max, frame, axes, thickness):
        a_min, a_max = self.bounding_box()
        a_min = list(a_min); a_max = list(a_max)
        for i in range(3):
            if a_min[i] == a_max[i]:
                a_min[i] = -1
                a_max[i] = 1

        # Rescale in each direction
        scale = [(xyz_max[i] - xyz_min[i]) / (a_max[i] - a_min[i]) for i in range(3)]
        X = self.scale(scale)
        a_min = [scale[i]*a_min[i] for i in range(3)]
        a_max = [scale[i]*a_max[i] for i in range(3)]

        # Translate so lower left corner of original bounding box
        # is in the right spot
        T = [xyz_min[i] - a_min[i] for i in range(3)]
        X = X.translate(T)
        if frame:
            from shapes2 import frame3d
            F = frame3d(xyz_min, xyz_max, opacity=0.5, color=(0,0,0), thickness=thickness)
            X += F

        if axes:
            # draw axes
            from shapes import Arrow
            A = (Arrow((min(0,a_min[0]),0, 0), (max(0,a_max[0]), 0,0),
                             thickness, color="blue"),
                 Arrow((0,min(0,a_min[1]), 0), (0, max(0,a_max[1]), 0),
                             thickness, color="blue"),
                 Arrow((0, 0, min(0,a_min[2])), (0, 0, max(0,a_max[2])),
                             thickness, color="blue"))
            X += sum(A).translate([-z for z in T])

        return X

    def show(self, viewer="jmol", filename=None, verbosity=0, figsize=5,
             aspect_ratio = None, zoom=1,
             frame=True, axes = False, **kwds):
        """
        INPUT:
            viewer -- string (default: 'jmol') which viewing system to use.
                      'jmol': an embedded non-OpenGL 3d java applet
                      'tachyon': an embedded ray tracer
                      'java3d': a popup OpenGL 3d java applet
            filename -- string (default: a temp file); file to save the image to
            verbosity -- display information about rendering the figure
            figsize -- (default: 5); x or pair [x,y] for numbers, e.g., [5,5]; controls
                       the size of the output figure.  E.g., with Tachyon the number of
                       pixels in each direction is 100 times figsize[0].
                       This is ignored for the jmol embedded renderer.
            **kwds -- other options, which make sense for particular rendering engines
        """
        import sage.misc.misc
        if filename is None:
            filename = sage.misc.misc.tmp_filename()
        if not isinstance(figsize, (list,tuple)):
            figsize = [figsize, figsize]
        from sage.plot.plot import EMBEDDED_MODE, DOCTEST_MODE
        ext = None

        # Tachyon resolution options
        if DOCTEST_MODE:
            opts = '-res 10 10'
            filename = sage.misc.misc.SAGE_TMP + "/tmp"
        elif EMBEDDED_MODE:
            opts = '-res %s %s'%(figsize[0]*100, figsize[1]*100)
            filename = sage.misc.misc.graphics_filename()[:-4]
        else:
            opts = '-res %s %s'%(figsize[0]*100, figsize[1]*100)

        if DOCTEST_MODE or viewer=='tachyon' or (viewer=='java3d' and EMBEDDED_MODE):
            T = self._prepare_for_tachyon(frame, axes, aspect_ratio, zoom)
            tachyon_rt(T.tachyon(**kwds), filename+".png", verbosity, True, opts)
            ext = "png"
            import sage.misc.viewer
            viewer_app = sage.misc.viewer.browser()

        if DOCTEST_MODE or viewer=='java3d':
            f = open(filename+".obj", "w")
            f.write("mtllib %s.mtl\n" % filename)
            f.write(self.obj())
            f.close()
            f = open(filename+".mtl", "w")
            f.write(self.mtl_str())
            f.close()
            ext = "obj"
            viewer_app = sage.misc.misc.SAGE_LOCAL + "/java/java3d/start_viewer"

        if DOCTEST_MODE or viewer=='jmol':
            # Temporary hack: encode the desired applet size in the end of the filename:
            # (This will be removed once we have dynamic resizing of applets in the browser.)
            base, ext = os.path.splitext(filename)
            fg = figsize[0]
            #if fg >= 2:
            #    fg = 2
            filename = '%s-size%s%s'%(base, fg*100, ext)
            ext = "jmol"
            archive_name = "%s.%s.zip" % (filename, ext)

            T = self._prepare_for_jmol(frame, axes, aspect_ratio, zoom)
            T.export_jmol(archive_name, force_reload=EMBEDDED_MODE, **kwds)
            viewer_app = sage.misc.misc.SAGE_LOCAL + "/java/jmol/jmol"

            # We need a script to load the file
            f = open(filename + '.jmol', 'w')
            f.write('set defaultdirectory "%s"\n' % archive_name)
            f.write('script SCRIPT\n')
            f.close()

        if ext is None:
            raise ValueError, "Unknown 3d plot type: %s" % viewer

        if not DOCTEST_MODE and not EMBEDDED_MODE:
            if verbosity:
                pipes = "2>&1"
            else:
                pipes = "2>/dev/null 1>/dev/null &"
            os.system('%s "%s.%s" %s' % (viewer_app, filename, ext, pipes))

class Graphics3dGroup(Graphics3d):
    def __init__(self, all=[]):
        self.all = all

    def bounding_box(self):
        # Box that contains the bounding boxes of
        # all the objects that make up self.
        v = [obj.bounding_box() for obj in self.all]
        return min3([a[0] for a in v]), max3([a[1] for a in v])

    def transform(self, **kwds):
        return TransformGroup(self.all, **kwds)

    def tachyon_repr(self, render_params):
        return [g.tachyon_repr(render_params) for g in self.all]

    def x3d_str(self):
        return "\n".join([g.x3d_str() for g in self.all])

    def obj_repr(self, render_params):
        return [g.obj_repr(render_params) for g in self.all]

    def jmol_repr(self, render_params):
        return [g.jmol_repr(render_params) for g in self.all]

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

    def bounding_box(self):
        try:
            return self._bounding_box
        except AttributeError:
            pass
        # Get the box before transformation
        a = Graphics3dGroup.bounding_box(self)
        # The corners of the box
        import sage.misc.mrange
        corners = []
        for f in sage.misc.mrange.cartesian_product_iterator([[0,1]]*3):
            corners.append([a[f[i]][i] for i in range(3)])
        # Transform the corners of the box
        T = self.get_transformation()
        w = [T.transform_point(p) for p in corners]
        # Figure out what the new bounding box is
        a_min = [min([z[i] for z in w]) for i in range(3)]
        a_max = [max([z[i] for z in w]) for i in range(3)]
        self._bounding_box = a_min, a_max
        return self._bounding_box

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

    def jmol_repr(self, render_params):
        render_params.push_transform(self.get_transformation())
        rep = [g.jmol_repr(render_params) for g in self.all]
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



cdef class PrimitiveObject(Graphics3d):
    def __init__(self, **kwds):
        try:
            self.texture = kwds['texture']
            if not is_Texture(self.texture):
                self.texture = Texture(self.texture)
        except KeyError:
            try:
                self.texture = kwds['color']
                if not is_Texture(self.texture):
                    if kwds.has_key('opacity'):
                        self.texture = Texture(self.texture, opacity=kwds['opacity'])
                    else:
                        self.texture = Texture(self.texture)
            except KeyError:
                self.texture = default_texture

    def set_texture(self, texture, **kwds):
        if not is_Texture(texture):
            texture = Texture(texture, **kwds)
        self.texture = texture

    def x3d_str(self):
        return "<Shape>" + self.x3d_geometry() + self.texture.x3d_str() + "</Shape>\n"

    def tachyon_repr(self, render_params):
        return self.triangulation().tachyon_repr(render_params)

    def obj_repr(self, render_params):
        return self.triangulation().obj_repr(render_params)

    def jmol_repr(self, render_params):
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
        # Use == not "other is 0" here, since e.g., Sage integer zero is not 0.
        if other == 0 or other is None:
            return self
        elif self == 0 or self is None:
            return other
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
        self.randomize_counter = 0
        self.output_file = sage.misc.misc.tmp_filename()
        self.obj_vertex_offset = 1
        self.transform_list = []
        self.transform = None
        self.ds = 1
        self.crease_threshold = .8
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
        if self.randomize_counter:
            self._uniq_counter = randint(1,1000000)
        else:
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


def min3(v):
    """
    Return the componentwise minimum of a list of 3-tuples.

    EXAMPLES:
        sage: from sage.plot.plot3d.base import min3, max3
        sage: min3([(-1,2,5), (-3, 4, 2)])
        (-3, 2, 2)
    """
    return tuple([min([a[i] for a in v]) for i in range(3)])

def max3(v):
    """
    Return the componentwise maximum of a list of 3-tuples.

    EXAMPLES:
        sage: from sage.plot.plot3d.base import min3, max3
        sage: max3([(-1,2,5), (-3, 4, 2)])
        (-1, 4, 5)
    """
    return tuple([max([a[i] for a in v]) for i in range(3)])
