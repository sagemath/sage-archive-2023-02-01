r"""
Base classes for 3D Graphics objects and plotting.

AUTHORS:

- Robert Bradshaw (2007-02): inital version

- Robert Bradshaw (2007-08): Cythonization, much optimization

- William Stein (2008)

TODO: - finish integrating tachyon - good default lights, camera
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

from texture import Texture, is_Texture
from transform cimport Transformation, point_c, face_c
include "point_c.pxi"

from sage.interfaces.tachyon import tachyon_rt

from sage.plot.plot import show_default


default_texture = Texture()
pi = RDF.pi()

cdef class Graphics3d(SageObject):

    def __repr__(self):
        if show_default():
            self.show()
            return ''
        else:
            return self.__str__()

    def __str__(self):
        return "Graphics3d Object"

    def __add__(left, right):
        # Use == not "other is 0" here, since e.g., Sage integer zero is not 0.
        if right == 0 or right is None:
            return left
        elif left == 0 or left is None:
            return right
        elif not isinstance(left, Graphics3d):
            left = left.plot3d()
        elif not isinstance(right, Graphics3d):
            right = right.plot3d()
        return Graphics3dGroup([left, right])

    def _set_extra_kwds(self,kwds):
        self._extra_kwds = kwds

    def aspect_ratio(self, v=None):
        if not v is None:
            if not isinstance(v, (tuple, list)):
                raise TypeError, "v must be a list or tuple of length 3"
            self._aspect_ratio = [float(a) for a in v]
        else:
            if self._aspect_ratio is None:
                self._aspect_ratio = [1.0,1.0,1.0]
            return self._aspect_ratio

    def frame_aspect_ratio(self, v=None):
        if not v is None:
            self._frame_aspect_ratio = v
        else:
            if self._frame_aspect_ratio is None:
                self._frame_aspect_ratio = [1,1,1]
            return self._frame_aspect_ratio

    def _determine_frame_aspect_ratio(self, aspect_ratio):
        a_min, a_max = self._safe_bounding_box()
        return [(a_max[i] - a_min[i])*aspect_ratio[i] for i in range(3)]

    def _safe_bounding_box(self):
        # bounding box but where no side length is 0
        a_min, a_max = self.bounding_box()
        a_min = list(a_min); a_max = list(a_max)
        for i in range(3):
            if a_min[i] == a_max[i]:
                a_min[i] = a_min[i] - 1
                a_max[i] = a_max[i] + 1
        return a_min, a_max


    def bounding_box(self):
        # default
        return ((0.0, 0.0, 0.0), (0.0, 0.0, 0.0))

    def transform(self, **kwds):
        return TransformGroup([self], **kwds)

    def translate(self, *x):
        if len(x)==1:
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

    def testing_render_params(self):
        params = RenderParams(ds=.075)
        params.output_archive = zipfile.ZipFile('/dev/null', 'w', zipfile.ZIP_STORED, True)
        return params

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
                    mesh = False, dots=False,
                    perspective_depth = True,
                    orientation = (-764,-346,-545,76.39), **ignored_kwds):
                    # orientation chosen to look same as tachyon
        render_params = self.default_render_params()
        render_params.output_file = filename
        render_params.force_reload = render_params.randomize_counter = force_reload
        render_params.output_archive = zipfile.ZipFile(filename, 'w', zipfile.ZIP_DEFLATED, True)
        # Render the data
        all = flatten_list([self.jmol_repr(render_params), ""])

        f = StringIO()

        if len(render_params.atom_list):
            # Load the atom model
            f.write('data "model list"\n')
            f.write('%s\nempty\n' % (len(render_params.atom_list) + 1))
            for atom in render_params.atom_list:
                f.write('Xx %s %s %s\n' % atom)
            f.write('Xx 5.5 5.5 5.5\n') # so the zoom fits the box
            f.write('end "model list"; show data\n')
            f.write('select *\n')
            f.write('wireframe off; spacefill off\n')
            f.write('set labelOffset 0 0\n')


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

        f.write('centerAt absolute {0 0 0}\n')
        f.write('zoom %s\n'%zoom)
        f.write('frank OFF\n') # jmol logo

        if perspective_depth:
            f.write('set perspectivedepth ON\n')
        else:
            f.write('set perspectivedepth OFF\n')

        # Put the rest of the object in
        f.write("\n".join(all))

        render_params.output_archive.writestr('SCRIPT', f.getvalue())
        render_params.output_archive.close()

    def jmol_repr(self, render_params):
        return []

    def tachyon_repr(self, render_params):
        return []

    def obj_repr(self, render_params):
        return []

    def texture_set(self):
        return set()

    def mtl_str(self):
        return "\n\n".join([t.mtl_str() for t in self.texture_set()]) + "\n"

    def flatten(self, T=None):
        if T is None:
            return self
        else:
            return self.transform(T=T)

    def _rescale_for_frame_aspect_ratio_and_zoom(self, b, frame_aspect_ratio, zoom):
        if frame_aspect_ratio is None:
            return (b*zoom,b*zoom,b*zoom), (-b*zoom,-b*zoom,-b*zoom)
        box = [b*w for w in frame_aspect_ratio]
        # Now take the maximum length in box and rescale to b.
        s = b / max(box)
        box_max = tuple([s*w*zoom for w in box])
        box_min = tuple([-w*zoom for w in box_max])
        return box_min, box_max

    def _prepare_for_jmol(self, frame, axes, frame_aspect_ratio, aspect_ratio, zoom):
        from sage.plot.plot import EMBEDDED_MODE
        if EMBEDDED_MODE:
            s = 6
        else:
            s = 3
        box_min, box_max = self._rescale_for_frame_aspect_ratio_and_zoom(s, frame_aspect_ratio, zoom)
        a_min, a_max = self._box_for_aspect_ratio(aspect_ratio, box_min, box_max)
        return self._transform_to_bounding_box(box_min, box_max, a_min, a_max, frame=frame,
                                               axes=axes, thickness=1,
                                               labels = True)   # jmol labels are implemented

    def _prepare_for_tachyon(self, frame, axes, frame_aspect_ratio, aspect_ratio, zoom):
        box_min, box_max = self._rescale_for_frame_aspect_ratio_and_zoom(1.0, frame_aspect_ratio, zoom)
        a_min, a_max = self._box_for_aspect_ratio(aspect_ratio, box_min, box_max)
        return self._transform_to_bounding_box(box_min, box_max, a_min, a_max,
                                               frame=frame, axes=axes, thickness=0.5,
                                               labels = False)  # no tachyon text implemented yet

    def _box_for_aspect_ratio(self, aspect_ratio, box_min, box_max):
        # 1. Find a box around self so that when self gets rescaled into the
        # box defined by box_min, box_max, it has the right aspect ratio
        a_min, a_max = self._safe_bounding_box()

        if aspect_ratio == "automatic":
            return a_min, a_max

        longest_side = 0; longest_length = a_max[0] - a_min[0]
        shortest_side = 0; shortest_length = a_max[0] - a_min[0]

        for i in range(3):
            s = a_max[i] - a_min[i]
            if s > longest_length:
                longest_length = s
                longest_side = i
            if s < shortest_length:
                shortest_length = s
                shortest_side = i

        # 2. Rescale aspect_ratio so the shortest side is 1.
        r = float(aspect_ratio[shortest_side])
        aspect_ratio = [a/r for a in aspect_ratio]

        # 3. Extend the bounding box of self by rescaling so the sides
        # have the same ratio as aspect_ratio, and without changing
        # the longest side.
        long_box_side = box_max[longest_side] - box_min[longest_side]
        sc = [1.0,1.0,1.0]
        for i in range(3):
            # compute the length we want:
            new_length = longest_length / aspect_ratio[i]
            # change the side length by a_min and a_max so
            # that a_max[i] - a_min[i] = new_length

            # We have to take into account the ratio of the
            # sides after transforming to the bounding box.
            z = long_box_side / (box_max[i] - box_min[i])
            w = new_length / ((a_max[i] - a_min[i]) * z)
            sc[i] = w

        w = min(sc)
        sc = [z/w for z in sc]
        for i in range(3):
            a_min[i] *= sc[i]
            a_max[i] *= sc[i]

        return a_min, a_max

    def _transform_to_bounding_box(self, xyz_min, xyz_max, a_min, a_max, frame, axes, thickness, labels):

        a_min_orig = a_min; a_max_orig = a_max

        # Rescale in each direction
        scale = [float(xyz_max[i] - xyz_min[i]) / (a_max[i] - a_min[i]) for i in range(3)]
        X = self.scale(scale)
        a_min = [scale[i]*a_min[i] for i in range(3)]
        a_max = [scale[i]*a_max[i] for i in range(3)]

        # Translate so lower left corner of original bounding box
        # is in the right spot
        T = [xyz_min[i] - a_min[i] for i in range(3)]
        X = X.translate(T)
        if frame:
            from shapes2 import frame3d, frame_labels
            F = frame3d(xyz_min, xyz_max, opacity=0.5, color=(0,0,0), thickness=thickness)
            if labels:
                F_text = frame_labels(xyz_min, xyz_max, a_min_orig, a_max_orig)

            X += F + F_text

        if axes:
            # draw axes
            from shapes import arrow3d
            A = (arrow3d((min(0,a_min[0]),0, 0), (max(0,a_max[0]), 0,0),
                             thickness, color="blue"),
                 arrow3d((0,min(0,a_min[1]), 0), (0, max(0,a_max[1]), 0),
                             thickness, color="blue"),
                 arrow3d((0, 0, min(0,a_min[2])), (0, 0, max(0,a_max[2])),
                             thickness, color="blue"))
            X += sum(A).translate([-z for z in T])

        return X

    def show(self, **kwds):
        """
        INPUT:


        -  ``viewer`` - string (default: 'jmol'), how to view
           the plot 'jmol': interactive 3d (java) 'tachyon': a static png
           image (ray traced) 'java3d': interactive opengl based 3d

        -  ``filename`` - string (default: a temp file); file
           to save the image to

        -  ``verbosity`` - display information about rendering
           the figure

        -  ``figsize`` - (default: 5); x or pair [x,y] for
           numbers, e.g., [5,5]; controls the size of the output figure. E.g.,
           with Tachyon the number of pixels in each direction is 100 times
           figsize[0]. This is ignored for the jmol embedded renderer.

        -  ``aspect_ratio`` - (default: "automatic") - aspect
           ratio of the coordinate system itself. Give [1,1,1] to make spheres
           look round.

        -  ``frame_aspect_ratio`` - (default: "automatic")
           aspect ratio of frame that contains the 3d scene.

        -  ``zoom`` - (default: 1) how zoomed in

        -  ``frame`` - (default: True) if True, draw a
           bounding frame with labels

        -  ``axes`` - (deault: False) if True, draw coordinate
           axes


        -  ``**kwds`` - other options, which make sense for particular
           rendering engines

        CHANGING DEFAULTS: Defaults can be uniformly changed by importing a
        dictionary and changing it. For example, here we change the default
        so images display without a frame instead of with one::

            sage: from sage.plot.plot3d.base import SHOW_DEFAULTS
            sage: SHOW_DEFAULTS['frame'] = False

        This sphere will not have a frame around it::

            sage: sphere((0,0,0))

        We change the default back::

            sage: SHOW_DEFAULTS['frame'] = True

        Now this sphere is enclosed in a frame::

            sage: sphere((0,0,0))

        EXAMPLES: We illustrate use of the aspect_ratio option::

            sage: x, y = var('x,y')
            sage: p = plot3d(2*sin(x*y), (x, -pi, pi), (y, -pi, pi))
            sage: p.show(aspect_ratio=[1,1,1])

        This looks flattened, but filled with the plot::

            sage: p.show(frame_aspect_ratio=[1,1,1/16])

        This looks flattened, but the plot is square and smaller::

            sage: p.show(aspect_ratio=[1,1,1], frame_aspect_ratio=[1,1,1/8])
        """
        ek = self._extra_kwds
        if ek is not None:
            for key in ek.keys():
                if not kwds.has_key(key):
                    kwds[key] = ek[key]

        for key in SHOW_DEFAULTS.keys():
            if not kwds.has_key(key):
                kwds[key] = SHOW_DEFAULTS[key]

        # must have one line for every named argument:
        if kwds.has_key('viewer'): viewer = kwds['viewer']; del kwds['viewer']
        if kwds.has_key('filename'): filename = kwds['filename']; del kwds['filename']
        if kwds.has_key('verbosity'): verbosity = kwds['verbosity']; del kwds['verbosity']
        if kwds.has_key('figsize'): figsize = kwds['figsize']; del kwds['figsize']
        if kwds.has_key('aspect_ratio'): aspect_ratio = kwds['aspect_ratio']; del kwds['aspect_ratio']
        if kwds.has_key('frame_aspect_ratio'): frame_aspect_ratio = kwds['frame_aspect_ratio']; del kwds['frame_aspect_ratio']
        if kwds.has_key('zoom'): zoom = kwds['zoom']; del kwds['zoom']
        if kwds.has_key('frame'): frame = kwds['frame']; del kwds['frame']
        if kwds.has_key('axes'): axes = kwds['axes']; del kwds['axes']

        if aspect_ratio == 1:
            aspect_ratio = (1, 1, 1)
        if not isinstance(aspect_ratio, (str, list, tuple)):
            raise TypeError, "aspect ratio must be a string, list, tuple, or 1"
        if aspect_ratio != "automatic" and frame_aspect_ratio == "automatic":
            # set the aspect_ratio of the frame to be the same as that of the
            # object we are rendering given the aspect_ratio we'll use for it.
            frame_aspect_ratio = self._determine_frame_aspect_ratio(aspect_ratio)
        elif frame_aspect_ratio == "automatic":
            frame_aspect_ratio = self.frame_aspect_ratio()

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
            T = self._prepare_for_tachyon(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
            tachyon_rt(T.tachyon(), filename+".png", verbosity, True, opts)
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
            if EMBEDDED_MODE:
                # jmol doesn't seem to correctly parse the ?params part of a URL
                archive_name = "%s-%s.%s.zip" % (filename, randint(0, 1 << 30), ext)

            T = self._prepare_for_jmol(frame, axes, frame_aspect_ratio, aspect_ratio, zoom)
            T.export_jmol(archive_name, force_reload=EMBEDDED_MODE, zoom=zoom*100, **kwds)
            viewer_app = "sage-native-execute " + sage.misc.misc.SAGE_LOCAL + "/java/jmol/jmol"

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

# if you add any default parameters you must update some code below
SHOW_DEFAULTS = {'viewer':'jmol',
                 'verbosity':0,
                 'figsize':5,
                 'aspect_ratio':"automatic",
                 'frame_aspect_ratio':"automatic",
                 'zoom':1,
                 'frame':True,
                 'axes':False}




class Graphics3dGroup(Graphics3d):
    def __init__(self, all=[]):
        self.all = all
        self.frame_aspect_ratio(optimal_aspect_ratios([a.frame_aspect_ratio() for a in all]))
        self.aspect_ratio(optimal_aspect_ratios([a.aspect_ratio() for a in all]))
        self._set_extra_kwds(optimal_extra_kwds([a._extra_kwds for a in all if a._extra_kwds is not None]))

    def bounding_box(self):
        # Box that contains the bounding boxes of
        # all the objects that make up self.
        v = [obj.bounding_box() for obj in self.all]
        return min3([a[0] for a in v]), max3([a[1] for a in v])

    def transform(self, **kwds):
        T = TransformGroup(self.all, **kwds)
        T._set_extra_kwds(self._extra_kwds)
        return T

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
        self.frame_aspect_ratio(optimal_aspect_ratios([a.frame_aspect_ratio() for a in all]))
        self.aspect_ratio(optimal_aspect_ratios([a.aspect_ratio() for a in all]))
        self._set_extra_kwds(optimal_extra_kwds([a._extra_kwds for a in all if a._extra_kwds is not None]))

    def bounding_box(self):
        try:
            return self._bounding_box
        except AttributeError:
            pass

        cdef Transformation T = self.get_transformation()
        w = sum([T.transform_bounding_box(obj.bounding_box()) for obj in self.all], ())
        self._bounding_box = point_list_bounding_box(w)
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
        if kwds.has_key('texture'):
            self.texture = kwds['texture']
            if not is_Texture(self.texture):
                self.texture = Texture(self.texture)
        else:
            self.texture = Texture(kwds)

    def set_texture(self, texture, **kwds):
        if not is_Texture(texture):
            texture = Texture(texture, **kwds)
        self.texture = texture

    def get_texture(self):
        return self.texture

    def x3d_str(self):
        return "<Shape>" + self.x3d_geometry() + self.texture.x3d_str() + "</Shape>\n"

    def tachyon_repr(self, render_params):
        return self.triangulation().tachyon_repr(render_params)

    def obj_repr(self, render_params):
        return self.triangulation().obj_repr(render_params)

    def jmol_repr(self, render_params):
        return self.triangulation().jmol_repr(render_params)

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
    This class is a container for all parameters that may be needed to
    render triangulate/render an object to a certain format. It can
    contain both cumulative and global parameters.
    """

    _uniq_counter = 0
    randomize_counter = 0
    force_reload = False
    mesh = False
    dots = False

    def __init__(self, **kwds):
        self.output_file = sage.misc.misc.tmp_filename()
        self.obj_vertex_offset = 1
        self.transform_list = []
        self.transform = None
        self.ds = 1
        self.crease_threshold = .8
        self.__dict__.update(kwds)
        # for jmol, some things (such as labels) must be attached to atoms
        self.atom_list = []

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
    This is an optimized routine to turn a list of lists (of lists ...)
    into a single list. We generate data in a non-flat format to avoid
    multiple data copying, and then concatenate it all at the end.

    This is NOT recursive, otherwise there would be a lot of redundant
    copying (which we are trying to avoid in the first place, though at
    least it would be just the pointers).
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

    EXAMPLES::

        sage: from sage.plot.plot3d.base import min3, max3
        sage: min3([(-1,2,5), (-3, 4, 2)])
        (-3, 2, 2)
    """
    return tuple([min([a[i] for a in v]) for i in range(3)])

def max3(v):
    """
    Return the componentwise maximum of a list of 3-tuples.

    EXAMPLES::

        sage: from sage.plot.plot3d.base import min3, max3
        sage: max3([(-1,2,5), (-3, 4, 2)])
        (-1, 4, 5)
    """
    return tuple([max([a[i] for a in v]) for i in range(3)])

def point_list_bounding_box(v):
    """
    EXAMPLES::

        sage: from sage.plot.plot3d.base import point_list_bounding_box
        sage: point_list_bounding_box([(1,2,3),(4,5,6),(-10,0,10)])
        ((-10.0, 0.0, 3.0), (4.0, 5.0, 10.0))
    """
    cdef point_c lower, upper, cur
    cur.x, cur.y, cur.z = v[0]
    upper = lower = cur
    for P in v:
        cur.x, cur.y, cur.z = P
        point_c_lower_bound(&lower, lower, cur)
        point_c_upper_bound(&upper, upper, cur)
    return (lower.x, lower.y, lower.z), (upper.x, upper.y, upper.z)

def optimal_aspect_ratios(ratios):
    # average the aspect ratios
    n = len(ratios)
    if n > 0:
        return [max([z[i] for z in ratios]) for i in range(3)]
    else:
        return [1.0,1.0,1.0]

def optimal_extra_kwds(v):
    """
    Given a list v of dictionaries, this function merges them such that
    later dictionaries have precedence.
    """
    if len(v) == 0:
        return {}
    a = dict(v[0])   # make a copy!
    for b in v[1:]:
        for k,w in b.iteritems():
            a[k] = w
    return a

