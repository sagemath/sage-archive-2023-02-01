r"""
The Tachyon 3D Ray Tracer

Given any 3D graphics object one can compute a raytraced
representation by typing ``show(viewer='tachyon')``.
For example, we draw two translucent spheres that contain a red
tube, and render the result using Tachyon.

::

    sage: S = sphere(opacity=0.8, aspect_ratio=[1,1,1])
    sage: L = line3d([(0,0,0),(2,0,0)], thickness=10, color='red')
    sage: M = S + S.translate((2,0,0)) + L
    sage: M.show(viewer='tachyon')

One can also directly control Tachyon, which gives a huge amount of
flexibility. For example, here we directly use Tachyon to draw 3
spheres on the coordinate axes. Notice that the result is
gorgeous::

    sage: t = Tachyon(xres=500,yres=500, camera_center=(2,0,0))
    sage: t.light((4,3,2), 0.2, (1,1,1))
    sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
    sage: t.texture('t3', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
    sage: t.texture('t4', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
    sage: t.sphere((0,0.5,0), 0.2, 't2')
    sage: t.sphere((0.5,0,0), 0.2, 't3')
    sage: t.sphere((0,0,0.5), 0.2, 't4')
    sage: t.show()

AUTHOR:

- John E. Stone (johns@megapixel.com): wrote tachyon ray tracer

- William Stein: sage-tachyon interface

- Joshua Kantor: 3d function plotting

- Tom Boothby: 3d function plotting n'stuff

- Leif Hille: key idea for bugfix for texfunc issue (trac #799)

- Marshall Hampton: improved doctests, rings, axis-aligned boxes.

TODO:

- clean up trianglefactory stuff
"""

from tri_plot import Triangle, SmoothTriangle, TriangleFactory, TrianglePlot


from sage.interfaces.tachyon import tachyon_rt

from sage.structure.sage_object import SageObject

from sage.misc.misc import SAGE_TMP
from sage.misc.temporary_file import tmp_filename, graphics_filename

#from sage.ext import fast_tachyon_routines

import os

from math import sqrt

class Tachyon(SageObject):
    r"""
    Create a scene the can be rendered using the Tachyon ray tracer.

    INPUT:

    - ``xres`` - (default 350)
    - ``yres`` - (default 350)
    - ``zoom`` - (default 1.0)
    - ``antialiasing`` - (default False)
    - ``aspectratio``  - (default 1.0)
    - ``raydepth`` - (default 5)
    - ``camera_center`` - (default (-3, 0, 0))
    - ``updir`` - (default (0, 0, 1))
    - ``look_at`` - (default (0,0,0))
    - ``viewdir`` - (default None)
    - ``projection`` - (default 'PERSPECTIVE')

    OUTPUT: A Tachyon 3d scene.

    Note that the coordinates are by default such that `z` is
    up, positive `y` is to the {left} and `x` is toward
    you. This is not oriented according to the right hand rule.

    EXAMPLES: Spheres along the twisted cubic.

    ::

        sage: t = Tachyon(xres=512,yres=512, camera_center=(3,0.3,0))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2,diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in srange(-1,1,0.05):
        ....:    k += 1
        ....:    t.sphere((i,i^2-0.5,i^3), 0.1, 't%s'%(k%3))
        sage: t.show()

    Another twisted cubic, but with a white background, got by putting
    infinite planes around the scene.

    ::

        sage: t = Tachyon(xres=512,yres=512, camera_center=(3,0.3,0), raydepth=8)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2,diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: t.texture('white', color=(1,1,1))
        sage: t.plane((0,0,-1), (0,0,1), 'white')
        sage: t.plane((0,-20,0), (0,1,0), 'white')
        sage: t.plane((-20,0,0), (1,0,0), 'white')

    ::

        sage: k=0
        sage: for i in srange(-1,1,0.05):
        ....:    k += 1
        ....:    t.sphere((i,i^2 - 0.5,i^3), 0.1, 't%s'%(k%3))
        ....:    t.cylinder((0,0,0), (0,0,1), 0.05,'t1')
        sage: t.show()

    Many random spheres::

        sage: t = Tachyon(xres=512,yres=512, camera_center=(2,0.5,0.5), look_at=(0.5,0.5,0.5), raydepth=4)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2, diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in range(100):
        ....:    k += 1
        ....:    t.sphere((random(),random(), random()), random()/10, 't%s'%(k%3))
        sage: t.show()

    Points on an elliptic curve, their height indicated by their height
    above the axis::

        sage: t = Tachyon(camera_center=(5,2,2), look_at=(0,1,0))
        sage: t.light((10,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
        sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
        sage: E = EllipticCurve('37a')
        sage: P = E([0,0])
        sage: Q = P
        sage: n = 100
        sage: for i in range(n):   # increase 20 for a better plot
        ....:    Q = Q + P
        ....:    t.sphere((Q[1], Q[0], ZZ(i)/n), 0.1, 't%s'%(i%3))
        sage: t.show()

    A beautiful picture of rational points on a rank 1 elliptic curve.

    ::

        sage: t = Tachyon(xres=1000, yres=800, camera_center=(2,7,4), look_at=(2,0,0), raydepth=4)
        sage: t.light((10,3,2), 1, (1,1,1))
        sage: t.light((10,-3,2), 1, (1,1,1))
        sage: t.texture('black', color=(0,0,0))
        sage: t.texture('red', color=(1,0,0))
        sage: t.texture('grey', color=(.9,.9,.9))
        sage: t.plane((0,0,0),(0,0,1),'grey')
        sage: t.cylinder((0,0,0),(1,0,0),.01,'black')
        sage: t.cylinder((0,0,0),(0,1,0),.01,'black')
        sage: E = EllipticCurve('37a')
        sage: P = E([0,0])
        sage: Q = P
        sage: n = 100
        sage: for i in range(n):
        ....:    Q = Q + P
        ....:    c = i/n + .1
        ....:    t.texture('r%s'%i,color=(float(i/n),0,0))
        ....:    t.sphere((Q[0], -Q[1], .01), .04, 'r%s'%i)
        sage: t.show()    # long time, e.g., 10-20 seconds

    A beautiful spiral.

    ::

        sage: t = Tachyon(xres=800,yres=800, camera_center=(2,5,2), look_at=(2.5,0,0))
        sage: t.light((0,0,100), 1, (1,1,1))
        sage: t.texture('r', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
        sage: for i in srange(0,50,0.1):
        ....:    t.sphere((i/10,sin(i),cos(i)), 0.05, 'r')
        sage: t.texture('white', color=(1,1,1), opacity=1, specular=1, diffuse=1)
        sage: t.plane((0,0,-100), (0,0,-100), 'white')
        sage: t.show()
    """
    def __init__(self,
                 xres=350, yres=350,
                 zoom = 1.0,
                 antialiasing = False,
                 aspectratio = 1.0,
                 raydepth = 8,
                 camera_center = (-3, 0, 0),
                 updir = (0, 0, 1),
                 look_at = (0,0,0),
                 viewdir = None,
                 projection = 'PERSPECTIVE'):
        r"""
        Creates an instance of the Tachyon class.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t._xres
            350
        """
        self._xres = xres
        self._yres = yres
        self._zoom = zoom
        self._aspectratio = aspectratio
        self._antialiasing = antialiasing
        self._raydepth = raydepth
        self._camera_center = camera_center
        self._updir = updir
        self._projection = projection
        self._objects = []
        if viewdir is None:
            self._viewdir = [look_at[i] - camera_center[i] for i in range(3)]
        else:
            self._viewdir = viewdir



    def __repr__(self):
        r"""
        Returns the string representation of the Tachyon object,
        which is just the scene string input to tachyon.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,1), 1,(1,1,1))
            sage: q.texture('s')
            sage: q.sphere((0,0,0),1,'s')
            sage: q.__repr__()[-20:]
            '  \n        end_scene'
        """
        return self.str()

    def save_image(self, filename=None, *args, **kwds):
        r"""
        Save an image representation of self.  The image type is
        determined by the extension of the filename.  For example,
        this could be ``.png``, ``.jpg``, ``.gif``, ``.pdf``,
        ``.svg``.  Currently this is implemented by calling the
        :meth:`save` method of self, passing along all arguments and
        keywords.

        .. Note::

            Not all image types are necessarily implemented for all
            graphics types.  See :meth:`save` for more details.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,11), 1,(1,1,1))
            sage: q.texture('s')
            sage: q.sphere((0,-1,1),1,'s')
            sage: tempname = tmp_filename()
            sage: q.save_image(tempname)

        TESTS:

        :meth:`save_image` is used for generating animations::

            sage: def tw_cubic(t):
            ....:     q = Tachyon()
            ....:     q.light((1,1,11), 1,(1,1,1))
            ....:     q.texture('s')
            ....:     for i in srange(-1,t,0.05):
            ....:         q.sphere((i,i^2-0.5,i^3), 0.1, 's')
            ....:     return q

            sage: a = animate([tw_cubic(t) for t in srange(-1,1,.3)])
            sage: a
            Animation with 7 frames
            sage: a.show() # optional -- ImageMagick
        """
        self.save(filename, *args, **kwds)

    def save(self, filename='sage.png', verbose=0, block=True, extra_opts=''):
        r"""
        INPUT:


        -  ``filename`` - (default: 'sage.png') output
           filename; the extension of the filename determines the type.
           Supported types include:

        -  ``tga`` - 24-bit (uncompressed)

        -  ``bmp`` - 24-bit Windows BMP (uncompressed)

        -  ``ppm`` - 24-bit PPM (uncompressed)

        -  ``rgb`` - 24-bit SGI RGB (uncompressed)

        -  ``png`` - 24-bit PNG (compressed, lossless)

        -  ``verbose`` - integer; (default: 0)

        -  ``0`` - silent

        -  ``1`` - some output

        -  ``2`` - very verbose output

        -  ``block`` - bool (default: True); if False, run the
           rendering command in the background.

        -  ``extra_opts`` - passed directly to tachyon command
           line. Use tachyon_rt.usage() to see some of the possibilities.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,11), 1,(1,1,1))
            sage: q.texture('s')
            sage: q.sphere((0,0,0),1,'s')
            sage: tempname = tmp_filename()
            sage: q.save(tempname)
            sage: os.system('rm ' + tempname)
            0
        """
        tachyon_rt(self.str(), filename, verbose, block, extra_opts)

    def show(self, verbose=0, extra_opts=''):
        r"""
        Creates a PNG file of the scene.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((-1,-1,10), 1,(1,1,1))
            sage: q.texture('s')
            sage: q.sphere((0,0,0),1,'s')
            sage: q.show(verbose = False)
        """
        import sage.plot.plot
        if sage.doctest.DOCTEST_MODE:
            filename = graphics_filename()
            self.save(os.path.join(SAGE_TMP, 'test.png'), verbose=verbose, extra_opts=extra_opts)
            return
        if sage.plot.plot.EMBEDDED_MODE:
            filename = graphics_filename()
            self.save(filename, verbose=verbose, extra_opts=extra_opts)
            return
        filename = tmp_filename(ext='.png')
        self.save(filename, verbose=verbose, extra_opts=extra_opts)
        os.system('%s %s 2>/dev/null 1>/dev/null &'%(sage.misc.viewer.png_viewer(), filename))

    def _res(self):
        r"""
        An internal function that writes the tachyon string for the
        resolution (x and y size of the image).

        EXAMPLES::

            sage: t = Tachyon(xres = 300, yres = 700)
            sage: t._res()
            '\nresolution 300 700\n'
        """
        return '\nresolution %s %s\n'%(self._xres, self._yres)

    def _camera(self):
        r"""
        An internal function that writes the tachyon string for the
        camera and other rendering information (ray depth, antialiasing).

        EXAMPLES::

            sage: t = Tachyon(raydepth = 16, zoom = 2, antialiasing = True)
            sage: t._camera().split()[3:10]
            ['aspectratio', '1.0', 'antialiasing', '1', 'raydepth', '16', 'center']
        """
        return r"""
           camera
              zoom %s
              aspectratio %s
              antialiasing %s
              raydepth %s
              center %s
              viewdir %s
              updir %s
           end_camera
        """%(float(self._zoom), float(self._aspectratio),
             int(self._antialiasing),
             int(self._raydepth),
             tostr(self._camera_center),
             tostr(self._viewdir),
             tostr(self._updir))

    def str(self):
        r"""
        Returns the complete tachyon scene file as a string.

        EXAMPLES::

            sage: t = Tachyon(xres=500,yres=500, camera_center=(2,0,0))
            sage: t.light((4,3,2), 0.2, (1,1,1))
            sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
            sage: t.texture('t3', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
            sage: t.texture('t4', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
            sage: t.sphere((0,0.5,0), 0.2, 't2')
            sage: t.sphere((0.5,0,0), 0.2, 't3')
            sage: t.sphere((0,0,0.5), 0.2, 't4')
            sage: t.str().find('PLASTIC')
            567
        """
        return r"""
        begin_scene
        %s
        %s
        %s
        end_scene"""%(
            self._res(),
            self._camera(),
            '\n'.join([x.str() for x in self._objects])
            )

    def light(self, center, radius, color):
        r"""
        Creates a light source of the given center, radius, and color.

        EXAMPLES::

            sage: q = Tachyon()
            sage: q.light((1,1,1),1.0,(.2,0,.8))
            sage: q.str().split('\n')[17]
            '        light center  1.0 1.0 1.0 '
        """
        self._objects.append(Light(center, radius, color))

    def texfunc(self, type=0, center=(0,0,0), rotate=(0,0,0), scale=(1,1,1)):
        r"""
        INPUT:

        -  ``type`` - (default: 0)

           0. No special texture, plain shading
           1. 3D checkerboard function, like a rubik's cube
           2. Grit Texture, randomized surface color
           3. 3D marble texture, uses object's base color
           4. 3D wood texture, light and dark brown, not very good yet
           5. 3D gradient noise function (can't remember what it looks
              like)
           6. Don't remember
           7. Cylindrical Image Map, requires ppm filename (don't know
              how to specify name in sage?!)
           8. Spherical Image Map, requires ppm filename (don't know
              how to specify name in sage?!)
           9. Planar Image Map, requires ppm filename (don't know how
              to specify name in sage?!)

        -  ``center`` - (default: (0,0,0))
        -  ``rotate`` - (default: (0,0,0))
        -  ``scale`` - (default: (1,1,1))


        EXAMPLES: We draw an infinite checkboard::

            sage: t = Tachyon(camera_center=(2,7,4), look_at=(2,0,0))
            sage: t.texture('black', color=(0,0,0), texfunc=1)
            sage: t.plane((0,0,0),(0,0,1),'black')
            sage: t.show()
        """
        type = int(type)
        if type < 0 or type > 9:
            raise ValueError, "type must be an integer between 0 and 9"
        return Texfunc(type,center,rotate,scale).str()

    def texture(self, name, ambient=0.2, diffuse=0.8,
                specular=0.0, opacity=1.0,
                color=(1.0,0.0, 0.5), texfunc=0, phong=0, phongsize=.5, phongtype="PLASTIC"):
        r"""
        INPUT:


        -  ``name`` - string; the name of the texture (to be
           used later)

        -  ``ambient`` - (default: 0.2)

        -  ``diffuse`` - (default: 0.8)

        -  ``specular`` - (default: 0.0)

        -  ``opacity`` - (default: 1.0)

        -  ``color`` - (default: (1.0,0.0,0.5))

        -  ``texfunc`` - (default: 0); a texture function; this
           is either the output of self.texfunc, or a number between 0 and 9,
           inclusive. See the docs for self.texfunc.

        -  ``phong`` - (default: 0)

        -  ``phongsize`` - (default: 0.5)

        -  ``phongtype`` - (default: "PLASTIC")

        EXAMPLES:

        We draw a scene with 4 spheres that illustrates various uses of
        the texture command::

            sage: t = Tachyon(camera_center=(2,5,4), look_at=(2,0,0), raydepth=6)
            sage: t.light((10,3,4), 1, (1,1,1))
            sage: t.texture('mirror', ambient=0.05, diffuse=0.05, specular=.9, opacity=0.9, color=(.8,.8,.8))
            sage: t.texture('grey', color=(.8,.8,.8), texfunc=3)
            sage: t.plane((0,0,0),(0,0,1),'grey')
            sage: t.sphere((4,-1,1), 1, 'mirror')
            sage: t.sphere((0,-1,1), 1, 'mirror')
            sage: t.sphere((2,-1,1), 0.5, 'mirror')
            sage: t.sphere((2,1,1), 0.5, 'mirror')
            sage: show(t)  # known bug (:trac:`7232`)
        """
        if texfunc and not isinstance(texfunc, Texfunc):
            texfunc = self.texfunc(int(texfunc))
        self._objects.append(Texture(name, ambient, diffuse,
                                     specular, opacity, color, texfunc,
                                     phong,phongsize,phongtype))

    def texture_recolor(self, name, colors):
        r"""
        Recolors default textures.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('s')
            sage: q = t.texture_recolor('s',[(0,0,1)])
            sage: t._objects[1]._color
            (0, 0, 1)
        """
        base_tex = None
        names = []
        ident = "SAGETEX%d"%len(self._objects) #don't collide with other texture names

        for o in self._objects:
            if isinstance(o, Texture) and o._name == name:
                base_tex = o
                break
        if base_tex is None:
            base_tex = Texture(name)

        for i in range(len(colors)):
            n = "%s_%d"%(ident,i)
            self._objects.append(base_tex.recolor(n, colors[i]))
            names.append(n)

        return names

    def sphere(self, center, radius, texture):
        r"""
        Creates the scene information for a sphere with the given
        center, radius, and texture.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('sphere_texture')
            sage: t.sphere((1,2,3), .1, 'sphere_texture')
            sage: t._objects[1].str()
            '\n        sphere center  1.0 2.0 3.0  rad 0.1 sphere_texture\n        '
        """
        self._objects.append(Sphere(center, radius, texture))

    def ring(self, center, normal, inner, outer, texture):
        r"""
        Creates the scene information for a ring with the given parameters.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.ring([0,0,0], [0,0,1], 1.0, 2.0, 's')
            sage: t._objects[0]._center
            [0, 0, 0]
        """
        self._objects.append(Ring(center, normal, inner, outer, texture))

    def cylinder(self, center, axis, radius, texture):
        r"""
        Creates the scene information for a infinite cylinder with the
        given center, axis direction, radius, and texture.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('c')
            sage: t.cylinder((0,0,0),(-1,-1,-1),.1,'c')
        """
        self._objects.append(Cylinder(center, axis, radius, texture))

    def plane(self, center, normal, texture):
        r"""
        Creates an infinite plane with the given center and normal.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.plane((0,0,0),(1,1,1),'s')
            sage: t.str()[338:380]
            'plane center  0.0 0.0 0.0  normal  1.0 1.0'
        """
        self._objects.append(Plane(center, normal, texture))

    def axis_aligned_box(self, min_p, max_p, texture):
        r"""
        Creates an axis-aligned box with minimal point ``min_p`` and
        maximum point ``max_p``.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.axis_aligned_box((0,0,0),(2,2,2),'s')
        """
        self._objects.append(Axis_aligned_box(min_p, max_p, texture))

    def fcylinder(self, base, apex, radius, texture):
        r"""
        Finite cylinders are almost the same as infinite ones, but the
        center and length of the axis determine the extents of the
        cylinder.  The finite cylinder is also really a shell, it
        doesn't have any caps. If you need to close off the ends of
        the cylinder, use two ring objects, with the inner radius set
        to 0.0 and the normal set to be the axis of the cylinder.
        Finite cylinders are built this way to enhance speed.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.fcylinder((1,1,1),(1,2,3),.01,'s')
            sage: len(t.str())
            423
        """
        self._objects.append(FCylinder(base, apex, radius, texture))

    def triangle(self, vertex_1, vertex_2, vertex_3, texture):
        r"""
        Creates a triangle with the given vertices and texture.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('s')
            sage: t.triangle([1,2,3],[4,5,6],[7,8,10],'s')
            sage: t._objects[1]
            [1, 2, 3] [4, 5, 6] [7, 8, 10] s

        """
        self._objects.append(TachyonTriangle(vertex_1,vertex_2,vertex_3,texture))

    def smooth_triangle(self, vertex_1, vertex_2, vertex_3, normal_1, normal_2, normal_3, texture):
        r"""
        Creates a triangle along with a normal vector for smoothing.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.light((1,1,1),.1,(1,1,1))
            sage: t.texture('s')
            sage: t.smooth_triangle([0,0,0],[0,0,1],[0,1,0],[0,1,1],[-1,1,2],[3,0,0],'s')
            sage: t._objects[2]
            [0, 0, 0] [0, 0, 1] [0, 1, 0] s [0, 1, 1] [-1, 1, 2] [3, 0, 0]
        """
        self._objects.append(TachyonSmoothTriangle(vertex_1, vertex_2, vertex_3, normal_1, normal_2, normal_3, texture))

    def fractal_landscape(self, res, scale, center, texture):
        r"""
        Axis-aligned fractal landscape.  Not very useful at the moment.

        EXAMPLES::

            sage: t = Tachyon()
            sage: t.texture('s')
            sage: t.fractal_landscape([30,30],[80,80],[0,0,0],'s')
            sage: len(t._objects)
            2
        """
        self._objects.append(FractalLandscape(res, scale, center, texture))

    def plot(self,f,(xmin,xmax),(ymin,ymax),texture,grad_f=None,
                  max_bend=.7,max_depth=5,initial_depth=3, num_colors=None):
        r"""
        INPUT:


        -  ``f`` - Function of two variables, which returns a
           float (or coercible to a float) (xmin,xmax)

        -  ``(ymin,ymax)`` - defines the rectangle to plot over
           texture: Name of texture to be used Optional arguments:

        -  ``grad_f`` - gradient function. If specified,
           smooth triangles will be used.

        -  ``max_bend`` - Cosine of the threshold angle
           between triangles used to determine whether or not to recurse after
           the minimum depth

        -  ``max_depth`` - maximum recursion depth. Maximum
           triangles plotted = `2^{2*max_depth}`

        -  ``initial_depth`` - minimum recursion depth. No
           error-tolerance checking is performed below this depth. Minimum
           triangles plotted: `2^{2*min_depth}`

        -  ``num_colors`` - Number of rainbow bands to color
           the plot with. Texture supplied will be cloned (with different
           colors) using the texture_recolor method of the Tachyon object.


        Plots a function by constructing a mesh with nonstandard sampling
        density without gaps. At very high resolutions (depths 10) it
        becomes very slow. Cython may help. Complexity is approx.
        `O(2^{2*maxdepth})`. This algorithm has been optimized for
        speed, not memory - values from f(x,y) are recycled rather than
        calling the function multiple times. At high recursion depth, this
        may cause problems for some machines.

        Flat Triangles::

            sage: t = Tachyon(xres=512,yres=512, camera_center=(4,-4,3),viewdir=(-4,4,-3), raydepth=4)
            sage: t.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return float(sin(x*y))
            sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,  opacity=1.0, color=(1.0,0,0))
            sage: t.plot(f,(-4,4),(-4,4),"t0",max_depth=5,initial_depth=3, num_colors=60)  # increase min_depth for better picture
            sage: t.show()

        Plotting with Smooth Triangles (requires explicit gradient
        function)::

            sage: t = Tachyon(xres=512,yres=512, camera_center=(4,-4,3),viewdir=(-4,4,-3), raydepth=4)
            sage: t.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return float(sin(x*y))
            sage: def g(x,y): return ( float(y*cos(x*y)), float(x*cos(x*y)), 1 )
            sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,  opacity=1.0, color=(1.0,0,0))
            sage: t.plot(f,(-4,4),(-4,4),"t0",max_depth=5,initial_depth=3, grad_f = g)  # increase min_depth for better picture
            sage: t.show()

        Preconditions: f is a scalar function of two variables, grad_f is
        None or a triple-valued function of two variables, min_x !=
        max_x, min_y != max_y

        ::

            sage: f = lambda x,y: x*y
            sage: t = Tachyon()
            sage: t.plot(f,(2.,2.),(-2.,2.),'')
            Traceback (most recent call last):
            ...
            ValueError: Plot rectangle is really a line.  Make sure min_x != max_x and min_y != max_y.
        """
        factory = TachyonTriangleFactory(self,texture)
        plot = TrianglePlot(factory, f, (xmin, xmax), (ymin, ymax), g = grad_f,
                             min_depth=initial_depth, max_depth=max_depth, max_bend=max_bend, num_colors = num_colors)
        self._objects.append(plot)


    def parametric_plot(self, f, t_0, t_f, tex, r=.1, cylinders = True, min_depth=4, max_depth=8, e_rel = .01, e_abs = .01):
        r"""
        Plots a space curve as a series of spheres and finite cylinders.
        Example (twisted cubic) ::

            sage: f = lambda t: (t,t^2,t^3)
            sage: t = Tachyon(camera_center=(5,0,4))
            sage: t.texture('t')
            sage: t.light((-20,-20,40), 0.2, (1,1,1))
            sage: t.parametric_plot(f,-5,5,'t',min_depth=6)
        """

        self._objects.append(ParametricPlot(f, t_0, t_f, tex, r=r, cylinders=cylinders,min_depth=min_depth,max_depth=max_depth,e_rel=.01,e_abs=.01))

#Doesn't seem to be used:
#    def collect(self, objects):
#        """
#        Add a set of objects to the scene from a collection.
#
#        EXAMPLES::
#
#            sage: t = Tachyon()
#            sage: t.texture('s')
#            sage: for i in range(10): t.sphere((0,0,i),i,'s')
#        """
#        self._objects.extend(objects)

class Light:
    r"""
    Represents lighting objects.

    EXAMPLES::

        sage: from sage.plot.plot3d.tachyon import Light
        sage: q = Light((1,1,1),1,(1,1,1))
        sage: q._center
        (1, 1, 1)
    """
    def __init__(self, center, radius, color):
        r"""
        Stores the center, radius and color.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Light
            sage: q = Light((1,1,1),1,(1,1,1))
            sage: q._color
            (1, 1, 1)
        """
        self._center = center
        self._radius = radius
        self._color = color

    def str(self):
        r"""
        Returns the tachyon string defining the light source.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Light
            sage: q = Light((1,1,1),1,(1,1,1))
            sage: q._radius
            1
        """
        return r"""
        light center %s
              rad %s
              color %s
        """%(tostr(self._center), float(self._radius),
             tostr(self._color))

class Texfunc:
    def __init__(self, type=0,center=(0,0,0), rotate=(0,0,0), scale=(1,1,1)):
        r"""
        Creates a texture function.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texfunc
            sage: t = Texfunc()
            sage: t._type
            0
        """
        self._type = type
        self._center = center
        self._rotate = rotate
        self._scale = scale

    def str(self):
        r"""
        Returns the scene string for this texture function.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texfunc
            sage: t = Texfunc()
            sage: t.str()
            '0 center  0.0 0.0 0.0  rotate  0.0 0.0 0.0  scale  1.0 1.0 1.0 '
        """
        if type == 0:
            return "0"
        return r"""%d center %s rotate %s scale %s"""%(self._type,
                                                      tostr(self._center),
                                                      tostr(self._rotate),
                                                      tostr(self._scale))

class Texture:
    def __init__(self, name, ambient=0.2, diffuse=0.8,
                 specular=0.0, opacity=1.0,
                 color=(1.0,0.0, 0.5), texfunc=0, phong=0, phongsize=0, phongtype="PLASTIC"):
        r"""
        Stores texture information.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texture
            sage: t = Texture('w')
            sage: t.str().split()[2:6]
            ['ambient', '0.2', 'diffuse', '0.8']
        """
        self._name = name
        self._ambient = ambient
        self._diffuse = diffuse
        self._specular = specular
        self._opacity = opacity
        self._color = color
        self._texfunc = texfunc
        self._phong = phong
        self._phongsize = phongsize
        self._phongtype = phongtype

    def recolor(self, name, color):
        r"""
        Returns a texture with the new given color.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texture
            sage: t2 = Texture('w')
            sage: t2w = t2.recolor('w2', (.1,.2,.3))
            sage: t2ws = t2w.str()
            sage: color_index = t2ws.find('color')
            sage: t2ws[color_index:color_index+20]
            'color  0.1 0.2 0.3  '
        """
        return Texture(name, self._ambient, self._diffuse, self._specular, self._opacity,
                             color, self._texfunc, self._phong, self._phongsize, self._phongtype)

    def str(self):
        r"""
        Returns the scene string for this texture.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Texture
            sage: t = Texture('w')
            sage: t.str().split()[2:6]
            ['ambient', '0.2', 'diffuse', '0.8']

        """
        return r"""
        texdef %s ambient %s diffuse %s specular %s opacity %s
        phong %s %s phong_size %s
        color %s texfunc %s
        """%(self._name,
             self._ambient,
             self._diffuse,
             self._specular,
             self._opacity,
             self._phongtype,
             self._phong,
             self._phongsize,
             tostr(self._color),
             self._texfunc)

class Sphere:
    r"""
    A class for creating spheres in tachyon.
    """
    def __init__(self, center, radius, texture):
        r"""
        Stores the center, radius, and texture information in a class.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Sphere
            sage: t.texture('r', color=(.8,0,0), ambient = .1)
            sage: s = Sphere((1,1,1),1,'r')
            sage: s._radius
            1
        """
        self._center = center
        self._radius = radius
        self._texture = texture

    def str(self):
        r"""
        Returns the scene string for the sphere.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Sphere
            sage: t.texture('r', color=(.8,0,0), ambient = .1)
            sage: s = Sphere((1,1,1),1,'r')
            sage: s.str()
            '\n        sphere center  1.0 1.0 1.0  rad 1.0 r\n        '
        """
        return r"""
        sphere center %s rad %s %s
        """%(tostr(self._center), float(self._radius), self._texture)

class Ring:
    r"""
    An annulus of zero thickness.
    """
    def __init__(self, center, normal, inner, outer, texture):
        r"""
        Creates a ring with the given center, normal, inner radius,
        outer radius, and texture.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Ring
            sage: r = Ring((1,1,1), (1,1,0), 1.0, 2.0, 's')
            sage: r._center
            (1, 1, 1)
        """
        self._center = center
        self._normal = normal
        self._inner = inner
        self._outer = outer
        self._texture = texture

    def str(self):
        r"""
        Returns the scene string of the ring.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Ring
            sage: r = Ring((0,0,0), (1,1,0), 1.0, 2.0, 's')
            sage: r.str()
            '\n        ring center  0.0 0.0 0.0  normal  1.0 1.0 0.0  inner 1.0 outer 2.0 s\n        '
        """
        return r"""
        ring center %s normal %s inner %s outer %s %s
        """%(tostr(self._center), tostr(self._normal), float(self._inner), float(self._outer), self._texture)

class FractalLandscape:
    r"""
    Axis-aligned fractal landscape.
    Does not seem very useful at the moment, but perhaps will be improved in the future.
    """
    def __init__(self, res, scale, center, texture):
        r"""
        Creates a fractal landscape in tachyon.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FractalLandscape
            sage: fl = FractalLandscape([20,20],[30,30],[1,2,3],'s')
            sage: fl._center
            [1, 2, 3]
        """
        self._res = res
        self._scale = scale
        self._center = center
        self._texture = texture

    def str(self):
        r"""
        Returns the scene string of the fractal landscape.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FractalLandscape
            sage: fl = FractalLandscape([20,20],[30,30],[1,2,3],'s')
            sage: fl.str()
            '\n        scape res  20 20  scale  30 30  center  1.0 2.0 3.0  s\n        '
        """
        return r"""
        scape res %s scale %s center %s %s
        """%(tostr(self._res, 2, int), tostr(self._scale, 2, int), tostr(self._center), self._texture)

class Cylinder:
    r"""
    An infinite cylinder.
    """
    def __init__(self, center, axis, radius, texture):
        r"""
        Creates a cylinder with the given parameters.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Cylinder
            sage: c = Cylinder((0,0,0),(1,1,1),.1,'s')
            sage: c.str()
            '\n        cylinder center  0.0 0.0 0.0  axis  1.0 1.0 1.0  rad 0.1 s\n        '
        """
        self._center = center
        self._axis = axis
        self._radius = radius
        self._texture = texture

    def str(self):
        r"""
        Returns the scene string of the cylinder.

        EXAMPLES::

            sage: t = Tachyon()
            sage: from sage.plot.plot3d.tachyon import Cylinder
            sage: c = Cylinder((0,0,0),(1,1,1),.1,'s')
            sage: c.str()
            '\n        cylinder center  0.0 0.0 0.0  axis  1.0 1.0 1.0  rad 0.1 s\n        '
            """
        return r"""
        cylinder center %s axis %s rad %s %s
        """%(tostr(self._center), tostr(self._axis), float(self._radius), self._texture)

class Plane:
    r"""
    An infinite plane.
    """
    def __init__(self, center, normal, texture):
        r"""
        Creates the plane object.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Plane
            sage: p = Plane((1,2,3),(1,2,4),'s')
            sage: p.str()
            '\n        plane center  1.0 2.0 3.0  normal  1.0 2.0 4.0  s\n        '
        """
        self._center = center
        self._normal = normal
        self._texture = texture

    def str(self):
        r"""
        Returns the scene string of the plane.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Plane
            sage: p = Plane((1,2,3),(1,2,4),'s')
            sage: p.str()
            '\n        plane center  1.0 2.0 3.0  normal  1.0 2.0 4.0  s\n        '
        """
        return r"""
        plane center %s normal %s %s
        """%(tostr(self._center), tostr(self._normal), self._texture)

class FCylinder:
    r"""
    A finite cylinder.
    """
    def __init__(self, base, apex, radius, texture):
        r"""
        Creates a finite cylinder object.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FCylinder
            sage: fc = FCylinder((0,0,0),(1,1,1),.1,'s')
            sage: fc.str()
            '\n        fcylinder base  0.0 0.0 0.0  apex  1.0 1.0 1.0  rad 0.1 s\n        '
        """
        self._center = base
        self._axis = apex
        self._radius = radius
        self._texture = texture

    def str(self):
        r"""
        Returns the scene string of the finite cylinder.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import FCylinder
            sage: fc = FCylinder((0,0,0),(1,1,1),.1,'s')
            sage: fc.str()
            '\n        fcylinder base  0.0 0.0 0.0  apex  1.0 1.0 1.0  rad 0.1 s\n        '
        """
        return r"""
        fcylinder base %s apex %s rad %s %s
        """%(tostr(self._center), tostr(self._axis), float(self._radius), self._texture)

class Axis_aligned_box():
    r"""
    Box with axis-aligned edges with the given min and max coordinates.
    """
    def __init__(self, min_p, max_p, texture):
        r"""
        Creates the axis-aligned box object.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Axis_aligned_box
            sage: aab = Axis_aligned_box((0,0,0),(1,1,1),'s')
            sage: aab.str()
            '\n        box min  0.0 0.0 0.0  max  1.0 1.0 1.0  s\n        '
        """
        self._min_p = min_p
        self._max_p = max_p
        self._texture = texture

    def str(self):
        r"""
        Returns the scene string of the axis-aligned box.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import Axis_aligned_box
            sage: aab = Axis_aligned_box((0,0,0),(1,1,1),'s')
            sage: aab.str()
            '\n        box min  0.0 0.0 0.0  max  1.0 1.0 1.0  s\n        '
        """
        return r"""
        box min %s max %s %s
        """%(tostr(self._min_p), tostr(self._max_p), self._texture)

class TachyonTriangle(Triangle):
    r"""
    Basic triangle class.
    """
    def str(self):
        r"""
        Returns the scene string for a triangle.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangle
            sage: t = TachyonTriangle([-1,-1,-1],[0,0,0],[1,2,3])
            sage: t.str()
            '\n        TRI V0  -1.0 -1.0 -1.0   V1  0.0 0.0 0.0    V2  1.0 2.0 3.0 \n            0\n        '
        """
        return r"""
        TRI V0 %s  V1 %s   V2 %s
            %s
        """%(tostr(self._a), tostr(self._b),tostr(self._c), self._color)

class TachyonSmoothTriangle(SmoothTriangle):
    r"""
    A triangle along with a normal vector, which is used for smoothing.
    """
    def str(self):
        r"""
        Returns the scene string for a smoothed triangle.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonSmoothTriangle
            sage: t = TachyonSmoothTriangle([-1,-1,-1],[0,0,0],[1,2,3],[1,0,0],[0,1,0],[0,0,1])
            sage: t.str()
            '\n        STRI V0  ...  1.0 0.0 0.0  N1  0.0 1.0 0.0   N2  0.0 0.0 1.0 \n             0\n        '
        """
        return r"""
        STRI V0 %s V1 %s  V2 %s
             N0 %s N1 %s  N2 %s
             %s
        """%(tostr(self._a),  tostr(self._b),  tostr(self._c),
             tostr(self._da), tostr(self._db), tostr(self._dc), self._color)



class TachyonTriangleFactory(TriangleFactory):
    r"""
    A class to produce triangles of various rendering types.
    """
    def __init__(self, tach, tex):
        r"""
        Initializes with tachyon instance and texture.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttf._texture
            's'
        """
        self._tachyon = tach
        self._texture = tex

    def triangle(self,a,b,c,color=None):
        r"""
        Creates a TachyonTriangle with vertices a, b, and c.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttft = ttf.triangle([1,2,3],[3,2,1],[0,2,1])
            sage: ttft.str()
            '\n        TRI V0  1.0 2.0 3.0   V1  3.0 2.0 1.0    V2  0.0 2.0 1.0 \n            s\n        '
        """
        if color is None:
            return TachyonTriangle(a,b,c,self._texture)
        else:
            return TachyonTriangle(a,b,c,color)

    def smooth_triangle(self,a,b,c,da,db,dc,color=None):
        r"""
        Creates a TachyonSmoothTriangle.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttfst = ttf.smooth_triangle([0,0,0],[1,0,0],[0,0,1],[1,1,1],[1,2,3],[-1,-1,2])
            sage: ttfst.str()
            '\n        STRI V0  0.0 0.0 0.0  ...'
        """
        if color is None:
            return TachyonSmoothTriangle(a,b,c,da,db,dc,self._texture)
        else:
            return TachyonSmoothTriangle(a,b,c,da,db,dc,color)

    def get_colors(self, list):
        r"""
        Returns a list of color labels.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import TachyonTriangleFactory
            sage: t = Tachyon()
            sage: t.texture('s')
            sage: ttf = TachyonTriangleFactory(t, 's')
            sage: ttf.get_colors([1])
            ['SAGETEX1_0']
        """
        return self._tachyon.texture_recolor(self._texture, list)

# following classes TachyonPlot and PlotBlock seems broken and not used anywhere, so commented out.  Please write to the sage-devel google-group if you are the author of these classes to comment.

#class TachyonPlot:
    #Recursively plots a function of two variables by building squares of 4 triangles, checking at
    # every stage whether or not each square should be split into four more squares.  This way,
    # more planar areas get fewer triangles, and areas with higher curvature get more triangles

#    def str(self):
#        return "".join([o.str() for o in self._objects])

#    def __init__(self, tachyon, f, (min_x, max_x), (min_y, max_y), tex, g = None,
#                              min_depth=4, max_depth=8, e_rel = .01, e_abs = .01, num_colors = None):
#        self._tachyon = tachyon
#        self._f = f
#        self._g = g
#        self._tex = tex
#        self._min_depth = min_depth
#        self._max_depth = max_depth
#        self._e_rel = e_rel
#        self._e_abs = e_abs
#        self._objects = []
#        self._eps = min(max_x - min_x, max_y - min_y)/(2**max_depth)
#        if self._eps == 0:
#            raise ValueError, 'Plot rectangle is really a line.  Make sure min_x != #max_x and min_y != max_y.'
#        self._num_colors = num_colors
#        if g is None:
#            def fcn(x,y):
#                return [self._f(x,y)]
#        else:
#            def fcn(x,y):
#                return [self._f(x,y), self._g(x,y)]

#        self._fcn = fcn


#        # generate the necessary data to kick-start the recursion
#        mid_x = (min_x + max_x)/2
#        mid_y = (min_y + max_y)/2
#        sw_z = fcn(min_x,min_y)
#        nw_z = fcn(min_x,max_y)
#        se_z = fcn(max_x,min_y)
#        ne_z = fcn(max_x,max_y)
#        mid_z = fcn(mid_x,mid_y)

#        self._min = min(sw_z[0], nw_z[0], se_z[0], ne_z[0], mid_z[0])
#        self._max = max(sw_z[0], nw_z[0], se_z[0], ne_z[0], mid_z[0])

#        # jump in and start building blocks
#        outer = self.plot_block(min_x, mid_x, max_x, min_y, mid_y, max_y, sw_z, nw_z, se_z, ne_z, mid_z, 0)
#
#        # build the boundary triangles
#        self.triangulate(outer.left, outer.left_c)
#        self.triangulate(outer.top, outer.top_c)
#        self.triangulate(outer.right, outer.right_c)
#        self.triangulate(outer.bottom, outer.bottom_c)

#        zrange = self._max - self._min
#        if num_colors is not None and zrange != 0:
#            colors = tachyon.texture_recolor(tex, [hue(float(i/num_colors)) for i in range(num_colors)])

#            for o in self._objects:
#                avg_z = (o._vertex_1[2] + o._vertex_2[2] + o._vertex_3[2])/3
#                o._texture = colors[int(num_colors * (avg_z - self._min) / zrange)]

#    def plot_block(self, min_x, mid_x, max_x, min_y, mid_y, max_y, sw_z, nw_z, se_z, ne_z, mid_z, depth):

#        if depth < self._max_depth:
#            # recursion is still an option -- step in one last level if we're within tolerance
#            # and just keep going if we're not.
#            # assumption: it's cheap to build triangles, so we might as well use all the data
#            # we calculate

#            # big square boundary midpoints
#            mid_w_z = self._fcn(min_x, mid_y)
#            mid_n_z = self._fcn(mid_x, max_y)
#            mid_e_z = self._fcn(max_x, mid_y)
#            mid_s_z = self._fcn(mid_x, min_y)

#            # midpoints locations of sub_squares
#            qtr1_x = (min_x + mid_x)/2
#            qtr1_y = (min_y + mid_y)/2
#            qtr3_x = (mid_x + max_x)/2
#           qtr3_y = (mid_y + max_y)/2

#            # function evaluated at these midpoints
#            mid_sw_z = self._fcn(qtr1_x,qtr1_y)
#            mid_nw_z = self._fcn(qtr1_x,qtr3_y)
#            mid_se_z = self._fcn(qtr3_x,qtr1_y)
#            mid_ne_z = self._fcn(qtr3_x,qtr3_y)

#            # linearization estimates of midpoints
#            est_sw_z = (mid_z[0] + sw_z[0])/2
#            est_nw_z = (mid_z[0] + nw_z[0])/2
#            est_se_z = (mid_z[0] + se_z[0])/2
#           est_ne_z = (mid_z[0] + ne_z[0])/2

#            self.extrema([mid_w_z[0], mid_n_z[0], mid_e_z[0], mid_s_z[0], mid_sw_z[0], mid_se_z[0], mid_nw_z[0], mid_sw_z[0]])

#            tol_check = [(est_sw_z, mid_sw_z[0]), (est_nw_z, mid_nw_z[0]), (est_se_z, mid_se_z[0]), (est_ne_z, mid_ne_z[0])]

#            if depth < self._min_depth or not self.tol_list(tol_check):
#                next_depth = depth + 1
#            else:
#                #lie about the depth to halt recursion
#                next_depth = self._max_depth

#            # recurse into the sub-squares
#            sw = self.plot_block(min_x, qtr1_x, mid_x, min_y, qtr1_y, mid_y, sw_z, mid_w_z, mid_s_z, mid_z, mid_sw_z, next_depth)
#            nw = self.plot_block(min_x, qtr1_x, mid_x, mid_y, qtr3_y, max_y, mid_w_z, nw_z, mid_z, mid_n_z, mid_nw_z, next_depth)
#            se = self.plot_block(mid_x, qtr3_x, max_x, min_y, qtr1_y, mid_y, mid_s_z, mid_z, se_z, mid_e_z, mid_se_z, next_depth)
#            ne = self.plot_block(mid_x, qtr3_x, max_x, mid_y, qtr3_y, max_y, mid_z, mid_n_z, mid_e_z, ne_z, mid_ne_z, next_depth)

#            # join the sub-squares
#            self.interface(1, sw.right, sw.right_c, se.left, se.left_c)
#            self.interface(1, nw.right, nw.right_c, ne.left, ne.left_c)
#            self.interface(0, sw.top, sw.top_c, nw.bottom, nw.bottom_c)
#            self.interface(0, se.top, se.top_c, ne.bottom, ne.bottom_c)

#            #get the boundary information about the subsquares
#            left     = sw.left     + nw.left[1:]
#            left_c   = sw.left_c   + nw.left_c
#            right    = se.right    + ne.right[1:]
#            right_c  = se.right_c  + ne.right_c
#            top      = nw.top      + ne.top[1:]
#            top_c    = nw.top_c    + ne.top_c
#            bottom   = sw.bottom   + se.bottom[1:]
#            bottom_c = sw.bottom_c + se.bottom_c

#        else:
#            # just build the square we're in
#            if self._g is None:
#                sw = [(min_x,min_y,sw_z[0])]
#                nw = [(min_x,max_y,nw_z[0])]
#                se = [(max_x,min_y,se_z[0])]
#                ne = [(max_x,max_y,ne_z[0])]
#                c  = [[(mid_x,mid_y,mid_z[0])]]
#            else:
#                sw = [(min_x,min_y,sw_z[0]),sw_z[1]]
#                nw = [(min_x,max_y,nw_z[0]),nw_z[1]]
#               se = [(max_x,min_y,se_z[0]),se_z[1]]
#               ne = [(max_x,max_y,ne_z[0]),ne_z[1]]
#               c  = [[(mid_x,mid_y,mid_z[0]),mid_z[1]]]


#            left = [sw,nw]
#            left_c = c
#            top = [nw,ne]
#            top_c = c
#            right = [se,ne]
#            right_c = c
#           bottom = [sw,se]
#            bottom_c = c

#        return PlotBlock(left, left_c, top, top_c, right, right_c, bottom, bottom_c)

#    def tol(self, (est, val)):
#        # Check relative, then absolute tolerance.  If both fail, return False
#        # This is a zero-safe error checker

#        if abs(est - val) < self._e_rel*abs(val):
#            return True
#        if abs(est - val) < self._e_abs:
#            return True
#        return False

#    def tol_list(self, l):
#        # Pass in a list of pairs of numbers, (est, val) to be passed to self.tol
#        # returns False if any pair does not fall within tolerance level

#        for p in l:
#            if not self.tol(p):
#                return False
#        return True

#    def interface(self, n, p, p_c, q, q_c):
#        # Takes a pair of lists of points, and compares the (n)th coordinate, and
#        # "zips" the lists together into one.  The "centers", supplied in p_c and
#        # q_c are matched up such that the lists describe triangles whose sides
#       # are "perfectly" aligned.  This algorithm assumes that p and q start and
#        # end at the same point, and are sorted smallest to largest.

#        m   = [p[0]] # a sorted union of p and q
#        mpc = [p_c[0]] # centers from p_c corresponding to m
#        mqc = [q_c[0]] # centers from q_c corresponding to m

#        i = 1
#        j = 1

#        while i < len(p_c) or j < len(q_c):
#            if abs(p[i][0][n] - q[j][0][n]) < self._eps:
#                m.append(p[i])
#                mpc.append(p_c[i])
#                mqc.append(q_c[j])
#               i += 1
#               j += 1
#            elif p[i][0][n] < q[j][0][n]:
#                m.append(p[i])
#                mpc.append(p_c[i])
#                mqc.append(mqc[-1])
#                i += 1
#            else:
#                m.append(q[j])
#               mpc.append(mpc[-1])
#                mqc.append(q_c[j])
#                j += 1

#        m.append(p[-1])

#        self.triangulate(m, mpc)
#        self.triangulate(m, mqc)


#    def triangulate(self, p, c):
#        # Pass in a list of edge points (p) and center points (c).
#        # Triangles will be rendered between consecutive edge points and the
#        # center point with the same index number as the earlier edge point.

#        if self._g is None:
#            for i in range(0,len(p)-1):
#                self._objects.append(Triangle(p[i][0], p[i+1][0], c[i][0], self._tex))
#        else:
#            for i in range(0,len(p)-1):
#                self._objects.append(SmoothTriangle(p[i][0], p[i+1][0], c[i][0],p[i][1], p[i+1][1], c[i][1], self._tex))


#    def extrema(self, list):
#        if self._num_colors is not None:
#            self._min = min(list+[self._min])
#            self._max = max(list+[self._max])


#class PlotBlock:
#   def __init__(self, left, left_c, top, top_c, right, right_c, bottom, bottom_c):
#       self.left = left
#       self.left_c = left_c
#       self.top = top
#      self.top_c = top_c
#       self.right = right
#       self.right_c = right_c
#       self.bottom = bottom
#       self.bottom_c = bottom_c

class ParametricPlot:
    r"""
    Parametric plotting routines.
    """
    def str(self):
        r"""
        Returns the tachyon string representation of the parameterized curve.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: t = var('t')
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q.str()[9:69]
            'sphere center  0.0 0.0 0.0  rad 0.1 s\n        \n        fcyli'
        """
        return "".join([o.str() for o in self._objects])

    def __init__(self, f, t_0, t_f, tex, r=.1, cylinders = True, min_depth=4, max_depth=8, e_rel = .01, e_abs = .01):
        r"""
        Creates the parametric plotting class.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: t = var('t')
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q._e_rel
            0.01
        """
        self._e_rel = e_rel
        self._e_abs = e_abs
        self._r = r
        self._f = f
        self._tex = tex
        self._cylinders = cylinders
        self._max_depth = max_depth
        self._min_depth = min_depth

        f_0 = f(t_0)
        f_f = f(t_f)
        self._objects = [Sphere(f_0, r, texture=tex) ]

        self._plot_step(0, t_0, t_f, f_0, f_f)

    def _plot_step(self, depth, t_0,t_f,f_0,f_f):
        r"""
        Recursively subdivides interval, eventually plotting with cylinders and spheres.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: t = var('t')
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q._plot_step(8,0,1,[0,0,0],[1,1,1])
            sage: len(q._objects)
            515
        """
        if depth < self._max_depth:
            t_mid = (t_f + t_0)/2
            f_mid = ((f_f[0] + f_0[0])/2, (f_f[1] + f_0[1])/2, (f_f[2] + f_0[2])/2)
            f_val = self._f(t_mid)
            if depth < self._min_depth or self.tol(f_mid, f_val):
                new_depth = depth + 1
            else:
                new_depth = self._max_depth

            self._plot_step(new_depth, t_0,t_mid, f_0, f_val)
            self._plot_step(new_depth, t_mid,t_f, f_val, f_f)
        else:
            if self._cylinders:
                self._objects.append(FCylinder(f_0,f_f,self._r,self._tex))
            self._objects.append(Sphere(f_f,self._r,self._tex))


    def tol(self, est, val):
        r"""
        Check relative, then absolute tolerance.  If both fail, return False.
        This is a zero-safe error checker.

        EXAMPLES::

            sage: from sage.plot.plot3d.tachyon import ParametricPlot
            sage: t = var('t')
            sage: f = lambda t: (t,t^2,t^3)
            sage: q = ParametricPlot(f,0,1,'s')
            sage: q.tol([0,0,0],[1,0,0])
            False
            sage: q.tol([0,0,0],[.0001,0,0])
            True
        """
        delta = sqrt((val[0]-est[0])**2 + (val[1]-est[1])**2 + (val[2]-val[2])**2)
        if delta < self._e_abs:
            return True

        r = sqrt(val[0]**2+val[1]**2+val[2]**2)
        if delta < self._e_rel*r:
            return True

        return False

def tostr(s, length = 3, out_type = float):
    r"""
    Converts vector information to a space-separated string.

    EXAMPLES::

        sage: from sage.plot.plot3d.tachyon import tostr
        sage: tostr((1,1,1))
        ' 1.0 1.0 1.0 '
        sage: tostr('2 3 2')
        '2 3 2'
    """
    if isinstance(s, str):
        return s
    output = ' '
    for an_item in s:
        output = output + str(out_type(an_item)) + ' '
    return output

