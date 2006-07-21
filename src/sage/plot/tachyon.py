r"""
Interface to the Tachyon Ray Tracer

AUTHOR:
    -- John E. Stone (johns@megapixel.com) -- wrote tachyon ray tracer
    -- William Stein -- write interface

TODO:
   -- currently only spheres, lights, and textures are wrapped.  need to add triangles, etc.
"""

from sage.interfaces.tachyon import tachyon_rt

from sage.ext.sage_object import SageObject

import os

class Tachyon(SageObject):
    """
    A scene the can be rendered using the Tachyon ray tracer.

    EXAMPLES:
    Two spheres.
        sage: t = Tachyon(xres=350,yres=350)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0.8,0.7,0.9))
        sage: t.sphere((0,-0.5,0), 0.2, 't1')
        sage: t.sphere((0,0.5,0), 0.4, 't1')
        sage: t.save()

    Sphere's along the twisted cubic.
        sage: t = Tachyon(xres=512,yres=512, center=(0,0,5))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2, diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in srange(-5,1.5,0.1):
        ...    k += 1
        ...    t.sphere((i,i^2-0.5,i^3), 0.1, 't%s'%(k%3))
        ...
        sage: t.save()

    Many random spheres:
        sage: t = Tachyon(xres=512,yres=512, center=(0.5,0.5,2), raydepth=4)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2, diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in range(100):
        ...    k += 1
        ...    t.sphere((random(),random(), random()), random()/10, 't%s'%(k%3))
        ...
        sage: t.save()         # long (several seconds)

    Points on an elliptic curve:
        sage: t = Tachyon(xres=512,yres=512, center=(0.5,1,10))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: E = EllipticCurve('37a')
        sage: P = E([0,0])
        sage: Q = P
        sage: def f(r):
        ...    h = log(log(r.height() + 1))
        ...    return h
        ...
        sage: for i in range(20):   # increase 20 for a better plot
        ...    Q = Q + P
        ...    t.sphere((Q[0], -f(Q), Q[1]), 0.1, 't0')
        ...
        sage: t.save()
    """
    def __init__(self, xres=350, yres=350,
                 zoom = 1.0,
                 antialiasing = False,
                 aspectratio = 1.0,
                 raydepth = 12,
                 center = (0, 0, 2),
                 viewdir = (0, 0, -1),
                 updir = (0, 1, 0),
                 projection = 'PERSPECTIVE'):
        self._xres = xres
        self._yres = yres
        self._zoom = zoom
        self._aspectratio = aspectratio
        self._antialiasing = antialiasing
        self._raydepth = raydepth
        self._center = center
        self._viewdir = viewdir
        self._updir = updir
        self._projection = projection
        self._objects = []

    def __repr__(self):
        return self.str()

    def save(self, filename='sage.png', verbose=0, block=True, extra_opts=''):
        """
            filename -- (default: 'sage.png')
                       output filename; the extension of
                       the filename determines the type.
                       Supported types include:
                         tga -- 24-bit (uncompressed)
                         bmp -- 24-bit Windows BMP (uncompressed)
                         ppm -- 24-bit PPM (uncompressed)
                         rgb -- 24-bit SGI RGB (uncompressed)
                         png -- 24-bit PNG (compressed, lossless)
            verbose -- integer; (default: 0)
                       0 -- silent
                       1 -- some output
                       2 -- very verbose output

            block -- bool (default: True); if False, run the rendering
                     command in the background.

            extra_opts -- passed directly to tachyon command line.
                     Use tachyon_rt.usage() to see some of the possibilities.
        """
        tachyon_rt(self.str(), filename, verbose, block, extra_opts)

    def show(self, verbose=0, extra_opts=''):
        import sage.server.support
        if sage.server.support.EMBEDDED_MODE:
            i = 0
            while os.path.exists('sage%s.png'%i):
                i += 1
            filename = 'sage%s.png'%i
            self.save(filename, verbose=verbose, extra_opts=extra_opts)
        else:
            raise NotImplementedError


    def _res(self):
        return '\nresolution %s %s\n'%(self._xres, self._yres)

    def _camera(self):
        return """
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
             tostr(self._center),
             tostr(self._viewdir),
             tostr(self._updir))

    def str(self):
        return """
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
        self._objects.append(Light(center, radius, color))

    def texture(self, name, ambient=0.2, diffuse=0.8,
                specular=0.0, opacity=1.0,
                color=(1.0,0.0, 0.5), texfunc=0):
        self._objects.append(Texture(name, ambient, diffuse,
                                     specular, opacity, color, texfunc))

    def sphere(self, center, radius, texture):
        self._objects.append(Sphere(center, radius, texture))


class Light:
    def __init__(self, center, radius, color):
        self._center = center
        self._radius = radius
        self._color = color

    def str(self):
        return """
        light center %s
              rad %s
              color %s
        """%(tostr(self._center), float(self._radius),
             tostr(self._color))

class Texture:
    def __init__(self, name, ambient=0.2, diffuse=0.8,
                 specular=0.0, opacity=1.0,
                 color=(1.0,0.0, 0.5), texfunc=0):
        self._name = name
        self._ambient = ambient
        self._diffuse = diffuse
        self._specular = specular
        self._opacity = opacity
        self._color = color
        self._texfunc = texfunc

    def str(self):
        return """
        texdef %s ambient %s diffuse %s specular %s opacity %s
        color %s texfunc %s
        """%(self._name,
             self._ambient,
             self._diffuse,
             self._specular,
             self._opacity,
             tostr(self._color),
             self._texfunc)

class Sphere:
    def __init__(self, center, radius, texture):
        self._center = center
        self._radius = radius
        self._texture = texture

    def str(self):
        return """
        sphere center %s rad %s %s
        """%(tostr(self._center), float(self._radius), self._texture)


def tostr(s):
    if isinstance(s, str):
        return s
    return ' %s %s %s '%(float(s[0]), float(s[1]), float(s[2]))





