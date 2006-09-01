r"""
Interface to the Tachyon Ray Tracer

AUTHOR:
    -- John E. Stone (johns@megapixel.com) -- wrote tachyon ray tracer
    -- William Stein -- sage-tachyon interface
    -- Joshua Kantor -- 3d function plotting
    -- Tom Boothby -- 3d function plotting n'stuff

TODO:
   -- currently only spheres, lights, and textures are wrapped.  need to add triangles, etc.
"""

from colorsys import hsv_to_rgb

from sage.interfaces.tachyon import tachyon_rt

from sage.ext.sage_object import SageObject

#from sage.ext import fast_tachyon_routines

import os

from math import modf,fabs

class Tachyon(SageObject):
    """
    Create a scene the can be rendered using the Tachyon ray tracer.

    INPUT:
                 xres=350, yres=350,
                 zoom = 1.0,
                 antialiasing = False,
                 aspectratio = 1.0,
                 raydepth = 12,
                 camera_center = (-3, 0, 0),
                 updir = (0, 0, 1),
                 look_at = (0,0,0),
                 viewdir = None,
                 projection = 'PERSPECTIVE'

    OUTPUT:
        A Tachyon 3d scene.

    Note that the coordinates are by default such that z is up,
    positive y is to the *left* and x is toward you.  This is
    not oriented according to the right hand rule.

    EXAMPLES:
    Three spheres on the coordinate axes:

        sage: t = Tachyon(xres=500,yres=500, camera_center=(2,0,0))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
        sage: t.texture('t3', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
        sage: t.texture('t4', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
        sage: t.sphere((0,0.5,0), 0.2, 't2')
        sage: t.sphere((0.5,0,0), 0.2, 't3')
        sage: t.sphere((0,0,0.5), 0.2, 't4')
        sage: t.save()

    Sphere's along the twisted cubic.
        sage: t = Tachyon(xres=512,yres=512, camera_center=(3,0.3,0))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2,diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in srange(-1,1,0.05):
        ...    k += 1
        ...    t.sphere((i,i^2-0.5,i^3), 0.1, 't%s'%(k%3))
        ...
        sage.: t.save()

    Another twisted cubic, but with a white background, got by putting
    infinite planes around the scene.
        sage: t = Tachyon(xres=512,yres=512, camera_center=(3,0.3,0), raydepth=8)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2,diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: t.texture('white', color=(1,1,1))
        sage: t.plane((0,0,-1), (0,0,1), 'white')
        sage: t.plane((0,-20,0), (0,1,0), 'white')
        sage: t.plane((-20,0,0), (1,0,0), 'white')
        sage:
        sage: k=0
        sage: for i in srange(-1,1,0.05):
        ...    k += 1
        ...    t.sphere((i,i^2 - 0.5,i^3), 0.1, 't%s'%(k%3))
        ...    t.cylinder((0,0,0), (0,0,1), 0.05,'t1')
        ...
        sage.: t.save()

    Many random spheres:
        sage: t = Tachyon(xres=512,yres=512, camera_center=(2,0.5,0.5), look_at=(0.5,0.5,0.5), raydepth=4)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2, diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in range(100):
        ...    k += 1
        ...    t.sphere((random(),random(), random()), random()/10, 't%s'%(k%3))
        ...
        sage.: t.save()


    Points on an elliptic curve, their height indicated by their height above the axis:
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
        ...    Q = Q + P
        ...    t.sphere((Q[1], Q[0], ZZ(i)/n), 0.1, 't%s'%(i%3))
        ...
        sage.: t.save()

    A beautiful picture of rational points on a rank 1 elliptic curve.
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
        ...    Q = Q + P
        ...    c = i/n + .1
        ...    t.texture('r%s'%i,color=(float(i/n),0,0))
        ...    t.sphere((Q[0], -Q[1], .01), .04, 'r%s'%i)
        ...
        ...
        sage.: t.save()    # 10-20 seconds
    """
    def __init__(self,
                 xres=350, yres=350,
                 zoom = 1.0,
                 antialiasing = False,
                 aspectratio = 1.0,
                 raydepth = 12,
                 camera_center = (-3, 0, 0),
                 updir = (0, 0, 1),
                 look_at = (0,0,0),
                 viewdir = None,
                 projection = 'PERSPECTIVE'):
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
             tostr(self._camera_center),
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
                color=(1.0,0.0, 0.5), texfunc=0, phong=0, phongsize=.5, phongtype="PLASTIC"):
        self._objects.append(Texture(name, ambient, diffuse,
                                     specular, opacity, color, texfunc,
                                     phong,phongsize,phongtype))

    def sphere(self, center, radius, texture):
        self._objects.append(Sphere(center, radius, texture))

    def cylinder(self, center, axis, radius, texture):
        self._objects.append(Cylinder(center, axis, radius, texture))

    def plane(self, center, normal, texture):
        self._objects.append(Plane(center, normal, texture))

    def fcylinder(self, base, apex, radius, texture):
        self._objects.append(FCylinder(base, apex, radius, texture))

    def triangle(self, vertex_1, vertex_2, vertex_3, texture):
	 self._objects.append(Triangle(vertex_1,vertex_2,vertex_3,texture))

    def plot(self,f,(xmin,xmax),(ymin,ymax),texture,max_var=.1,max_depth=5,initial_depth=3):
        r"""
        Plots a function by constructing a mesh with nonstandard sampling density
        without gaps. At very high resolutions (depths > 10) it becomes very
        slow.  Pyrex may help.  Complexity is approx.
        $O(2^{2*maxdepth})$.  This
        algorithm has been optimized for speed, not memory -- values from f(x,y) are
        recycled rather than calling the function multiple times.  At high recursion
        depth, this may cause problems for some machines.

        EXAMPLE:
            sage: t = Tachyon(xres=512,yres=512, camera_center=(4,-4,3),viewdir=(-4,4,-3), raydepth=4)
            sage: t.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return(float(math.sin(x*y)))
            sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,  opacity=1.0, color=(1.0,0,0))
            sage: t.plot(f,(-4,4),(-4,4),"t0",max_depth=8,initial_depth=6)  # increase max_depth for better picture
            sage.: t.show()
        """
        TachyonPlot(self, f, (xmin, xmax), (ymin, ymax), texture, min_depth=initial_depth, max_depth=max_depth, e_rel = .01, e_abs = max_var)

    def collect(self, objects):
        """Add a set of objects to the scene from a collection"""
        self._objects.extend(objects)





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
                 color=(1.0,0.0, 0.5), texfunc=0, phong=0, phongsize=0, phongtype="PLASTIC"):
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

    def str(self):
        return """
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
    def __init__(self, center, radius, texture):
        self._center = center
        self._radius = radius
        self._texture = texture

    def str(self):
        return """
        sphere center %s rad %s %s
        """%(tostr(self._center), float(self._radius), self._texture)


class Cylinder:
    def __init__(self, center, axis, radius, texture):
        self._center = center
        self._axis = axis
        self._radius = radius
        self._texture = texture

    def str(self):
        return """
        cylinder center %s axis %s rad %s %s
        """%(tostr(self._center), tostr(self._axis), float(self._radius), self._texture)

class Plane:
    def __init__(self, center, normal, texture):
        self._center = center
        self._normal = normal
        self._texture = texture

    def str(self):
        return """
        plane center %s normal %s %s
        """%(tostr(self._center), tostr(self._normal), self._texture)

class FCylinder:
    def __init__(self, base, apex, radius, texture):
        self._center = base
        self._axis = apex
        self._radius = radius
        self._texture = texture

    def str(self):
        return """
        fcylinder base %s apex %s rad %s %s
        """%(tostr(self._center), tostr(self._axis), float(self._radius), self._texture)

class Triangle:
      def __init__(self,vertex_1,vertex_2,vertex_3,texture):
	  self._vertex_1 = vertex_1
	  self._vertex_2 = vertex_2
	  self._vertex_3 = vertex_3
	  self._texture = texture


      def str(self):
	  return """
	  TRI
	  V0 %s
	  V1 %s
	  V2 %s
	  %s
	  """%(tostr(self._vertex_1), tostr(self._vertex_2),tostr(self._vertex_3), self._texture)




class TachyonPlot:
    #Recursively plots a function of two variables by building squares of 4 triangles, checking at
    # every stage whether or not each square should be split into four more squares.  This way,
    # more planar areas get fewer triangles, and areas with higher curvature get more trianges

    def __init__(self, tachyon, f, (min_x, max_x), (min_y, max_y), tex, min_depth=4, max_depth=8, e_rel = .01, e_abs = .01):
        self._tachyon = tachyon
        self._f = f
        self._tex = tex
        self._min_depth = min_depth
        self._max_depth = max_depth
        self._e_rel = e_rel
        self._e_abs = e_abs
        self._trianglist = []
        self._eps = min(max_x - min_x, max_y - min_y)/(2**max_depth)

        # generate the necessary data to kick-start the recursion
        mid_x = (min_x + max_x)/2
        mid_y = (min_y + max_y)/2
        sw_z = f(min_x,min_y)
        nw_z = f(min_x,max_y)
        se_z = f(max_x,min_y)
        ne_z = f(max_x,max_y)
        mid_z = f(mid_x,mid_y)

        # jump in and start building blocks
        outer = self.plot_block(min_x, mid_x, max_x, min_y, mid_y, max_y, sw_z, nw_z, se_z, ne_z, mid_z, 0)

        # build the boundary triangles
        self.triangulate(outer.left, outer.left_c)
        self.triangulate(outer.top, outer.top_c)
        self.triangulate(outer.right, outer.right_c)
        self.triangulate(outer.bottom, outer.bottom_c)

    def plot_block(self, min_x, mid_x, max_x, min_y, mid_y, max_y, sw_z, nw_z, se_z, ne_z, mid_z, depth):

        if depth < self._max_depth:
            # recursion is still an option -- step in one last level if we're within tolerance
            # and just keep going if we're not.
            # assumption: it's cheap to build triangles, so we might as well use all the data
            # we calculate

            # big square boundary midpoints
            mid_w_z = self._f(min_x, mid_y)
            mid_n_z = self._f(mid_x, max_y)
            mid_e_z = self._f(max_x, mid_y)
            mid_s_z = self._f(mid_x, min_y)

            # midpoints locations of sub_squares
            qtr1_x = (min_x + mid_x)/2
            qtr1_y = (min_y + mid_y)/2
            qtr3_x = (mid_x + max_x)/2
            qtr3_y = (mid_y + max_y)/2

            # function evaluated at these midpoints
            mid_sw_z = self._f(qtr1_x,qtr1_y)
            mid_nw_z = self._f(qtr1_x,qtr3_y)
            mid_se_z = self._f(qtr3_x,qtr1_y)
            mid_ne_z = self._f(qtr3_x,qtr3_y)

            # linearization estimates of midpoints
            est_sw_z = (mid_z + sw_z)/2
            est_nw_z = (mid_z + nw_z)/2
            est_se_z = (mid_z + se_z)/2
            est_ne_z = (mid_z + ne_z)/2

            tol_check = [(est_sw_z, mid_sw_z), (est_nw_z, mid_nw_z), (est_se_z, mid_se_z), (est_ne_z, mid_ne_z)]

            if depth < self._min_depth or not self.tol_list(tol_check):
                next_depth = depth + 1
            else:
                #lie about the depth to halt recursion
                next_depth = self._max_depth

            # recurse into the sub-squares
            sw = self.plot_block(min_x, qtr1_x, mid_x, min_y, qtr1_y, mid_y, sw_z, mid_w_z, mid_s_z, mid_z, mid_sw_z, next_depth)
            nw = self.plot_block(min_x, qtr1_x, mid_x, mid_y, qtr3_y, max_y, mid_w_z, nw_z, mid_z, mid_n_z, mid_nw_z, next_depth)
            se = self.plot_block(mid_x, qtr3_x, max_x, min_y, qtr1_y, mid_y, mid_s_z, mid_z, se_z, mid_e_z, mid_se_z, next_depth)
            ne = self.plot_block(mid_x, qtr3_x, max_x, mid_y, qtr3_y, max_y, mid_z, mid_n_z, mid_e_z, ne_z, mid_ne_z, next_depth)

            # join the sub-squares
            self.interface(1, sw.right, sw.right_c, se.left, se.left_c)
            self.interface(1, nw.right, nw.right_c, ne.left, ne.left_c)
            self.interface(0, sw.top, sw.top_c, nw.bottom, nw.bottom_c)
            self.interface(0, se.top, se.top_c, ne.bottom, ne.bottom_c)

            #get the boundary information about the subsquares
            left     = sw.left     + nw.left[1:]
            left_c   = sw.left_c   + nw.left_c
            right    = se.right    + ne.right[1:]
            right_c  = se.right_c  + ne.right_c
            top      = nw.top      + ne.top[1:]
            top_c    = nw.top_c    + ne.top_c
            bottom   = sw.bottom   + se.bottom[1:]
            bottom_c = sw.bottom_c + se.bottom_c

        else:
            # just build the square we're in
            sw = (min_x,min_y,sw_z)
            nw = (min_x,max_y,nw_z)
            se = (max_x,min_y,se_z)
            ne = (max_x,max_y,ne_z)
            c  = [(mid_x,mid_y,mid_z)]

            left = [sw,nw]
            left_c = c
            top = [nw,ne]
            top_c = c
            right = [se,ne]
            right_c = c
            bottom = [sw,se]
            bottom_c = c

        return PlotBlock(left, left_c, top, top_c, right, right_c, bottom, bottom_c)

    def tol(self, (est, val)):
        # Check relative, then absolute tolerance.  If both fail, return False
        # This is a zero-safe error checker

        if abs(est - val) < self._e_rel*abs(val):
            return True
        if abs(est - val) < self._e_abs:
            return True
        return False

    def tol_list(self, l):
        # Pass in a list of pairs of numbers, (est, val) to be passed to self.tol
        # returns False if any pair does not fall within tolerance level

        for p in l:
            if not self.tol(p):
                return False
        return True

    def interface(self, n, p, p_c, q, q_c):
        # Takes a pair of lists of points, and compares the (n)th coordinate, and
        # "zips" the lists together into one.  The "centers", supplied in p_c and
        # q_c are matched up such that the lists describe triangles whose sides
        # are "perfectly" aligned.  This algorithm assumes that p and q start and
        # end at the same point, and are sorted smallest to largest.

        m   = [p[0]] # a sorted union of p and q
        mpc = [p_c[0]] # centers from p_c corresponding to m
        mqc = [q_c[0]] # centers from q_c corresponding to m

        i = 1
        j = 1

        while i < len(p_c) or j < len(q_c):
            if abs(p[i][n] - q[j][n]) < self._eps:
                m.append(p[i])
                mpc.append(p_c[i])
                mqc.append(q_c[j])
                i += 1
                j += 1
            elif p[i][n] < q[j][n]:
                m.append(p[i])
                mpc.append(p_c[i])
                mqc.append(mqc[-1])
                i += 1
            else:
                m.append(q[j])
                mpc.append(mpc[-1])
                mqc.append(q_c[j])
                j += 1

        m.append(p[-1])

        self.triangulate(m, mpc)
        self.triangulate(m, mqc)


    def triangulate(self, p, c):
        # Pass in a list of edge points (p) and center points (c).
        # Triangles will be rendered between consecutive edge points and the
        # center point with the same index number as the earlier edge point.

        for i in range(0,len(p)-1):
            self._tachyon.triangle(p[i], p[i+1], c[i], self._tex)


class PlotBlock:
   def __init__(self, left, left_c, top, top_c, right, right_c, bottom, bottom_c):
       self.left = left
       self.left_c = left_c
       self.top = top
       self.top_c = top_c
       self.right = right
       self.right_c = right_c
       self.bottom = bottom
       self.bottom_c = bottom_c



def tostr(s):
    if isinstance(s, str):
        return s
    return ' %s %s %s '%(float(s[0]), float(s[1]), float(s[2]))



def hue(h, s=1, v=1):
    """
      hue(h,s=1,v=1) where 'h' stands for hue,
      's' stands for saturation, 'v' stands for value.
      hue returns a list of rgb intensities (r, g, b)
      All values are in range 0 to 1.

      INPUT:
         h, s, v -- real numbers between 0 and 1.  Note that
                    if any are not in this range they are automatically
                    normalized to be in this range by reducing them
                    modulo 1.
      OUTPUT:
         A valid RGB tuple.

      EXAMPLES:
        sage: hue(0.6)
        (0.0, 0.40000000000000036, 1.0)

        hue is an easy way of getting a broader
        range of colors for graphics

        sage: p = plot(sin, -2, 2, rgbcolor=hue(0.6))

    """
    h = float(h); s = float(s); v = float(v)
    if h != 1:
        h = modf(h)[0]
        if h < 0:
            h += 1
    if s != 1:
        s = modf(s)[0]
        if s < 0:
            s += 1
    if v != 1:
        v = modf(v)[0]
        if v < 0:
            v += 1
    c = hsv_to_rgb(h, s, v)
    return (float(c[0]), float(c[1]), float(c[2]))
