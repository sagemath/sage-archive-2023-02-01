"""
Color and texture mapping on surfaces

The class below allows for a texture attribute of a surface to be varied
gradually. The precise way in which the attribute varies can be
specified by means of a scalar function, the most straightforward ones
being projection onto the x, y, or z-axis.  In this way, it is
possible to apply a colormap to a surface, to create partially
translucent objects, or to color in surfaces in a variety of ways (for
instance, to visualize Gaussian curvature).

Please note that similar colorings of surfaces can also be obtained
using the keyword ``color_data`` of the functions
:func:`sage.plot.plot3d.parametric_plot3d.parametric_plot3d` and
:func:`sage.plot.plot3d.implicit_plot3d.implicit_plot3d`.

Let us first load the methods and define some colormaps::

    sage: from sage.plot.plot3d.texture_map import ColormapTransform,OpacityTransform
    sage: cmsel1 = colormaps.autumn
    sage: cmsel2 = colormaps.spectral

Now for a first example, a wavy surface and a sphere::

    sage: r, v = var('r,v')
    sage: c = ColormapTransform(cmsel1)
    sage: p1 = plot3d(0.2*(r**2 + v**2) + cos(2*r)*sin(2*v),(r,-2,2), (v,-2,2), opacity=0.9)
    sage: p2 = sphere((0,0,0),1,color='black',opacity=0.5)
    sage: p3 = c.apply(p1)
    sage: (p2+p3).show(aspect_ratio=(1,1,1))

An example of an implicit surface with a color map::

    sage: c = ColormapTransform(cmsel1)
    sage: x, y, z = var('x,y,z')
    sage: p1 = implicit_plot3d(x^2+y^2+z^2==4, (x, -3, 3), (y, -3,3), (z, -3,3))
    sage: p2 = c.apply(p1)
    sage: p2.show()

By default, the color values are selected according to the z-values of
the points on the surface. Other functions are possible, too.

One can for instance use a uniform gradient along the x-axis::

    sage: r, v = var('r,v')
    sage: p = plot3d(0.2*(r**2 + v**2) + cos(2*r)*sin(2*v),(r,-2,2), (v,-2,2))
    sage: c = ColormapTransform(cmsel1, 'height_x')
    sage: c.apply(p).show(aspect_ratio=(1,1,1))

A non-uniform coloring, where a function is specified::

    sage: x, y, z = var('x,y,z')
    sage: p = plot3d(0.2*(r**2 + v**2) + cos(2*r)*sin(2*v),(r,-2,2), (v,-2,2))
    sage: c = ColormapTransform(fun=lambda x, y, z: x+y, bounds=(-4, 4), cmap=cmsel1)
    sage: c.apply(p)

Another non-uniform coloring, this time with a nonlinear function::

    sage: r, v = var('r,v')
    sage: p = implicit_plot3d(x^2+y^2+z^2==4, (x, -3, 3), (y, -3,3), (z, -3,3))
    sage: c = ColormapTransform(cmsel1, fun=lambda x, y, z: cos(2*pi*z/3-pi/3), bounds=(-1, 1))
    sage: c.apply(p).show(aspect_ratio=(1,1,1))

A more complicated surface::

    sage: u, v = var('u,v')
    sage: fx = u -u^3/3  + u*v^2
    sage: fy = v -v^3/3  + v*u^2
    sage: fz = u^2 - v^2
    sage: p = parametric_plot3d([fx, fy, fz], (u, -2, 2), (v, -2, 2))
    sage: c = ColormapTransform(cmsel2)
    sage: c.apply(p)

We can also gradually vary other texture attributes, such as opacity::

    sage: p1 = sphere((0,0,0), .5, color='red')
    sage: p2 = sphere((0,0,0),  1)
    sage: opacity_values = sxrange(0, 1, .05)
    sage: c = OpacityTransform(opacity_values)
    sage: p3 = c.apply(p2)
    sage: (p1 + p3).show()

The next three examples show implicit surfaces colored by Gaussian curvature.

Let us first define a function computing Gaussian curvature::

    sage: def gaussian_curvature_implicit_surface(f):
    ....:     fx = diff(f, x); fy = diff(f, y); fz = diff(f, z)
    ....:     fxx = diff(fx, x); fxy = diff(fx, y); fxz = diff(fx, z)
    ....:     fyy = diff(fy, y); fyz = diff(fy, z); fzz = diff(fz, z)
    ....:     T1 = (fz*(fxx*fz - 2*fx*fxz) + fx^2*fzz).simplify_full()
    ....:     T2 = (fz*(fyy*fz - 2*fy*fyz) + fy^2*fzz).simplify_full()
    ....:     T3 = (fz*(-fx*fyz + fxy*fz - fxz*fy) + fx*fy*fzz).simplify_full()
    ....:     T4 = (fz^2*(fx^2 + fy^2 + fz^2)^2).simplify_full()
    ....:     K = ((T1*T2 - T3^2)/T4).simplify_full()
    ....:     return fast_float(K, 'x', 'y', 'z')

Our first example is an ellipsoid (positive curvature)::

    sage: f = x^2/4 + y^2 + z^2 - 1
    sage: p = implicit_plot3d(f, (x, -2, 2), (y, -1, 1), (z, -1, 1))
    sage: K = gaussian_curvature_implicit_surface(f)
    sage: c = ColormapTransform(cmsel1, fun=K)
    sage: c.apply(p).show(aspect_ratio=(1, 1, 1))

Our second example is a saddle surface (negative curvature at the origin)::

    sage: f = z+x^2-y^2
    sage: p = implicit_plot3d(f, (x, -2, 2), (y, -2, 2), (z, -2, 2), plot_points=60, contour=0)
    sage: K = gaussian_curvature_implicit_surface(f)
    sage: c = ColormapTransform(colormaps.gist_rainbow, fun=K, bounds=(-4, 0))
    sage: c.apply(p)

Now a more complicated example::

    sage: T = RDF(golden_ratio)
    sage: f = 2 - (cos(x + T*y) + cos(x - T*y) + cos(y + T*z)
    ....:   + cos(y - T*z) + cos(z - T*x) + cos(z + T*x))
    sage: r = 4.77
    sage: p = implicit_plot3d(f, (x, -r, r), (y, -r, r), (z, -r, r),
    ....:  plot_points=60); p
    sage: K = gaussian_curvature_implicit_surface(f) # long time
    sage: c = ColormapTransform(cmsel1, fun=K)  # long time
    sage: c.apply(p)   # long time

AUTHOR:

- Joris Vankerschaver
"""
from sage.plot.plot3d.base import Graphics3dGroup
from sage.structure.sage_object import SageObject
from matplotlib.colors import LinearSegmentedColormap
from sage.misc.misc import sxrange


class GradualTextureTransform(SageObject):
    """
    General class to apply texture transformation on graphic objects
    """
    def __init__(self, attribute, values, fun=None, bounds=None):
        """
        INPUT:

        - ``attribute`` -- texture attribute to transform.  If this is
            not one of ``color``, ``opacity``, ``ambient``,
            ``diffuse``, ``specular`` or ``shininess``, the
            transformation has no effect.

        - ``values`` -- range of values over which the attribute
            should range.  In the case where the color attribute is
            transformed, this can be a colormap.

        - ``fun`` -- scalar-valued function specifying the change in
            attribute.  If this keyword is omitted, the projection x,
            y, z -> z will be used.  This could be either a symbolic
            expression, or one of the following: 'height_x',
            'height_y', 'height_z', denoting respectively the height
            above the x-, y-, or z-plane.

        - ``bounds`` -- tuple representing the range for the scalar
            function ``fun``.  If this keyword is omitted, some preset
            defaults will be used.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture_map import ColormapTransform
            sage: cmsel1 = colormaps.autumn
            sage: c = ColormapTransform(cmsel1, 'height_x')
            sage: x,y,z = var('x,y,z')
            sage: p1 = implicit_plot3d(x^2+y^2+z^2==4, (x, -3, 3), (y, -3,3), (z, -3,3))
            sage: p2 = c.apply(p1)
        """
        self.attribute = attribute
        self.fun = fun
        self.bounds = bounds
        if isinstance(values, LinearSegmentedColormap):
            self.values = [values(i) for i in sxrange(0, 1, 0.05)]
        else:
            self.values = list(values)

        self.projection_axis = 2

        if fun is None:
            self.fun = lambda x, y, z: z
            self.projection_axis = 2
        elif isinstance(fun, str):
            if fun is "height_x":
                self.fun = lambda x, y, z: x
                self.projection_axis = 0
            elif fun is "height_y":
                self.fun = lambda x, y, z: y
                self.projection_axis = 1
            elif fun is "height_z":
                self.fun = lambda x, y, z: z
                self.projection_axis = 2
            else:
                raise ValueError("Not a valid texture function")

    def apply(self, group):
        """
        Apply this texture transformation to a graphical object.

        INPUT:

        - ``group`` -- a graphical object (``Graphics3d``) or group of
          objects (``Graphics3dGroup``).

        OUTPUT:

        A ``Graphics3dGroup`` representing the transformed object.

        EXAMPLES::

            sage: from sage.plot.plot3d.texture_map import ColormapTransform
            sage: cmsel1 = colormaps.autumn
            sage: c = ColormapTransform(cmsel1)
            sage: x, y, z = var('x,y,z')
            sage: p1 = implicit_plot3d(x^2+y^2+z^2==4, (x, -3, 3), (y, -3,3), (z, -3,3))
            sage: p2 = c.apply(p1)
            sage: p2.show()
        """
        if not isinstance(group, Graphics3dGroup):
            group = Graphics3dGroup([group])

        # Code below adapted from plot3d_adaptive
        if self.bounds is None:
            # Determine limits for function from bounding box
            bounds = group.bounding_box()

            min_val = bounds[0][self.projection_axis]
            max_val = bounds[1][self.projection_axis]
        else:
            min_val = self.bounds[0]
            max_val = self.bounds[1]

        if max_val == min_val:
            span = 0
        else:
            span = (len(self.values) - 1) / (max_val - min_val)

        def normalized_fun(x, y, z):
            val = self.fun(x, y, z)
            if val > max_val:
                val = max_val
            if val < min_val:
                val = min_val
            return int((val - min_val) * span)

        new_all = []
        kwds = {}
        for obj in group.all:
            try:
                obj.triangulate()
                # Faces have to exist before object can be partitioned
            except AttributeError:
                # assume object is already triangulated
                # e.g. in the case obj is an IndexFaceSet
                pass
            parts = obj.partition(normalized_fun)

            for k, G in parts.iteritems():
                kwds[self.attribute] = self.values[k]
                texture_dict = (obj.get_texture().__dict__).copy()
                _ = texture_dict.pop('id',None) # want new uniquely-generated id
                texture_dict['name'] = None # forget name
                texture_dict.update(kwds)
                G.set_texture(texture_dict)
                new_all.append(G)

        return Graphics3dGroup(new_all)


class ColormapTransform(GradualTextureTransform):
    """
    Transformation to apply a colormap to a given surface.

    See :class:`GradualTextureTransform` for more details.

    EXAMPLES::

        sage: from sage.plot.plot3d.texture_map import ColormapTransform
        sage: cmsel1 = colormaps.autumn
        sage: c = ColormapTransform(cmsel1)
        sage: x, y, z = var('x,y,z')
        sage: p1 = implicit_plot3d(x^2+y^2+z^2==4, (x, -3, 3), (y, -3,3), (z, -3,3))
        sage: p2 = c.apply(p1)
        sage: p2.show()
    """
    def __init__(self, cmap, fun=None, bounds=None):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.texture_map import ColormapTransform
            sage: cms = colormaps.summer
            sage: c = ColormapTransform(cms)  # indirect doctest
            sage: x, y, z = var('x,y,z')
            sage: p1 = implicit_plot3d(x^2+2*y^2+z^2==4, (x, -3, 3), (y, -3,3), (z, -3,3))
            sage: p2 = c.apply(p1)
            sage: p2.show()
        """
        GradualTextureTransform.__init__(self, 'color', cmap, fun, bounds)


class OpacityTransform(GradualTextureTransform):
    """
    Transformation to change the opacity of given surface in a non-uniform way.

    See :class:`GradualTextureTransform` for more details.

    EXAMPLES::

        sage: from sage.plot.plot3d.texture_map import OpacityTransform
        sage: p1 = sphere((0,0,0), .5, color='red')
        sage: p2 = sphere((0,0,0),  1)
        sage: opacity_values = sxrange(0, 1, .05)
        sage: c = OpacityTransform(opacity_values)
        sage: p3 = c.apply(p2)
        sage: (p1 + p3).show()
    """
    def __init__(self, values, fun=None, bounds=None):
        """
        EXAMPLES::

            sage: from sage.plot.plot3d.texture_map import OpacityTransform
            sage: p1 = sphere((0,0.6,0), .5, color='yellow')
            sage: p2 = sphere((0,0,0),  1)
            sage: opacity_values = sxrange(0.1, 0.9, .05)
            sage: c = OpacityTransform(opacity_values)  # indirect doctest
            sage: p3 = c.apply(p2)
            sage: (p1 + p3).show()
        """
        GradualTextureTransform.__init__(self, 'opacity', values, fun, bounds)
