.. -*- coding: utf-8 -*-

.. linkall

.. _vector_calc_cylindrical:


Vector calculus 3: Cylindrical coordinates
==========================================

This tutorial introduces some vector calculus functionalities of SageMath
within the 3-dimensional Euclidean space.
The corresponding tools have been developed via the
`SageManifolds <https://sagemanifolds.obspm.fr>`__ project.

The tutorial is also available as a Jupyter notebook, either
`passive <https://nbviewer.jupyter.org/github/sagemanifolds/SageManifolds/blob/master/Notebooks/SM_vector_calc_cylindrical.ipynb>`__ (``nbviewer``)
or `interactive <https://mybinder.org/v2/gh/sagemanifolds/SageManifolds/master?filepath=Notebooks/SM_vector_calc_cylindrical.ipynb>`__ (``binder``).


The 3-dimensional Euclidean space
---------------------------------

We start by declaring the 3-dimensional Euclidean space
:math:`\mathbb{E}^3`, with :math:`(\rho,\phi,z)` as cylindrical
coordinates:

::

    sage: E.<rh,ph,z> = EuclideanSpace(coordinates='cylindrical')
    sage: E
    Euclidean space E^3

:math:`\mathbb{E}^3` is endowed with the chart of cylindrical
coordinates:

::

    sage: E.atlas()
    [Chart (E^3, (rh, ph, z))]

as well as with the associated orthonormal vector frame
:math:`(e_\rho, e_\phi, e_z)`:

::

    sage: E.frames()
    [Coordinate frame (E^3, (d/drh,d/dph,d/dz)),
     Vector frame (E^3, (e_rh,e_ph,e_z))]

In the above output, ``(d/drh,d/dph,d/dz)`` =
:math:`\left(\frac{\partial}{\partial\rho}, \frac{\partial}{\partial\phi}, \frac{\partial}{\partial z}\right)`
is the coordinate frame associated with :math:`(\rho,\phi,z)`; it is not
an orthonormal frame and will not be used below.

Vector fields
-------------

We define a vector field on :math:`\mathbb{E}^3` from its components in
the orthonormal vector frame :math:`(e_\rho,e_\phi,e_z)`:

::

    sage: v = E.vector_field(rh*(1+sin(2*ph)), 2*rh*cos(ph)^2, z,
    ....:                    name='v')
    sage: v.display()
    v = rh*(sin(2*ph) + 1) e_rh + 2*rh*cos(ph)^2 e_ph + z e_z

We can access to the components of :math:`v` via the square bracket
operator:

::

    sage: v[1]
    rh*(sin(2*ph) + 1)
    sage: v[:]
    [rh*(sin(2*ph) + 1), 2*rh*cos(ph)^2, z]

A vector field can evaluated at any point of :math:`\mathbb{E}^3`:

::

    sage: p = E((1, pi, 0), name='p')
    sage: p
    Point p on the Euclidean space E^3
    sage: p.coordinates()
    (1, pi, 0)
    sage: vp = v.at(p)
    sage: vp
    Vector v at Point p on the Euclidean space E^3
    sage: vp.display()
    v = e_rh + 2 e_ph

We may define a vector field with generic components:

::

    sage: u = E.vector_field(function('u_rho')(rh,ph,z),
    ....:                    function('u_phi')(rh,ph,z),
    ....:                    function('u_z')(rh,ph,z),
    ....:                    name='u')
    sage: u.display()
    u = u_rho(rh, ph, z) e_rh + u_phi(rh, ph, z) e_ph + u_z(rh, ph, z) e_z
    sage: u[:]
    [u_rho(rh, ph, z), u_phi(rh, ph, z), u_z(rh, ph, z)]
    sage: up = u.at(p)
    sage: up.display()
    u = u_rho(1, pi, 0) e_rh + u_phi(1, pi, 0) e_ph + u_z(1, pi, 0) e_z


Algebraic operations on vector fields
-------------------------------------

Dot product
~~~~~~~~~~~

The dot (or scalar) product of the vector fields :math:`u` and :math:`v`
is obtained by the method ``dot_product``, which admits ``dot`` as a
shortcut alias:

::

    sage: s = u.dot(v)
    sage: s
    Scalar field u.v on the Euclidean space E^3

:math:`s= u\cdot v` is a *scalar field*, i.e. a map
:math:`\mathbb{E}^3 \rightarrow \mathbb{R}`:

::

    sage: s.display()
    u.v: E^3 --> R
       (rh, ph, z) |--> 2*rh*cos(ph)^2*u_phi(rh, ph, z)
        + (2*cos(ph)*sin(ph) + 1)*rh*u_rho(rh, ph, z) + z*u_z(rh, ph, z)

It maps points of :math:`\mathbb{E}^3` to real numbers:

::

    sage: s(p)
    2*u_phi(1, pi, 0) + u_rho(1, pi, 0)

Its coordinate expression is

::

    sage: s.expr()
    2*rh*cos(ph)^2*u_phi(rh, ph, z) + (2*cos(ph)*sin(ph) + 1)*rh*u_rho(rh, ph, z)
     + z*u_z(rh, ph, z)


Norm
~~~~

The norm of a vector field is

::

    sage: s = norm(u)
    sage: s
    Scalar field |u| on the Euclidean space E^3
    sage: s.display()
    |u|: E^3 --> R
       (rh, ph, z) |--> sqrt(u_phi(rh, ph, z)^2 + u_rho(rh, ph, z)^2 + u_z(rh, ph, z)^2)
    sage: s.expr()
    sqrt(u_phi(rh, ph, z)^2 + u_rho(rh, ph, z)^2 + u_z(rh, ph, z)^2)

The norm is related to the dot product by :math:`\|u\|^2 = u\cdot u`, as
we can check:

::

    sage: norm(u)^2 == u.dot(u)
    True

For :math:`v`, we have

::

    sage: norm(v).expr()
    sqrt((4*cos(ph)^2 + 4*cos(ph)*sin(ph) + 1)*rh^2 + z^2)


Cross product
~~~~~~~~~~~~~

The cross product of :math:`u` by :math:`v` is obtained by the method
``cross_product``, which admits ``cross`` as a shortcut alias:

::

    sage: s = u.cross(v)
    sage: s
    Vector field u x v on the Euclidean space E^3
    sage: s.display()
    u x v = (-2*rh*cos(ph)^2*u_z(rh, ph, z) + z*u_phi(rh, ph, z)) e_rh
     + ((2*cos(ph)*sin(ph) + 1)*rh*u_z(rh, ph, z) - z*u_rho(rh, ph, z)) e_ph
     + (2*rh*cos(ph)^2*u_rho(rh, ph, z) - (2*cos(ph)*sin(ph) + 1)*rh*u_phi(rh, ph, z)) e_z


Scalar triple product
~~~~~~~~~~~~~~~~~~~~~

Let us introduce a third vector field. As a example, we do not pass the
components as arguments of ``vector_field``, as we did for :math:`u` and
:math:`v`; instead, we set them in a second stage, via the square
bracket operator, any unset component being assumed to be zero:

::

    sage: w = E.vector_field(name='w')
    sage: w[1] = rh
    sage: w[3] = z
    sage: w.display()
    w = rh e_rh + z e_z

The scalar triple product of the vector fields :math:`u`, :math:`v` and
:math:`w` is obtained as follows:

::

    sage: triple_product = E.scalar_triple_product()
    sage: s = triple_product(u, v, w)
    sage: s
    Scalar field epsilon(u,v,w) on the Euclidean space E^3
    sage: s.expr()
    -2*rh^2*cos(ph)^2*u_z(rh, ph, z) - 2*(rh*cos(ph)*sin(ph)*u_phi(rh, ph, z)
     - rh*cos(ph)^2*u_rho(rh, ph, z))*z

Let us check that the scalar triple product of :math:`u`, :math:`v` and
:math:`w` is :math:`u\cdot(v\times w)`:

::

    sage: s == u.dot(v.cross(w))
    True


Differential operators
----------------------

While the standard operators :math:`\mathrm{grad}`,
:math:`\mathrm{div}`, :math:`\mathrm{curl}`, etc. involved in vector
calculus are accessible via the dot notation (e.g. ``v.div()``), let us
import functions ``grad``, ``div``, ``curl``, etc. that allow for using
standard mathematical notations (e.g. ``div(v)``):

::

    sage: from sage.manifolds.operators import *


Gradient of a scalar field
~~~~~~~~~~~~~~~~~~~~~~~~~~

We first introduce a scalar field, via its expression in terms of
Cartesian coordinates; in this example, we consider a unspecified
function of :math:`(\rho,\phi,z)`:

::

    sage: F = E.scalar_field(function('f')(rh,ph,z), name='F')
    sage: F.display()
    F: E^3 --> R
       (rh, ph, z) |--> f(rh, ph, z)

The value of :math:`F` at a point:

::

    sage: F(p)
    f(1, pi, 0)

The gradient of :math:`F`:

::

    sage: grad(F)
    Vector field grad(F) on the Euclidean space E^3
    sage: grad(F).display()
    grad(F) = d(f)/drh e_rh + d(f)/dph/rh e_ph + d(f)/dz e_z
    sage: norm(grad(F)).display()
    |grad(F)|: E^3 --> R
       (rh, ph, z) |--> sqrt(rh^2*(d(f)/drh)^2 + rh^2*(d(f)/dz)^2 + (d(f)/dph)^2)/rh


Divergence
~~~~~~~~~~

The divergence of a vector field:

::

    sage: s = div(u)
    sage: s.display()
    div(u): E^3 --> R
       (rh, ph, z) |--> (rh*d(u_rho)/drh + rh*d(u_z)/dz + u_rho(rh, ph, z) + d(u_phi)/dph)/rh
    sage: s.expr().expand()
    u_rho(rh, ph, z)/rh + diff(u_phi(rh, ph, z), ph)/rh + diff(u_rho(rh, ph, z), rh)
     + diff(u_z(rh, ph, z), z)

For :math:`v` and :math:`w`, we have

::

    sage: div(v).expr()
    3
    sage: div(w).expr()
    3

An identity valid for any scalar field :math:`F` and any vector field
:math:`u`:

::

    sage: div(F*u) == F*div(u) + u.dot(grad(F))
    True


Curl
~~~~

The curl of a vector field:

::

    sage: s = curl(u)
    sage: s
    Vector field curl(u) on the Euclidean space E^3
    sage: s.display()
    curl(u) = -(rh*d(u_phi)/dz - d(u_z)/dph)/rh e_rh + (d(u_rho)/dz - d(u_z)/drh) e_ph
     + (rh*d(u_phi)/drh + u_phi(rh, ph, z) - d(u_rho)/dph)/rh e_z

To use the notation ``rot`` instead of ``curl``, simply do

::

    sage: rot = curl

An alternative is

::

    sage: from sage.manifolds.operators import curl as rot

We have then

::

    sage: rot(u).display()
    curl(u) = -(rh*d(u_phi)/dz - d(u_z)/dph)/rh e_rh + (d(u_rho)/dz - d(u_z)/drh) e_ph
     + (rh*d(u_phi)/drh + u_phi(rh, ph, z) - d(u_rho)/dph)/rh e_z
    sage: rot(u) == curl(u)
    True

For :math:`v` and :math:`w`, we have

::

    sage: curl(v).display()
    curl(v) = 2 e_z
    sage: curl(w).display()
    curl(w) = 0

The curl of a gradient is always zero:

::

    sage: curl(grad(F)).display()
    curl(grad(F)) = 0

The divergence of a curl is always zero:

::

    sage: div(curl(u)).display()
    div(curl(u)): E^3 --> R
       (rh, ph, z) |--> 0

An identity valid for any scalar field :math:`F` and any vector field
:math:`u`:

::

    sage: curl(F*u) == grad(F).cross(u) + F*curl(u)
    True


Laplacian
~~~~~~~~~

The Laplacian of a scalar field:

::

    sage: s = laplacian(F)
    sage: s.display()
    Delta(F): E^3 --> R
       (rh, ph, z) |--> (rh^2*d^2(f)/drh^2 + rh^2*d^2(f)/dz^2 + rh*d(f)/drh
        + d^2(f)/dph^2)/rh^2
    sage: s.expr().expand()
    diff(f(rh, ph, z), rh)/rh + diff(f(rh, ph, z), ph, ph)/rh^2
     + diff(f(rh, ph, z), rh, rh) + diff(f(rh, ph, z), z, z)

For a scalar field, the Laplacian is nothing but the divergence of the
gradient:

::

    sage: laplacian(F) == div(grad(F))
    True

The Laplacian of a vector field:

::

    sage: Du = laplacian(u)
    sage: Du.display()
    Delta(u) = (rh^2*d^2(u_rho)/drh^2 + rh^2*d^2(u_rho)/dz^2 + rh*d(u_rho)/drh
     - u_rho(rh, ph, z) - 2*d(u_phi)/dph + d^2(u_rho)/dph^2)/rh^2 e_rh
     + (rh^2*d^2(u_phi)/drh^2 + rh^2*d^2(u_phi)/dz^2 + rh*d(u_phi)/drh
     - u_phi(rh, ph, z) + d^2(u_phi)/dph^2 + 2*d(u_rho)/dph)/rh^2 e_ph
     + (rh^2*d^2(u_z)/drh^2 + rh^2*d^2(u_z)/dz^2 + rh*d(u_z)/drh
     + d^2(u_z)/dph^2)/rh^2 e_z

Since this expression is quite lengthy, we may ask for a display
component by component:

::

    sage: Du.display_comp()
    Delta(u)^1 = (rh^2*d^2(u_rho)/drh^2 + rh^2*d^2(u_rho)/dz^2 + rh*d(u_rho)/drh
     - u_rho(rh, ph, z) - 2*d(u_phi)/dph + d^2(u_rho)/dph^2)/rh^2
    Delta(u)^2 = (rh^2*d^2(u_phi)/drh^2 + rh^2*d^2(u_phi)/dz^2 + rh*d(u_phi)/drh
     - u_phi(rh, ph, z) + d^2(u_phi)/dph^2 + 2*d(u_rho)/dph)/rh^2
    Delta(u)^3 = (rh^2*d^2(u_z)/drh^2 + rh^2*d^2(u_z)/dz^2 + rh*d(u_z)/drh
     + d^2(u_z)/dph^2)/rh^2

We may expand each component:

::

    sage: for i in E.irange():
    ....:     s = Du[i].expand()
    sage: Du.display_comp()
    Delta(u)^1 = d(u_rho)/drh/rh - u_rho(rh, ph, z)/rh^2 - 2*d(u_phi)/dph/rh^2
     + d^2(u_rho)/dph^2/rh^2 + d^2(u_rho)/drh^2 + d^2(u_rho)/dz^2
    Delta(u)^2 = d(u_phi)/drh/rh - u_phi(rh, ph, z)/rh^2 + d^2(u_phi)/dph^2/rh^2
     + 2*d(u_rho)/dph/rh^2 + d^2(u_phi)/drh^2 + d^2(u_phi)/dz^2
    Delta(u)^3 = d(u_z)/drh/rh + d^2(u_z)/dph^2/rh^2 + d^2(u_z)/drh^2 + d^2(u_z)/dz^2

::

    sage: Du[1]
    d(u_rho)/drh/rh - u_rho(rh, ph, z)/rh^2 - 2*d(u_phi)/dph/rh^2
     + d^2(u_rho)/dph^2/rh^2 + d^2(u_rho)/drh^2 + d^2(u_rho)/dz^2
    sage: Du[2]
    d(u_phi)/drh/rh - u_phi(rh, ph, z)/rh^2 + d^2(u_phi)/dph^2/rh^2
     + 2*d(u_rho)/dph/rh^2 + d^2(u_phi)/drh^2 + d^2(u_phi)/dz^2
    sage: Du[3]
    d(u_z)/drh/rh + d^2(u_z)/dph^2/rh^2 + d^2(u_z)/drh^2 + d^2(u_z)/dz^2

As a test, we may check that these formulas coincide with those of
Wikipedia's article `*Del in cylindrical and spherical
coordinates* <https://en.wikipedia.org/wiki/Del_in_cylindrical_and_spherical_coordinates#Del_formula>`__.

For :math:`v` and :math:`w`, we have

::

    sage: laplacian(v).display()
    Delta(v) = 0
    sage: laplacian(w).display()
    Delta(w) = 0

We have

::

    sage: curl(curl(u)).display()
    curl(curl(u)) = -(rh^2*d^2(u_rho)/dz^2 - rh^2*d^2(u_z)/drhdz
     - rh*d^2(u_phi)/drhdph - d(u_phi)/dph + d^2(u_rho)/dph^2)/rh^2 e_rh
     - (rh^2*d^2(u_phi)/drh^2 + rh^2*d^2(u_phi)/dz^2 + rh*d(u_phi)/drh
     - rh*d^2(u_rho)/drhdph - rh*d^2(u_z)/dphdz - u_phi(rh, ph, z)
     + d(u_rho)/dph)/rh^2 e_ph + (rh^2*d^2(u_rho)/drhdz - rh^2*d^2(u_z)/drh^2
     + rh*d^2(u_phi)/dphdz + rh*d(u_rho)/dz - rh*d(u_z)/drh - d^2(u_z)/dph^2)/rh^2 e_z
    sage: grad(div(u)).display()
    grad(div(u)) = (rh^2*d^2(u_rho)/drh^2 + rh^2*d^2(u_z)/drhdz + rh*d^2(u_phi)/drhdph
     + rh*d(u_rho)/drh - u_rho(rh, ph, z) - d(u_phi)/dph)/rh^2 e_rh
     + (rh*d^2(u_rho)/drhdph + rh*d^2(u_z)/dphdz + d^2(u_phi)/dph^2
     + d(u_rho)/dph)/rh^2 e_ph + (rh*d^2(u_rho)/drhdz + rh*d^2(u_z)/dz^2
     + d^2(u_phi)/dphdz + d(u_rho)/dz)/rh e_z

and we may check a famous identity:

::

    sage: curl(curl(u)) == grad(div(u)) - laplacian(u)
    True


Customizations
--------------

Customizing the symbols of the orthonormal frame vectors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the vectors of the orthonormal frame associated with
cylindrical coordinates are denoted :math:`(e_\rho,e_\phi,z)`:

::

    sage: frame = E.cylindrical_frame()
    sage: frame
    Vector frame (E^3, (e_rh,e_ph,e_z))

But this can be changed, thanks to the method ``set_name``:

::

    sage: frame.set_name('a', indices=('rh', 'ph', 'z'),
    ....:                latex_indices=(r'\rho', r'\phi', 'z'))
    sage: frame
    Vector frame (E^3, (a_rh,a_ph,a_z))
    sage: v.display()
    v = rh*(sin(2*ph) + 1) a_rh + 2*rh*cos(ph)^2 a_ph + z a_z

::

    sage: frame.set_name(('hrh', 'hph', 'hz'),
    ....:                latex_symbol=(r'\hat{\rho}', r'\hat{\phi}', r'\hat{z}'))
    sage: frame
    Vector frame (E^3, (hrh,hph,hz))
    sage: v.display()
    v = rh*(sin(2*ph) + 1) hrh + 2*rh*cos(ph)^2 hph + z hz


Customizing the coordinate symbols
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The coordinates symbols are defined within the angle brackets ``<...>``
at the construction of the Euclidean space. Above we did

::

    sage: E.<rh,ph,z> = EuclideanSpace(coordinates='cylindrical')

which resulted in the coordinate symbols :math:`(\rho,\phi,z)` and in
the corresponding Python variables ``rh``, ``ph`` and ``z`` (SageMath
symbolic expressions). Using other symbols, for instance
:math:`(R,\Phi,Z)`, is possible through the optional argument
``symbols`` of the function ``EuclideanSpace``. It has to be a string,
usually prefixed by ``r`` (for *raw* string, in order to allow for the
backslash character of LaTeX expressions). This string contains the
coordinate fields separated by a blank space; each field contains the
coordinate’s text symbol and possibly the coordinate’s LaTeX symbol
(when the latter is different from the text symbol), both symbols being
separated by a colon (``:``):

::

    sage: E.<R,Ph,Z> = EuclideanSpace(coordinates='cylindrical', symbols=r'R Ph:\Phi Z')

We have then

::

    sage: E.atlas()
    [Chart (E^3, (R, Ph, Z))]
    sage: E.frames()
    [Coordinate frame (E^3, (d/dR,d/dPh,d/dZ)), Vector frame (E^3, (e_R,e_Ph,e_Z))]
    sage: E.cylindrical_frame()
    Vector frame (E^3, (e_R,e_Ph,e_Z))
    sage: v = E.vector_field(R*(1+sin(2*Ph)), 2*R*cos(Ph)^2, Z, name='v')
    sage: v.display()
    v = R*(sin(2*Ph) + 1) e_R + 2*R*cos(Ph)^2 e_Ph + Z e_Z
