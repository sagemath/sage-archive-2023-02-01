.. -*- coding: utf-8 -*-

.. linkall

.. _vector_calc_spherical:


Vector calculus 2: spherical coordinates
========================================

This notebook illustrates some vector calculus capabilities of SageMath
within the 3-dimensional Euclidean space. The corresponding tools have
been developed within the
`SageManifolds <https://sagemanifolds.obspm.fr>`__ project.

Click
`here <https://raw.githubusercontent.com/sagemanifolds/SageManifolds/master/Worksheets/v1.3/SM_vector_calc_spherical.ipynb>`__
to download the corresponding Jupyter notebook file.


The 3-dimensional Euclidean space
---------------------------------

We start by declaring the 3-dimensional Euclidean space
:math:`\mathbb{E}^3`, with :math:`(r,\theta,\phi)` as spherical
coordinates:

::

    sage: E.<r,th,ph> = EuclideanSpace(coordinates='spherical')
    sage: E
    Euclidean space E^3

:math:`\mathbb{E}^3` is endowed with the chart of spherical coordinates:

::

    sage: E.atlas()
    [Chart (E^3, (r, th, ph))]

as well as with the associated orthonormal vector frame
:math:`(e_r, e_\theta, e_\phi)`:

::

    sage: E.frames()
    [Coordinate frame (E^3, (d/dr,d/dth,d/dph)),
     Vector frame (E^3, (e_r,e_th,e_ph))]

In the above output, ``(d/dr,d/dth,d/dph)`` =
:math:`\left(\frac{\partial}{\partial r}, \frac{\partial}{\partial\theta}, \frac{\partial}{\partial \phi}\right)`
is the coordinate frame associated with :math:`(r,\theta,\phi)`; it is
not an orthonormal frame and will not be used below.

Vector fields
-------------

We define a vector field on :math:`\mathbb{E}^3` from its components in
the orthonormal vector frame :math:`(e_r,e_\theta,e_\phi)`:

::

    sage: v = E.vector_field(r*sin(2*ph)*sin(th)^2 + r,
    ....:                    r*sin(2*ph)*sin(th)*cos(th),
    ....:                    2*r*cos(ph)^2*sin(th), name='v')
    sage: v.display()
    v = (r*sin(2*ph)*sin(th)^2 + r) e_r + r*cos(th)*sin(2*ph)*sin(th) e_th
     + 2*r*cos(ph)^2*sin(th) e_ph

We can access to the components of :math:`v` via the square bracket
operator:

::

    sage: v[1]
    r*sin(2*ph)*sin(th)^2 + r
    sage: v[:]
    [r*sin(2*ph)*sin(th)^2 + r, r*cos(th)*sin(2*ph)*sin(th), 2*r*cos(ph)^2*sin(th)]

A vector field can evaluated at any point of :math:`\mathbb{E}^3`:

::

    sage: p = E((1, pi/2, pi), name='p')
    sage: p
    Point p on the Euclidean space E^3
    sage: p.coordinates()
    (1, 1/2*pi, pi)
    sage: vp = v.at(p)
    sage: vp
    Vector v at Point p on the Euclidean space E^3
    sage: vp.display()
    v = e_r + 2 e_ph

We may define a vector field with generic components:

::

    sage: u = E.vector_field(function('u_r')(r,th,ph),
    ....:                    function('u_theta')(r,th,ph),
    ....:                    function('u_phi')(r,th,ph),
    ....:                    name='u')
    sage: u.display()
    u = u_r(r, th, ph) e_r + u_theta(r, th, ph) e_th + u_phi(r, th, ph) e_ph
    sage: u[:]
    [u_r(r, th, ph), u_theta(r, th, ph), u_phi(r, th, ph)]
    sage: up = u.at(p)
    sage: up.display()
    u = u_r(1, 1/2*pi, pi) e_r + u_theta(1, 1/2*pi, pi) e_th
     + u_phi(1, 1/2*pi, pi) e_ph


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
       (r, th, ph) |--> 2*r*cos(ph)*sin(ph)*sin(th)^2*u_r(r, th, ph) +
        2*(r*cos(ph)*cos(th)*sin(ph)*u_theta(r, th, ph)
        + r*cos(ph)^2*u_phi(r, th, ph))*sin(th) + r*u_r(r, th, ph)

It maps points of :math:`\mathbb{E}^3` to real numbers:

::

    sage: s(p)
    2*u_phi(1, 1/2*pi, pi) + u_r(1, 1/2*pi, pi)

Its coordinate expression is

::

    sage: s.expr()
    2*r*cos(ph)*sin(ph)*sin(th)^2*u_r(r, th, ph)
     + 2*(r*cos(ph)*cos(th)*sin(ph)*u_theta(r, th, ph)
     + r*cos(ph)^2*u_phi(r, th, ph))*sin(th) + r*u_r(r, th, ph)


Norm
~~~~

The norm of a vector field is

::

    sage: s = norm(u)
    sage: s
    Scalar field |u| on the Euclidean space E^3
    sage: s.display()
    |u|: E^3 --> R
       (r, th, ph) |--> sqrt(u_phi(r, th, ph)^2 + u_r(r, th, ph)^2 + u_theta(r, th, ph)^2)
    sage: s.expr()
    sqrt(u_phi(r, th, ph)^2 + u_r(r, th, ph)^2 + u_theta(r, th, ph)^2)

The norm is related to the dot product by :math:`\|u\|^2 = u\cdot u`, as
we can check:

::

    sage: norm(u)^2 == u.dot(u)
    True

For :math:`v`, we have

::

    sage: norm(v).expr()
    sqrt(4*(cos(ph)^2 + cos(ph)*sin(ph))*sin(th)^2 + 1)*r


Cross product
~~~~~~~~~~~~~

The cross product of :math:`u` by :math:`v` is obtained by the method
``cross_product``, which admits ``cross`` as a shortcut alias:

::

    sage: s = u.cross(v)
    sage: s
    Vector field u x v on the Euclidean space E^3
    sage: s.display()
    u x v = -2*(r*cos(ph)*cos(th)*sin(ph)*u_phi(r, th, ph)
     - r*cos(ph)^2*u_theta(r, th, ph))*sin(th) e_r
     + (2*r*cos(ph)*sin(ph)*sin(th)^2*u_phi(r, th, ph)
     - 2*r*cos(ph)^2*sin(th)*u_r(r, th, ph) + r*u_phi(r, th, ph)) e_th
     + (2*r*cos(ph)*cos(th)*sin(ph)*sin(th)*u_r(r, th, ph)
     - 2*r*cos(ph)*sin(ph)*sin(th)^2*u_theta(r, th, ph) - r*u_theta(r, th, ph)) e_ph


Scalar triple product
~~~~~~~~~~~~~~~~~~~~~

Let us introduce a third vector field. As a example, we do not pass the
components as arguments of ``vector_field``, as we did for :math:`u` and
:math:`v`; instead, we set them in a second stage, via the square
bracket operator, any unset component being assumed to be zero:

::

    sage: w = E.vector_field(name='w')
    sage: w[1] = r
    sage: w.display()
    w = r e_r

The scalar triple product of the vector fields :math:`u`, :math:`v` and
:math:`w` is obtained as follows:

::

    sage: triple_product = E.scalar_triple_product()
    sage: s = triple_product(u, v, w)
    sage: s
    Scalar field epsilon(u,v,w) on the Euclidean space E^3
    sage: s.expr()
    -2*(r^2*cos(ph)*cos(th)*sin(ph)*u_phi(r, th, ph)
     - r^2*cos(ph)^2*u_theta(r, th, ph))*sin(th)

Let us check that the scalar triple product of :math:`u`, :math:`v` and
:math:`w` is :math:`u\cdot(v\times w)`:

::

    sage: s == u.dot(v.cross(w))
    True


Differential operators
----------------------

While the standard operators :math:`\mathrm{grad}`, :math:`\mathrm{div}`,
:math:`\mathrm{curl}`, etc. involved in vector calculus are accessible via
the dot notation (e.g. ``v.div()``), let us import functions ``grad``,
``div``, ``curl``, etc. that allow for using standard mathematical notations
(e.g. ``div(v)``):

::

    sage: from sage.manifolds.operators import *


Gradient of a scalar field
~~~~~~~~~~~~~~~~~~~~~~~~~~

We first introduce a scalar field, via its expression in terms of
Cartesian coordinates; in this example, we consider a unspecified
function of :math:`(r,\theta,\phi)`:

::

    sage: F = E.scalar_field(function('f')(r,th,ph), name='F')
    sage: F.display()
    F: E^3 --> R
       (r, th, ph) |--> f(r, th, ph)

The value of :math:`F` at a point:

::

    sage: F(p)
    f(1, 1/2*pi, pi)

The gradient of :math:`F`:

::

    sage: grad(F)
    Vector field grad(F) on the Euclidean space E^3
    sage: grad(F).display()
    grad(F) = d(f)/dr e_r + d(f)/dth/r e_th + d(f)/dph/(r*sin(th)) e_ph
    sage: norm(grad(F)).display()
    |grad(F)|: E^3 --> R
       (r, th, ph) |--> sqrt((r^2*(d(f)/dr)^2 + (d(f)/dth)^2)*sin(th)^2
        + (d(f)/dph)^2)/(r*sin(th))


Divergence
~~~~~~~~~~

The divergence of a vector field:

::

    sage: s = div(u)
    sage: s.display()
    div(u): E^3 --> R
       (r, th, ph) |--> ((r*d(u_r)/dr + 2*u_r(r, th, ph)
        + d(u_theta)/dth)*sin(th) + cos(th)*u_theta(r, th, ph)
        + d(u_phi)/dph)/(r*sin(th))
    sage: s.expr().expand()
    2*u_r(r, th, ph)/r + cos(th)*u_theta(r, th, ph)/(r*sin(th))
     + diff(u_theta(r, th, ph), th)/r + diff(u_phi(r, th, ph), ph)/(r*sin(th))
     + diff(u_r(r, th, ph), r)

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

::

    sage: s.display()
    curl(u) = (cos(th)*u_phi(r, th, ph) + sin(th)*d(u_phi)/dth
     - d(u_theta)/dph)/(r*sin(th)) e_r - ((r*d(u_phi)/dr + u_phi(r, th, ph))*sin(th)
     - d(u_r)/dph)/(r*sin(th)) e_th + (r*d(u_theta)/dr + u_theta(r, th, ph)
     - d(u_r)/dth)/r e_ph


To use the notation ``rot`` instead of ``curl``, simply do

::

    sage: rot = curl

An alternative is

::

    sage: from sage.manifolds.operators import curl as rot

We have then

::

    sage: rot(u).display()
    curl(u) = (cos(th)*u_phi(r, th, ph) + sin(th)*d(u_phi)/dth
     - d(u_theta)/dph)/(r*sin(th)) e_r - ((r*d(u_phi)/dr + u_phi(r, th, ph))*sin(th)
     - d(u_r)/dph)/(r*sin(th)) e_th + (r*d(u_theta)/dr + u_theta(r, th, ph)
     - d(u_r)/dth)/r e_ph
    sage: rot(u) == curl(u)
    True

For :math:`v` and :math:`w`, we have

::

    sage: curl(v).display()
    curl(v) = 2*cos(th) e_r - 2*sin(th) e_th
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
       (r, th, ph) |--> 0

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
       (r, th, ph) |--> ((r^2*d^2(f)/dr^2 + 2*r*d(f)/dr
        + d^2(f)/dth^2)*sin(th)^2 + cos(th)*sin(th)*d(f)/dth
        + d^2(f)/dph^2)/(r^2*sin(th)^2)
    sage: s.expr().expand()
    2*diff(f(r, th, ph), r)/r + cos(th)*diff(f(r, th, ph), th)/(r^2*sin(th))
     + diff(f(r, th, ph), th, th)/r^2 + diff(f(r, th, ph), ph, ph)/(r^2*sin(th)^2)
     + diff(f(r, th, ph), r, r)

For a scalar field, the Laplacian is nothing but the divergence of the
gradient:

::

    sage: laplacian(F) == div(grad(F))
    True

The Laplacian of a vector field:

::

    sage: Du = laplacian(u)
    sage: Du.display()
    Delta(u) = ((r^2*d^2(u_r)/dr^2 + 2*r*d(u_r)/dr - 2*u_r(r, th, ph)
     + d^2(u_r)/dth^2 - 2*d(u_theta)/dth)*sin(th)^2 - ((2*u_theta(r, th, ph)
     - d(u_r)/dth)*cos(th) + 2*d(u_phi)/dph)*sin(th) + d^2(u_r)/dph^2)/(r^2*sin(th)^2) e_r
     + ((r^2*d^2(u_theta)/dr^2 + 2*r*d(u_theta)/dr + 2*d(u_r)/dth + d^2(u_theta)/dth^2)*sin(th)^2
     + cos(th)*sin(th)*d(u_theta)/dth - 2*cos(th)*d(u_phi)/dph - u_theta(r, th, ph)
     + d^2(u_theta)/dph^2)/(r^2*sin(th)^2) e_th
     + ((r^2*d^2(u_phi)/dr^2 + 2*r*d(u_phi)/dr
     + d^2(u_phi)/dth^2)*sin(th)^2 + (cos(th)*d(u_phi)/dth + 2*d(u_r)/dph)*sin(th)
     + 2*cos(th)*d(u_theta)/dph - u_phi(r, th, ph) + d^2(u_phi)/dph^2)/(r^2*sin(th)^2) e_ph

Since this expression is quite lengthy, we may ask for a display
component by component:

::

    sage: Du.display_comp()
    Delta(u)^1 = ((r^2*d^2(u_r)/dr^2 + 2*r*d(u_r)/dr - 2*u_r(r, th, ph) + d^2(u_r)/dth^2
     - 2*d(u_theta)/dth)*sin(th)^2 - ((2*u_theta(r, th, ph) - d(u_r)/dth)*cos(th)
     + 2*d(u_phi)/dph)*sin(th) + d^2(u_r)/dph^2)/(r^2*sin(th)^2)
    Delta(u)^2 = ((r^2*d^2(u_theta)/dr^2 + 2*r*d(u_theta)/dr + 2*d(u_r)/dth
     + d^2(u_theta)/dth^2)*sin(th)^2 + cos(th)*sin(th)*d(u_theta)/dth
     - 2*cos(th)*d(u_phi)/dph - u_theta(r, th, ph) + d^2(u_theta)/dph^2)/(r^2*sin(th)^2)
    Delta(u)^3 = ((r^2*d^2(u_phi)/dr^2 + 2*r*d(u_phi)/dr + d^2(u_phi)/dth^2)*sin(th)^2
     + (cos(th)*d(u_phi)/dth + 2*d(u_r)/dph)*sin(th) + 2*cos(th)*d(u_theta)/dph
     - u_phi(r, th, ph) + d^2(u_phi)/dph^2)/(r^2*sin(th)^2)

We may expand each component:

::

    sage: for i in E.irange():
    ....:     s = Du[i].expand()
    sage: Du.display_comp()
    Delta(u)^1 = 2*d(u_r)/dr/r - 2*u_r(r, th, ph)/r^2
     - 2*cos(th)*u_theta(r, th, ph)/(r^2*sin(th)) + cos(th)*d(u_r)/dth/(r^2*sin(th))
     + d^2(u_r)/dth^2/r^2 - 2*d(u_theta)/dth/r^2 - 2*d(u_phi)/dph/(r^2*sin(th))
     + d^2(u_r)/dph^2/(r^2*sin(th)^2) + d^2(u_r)/dr^2
    Delta(u)^2 = 2*d(u_theta)/dr/r + 2*d(u_r)/dth/r^2 + cos(th)*d(u_theta)/dth/(r^2*sin(th))
     + d^2(u_theta)/dth^2/r^2 - 2*cos(th)*d(u_phi)/dph/(r^2*sin(th)^2)
     - u_theta(r, th, ph)/(r^2*sin(th)^2) + d^2(u_theta)/dph^2/(r^2*sin(th)^2)
     + d^2(u_theta)/dr^2
    Delta(u)^3 = 2*d(u_phi)/dr/r + cos(th)*d(u_phi)/dth/(r^2*sin(th))
     + d^2(u_phi)/dth^2/r^2 + 2*d(u_r)/dph/(r^2*sin(th))
     + 2*cos(th)*d(u_theta)/dph/(r^2*sin(th)^2) - u_phi(r, th, ph)/(r^2*sin(th)^2)
     + d^2(u_phi)/dph^2/(r^2*sin(th)^2) + d^2(u_phi)/dr^2

::

    sage: Du[1]
    2*d(u_r)/dr/r - 2*u_r(r, th, ph)/r^2 - 2*cos(th)*u_theta(r, th, ph)/(r^2*sin(th))
     + cos(th)*d(u_r)/dth/(r^2*sin(th)) + d^2(u_r)/dth^2/r^2 - 2*d(u_theta)/dth/r^2
     - 2*d(u_phi)/dph/(r^2*sin(th)) + d^2(u_r)/dph^2/(r^2*sin(th)^2) + d^2(u_r)/dr^2
    sage: Du[2]
    2*d(u_theta)/dr/r + 2*d(u_r)/dth/r^2 + cos(th)*d(u_theta)/dth/(r^2*sin(th))
     + d^2(u_theta)/dth^2/r^2 - 2*cos(th)*d(u_phi)/dph/(r^2*sin(th)^2)
     - u_theta(r, th, ph)/(r^2*sin(th)^2) + d^2(u_theta)/dph^2/(r^2*sin(th)^2)
     + d^2(u_theta)/dr^2
    sage: Du[3]
    2*d(u_phi)/dr/r + cos(th)*d(u_phi)/dth/(r^2*sin(th)) + d^2(u_phi)/dth^2/r^2
     + 2*d(u_r)/dph/(r^2*sin(th)) + 2*cos(th)*d(u_theta)/dph/(r^2*sin(th)^2)
     - u_phi(r, th, ph)/(r^2*sin(th)^2) + d^2(u_phi)/dph^2/(r^2*sin(th)^2) + d^2(u_phi)/dr^2

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
    curl(curl(u)) = ((r*d^2(u_theta)/drdth - d^2(u_r)/dth^2 + d(u_theta)/dth)*sin(th)^2
     + ((r*d(u_theta)/dr + u_theta(r, th, ph) - d(u_r)/dth)*cos(th) + r*d^2(u_phi)/drdph
     + d(u_phi)/dph)*sin(th) - d^2(u_r)/dph^2)/(r^2*sin(th)^2) e_r - ((r^2*d^2(u_theta)/dr^2
    - r*d^2(u_r)/drdth + 2*r*d(u_theta)/dr)*sin(th)^2 - sin(th)*d^2(u_phi)/dthdph
    - cos(th)*d(u_phi)/dph + d^2(u_theta)/dph^2)/(r^2*sin(th)^2) e_th - ((r^2*d^2(u_phi)/dr^2
    + 2*r*d(u_phi)/dr + d^2(u_phi)/dth^2)*sin(th)^2 + (cos(th)*d(u_phi)/dth - r*d^2(u_r)/drdph
    - d^2(u_theta)/dthdph)*sin(th) + cos(th)*d(u_theta)/dph - u_phi(r, th, ph))/(r^2*sin(th)^2) e_ph
    sage: grad(div(u)).display()
    grad(div(u)) = ((r*d(u_theta)/dr - u_theta(r, th, ph))*cos(th)
     + (r^2*d^2(u_r)/dr^2 + 2*r*d(u_r)/dr + r*d^2(u_theta)/drdth - 2*u_r(r, th, ph)
     - d(u_theta)/dth)*sin(th) + r*d^2(u_phi)/drdph - d(u_phi)/dph)/(r^2*sin(th)) e_r
     + ((r*d^2(u_r)/drdth + 2*d(u_r)/dth + d^2(u_theta)/dth^2)*sin(th)^2 + (cos(th)*d(u_theta)/dth
     + d^2(u_phi)/dthdph)*sin(th) - cos(th)*d(u_phi)/dph - u_theta(r, th, ph))/(r^2*sin(th)^2) e_th
     + ((r*d^2(u_r)/drdph + 2*d(u_r)/dph + d^2(u_theta)/dthdph)*sin(th) + cos(th)*d(u_theta)/dph
     + d^2(u_phi)/dph^2)/(r^2*sin(th)^2) e_ph

and we may check a famous identity:

::

    sage: curl(curl(u)) == grad(div(u)) - laplacian(u)
    True


Customizations
--------------

Customizing the symbols of the orthonormal frame vectors
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

By default, the vectors of the orthonormal frame associated with
spherical coordinates are denoted :math:`(e_r,e_\theta,e_\phi)`:

::

    sage: frame = E.spherical_frame()
    sage: frame
    Vector frame (E^3, (e_r,e_th,e_ph))

But this can be changed, thanks to the method ``set_name``:

::

    sage: frame.set_name('a', indices=('r', 'th', 'ph'),
    ....:                latex_indices=('r', r'\theta', r'\phi'))
    sage: frame
    Vector frame (E^3, (a_r,a_th,a_ph))
    sage: v.display()
    v = (r*sin(2*ph)*sin(th)^2 + r) a_r + r*cos(th)*sin(2*ph)*sin(th) a_th
     + 2*r*cos(ph)^2*sin(th) a_ph

::

    sage: frame.set_name(('hr', 'hth', 'hph'),
    ....:                latex_symbol=(r'\hat{r}', r'\hat{\theta}', r'\hat{\phi}'))
    sage: frame
    Vector frame (E^3, (hr,hth,hph))
    sage: v.display()
    v = (r*sin(2*ph)*sin(th)^2 + r) hr + r*cos(th)*sin(2*ph)*sin(th) hth
     + 2*r*cos(ph)^2*sin(th) hph


Customizing the coordinate symbols
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The coordinates symbols are defined within the angle brackets ``<...>``
at the construction of the Euclidean space. Above we did

::

    sage: E.<r,th,ph> = EuclideanSpace(coordinates='spherical')

which resulted in the coordinate symbols :math:`(r,\theta,\phi)` and in
the corresponding Python variables ``r``, ``th`` and ``ph`` (SageMath
symbolic expressions). Using other symbols, for instance
:math:`(R,\Theta,\Phi)`, is possible through the optional argument
``symbols`` of the function ``EuclideanSpace``. It has to be a string,
usually prefixed by ``r`` (for raw string, in order to allow for the
backslash character of LaTeX expressions). This string contains the
coordinate fields separated by a blank space; each field contains the
coordinate’s text symbol and possibly the coordinate’s LaTeX symbol
(when the latter is different from the text symbol), both symbols being
separated by a colon (``:``):

::

    sage: E.<R,Th,Ph> = EuclideanSpace(coordinates='spherical', symbols=r'R Th:\Theta Ph:\Phi')

We have then

::

    sage: E.atlas()
    [Chart (E^3, (R, Th, Ph))]
    sage: E.frames()
    [Coordinate frame (E^3, (d/dR,d/dTh,d/dPh)),
     Vector frame (E^3, (e_R,e_Th,e_Ph))]
    sage: E.spherical_frame()
    Vector frame (E^3, (e_R,e_Th,e_Ph))
    sage: v = E.vector_field(R*sin(2*Ph)*sin(Th)^2 + R,
    ....:                    R*sin(2*Ph)*sin(Th)*cos(Th),
    ....:                    2*R*cos(Ph)^2*sin(Th), name='v')
    sage: v.display()
    v = (R*sin(2*Ph)*sin(Th)^2 + R) e_R + R*cos(Th)*sin(2*Ph)*sin(Th) e_Th
     + 2*R*cos(Ph)^2*sin(Th) e_Ph
