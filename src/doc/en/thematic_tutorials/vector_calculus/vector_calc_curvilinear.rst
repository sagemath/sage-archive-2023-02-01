.. -*- coding: utf-8 -*-

.. linkall

How to use curvilinear coordinates for vector calculus?
=======================================================

This tutorial introduces some vector calculus functionalities of SageMath
within the 3-dimensional Euclidean space.
The corresponding tools have been developed via the
`SageManifolds <https://sagemanifolds.obspm.fr>`__ project.

The tutorial is also available as a Jupyter notebook, either
`passive <https://nbviewer.jupyter.org/github/sagemanifolds/SageManifolds/blob/master/Notebooks/SM_vector_calc_spherical.ipynb>`__ (``nbviewer``)
or `interactive <https://mybinder.org/v2/gh/sagemanifolds/SageManifolds/master?filepath=Notebooks/SM_vector_calc_spherical.ipynb>`__ (``binder``).

Using spherical coordinates within the Euclidean 3-space
--------------------------------------------------------

Let us introduce the 3-dimensional Euclidean space
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


