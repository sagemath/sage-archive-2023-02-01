r"""
Bundle Connections

Let `E \to M` be a smooth vector bundle of rank `n` over a smooth manifold `M`
and over a non-discrete topological field `K` (typically `K=\RR` or `K=\CC`). A
*bundle connection* on this vector bundle is a `K`-linear map

.. MATH::

    \nabla : C^\infty(M;E) \to C^\infty(M;E \otimes T^*M)

such that the Leibniz rule applies for each scalar field `f \in C^\infty(M)` and
section `s \in C^\infty(M;E)`:

.. MATH::

    \nabla(f \, s) = f \cdot \nabla s + s \otimes \mathrm{d}f .

If `e` is a local frame on `E`, we have

.. MATH::

    \nabla e_i = \sum^n_{j=1} e_j \otimes \omega^j_i ,

and the corresponding `n \times n`-matrix `(\omega^j_i)_{i,j}` consisting of
one forms is called *connection matrix of* `\nabla` *with respect to* `e`.

AUTHORS:

- Michael Jung (2019) : initial version

"""
# ******************************************************************************
#       Copyright (C) 2019 Michael Jung <micjung@uni-potsdam.de>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  https://www.gnu.org/licenses/
# ******************************************************************************

from sage.structure.sage_object import SageObject
from sage.rings.integer import Integer
from sage.manifolds.differentiable.vector_bundle import \
    DifferentiableVectorBundle


class BundleConnection(SageObject):
    r"""
    An instance of this class represents a bundle connection `\nabla` on a
    smooth vector bundle `E \to M`.

    INPUT:

    - ``vbundle`` -- the vector bundle on which the connection is defined
      (must be an instance of class
      :class:`~sage.manifolds.differentiable.vector_bundle.DifferentiableVectorBundle`)
    - ``name`` -- name given to the bundle connection
    - ``latex_name`` -- (default: ``None``) LaTeX symbol to denote the bundle
      connection; if ``None``, it is set to ``name``

    EXAMPLES:

    Define a bundle connection on a rank 2 vector bundle over some
    3-dimensional smooth manifold::

        sage: M = Manifold(3, 'M', start_index=1)
        sage: X.<x,y,z> = M.chart()
        sage: E = M.vector_bundle(2, 'E')
        sage: e = E.local_frame('e') # standard frame for E
        sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla'); nab
        Bundle connection nabla on the Differentiable real vector bundle E -> M
         of rank 2 over the base space 3-dimensional differentiable manifold M

    First, let us initialize all connection 1-forms w.r.t. the frame ``e`` to
    zero::

        sage: nab[e, :] = [[0, 0], [0, 0]]

    This line can be shortened by the following::

        sage: nab[e, :] = 0  # initialize to zero

    Now, we want to specify some non-zero entries::

        sage: nab[e, 1, 2] = [x*z, y*z, z^2]
        sage: nab[e, 2, 1] = [x, x^2, x^3]
        sage: nab[e, 1, 2].display()
        connection (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x*z dx + y*z dy + z^2 dz
        sage: nab[e, 2, 1].display()
        connection (2,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x dx + x^2 dy + x^3 dz

    The other entries remain zero::

        sage: nab[e, 1, 1].display()
        connection (1,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 0
        sage: nab[e, 2, 2].display()
        connection (2,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 0

    Notice, when we omit the frame, the default frame of the vector bundle is
    assumed (in this case ``e``)::

        sage: nab[2, 2].display()
        connection (2,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 0

    The same holds for the assignment::

        sage: nab[:] = 0
        sage: nab[e, 1, 2].display()
        connection (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 0

    We can also use :meth:`set_connection_form` to specify the connection
    1-forms::

        sage: nab[:] = 0  # re-initialize to zero
        sage: nab.set_connection_form(1, 2)[:] = [x*z, y*z, z^2]
        sage: nab.set_connection_form(2, 1)[:] = [x, x^2, x^3]
        sage: nab[1, 2].display()
        connection (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x*z dx + y*z dy + z^2 dz
        sage: nab[2, 1].display()
        connection (2,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = x dx + x^2 dy + x^3 dz

    .. NOTE::

        Notice that list assignments and :meth:`set_connection_form` delete
        the connection 1-forms w.r.t. other frames for consistency reasons. To
        avoid this behavior, :meth:`add_connection_form` must be used instead.

    After the connection has been specified, the curvature 2-forms can be
    derived::

        sage: Omega = nab.curvature_form
        sage: for i in E.irange():
        ....:     for j in E.irange():
        ....:         print(Omega(i ,j, e).display())
        curvature (1,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = -(x^3 - x*y)*z dx/\dy + (-x^4*z + x*z^2) dx/\dz +
         (-x^3*y*z + x^2*z^2) dy/\dz
         curvature (1,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = -x dx/\dz - y dy/\dz
         curvature (2,1) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = 2*x dx/\dy + 3*x^2 dx/\dz
         curvature (2,2) of bundle connection nabla w.r.t. Local frame
         (E|_M, (e_1,e_2)) = (x^3 - x*y)*z dx/\dy + (x^4*z - x*z^2) dx/\dz +
         (x^3*y*z - x^2*z^2) dy/\dz

    The derived forms certainly obey the structure equations, see
    :meth:`curvature_form` for details::

        sage: omega = nab.connection_form
        sage: check = []
        sage: for i in E.irange():  # long time
        ....:     for j in E.irange():
        ....:         check.append(Omega(i,j,e) == \
        ....:                       omega(i,j,e).exterior_derivative() + \
        ....:         sum(omega(k,j,e).wedge(omega(i,k,e))
        ....:             for k in E.irange()))
        sage: check  # long time
        [True, True, True, True]

    """

    def __init__(self, vbundle, name, latex_name=None):
        r"""
        Construct a bundle connection.

        TESTS::

            sage: M = Manifold(3, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: from sage.manifolds.differentiable.bundle_connection \
            ....:                                       import BundleConnection
            sage: nab = BundleConnection(E, 'nabla', latex_name=r'\nabla')
            sage: nab
            Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 2 over the base space 3-dimensional differentiable
             manifold M
            sage: X.<x,y,z> = M.chart()
            sage: e = E.local_frame('e')
            sage: nab[:] = 0
            sage: nab.set_connection_form(1, 0)[:] = [x*z, y*z, z^2]
            sage: TestSuite(nab).run()

        """
        if not isinstance(vbundle, DifferentiableVectorBundle):
            raise TypeError("the first argument must be a differentiable " +
                            "vector bundle")
        self._vbundle = vbundle
        self._base_space = vbundle.base_space()
        self._name = name
        if latex_name is None:
            self._latex_name = self._name
        else:
            self._latex_name = latex_name
        self._connection_forms = {}  # dict. of con. forms, with frames as keys
        self._coefficients = self._connection_forms
        # Initialization of derived quantities:
        self._init_derived()

    def _repr_(self):
        r"""
        String representation of the object.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: E = M.vector_bundle(3, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._repr_()
            'Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 3 over the base space 5-dimensional differentiable
             manifold M'
            sage: repr(nab)  # indirect doctest
            'Bundle connection nabla on the Differentiable real vector bundle
             E -> M of rank 3 over the base space 5-dimensional differentiable
             manifold M'

        """
        description = "Bundle connection"
        if self._name is not None:
            description += " " + self._name
        description += " on the {}".format(self._vbundle)
        return description

    def _latex_(self):
        r"""
        LaTeX representation of the object.

        TESTS::

            sage: M = Manifold(5, 'M')
            sage: E = M.vector_bundle(3, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._latex_()
            '\\nabla'
            sage: latex(nab)  # indirect doctest
            \nabla
            sage: nab = E.bundle_connection('D')
            sage: nab._latex_()
            'D'
            sage: latex(nab)  # indirect doctest
            D

        """
        return self._latex_name

    def _init_derived(self):
        r"""
        Initialize the derived quantities.

        TESTS::

            sage: M = Manifold(4, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._init_derived()

        """
        self._curvature_forms = {}  # dict. of dict. of curvature forms
        # (key: local frame)
        self._hash = -1

    def _del_derived(self):
        r"""
        Delete the derived quantities.

        TESTS::

            sage: M = Manifold(4, 'M')
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._del_derived()

        """
        self._curvature_forms.clear()

    def __eq__(self, other):
        r"""
        Comparison (equality) operator.

        INPUT:

        - ``other`` -- a bundle connection

        OUTPUT:

        - ``True`` if ``self`` is equal to ``other`` and ``False`` otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')  # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab[:] = 0
            sage: nab[0, 1][:] = [x^2, x]
            sage: nab[1, 0][:] = [y^2, y]
            sage: nab1 = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab1[:] = 0
            sage: (nab1 == nab) or (nab == nab1)
            False
            sage: nab1[0, 1][:] = [x, x^2]
            sage: nab1[1, 0][:] = [y, y^2]
            sage: (nab1 == nab) or (nab == nab1)
            False
            sage: nab1[0, 1][:] = [x^2, x]
            sage: nab1[1, 0][:] = [y^2, y]
            sage: (nab1 == nab) and (nab == nab1)
            True

        """
        if other is self:
            return True
        if not isinstance(other, BundleConnection):
            return False
        if other._base_space != self._base_space:
            return False
        if self._connection_forms == {}:
            return False
        for frame in self._connection_forms:
            if frame not in other._connection_forms:
                return False
            for ind in self._connection_forms[frame]:
                if (other._connection_forms[frame][ind] !=
                        self._connection_forms[frame][ind]):
                    return False
        return True

    def __ne__(self, other):
        r"""
        Inequality operator.

        INPUT:

        - ``other`` -- an affine connection

        OUTPUT:

        - ``True`` if ``self`` is different from ``other`` and ``False``
          otherwise

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')  # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab[:] = 0
            sage: nab[0, 1][:] = [x^2, x]
            sage: nab[1, 0][:] = [y^2, y]
            sage: nab1 = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab1[:] = 0
            sage: (nab1 != nab) and (nab != nab1)
            True
            sage: nab1[0, 1][:] = [x, x^2]
            sage: nab1[1, 0][:] = [y, y^2]
            sage: (nab1 != nab) and (nab != nab1)
            True
            sage: nab1[0, 1][:] = [x^2, x]
            sage: nab1[1, 0][:] = [y^2, y]
            sage: (nab1 != nab) or (nab != nab1)
            False

        """
        return not (self == other)

    def vector_bundle(self):
        r"""
        Return the vector bundle on which the bundle connection is defined.

        OUTPUT:

        - instance of class
          :class:`~sage.manifolds.differentiable.vector_bundle.DifferentiableVectorBundle`
          representing the vector bundle on which ``self`` is defined.

        EXAMPLES::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', r'\nabla')
            sage: nab.vector_bundle()
            Differentiable real vector bundle E -> M of rank 2 over the base
             space 3-dimensional differentiable manifold M

        """
        return self._vbundle

    def _new_forms(self, frame):
        r"""
        Create the connection forms w.r.t. the given frame.

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab._new_forms(e)  # random
            dict of forms

        """
        dom = frame._domain
        forms_dict = {}
        for i in self._vbundle.irange():
            for j in self._vbundle.irange():
                # set names:
                name = "connection ({},{}) of bundle ".format(i, j)
                name += "connection " + self._name + " w.r.t. {}".format(frame)
                latex_name = r"\omega^" + str(j) + r"_{\ \, " + str(i) + "}"
                form = dom.diff_form(1, name=name, latex_name=latex_name)
                forms_dict[(i, j)] = form
        return forms_dict

    def connection_forms(self, frame=None):
        r"""
        Return the connection forms relative to the given frame.

        If `e` is a local frame on `E`, we have

        .. MATH::

            \nabla e_i = \sum^n_{j=1} e_j \otimes \omega^j_i ,

        and the corresponding `n \times n`-matrix `(\omega^j_i)_{i,j}`
        consisting of one forms is called *connection matrix of* `\nabla` *with
        respect to* `e`.

        INPUT:

        - ``frame`` -- (default: ``None``) local frame relative to which the
          connection forms are required; if none is provided, the
          vector bundle's default frame is assumed

        OUTPUT:

        - connection forms relative to the frame ``frame``, as a dictionary
          with tuples `(i, j)` as key and one forms as instances of
          :class:`~sage.manifolds.differentiable.diff_form` as value
          representing the matrix entries.

        EXAMPLES:

        Connection forms of a bundle connection on a rank 2 vector bundle
        over a 3-dimensional manifold::

            sage: M = Manifold(3, 'M', start_index=1)
            sage: c_xyz.<x,y,z> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e')
            sage: nab = E.bundle_connection('nabla', r'\nabla')
            sage: nab[:] = 0  # initialize curvature forms
            sage: nab.connection_forms() # random
            {(1, 1): 1-form a on the 3-dimensional differentiable manifold M,
             (1, 2): 1-form zero on the 3-dimensional differentiable
             manifold M,
             (2, 1): 1-form zero on the 3-dimensional differentiable
             manifold M,
             (2, 2): 1-form b on the 3-dimensional differentiable manifold M}

        """
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._base_space)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided!")
        if frame not in self._connection_forms:
            # the connection forms must be computed
            #
            # Check whether frame is a subframe of a frame in which the
            # forms are already known:
            for oframe in self._connection_forms:
                if frame in oframe._subframes:
                    self._connection_forms[frame] = self._new_forms(frame)
                    comp_store = self._connection_forms[frame]
                    ocomp_store = self._connection_forms[oframe]
                    for ind, value in ocomp_store.items():
                        comp_store[ind] = value.restrict(frame._domain)
                    break
            else:
                # TODO: Compute coefficients out of known ones
                pass
        return self._connection_forms[frame]

    def connection_form(self, i, j, frame=None):
        r"""
        Return the connection 1-form corresponding to the given index and
        local frame.

        INPUT:

        - ``i``, ``j`` -- indices identifying the 1-form `\omega^j_i`
        - ``frame`` -- (default: ``None``) local frame relative to which the
          connection 1-forms are defined; if ``None``, the default frame of the
          vector bundle's corresponding section module is assumed.

        OUTPUT:

        - the 1-form `\omega^j_i`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab.set_connection_form(0, 1)[:] = [x^2, x]
            sage: nab.set_connection_form(1, 0)[:] = [y^2, y]
            sage: nab.connection_form(0, 1).display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy
            sage: nab.connection_form(1, 0).display()
            connection (1,0) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = y^2 dx + y dy

        """
        return self.connection_forms(frame)[(i, j)]

    def add_connection_form(self, i, j, form=None, frame=None):
        r"""
        Return the connection form `\omega^j_i` in a given frame for
        assignment.

        See method :meth:`connection_forms` for details about the definition of
        the connection forms.

        To delete the connection forms in other frames, use the method
        :meth:`set_connection_form` instead.

        INPUT:

        - ``i``, ``j`` -- indices identifying the 1-form `\omega^j_i`
        - ``frame`` -- (default: ``None``) local frame in which the connection
          1-form is defined; if ``None``, the default frame of the vector
          bundle is assumed.

        .. WARNING::

            If the connection has already forms in other frames, it is the
            user's responsibility to make sure that the 1-forms to be added
            are consistent with them.

        OUTPUT:

        - connection 1-form `\omega^j_i` in the given frame, as an instance of
          the class :class:`~sage.manifolds.differentiable.diff_form.DiffForm`;
          if such connection 1-form did not exist previously, it is created.
          See method :meth:`connection_forms` for the storage convention of the
          connection 1-forms.

        EXAMPLES::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab.add_connection_form(0, 1, frame=e)[:] = [x^2, x]
            sage: nab[e, 0, 1].display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy

        Since ``e`` is the vector bundle's default local frame, its mention may
        be omitted::

            sage: nab.add_connection_form(1, 0)[:] = [y^2, y]
            sage: nab[1, 0].display()
            connection (1,0) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = y^2 dx + y dy

        Adding connection 1-forms w.r.t. to another local frame::

            sage: f = E.local_frame('f')
            sage: nab.add_connection_form(1, 1, frame=f)[:] = [x, y]
            sage: nab[f, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (f_0,f_1)) = x dx + y dy

        The forms w.r.t. the frame ``e`` have been kept::

            sage: nab[e, 0, 1].display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy

        To delete them, use the method :meth:`set_connection_form` instead.

        """
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._base_space)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided!")
        # Are the components already known?
        if frame not in self._connection_forms:
            if frame not in self._vbundle._frames:
                raise ValueError("the {} is not".format(frame) +
                                 " a frame on the {}".format(self._base_space))
            self._connection_forms[frame] = self._new_forms(frame)
        self._del_derived()  # deletes the derived quantities
        if form:
            # TODO: Remove input `form` in Sage 9.3
            from sage.misc.superseded import deprecation
            msg = "the input 'form' is outdated and will be removed in a "
            msg += "future version of Sage"
            deprecation(30208, msg)
            self._connection_forms[frame][(i, j)].copy_from(form)
        return self._connection_forms[frame][(i, j)]

    def set_connection_form(self, i, j, form=None, frame=None):
        r"""
        Return the connection form `\omega^j_i` in a given frame for
        assignment.

        See method :meth:`connection_forms` for details about the definition of
        the connection forms.

        The connection forms with respect to other frames are deleted,
        in order to avoid any inconsistency. To keep them, use the method
        :meth:`add_connection_form` instead.

        INPUT:

        - ``i``, ``j`` -- indices identifying the 1-form `\omega^j_i`
        - ``frame`` -- (default: ``None``) local frame in which the connection
          1-form is defined; if ``None``, the default frame of the vector
          bundle is assumed.

        OUTPUT:

        - connection 1-form `\omega^j_i` in the given frame, as an instance of
          the class :class:`~sage.manifolds.differentiable.diff_form.DiffForm`;
          if such connection 1-form did not exist previously, it is created.
          See method :meth:`connection_forms` for the storage convention of the
          connection 1-forms.

        EXAMPLES:

        Setting the connection forms of a bundle connection w.r.t. some local
        frame::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: e = E.local_frame('e') # standard frame for E
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: nab.set_connection_form(0, 1)[:] = [x^2, x]
            sage: nab[0, 1].display()
            connection (0,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = x^2 dx + x dy

        Since ``e`` is the vector bundle's default local frame, its mention may
        be omitted::

            sage: nab.set_connection_form(1, 0)[:] = [y^2, y]
            sage: nab[1, 0].display()
            connection (1,0) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_0,e_1)) = y^2 dx + y dy

        Setting connection 1-forms w.r.t. to another local frame::

            sage: f = E.local_frame('f')
            sage: nab.set_connection_form(1, 1, frame=f)[:] = [x, y]
            sage: nab[f, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (f_0,f_1)) = x dx + y dy

        The forms w.r.t. the frame ``e`` have been deleted::

            sage: nab[e, 0, 1].display()
            Traceback (most recent call last):
            ...
            KeyError: Local frame (E|_M, (e_0,e_1))

        To keep them, use the method :meth:`add_connection_form` instead.

        """
        if form:
            # TODO: Remove input `form` in Sage 9.3
            from sage.misc.superseded import deprecation
            msg = "the input 'form' is outdated and will be removed in a "
            msg += "future version of Sage"
            deprecation(30208, msg)
        omega = self.add_connection_form(i, j, form=form, frame=frame)
        self.del_other_forms(frame)
        return omega

    def del_other_forms(self, frame=None):
        r"""
        Delete all the connection forms but those corresponding to ``frame``.

        INPUT:

        - ``frame`` -- (default: ``None``) local frame, the connection forms
          w.r.t. which are to be kept; if ``None``, the default frame of the
          vector bundle is assumed.

        EXAMPLES:

        We first create two sets of connection forms::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: e = E.local_frame('e')
            sage: nab.set_connection_form(1, 1, frame=e)[:] = [x^2, x]
            sage: f = E.local_frame('f')
            sage: nab.add_connection_form(1, 1, frame=f)[:] = [y^2, y]
            sage: nab[e, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[f, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (f_1,f_2)) = y^2 dx + y dy

        Let us delete the connection forms w.r.t. all frames except for
        frame ``e``::

            sage: nab.del_other_forms(e)
            sage: nab[e, 1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy

        The connection forms w.r.t. frame ``e`` have indeed been
        deleted::

            sage: nab[f, :]
            Traceback (most recent call last):
            ...
            KeyError: Local frame (E|_M, (f_1,f_2))

        """
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._base_space)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided!")
        if frame not in self._connection_forms:
            raise ValueError("the coefficients w.r.t. {}".format(frame) +
                             " have not been defined")
        to_be_deleted = []
        for other_frame in self._connection_forms:
            if other_frame != frame:
                to_be_deleted.append(other_frame)
        for other_frame in to_be_deleted:
            del self._connection_forms[other_frame]

    def curvature_form(self, i, j, frame=None):
        r"""
        Return the curvature 2-form corresponding to the given index and local
        frame.

        The *curvature 2-forms* with respect to the frame `e` are the 2-forms
        `\Omega^j_i` given by the formula

        .. MATH::

            \Omega^j_i = \mathrm{d} \omega^j_i + \sum^n_{k=1} \omega^j_k
                \wedge \omega^k_i

        INPUT:

        - ``i``, ``j`` -- indices identifying the 2-form `\Omega^j_i`
        - ``frame`` -- (default: ``None``) local frame relative to which the
          curvature 2-forms are defined; if ``None``, the default frame
          of the vector bundle is assumed.

        OUTPUT:

        - the 2-form `\Omega^j_i`, as an instance of
          :class:`~sage.manifolds.differentiable.diff_form.DiffForm`

        EXAMPLES::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(1, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: e = E.local_frame('e')
            sage: nab.set_connection_form(1, 1)[:] = [x^2, x]
            sage: curv = nab.curvature_form(1, 1); curv
            2-form curvature (1,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1)) on the 2-dimensional differentiable manifold M
            sage: curv.display()
            curvature (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1)) = dx/\dy

        """
        if frame is None:
            smodule = self._vbundle.section_module(domain=self._base_space)
            frame = smodule.default_frame()
            if frame is None:
                raise ValueError("a frame must be provided!")
        if frame not in self._curvature_forms:
            self._curvature_forms[frame] = {}
        if (i, j) not in self._curvature_forms[frame]:
            name = "curvature ({},{}) of bundle connection ".format(i, j) + \
                   self._name + " w.r.t. {}".format(frame)
            latex_name = r"\Omega^" + str(i) + r"_{\ \, " + \
                         str(j) + "}"
            omega = self.connection_form
            curv_form = omega(i, j, frame).exterior_derivative()
            curv_form += sum(omega(k, j, frame).wedge(omega(i, k, frame))
                             for k in self._vbundle.irange())
            curv_form.set_name(name=name, latex_name=latex_name)
            self._curvature_forms[frame][(i, j)] = curv_form
        return self._curvature_forms[frame][(i, j)]

    def __hash__(self):
        r"""
        Hash function.

        TESTS::

            sage: M = Manifold(2, 'M')
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla', latex_name=r'\nabla')
            sage: hash(nab) == nab.__hash__()
            True

        Let us check that ``nab`` can be used as a dictionary key::

            sage: {nab: 1}[nab]
            1

        """
        if self._hash == -1:
            self._hash = hash(repr(self))
        return self._hash

    def __getitem__(self, args):
        r"""
        Return a component of ``self`` with respect to some frame.

        INPUT:

        - ``args`` -- list of indices defining the component; if ``[:]`` is
          provided, all the components are returned

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla')
            sage: e = E.local_frame('e')
            sage: nab[:] = 0
            sage: nab[1, 2][:] = [x^2, x]
            sage: nab[1, 1][:] = [y, y]
            sage: nab[1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = y dx + y dy
            sage: nab[e, 1, 2].display()
            connection (1,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[e, :]
            [[1-form connection (1,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (1,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M],
             [1-form connection (2,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (2,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M]]
            sage: nab[:]
            [[1-form connection (1,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (1,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M],
             [1-form connection (2,1) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M,
             1-form connection (2,2) of bundle connection nabla w.r.t. Local
             frame (E|_M, (e_1,e_2)) on the 2-dimensional differentiable
             manifold M]]

        """
        # extract frame from first index:
        vb = self._vbundle
        if isinstance(args, (int, Integer, slice)):
            smodule = vb.section_module(domain=self._base_space)
            frame = smodule.default_frame()
        elif not isinstance(args[0], (int, Integer, slice)):
            frame = args[0]
            args = args[1:]
        else:
            smodule = vb.section_module(domain=self._base_space)
            frame = smodule.default_frame()
        # indexing:
        if isinstance(args, slice):
            indices = args
        elif isinstance(args[0], slice):
            indices = args[0]
        else:
            indices = args
        if isinstance(indices, slice):
            if indices.start is None and indices.stop is None:
                return [[self.connection_form(i, j, frame=frame)
                         for j in vb.irange()] for i in vb.irange()]
            else:
                raise NotImplementedError("[start:stop] syntax not "
                                          "implemented")
        if len(indices) != 2:
            raise ValueError("index must be a pair of integers")
        (i, j) = indices
        if i not in vb.irange() or j not in vb.irange():
            raise IndexError("index out of range")
        form = self.connection_form(indices[0], indices[1], frame=frame)
        return form

    def __setitem__(self, args, value):
        r"""
        Sets the components of ``self`` corresponding to the given indices.

        INPUT:

        - ``args`` -- list of indices (usually a pair of integers); if ``[:]``
          is provided, all the components are set
        - ``value`` -- the value to be set or a list of values if
          ``args = [:]``

        TESTS::

            sage: M = Manifold(2, 'M', start_index=1)
            sage: X.<x,y> = M.chart()
            sage: E = M.vector_bundle(2, 'E')
            sage: nab = E.bundle_connection('nabla')
            sage: e = E.local_frame('e')
            sage: a = M.one_form([x^2, x], name='a')
            sage: a.display()
            a = x^2 dx + x dy
            sage: nab[:] = 0
            sage: nab[1, 1] = a
            sage: nab[1, 1].display()
            connection (1,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[e, 2, 2] = a
            sage: nab[e, 2, 2].display()
            connection (2,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[:] = [[0, 0], [a, a]]
            sage: nab[e, 2, 1].display()
            connection (2,1) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy
            sage: nab[e, :] = [[a, a], [0, 0]]
            sage: nab[e, 1, 2].display()
            connection (1,2) of bundle connection nabla w.r.t. Local frame
             (E|_M, (e_1,e_2)) = x^2 dx + x dy

        """
        # extract frame from first index:
        vb = self._vbundle
        if isinstance(args, (int, Integer, slice)):
            smodule = vb.section_module(domain=self._base_space)
            frame = smodule.default_frame()
        elif not isinstance(args[0], (int, Integer, slice)):
            frame = args[0]
            args = args[1:]
        else:
            smodule = vb.section_module(domain=self._base_space)
            frame = smodule.default_frame()
        # determine indices:
        if isinstance(args, slice):
            indices = args
        elif isinstance(args[0], slice):
            indices = args[0]
        else:
            indices = args
        # set values:
        if isinstance(indices, (list, tuple)):
            if len(indices) != 2:
                raise IndexError("a tuple of integers must be provided")
            (i, j) = indices
            if i not in vb.irange() or j not in vb.irange():
                raise IndexError("index out of range")
            omega = self.set_connection_form(i, j, frame=frame)
            dom = frame.domain()
            dmodule = dom.diff_form_module(1)
            try:
                form = dmodule(value)
            except (TypeError, ValueError):
                msg = "{} must be convertible ".format(value)
                msg += "into an element of {}".format(dmodule)
                raise TypeError(msg)
            omega.copy_from(form)  # copy all components from value
        if isinstance(indices, slice):
            if indices.start is None and indices.stop is None:
                # check types:
                if isinstance(value, (int, Integer)) and value == 0:
                    for i in vb.irange():
                        for j in vb.irange():
                            self[frame, i, j] = 0
                elif not isinstance(value, (list, tuple)):
                    raise TypeError("in case of [:] syntax, zero or a "
                                    "list/tuple as value should be provided")
                elif any(not isinstance(row, (list, tuple)) for row in value):
                    raise TypeError("in case of [:] syntax, the list/tuple "
                                    "of value must contain lists/tuples")
                else:
                    # check lenghts:
                    rk = vb._rank
                    if len(value) != rk:
                        raise ValueError("value must have "
                                         "length {}".format(rk))
                    if any(len(row) != rk for row in value):
                        raise ValueError("lists in value must have length "
                                         "{}".format(rk))
                    # perform designation:
                    sind = vb._base_space._sindex
                    for i in vb.irange():
                        for j in vb.irange():
                            self[frame, i, j] = value[i - sind][j - sind]
            else:
                raise NotImplementedError("[start:stop] syntax not "
                                          "implemented")
