cpdef call_registered_function(unsigned int serial,
                           int nargs,
                           list args,
                           bint hold,
                           bint allow_numeric_result,
                           result_parent):
    r"""
    Call a ginac function.

    INPUT:

    - ``serial`` - serial number of the function

    - ``nargs`` - declared number of args (0 is variadic)

    - ``args`` - a list of :class:`Expression`s

    - ``allow_numeric_result`` - if ``True``, keep numeric results numeric;
      if ``False``, make all results symbolic expressions

    - ``result_parent`` - an instance of :class:`SymbolicRing`
    """
    cdef Py_ssize_t i
    cdef GEx res
    cdef GExVector vec
    if nargs == 0 or nargs > 3:
        for i from 0 <= i < len(args):
            vec.push_back((<Expression>args[i])._gobj)
        res = g_function_evalv(serial, vec, hold)
    elif nargs == 1:
        res = g_function_eval1(serial,
                (<Expression>args[0])._gobj, hold)
    elif nargs == 2:
        res = g_function_eval2(serial, (<Expression>args[0])._gobj,
                (<Expression>args[1])._gobj, hold)
    elif nargs == 3:
        res = g_function_eval3(serial,
                (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                (<Expression>args[2])._gobj, hold)

    if allow_numeric_result and is_a_numeric(res):
        return py_object_from_numeric(res)

    return new_Expression_from_GEx(result_parent, res)


cpdef unsigned int find_registered_function(name, int nargs) except -1:
    r"""
    Look up a function registered with ginac.

    Raise a ``ValueError`` if the function is not registered.

    OUTPUT:

    - serial number of the function
    """
    try:
        return find_function(str_to_bytes(name), nargs)
    except RuntimeError as err:
        raise ValueError("cannot find GiNaC function with name %s and %s arguments" % (name, nargs))


cpdef unsigned int register_or_update_function(self, name, latex_name, int nargs,
                                               evalf_params_first, bint update) except -1:
    r"""
    Register the function ``self`` with ginac.

    OUTPUT:

    - serial number of the function
    """
    cdef GFunctionOpt opt
    cdef unsigned int serial

    if update:
        serial = find_registered_function(name, nargs)
        opt = g_registered_functions().index(serial)
    else:
        opt = g_function_options_args(str_to_bytes(name), nargs)

    if hasattr(self, '_eval_'):
        opt.eval_func(self)

    if not evalf_params_first:
        opt.do_not_evalf_params()

    if hasattr(self, '_subs_'):
        opt.subs_func(self)

    if hasattr(self, '_evalf_'):
        opt.evalf_func(self)

    if hasattr(self, '_conjugate_'):
        opt.conjugate_func(self)

    if hasattr(self, '_real_part_'):
        opt.real_part_func(self)

    if hasattr(self, '_imag_part_'):
        opt.imag_part_func(self)

    if hasattr(self, '_derivative_'):
        opt.derivative_func(self)

    if hasattr(self, '_tderivative_'):
        opt.do_not_apply_chain_rule()
        opt.derivative_func(self)

    if hasattr(self, '_power_'):
        opt.power_func(self)

    if hasattr(self, '_series_'):
        opt.series_func(self)

    # custom print functions are called from python
    # so we don't register them with the ginac function_options object

    if latex_name:
        opt.latex_name(str_to_bytes(latex_name))

    if not update:
        serial = g_register_new(opt)

    g_foptions_assign(g_registered_functions().index(serial), opt)

    return serial
