"""
QQbar decorators

Python decorators for use with the algebraic field QQbar.

AUTHORS:

- Brent Baccala (7 Jun 2018) -- handle_AA_and_QQbar

Decorators
==========
"""

def handle_AA_and_QQbar(func):
    r"""
    Decorator to call a function that only accepts arguments in number fields.

    The argument list is scanned for ideals and/or polynomials over algebraic
    fields (``QQbar`` or ``AA``).  If any exist, they are converted to a common
    number field before calling the function, and the results are converted back.
    Lists, sets, and tuples are converted recursively.

    Keyword arguments are currently not converted.

    Only works for functions that don't depend on polynomial factorization.
    Thus, :meth:`intersection` is suitable, but :meth:`associated_primes`
    is not.
    """

    from sage.misc.decorators import sage_wraps

    @sage_wraps(func)
    def wrapper(*args, **kwds):

        from sage.misc.flatten import flatten
        from sage.rings.polynomial.multi_polynomial import MPolynomial
        from sage.rings.ideal import Ideal, Ideal_generic
        from sage.rings.qqbar import AlgebraicField_common, number_field_elements_from_algebraics

        if not any([isinstance(a, (MPolynomial, Ideal_generic))
                    and isinstance(a.base_ring(), AlgebraicField_common) for a in args]):
            return func(*args, **kwds)

        orig_elems = []
        for a in flatten(args, ltypes=(list, tuple, set)):
            if isinstance(a, Ideal_generic):
                for g in a.gens():
                    orig_elems.extend(g.coefficients())
            elif isinstance(a, MPolynomial):
                orig_elems.extend(a.coefficients())

        # We need minimal=True if these elements are over AA, because
        # same_field=True might trigger an exception otherwise.

        numfield, new_elems, morphism = number_field_elements_from_algebraics(orig_elems, same_field=True, minimal=True)

        elem_dict = dict(zip(orig_elems, new_elems))

        def polynomial_map(p):
            return p.map_coefficients(elem_dict.__getitem__, new_base_ring=numfield)

        def ideal_map(ideal):
            return Ideal([polynomial_map(g) for g in ideal.gens()])

        def generic_map(item):
            if isinstance(item, Ideal_generic):
                return ideal_map(item)
            elif isinstance(item, MPolynomial):
                return polynomial_map(item)
            elif isinstance(item, (list, tuple)):
                return map(generic_map, item)
            elif isinstance(item, set):
                return set(map(generic_map, list(item)))
            else:
                return item

        def inverse_polynomial_map(p):
            return p.map_coefficients(morphism)

        def inverse_ideal_map(ideal):
            return Ideal([inverse_polynomial_map(g) for g in ideal.gens()])

        def inverse_generic_map(item):
            if isinstance(item, Ideal_generic):
                return inverse_ideal_map(item)
            elif isinstance(item, MPolynomial):
                return inverse_polynomial_map(item)
            elif isinstance(item, (list, tuple)):
                return map(inverse_generic_map, item)
            elif isinstance(item, set):
                return set(map(inverse_generic_map, list(item)))
            else:
                return item

        args = map(generic_map, args)

        return inverse_generic_map(func(*args, **kwds))

    return wrapper
