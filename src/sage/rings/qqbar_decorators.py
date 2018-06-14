"""
QQbar decorators

Python decorators for use with the algebraic field QQbar.

AUTHORS:

- Brent Baccala (7 Jun 2018) -- handle_AA_and_QQbar

Decorators
==========
"""

from sage.misc.decorators import decorator_keywords, sage_wraps

@decorator_keywords
def handle_AA_and_QQbar(func, factor_field=False):
    r"""
    Decorator to call a function that only accepts arguments in number fields.

    The argument list is scanned for ideals and/or polynomials over algebraic
    fields (``QQbar`` or ``AA``).  If any exist, they are converted to a common
    number field before calling the function, and the results are converted back.
    Lists, sets, and tuples are converted recursively.

    Keyword arguments are currently not converted.

    If the decorator is used with the optional argument `factor_field=True`,
    the polynomials and/or ideals are converted to a larger number field in
    which they can be factored.  While slower, this option is required for
    methods that depend the ability to factor their arguments.
    """

    @sage_wraps(func)
    def wrapper(*args, **kwds):

        from sage.misc.flatten import flatten
        from sage.rings.polynomial.multi_polynomial import MPolynomial
        from sage.rings.ideal import Ideal, Ideal_generic
        from sage.rings.qqbar import AlgebraicField_common, number_field_elements_from_algebraics

        if not any([isinstance(a, (MPolynomial, Ideal_generic))
                    and isinstance(a.base_ring(), AlgebraicField_common) for a in args]):
            return func(*args, **kwds)

        polynomials = []

        for a in flatten(args, ltypes=(list, tuple, set)):
            if isinstance(a, Ideal_generic):
                polynomials.extend(a.gens())
            elif isinstance(a, MPolynomial):
                polynomials.append(a)

        orig_elems = flatten([p.coefficients() for p in polynomials])

        if factor_field is True:
            for p in polynomials:
                orig_elems.extend(flatten([f[0].coefficients() for f in p.factor()]))

        # We need minimal=True if these elements are over AA, because
        # same_field=True might trigger an exception otherwise.

        numfield, new_elems, morphism = number_field_elements_from_algebraics(orig_elems, same_field=True, minimal=True)

        elem_dict = dict(zip(orig_elems, new_elems))

        def forward_map(item):
            if isinstance(item, Ideal_generic):
                return Ideal([forward_map(g) for g in item.gens()])
            elif isinstance(item, MPolynomial):
                return item.map_coefficients(elem_dict.__getitem__, new_base_ring=numfield)
            elif isinstance(item, list):
                return map(forward_map, item)
            elif isinstance(item, tuple):
                return tuple(map(forward_map, item))
            elif isinstance(item, set):
                return set(map(forward_map, list(item)))
            else:
                return item

        def reverse_map(item):
            if isinstance(item, Ideal_generic):
                return Ideal([reverse_map(g) for g in item.gens()])
            elif isinstance(item, MPolynomial):
                return item.map_coefficients(morphism)
            elif isinstance(item, list):
                return map(reverse_map, item)
            elif isinstance(item, tuple):
                return tuple(map(reverse_map, item))
            elif isinstance(item, set):
                return set(map(reverse_map, list(item)))
            else:
                return item

        args = map(forward_map, args)

        return reverse_map(func(*args, **kwds))

    return wrapper
