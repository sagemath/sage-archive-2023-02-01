import inspect

def class_graph(top_module, depth=5, name_filter=None, classes=None):
    """
    Return a dictionary with class names as keys and lists of class
    names as values. The values represent direct superclasses of the
    keys.

    This object is understood by e.g. the Graph and DiGraph constructors and
    can thus be used to create an inheritance graph for the given module.

    INPUT:
        top_module -- the module to start in (e.g. sage)
        depth -- maximal recursion depth (default: 5)
        name_filter -- e.g. 'sage.rings' to only consider classes in sage.rings
        classes -- optional dictionary to be filled in (it is also returned)

    EXAMPLE:
        sage: C = class_graph(sage.rings.polynomial.padics, depth=2)
        sage: C
        {'Polynomial_padic_capped_relative_dense': ['Polynomial_generic_domain'],
        'Polynomial_padic_flat': ['Polynomial_generic_dense']}
        sage: Graph(C)
        Graph on 4 vertices
        sage: DiGraph(C)
        DiGraph on 4 vertices
    """
    # termination
    if depth == 0: return

    if classes is None: classes = dict()

    if inspect.ismodule(top_module):
        if name_filter is None: name_filter = top_module.__name__
        top_module = top_module.__dict__.values()

    for X in top_module:
        if inspect.ismodule(X):
            if '.all' in X.__name__:
                continue
            class_graph(X.__dict__.values(), depth=depth-1, name_filter=name_filter, classes=classes)
        if inspect.isclass(X):
            B = X.__bases__
            if len(B) > 0:
                if X.__module__.startswith(name_filter):
                    classes[X.__name__] = [e.__name__ for e in B]
                class_graph(B, depth=depth, name_filter=name_filter, classes=classes)
    return classes
