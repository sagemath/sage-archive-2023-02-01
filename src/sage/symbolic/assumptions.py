from sage.structure.sage_object import SageObject
from sage.symbolic.ring import is_SymbolicVariable
_assumptions = []

class GenericDeclaration(SageObject):

    def __init__(self, var, assumption):
        """
        This class represents generic assumptions, such as a variable being
        an integer or a function being increasing. It passes such
        information to maxima's declare (wrapped in a context so it is able
        to forget).

        INPUT:


        -  ``var`` - the variable about which assumptions are
           being made

        -  ``assumption`` - a Maxima feature, either user
           defined or in the list given by maxima('features')


        EXAMPLES::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'integer')
            sage: decl.assume()
            sage: sin(x*pi)
            sin(pi*x)
            sage: sin(x*pi).simplify()
            0
            sage: decl.forget()
            sage: sin(x*pi)
            sin(pi*x)

        Here is the list of acceptable features::

            sage: maxima('features')
            [integer,noninteger,even,odd,rational,irrational,real,imaginary,complex,analytic,increasing,decreasing,oddfun,evenfun,posfun,constant,commutative,lassociative,rassociative,symmetric,antisymmetric,integervalued]
        """
        self._var = var
        self._assumption = assumption
        self._context = None

    def __repr__(self):
        """
        EXAMPLES::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: GenericDeclaration(x, 'foo')
            x is foo
        """
        return "%s is %s" % (self._var, self._assumption)

    def __cmp__(self, other):
        """
        TESTS::

            sage: from sage.symbolic.assumptions import GenericDeclaration as GDecl
            sage: var('y')
            y
            sage: GDecl(x, 'integer') == GDecl(x, 'integer')
            True
            sage: GDecl(x, 'integer') == GDecl(x, 'rational')
            False
            sage: GDecl(x, 'integer') == GDecl(y, 'integer')
            False
        """
        if isinstance(self, GenericDeclaration) and isinstance(other, GenericDeclaration):
            return cmp( (self._var, self._assumption),
                        (other._var, other._assumption) )
        else:
            return cmp(type(self), type(other))

    def assume(self):
        """
        TEST::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'even')
            sage: decl.assume()
            sage: cos(x*pi).simplify()
            1
            sage: decl.forget()
        """
        from sage.calculus.calculus import maxima
        if self._context is None:
            # We get the list here because features may be added with time.
            valid_features = list(maxima("features"))
            if self._assumption not in [repr(x).strip() for x in list(valid_features)]:
                raise ValueError, "%s not a valid assumption, must be one of %s" % (self._assumption, valid_features)
            cur = maxima.get("context")
            self._context = maxima.newcontext('context' + maxima._next_var_name())
            maxima.eval("declare(%s, %s)" % (repr(self._var), self._assumption))
            maxima.set("context", cur)

        if not self in _assumptions:
            maxima.activate(self._context)
            _assumptions.append(self)

    def forget(self):
        """
        TEST::

            sage: from sage.symbolic.assumptions import GenericDeclaration
            sage: decl = GenericDeclaration(x, 'odd')
            sage: decl.assume()
            sage: cos(pi*x)
            cos(pi*x)
            sage: cos(pi*x).simplify()
            -1
            sage: decl.forget()
            sage: cos(x*pi).simplify()
            cos(pi*x)
        """
        from sage.calculus.calculus import maxima
        try:
            _assumptions.remove(self)
        except ValueError:
            return
        if self._context is not None:
            maxima.deactivate(self._context)

def preprocess_assumptions(args):
    """
    Turns a list of the form (var1, var2, ..., 'property') into a
    sequence of declarations (var1 is property), (var2 is property),
    ...

    EXAMPLES::

        sage: from sage.symbolic.assumptions import preprocess_assumptions
        sage: preprocess_assumptions([x, 'integer', x > 4])
        [x is integer, x > 4]
        sage: var('x,y')
        (x, y)
        sage: preprocess_assumptions([x, y, 'integer', x > 4, y, 'even'])
        [x is integer, y is integer, x > 4, y is even]
    """
    args = list(args)
    last = None
    for i, x in reversed(list(enumerate(args))):
        if isinstance(x, str):
            del args[i]
            last = x
        elif ((not hasattr(x, 'assume') or is_SymbolicVariable(x))
              and last is not None):
            args[i] = GenericDeclaration(x, last)
        else:
            last = None
    return args

def assume(*args):
    """
    Make the given assumptions.

    INPUT:

    -  ``*args`` - assumptions

    EXAMPLES::

        sage: assume(x > 0)
        sage: bool(sqrt(x^2) == x)
        True
        sage: forget()
        sage: bool(sqrt(x^2) == x)
        False

    An integer constraint::

        sage: var('n, P, r, r2')
        (n, P, r, r2)
        sage: assume(n, 'integer')
        sage: c = P*e^(r*n)
        sage: d = P*(1+r2)^n
        sage: solve(c==d,r2)
        [r2 == e^r - 1]

    ::

        sage: sin(n*pi)
        sin(pi*n)
        sage: sin(n*pi).simplify()
        0
        sage: forget()
        sage: sin(n*pi).simplify()
        sin(pi*n)

    TESTS:

    Test that you can do two non-relational
    declarations at once (fixing Trac ticket 7084)::

        sage: var('m')
        m
        sage: assume(n, 'integer'); assume(m, 'integer')
        sage: sin(n*pi).simplify()
        0
        sage: sin(m*pi).simplify()
        0
        sage: forget()
    """
    for x in preprocess_assumptions(args):
        if isinstance(x, (tuple, list)):
            assume(*x)
        else:
            try:
                x.assume()
            except KeyError:
                raise TypeError, "assume not defined for objects of type '%s'"%type(x)

def forget(*args):
    """
    Forget the given assumption, or call with no arguments to forget
    all assumptions.

    Here an assumption is some sort of symbolic constraint.

    INPUT:


    -  ``*args`` - assumptions (default: forget all
       assumptions)


    EXAMPLES: We define and forget multiple assumptions::

        sage: var('x,y,z')
        (x, y, z)
        sage: assume(x>0, y>0, z == 1, y>0)
        sage: list(sorted(assumptions(), lambda x,y:cmp(str(x),str(y))))
        [x > 0, y > 0, z == 1]
        sage: forget(x>0, z==1)
        sage: assumptions()
        [y > 0]
        sage: assume(y, 'even')
        sage: assumptions()
        [y > 0, y is even]
        sage: cos(y*pi).simplify()
        1
        sage: forget()
        sage: cos(y*pi).simplify()
        cos(pi*y)
        sage: assumptions()
        []
    """
    if len(args) == 0:
        _forget_all()
        return
    for x in preprocess_assumptions(args):
        if isinstance(x, (tuple, list)):
            assume(*x)
        else:
            try:
                x.forget()
            except KeyError:
                raise TypeError, "forget not defined for objects of type '%s'"%type(x)

def assumptions():
    """
    List all current symbolic assumptions.

    EXAMPLES::

        sage: var('x,y,z, w')
        (x, y, z, w)
        sage: forget()
        sage: assume(x^2+y^2 > 0)
        sage: assumptions()
        [x^2 + y^2 > 0]
        sage: forget(x^2+y^2 > 0)
        sage: assumptions()
        []
        sage: assume(x > y)
        sage: assume(z > w)
        sage: list(sorted(assumptions(), lambda x,y:cmp(str(x),str(y))))
        [x > y, z > w]
        sage: forget()
        sage: assumptions()
        []
    """
    return list(_assumptions)

def _forget_all():
    """
    Forget all symbolic assumptions.

    This is called by ``forget()``.

    EXAMPLES::

        sage: var('x,y')
        (x, y)
        sage: assume(x > 0, y < 0)
        sage: bool(x*y < 0)      # means definitely true
        True
        sage: bool(x*y > 0)      # might not be true
        False
        sage: forget()    # implicitly calls _forget_all
        sage: bool(x*y < 0)      # might not be true
        False
        sage: bool(x*y > 0)      # might not be true
        False
    """
    from sage.calculus.calculus import maxima
    global _assumptions
    if len(_assumptions) == 0:
        return
    try:
        maxima._eval_line('forget(facts());')
    except TypeError:
        pass
    #maxima._eval_line('forget([%s]);'%(','.join([x._maxima_init_() for x in _assumptions])))
    for x in _assumptions:
        if isinstance(x, GenericDeclaration):
            # these don't show up in facts()
            x.forget()
    _assumptions = []
