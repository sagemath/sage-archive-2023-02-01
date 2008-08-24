r"""

Support for symbolic functions.



"""
include "../ext/interrupt.pxi"
include "../ext/stdsage.pxi"
include "../ext/cdefs.pxi"
include "../libs/ginac/decl.pxi"


from expression cimport new_Expression_from_GEx, Expression
from ring import NSR as SR

cdef class SFunction:
    """

    EXAMPLES:
        sage: from sage.symbolic.function import function as nfunction
        sage: foo = nfunction("foo", 2)
        sage: x,y,z = var("x y z", ns=1)
        sage: foo(x,y) + foo(y,z)^2
        foo(x,y) + foo(y,z)^2

    """
    cdef unsigned int serial
    cdef object name
    cdef int nargs

    def __init__(self, name, nargs=0):
        try:
            self.serial = find_function(name, nargs)
        except ValueError, err:
            self.serial = g_register_new(g_function_options_args(name, nargs))

        self.name = name
        if nargs:
            self.nargs = nargs
        else:
            self.nargs = 0

    def __repr__(self):
        """
        EXAMPLES:
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2); foo
            foo
        """
        return self.name

    def __call__(self, *args, coerce=True):
        """

        EXAMPLES:
            sage: from sage.symbolic.function import function as nfunction
            sage: foo = nfunction("foo", 2)
            sage: x,y,z = var("x y z", ns=1)
            sage: foo(x,y)
            foo(x,y)

            sage: foo(y)
            Traceback (most recent call last):
            ...
            ValueError: Symbolic function foo expects 2 arguments, got 1.

            sage: bar = nfunction("bar")
            sage: bar(x)
            bar(x)
            sage: bar(x,y)
            bar(x,y)

        """
        cdef GExVector vec
        if self.nargs > 0 and len(args) != self.nargs:
            raise ValueError, "Symbolic function %s expects %s arguments, got %s."%(self.name, self.nargs, len(args))

        if coerce:
            args = map(SR, args)

        if self.nargs == 0:
            for i from 0 <= i < len(args):
                vec.push_back((<Expression>args[i])._gobj)
            exp = new_Expression_from_GEx(g_function_evalv(self.serial, vec))
            return exp

        elif self.nargs == 1:
            return new_Expression_from_GEx(g_function_eval1(self.serial,
                    (<Expression>args[0])._gobj))
        elif self.nargs == 2:
            return new_Expression_from_GEx(g_function_eval2(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj))
        elif self.nargs == 3:
            return new_Expression_from_GEx(g_function_eval3(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj))
        elif self.nargs == 4:
            return new_Expression_from_GEx(g_function_eval4(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj))
        elif self.nargs == 5:
            return new_Expression_from_GEx(g_function_eval5(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj))
        elif self.nargs == 6:
            return new_Expression_from_GEx(g_function_eval6(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj))
        elif self.nargs == 7:
            return new_Expression_from_GEx(g_function_eval7(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj))
        elif self.nargs == 8:
            return new_Expression_from_GEx(g_function_eval8(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj))
        elif self.nargs == 9:
            return new_Expression_from_GEx(g_function_eval9(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj))
        elif self.nargs == 10:
            return new_Expression_from_GEx(g_function_eval10(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj))
        elif self.nargs == 11:
            return new_Expression_from_GEx(g_function_eval11(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj))
        elif self.nargs == 12:
            return new_Expression_from_GEx(g_function_eval12(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj))
        elif self.nargs == 13:
            return new_Expression_from_GEx(g_function_eval13(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj,
                    (<Expression>args[12])._gobj))
        elif self.nargs == 14:
            return new_Expression_from_GEx(g_function_eval14(self.serial,
                    (<Expression>args[0])._gobj, (<Expression>args[1])._gobj,
                    (<Expression>args[2])._gobj, (<Expression>args[3])._gobj,
                    (<Expression>args[4])._gobj, (<Expression>args[5])._gobj,
                    (<Expression>args[6])._gobj, (<Expression>args[7])._gobj,
                    (<Expression>args[8])._gobj, (<Expression>args[9])._gobj,
                    (<Expression>args[10])._gobj, (<Expression>args[11])._gobj,
                    (<Expression>args[12])._gobj, (<Expression>args[13])._gobj))

function = SFunction

cdef new_SFunction_from_serial(int serial, char* name, int nargs):
    cdef SFunction f = PY_NEW(SFunction)
    f.serial = serial
    f.name = name
    f.nargs = nargs
    return f

cos = new_SFunction_from_serial(cos_serial, "cos", 1)
