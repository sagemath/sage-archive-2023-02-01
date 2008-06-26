r"""
Sage Input Formatting

AUTHORS:
    -- Carl Witty (2008-04): new file

This module provides a function \function{sage_input}.  This function
takes an arbitrary \sage value and produces a sequence of commands
that, if typed at the \code{sage:} prompt, will recreate the value.
(If this is not implemented for a particular value, then an exception
is raised instead.)  This might be useful in understanding a part of
Sage, or for debugging.  (For instance, if you have a value produced
in a complicated way in the middle of a debugging session, you could
use \function{sage_input} to find a simple way to produce the same
value.)  We attempt to produce commands that are readable and
idiomatic.

sage: sage_input(3)
3
sage: sage_input((polygen(RR) + RR(pi))^2, verify=True)
# Verified
R.<x> = RR[]
x^2 + 6.2831853071795862*x + 9.8696044010893580

With \code{verify=True}, \function{sage_input} also verifies the
results, by calling sage_eval on the result and verifying that it is
equal to the input.

sage: sage_input(GF(2)(1), verify=True)
# Verified
GF(2)(1)

We can generate code that works without the preparser, with
\code{preparse=False}; or we can generate code that will work whether
or not the preparser is enabled, with \code{preparse=None}.
Generating code with \code{preparse=False} may be useful to see how to
create a certain value in a Python or Cython source file.

sage: sage_input(5, verify=True)
# Verified
5
sage: sage_input(5, preparse=False)
ZZ(5)
sage: sage_input(5, preparse=None)
ZZ(5)
sage: sage_input(5r, verify=True)
# Verified
5r
sage: sage_input(5r, preparse=False)
5
sage: sage_input(5r, preparse=None)
int(5)

Adding \function{sage_input} support to your own classes is
straightforward.  You need to add a \method{_sage_input_} method which
returns a \class{SageInputExpression} (henceforth abbreviated as SIE)
which will reconstruct this instance of your class.

A \method{_sage_input_} method takes two parameters, conventionally
named \var{sib} and \var{coerced}.  The first argument is a
\class{SageInputBuilder}; it has methods to build SIEs.  The second
argument, \var{coerced}, is a boolean.  This is only useful if your
class is a subclass of \class{Element} (although it is always
present).  If \var{coerced} is \code{False}, then your method must
generate an expression which will evaluate to a value of the correct
type with the correct parent.  If \var{coerced} is \code{True}, then
your method may generate an expression of a type that has a canonical
coercion to your type.

Let's work through some examples.  We'll build a sequence of functions
that would be acceptable as \method{_sage_input_} methods for the
\class{Rational} class.

Here's the first and simplest version.

sage: def qq_sage_input_v1(self, sib, coerced):
...       return sib(self.numerator())/sib(self.denominator())

We see that given a \class{SageInputBuilder} \var{sib}, you can
construct a SIE for a value \var{v} simply with \code{sib(v)}, and you
can construct a SIE for a quotient with the division operator.  Of course,
the other operators also work, and so do function calls, method calls,
subscripts, etc.

We'll test with the following code, which you don't need to understand.
(It produces a list of 8 results, showing the formatted versions of
-5/7 and 3, with the preparser either enabled or disabled and either
with or without an automatic coercion to QQ.)

sage: from sage.misc.sage_input import SageInputBuilder
sage: def test_qq_formatter(fmt):
...       results = []
...       for v in [-5/7, QQ(3)]:
...           for pp in [False, True]:
...               for coerced in [False, True]:
...                   sib = SageInputBuilder(preparse=pp)
...                   results.append(sib.result(fmt(v, sib, coerced)))
...       return results

sage: test_qq_formatter(qq_sage_input_v1)
[-ZZ(5)/ZZ(7), -ZZ(5)/ZZ(7), -5/7, -5/7, ZZ(3)/ZZ(1), ZZ(3)/ZZ(1), 3/1, 3/1]

Let's try for some shorter, perhaps nicer-looking output.  We'll start
by getting rid of the \code{ZZ} in the denominators; even without the
preparser, \code{-ZZ(5)/7 == -ZZ(5)/ZZ(7)}.

sage: def qq_sage_input_v2(self, sib, coerced):
...       return sib(self.numerator())/sib.int(self.denominator())

The \method{int} method on \class{SageInputBuilder} returns a SIE for
an integer that is always represented in the simple way, without
coercions.  (So, depending on the preparser mode, it might read in as an
\class{Integer}, an \class{int}, or a \class{long}.)

sage: test_qq_formatter(qq_sage_input_v2)
[-ZZ(5)/7, -ZZ(5)/7, -5/7, -5/7, ZZ(3)/1, ZZ(3)/1, 3/1, 3/1]

Next let's get rid of the divisions by 1.  These are more complicated,
since if we're not careful we'll get results in \ZZ instead of \QQ.

sage: def qq_sage_input_v3(self, sib, coerced):
...       if self.denominator() == 1:
...           if coerced:
...               return sib.int(self.numerator())
...           else:
...               return sib.name('QQ')(sib.int(self.numerator()))
...       return sib(self.numerator())/sib.int(self.denominator())

We see that the \method{name} method gives an SIE representing a \sage
constant or function.

sage: test_qq_formatter(qq_sage_input_v3)
[-ZZ(5)/7, -ZZ(5)/7, -5/7, -5/7, QQ(3), 3, QQ(3), 3]

This is the prettiest output we're going to get, but let's make one
further refinement.  Other \class{_sage_input_} methods, like the one
for polynomials, analyze the structure of SIEs; they work better (give
prettier output) if negations are at the outside.

sage: def qq_sage_input_v4(self, sib, coerced):
...       num = self.numerator()
...       neg = (num < 0)
...       if neg: num = -num
...       if self.denominator() == 1:
...           if coerced:
...               v = sib.int(num)
...           else:
...               v = sib.name('QQ')(sib.int(num))
...       else:
...           v = sib(num)/sib.int(self.denominator())
...       if neg: v = -v
...       return v

sage: test_qq_formatter(qq_sage_input_v4)
[-ZZ(5)/7, -ZZ(5)/7, -5/7, -5/7, QQ(3), 3, QQ(3), 3]

"""

from sage.misc.functional import parent

##########################################################################
#
#       Copyright (C) 2008 Carl Witty <Carl.Witty@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#
##########################################################################

def sage_input(x, preparse=True, verify=False, allow_locals=False):
    r"""
    INPUTS:
        x -- the value we want to find an ``input form'' for
        preparse -- (default \code{True}) Whether to generate code
            that requires the preparser.  With \code{True}, generated
            code requires the preparser.  With \code{False}, generated
            code requires that the preparser not be used.  With \code{None},
            generated code will work whether or not the preparser is used.
        verify -- (default \code{False}) If \code{True}, then the
            answer will be evaluated with \function{sage_eval}, and
            an exception will be raised if the result is not equal to
            the original value.  (In fact, for \code{verify=True},
            \function{sage_input} is effectively run three times,
            with \var{preparse} set to \code{True}, \code{False}, and
            \code{None}, and all three results are checked.)  This is
            particularly useful for doctests.
        allow_locals -- (default \code{False}) If \code{True}, then
            values that \function{sage_input} cannot handle are returned
            in a dictionary, and the returned code assumes that this
            dictionary is passed as the \var{locals} parameter of
            \function{sage_eval}.  (Otherwise, if \function{sage_input}
            cannot handle a value, an exception is raised.)

    EXAMPLES:
        sage: sage_input(GF(2)(1))
        GF(2)(1)
        sage: sage_input((GF(2)(0), GF(2)(1)), verify=True)
        # Verified
        GF_2 = GF(2)
        (GF_2(0), GF_2(1))

    When the preparser is enabled, we use the \sage generator syntax.

        sage: K.<x> = GF(5)[]
        sage: sage_input(x^3 + 2*x, verify=True)
        # Verified
        R.<x> = GF(5)[]
        x^3 + 2*x
        sage: sage_input(x^3 + 2*x, preparse=False)
        R = GF(5)['x']
        x = R.gen()
        x**3 + 2*x

    The result of \function{sage_input} is actually a pair of strings with
    a special \method{__repr__} method to print nicely.

        sage: r = sage_input(RealField(20)(pi), verify=True)
        sage: r
        # Verified
        RealField(20)(3.1415939)
        sage: isinstance(r, tuple)
        True
        sage: len(r)
        2
        sage: tuple(r)
        ('# Verified\n', 'RealField(20)(3.1415939)')

    We cannot find an input form for a function.

        sage: sage_input((3, lambda x: x))
        Traceback (most recent call last):
        ...
        ValueError: Can't convert <function <lambda> at 0x...> to sage_input form

    But we can have \function{sage_input} continue anyway, and return
    an input form for the rest of the expression, with
    \code{allow_locals=True}.

        sage: r = sage_input((3, lambda x: x), verify=True, allow_locals=True)
        sage: r
        LOCALS:
            _sil1: <function <lambda> at 0x...>
        # Verified
        (3, _sil1)
        sage: tuple(r)
        ('# Verified\n', '(3, _sil1)', {'_sil1': <function <lambda> at 0x...>})
    """
    if not verify:
        sib = SageInputBuilder(allow_locals=allow_locals, preparse=preparse)
        return sib.result(sib(x))

    # In verify mode, we actually compute and verify the answer with
    # all three settings of preparse.

    for pp in (True, False, None):
        sib = SageInputBuilder(allow_locals=allow_locals, preparse=pp)
        ans = sib.result(sib(x))
        verify_si_answer(x, ans, pp)
        if pp == preparse:
            ans_l = list(ans)
            ans_l[0] = '# Verified\n' + ans_l[0]
            final_answer = SageInputAnswer(*ans_l)

    return final_answer

class SageInputBuilder:
    r"""
    An instance of this class is passed to \method{_sage_input_} methods.
    It keeps track of the current state of the \method{_sage_input_} process,
    and contains many utility methods for building \class{SageInputExpression}
    objects.

    In normal use, instances of \class{SageInputBuilder} are created
    internally by \function{sage_input}, but it may be useful to create
    an instance directly for testing or doctesting.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

    We can create a \class{SageInputBuilder}, use it to create some
    \class{SageInputExpression}s, and get a result.  (As mentioned
    above, this is only useful for testing or doctesting; normally
    you would just use \function{sage_input}.)

        sage: sib = SageInputBuilder()
        sage: sib.result((sib(3) + sib(4)) * (sib(5) + sib(6)))
        (3 + 4)*(5 + 6)
    """

    def __init__(self, allow_locals=False, preparse=True):
        r"""
        Initialize an instance of \class{SageInputBuilder}.

        In normal use, instances of \class{SageInputBuilder} are created
        internally by \function{sage_input}, but it may be useful to create
        an instance directly for testing or doctesting.

        INPUTS:
            allow_locals -- (default \code{False}) If true, then values
                that cannot be converted to input form will be stored in
                a dictionary, which must be passed as the \var{locals}
                when evaluating the result.
            preparse -- (default \code{True}) If true, then the result
                will assume that the preparser is enabled.  If false, then
                the result will assume that the preparser is disabled.
                If \code{None}, then the result will work whether or
                not the preparser is enabled.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: SageInputBuilder().preparse()
            True
            sage: SageInputBuilder(preparse=False).preparse()
            False
        """
        self._allow_locals = allow_locals
        self._preparse = preparse
        self._cached_types = set()
        self._cache = {}
        self._parent_gens = {}
        self._next_local = 1
        self._locals = {}

    def __call__(self, x, coerced=False):
        r"""
        Tries to convert an arbitrary value \var{x} into a
        \class{SageInputExpression} (an SIE).

        We first check to see if an SIE has been cached for \var{x};
        if so, we return it.  If \var{x} is already an SIE, we return
        it unchanged.

        If \var{x} has a \method{_sage_input_} method, we call that
        method.

        Otherwise, if \var{x} is a value of some Python type that
        we know how to deal with, we convert it directly.

        Finally, for values we don't know how to convert, if
        \code{self._allow_locals} is true, we add it to a
        ``locals'' dictionary.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib(sib(3)))
            3

            sage: sib = SageInputBuilder()
            sage: sib.result(sib(GF(17)(5)))
            GF(17)(5)

        The argument \code{coerced=True} will get passed to the
        \method{_sage_input_} method of the argument.
            sage: sib = SageInputBuilder()
            sage: sib.result(sib(GF(17)(5), True))
            5

        Since \function{sage_input} directly calls this method, all
        of the following are indirect doctests.
            sage: sage_input(True)
            True
            sage: sage_input(-5r, verify=True)
            # Verified
            -5r
            sage: sage_input(7r, preparse=False, verify=True)
            # Verified
            7
            sage: sage_input(-11r, preparse=None, verify=True)
            # Verified
            -int(11)
            sage: sage_input(long(-5), verify=True)
            # Verified
            -long(5)
            sage: sage_input(long(-7), preparse=False, verify=True)
            # Verified
            -7L
            sage: sage_input(long(11), preparse=None, verify=True)
            # Verified
            long(11)
            sage: sage_input(long(2^70), verify=True)
            # Verified
            1180591620717411303424r
            sage: sage_input(-long(2^80), preparse=False, verify=True)
            # Verified
            -1208925819614629174706176
            sage: sage_input(long(2^75), preparse=None, verify=True)
            # Verified
            long(37778931862957161709568)
            sage: sage_input("Hello, world\n", verify=True)
            # Verified
            'Hello, world\n'
            sage: sage_input("'", verify=True)
            # Verified
            "'"
            sage: sage_input('"', verify=True)
            # Verified
            '"'
            sage: sage_input(''' "'Hi,' she said." ''', verify=True)
            # Verified
            ' "\'Hi,\' she said." '
            sage: sage_input('Icky chars: \0\n\t\b\'\"\200\300\234', verify=True)
            # Verified
            'Icky chars: \x00\n\t\x08\'"\x80\xc0\x9c'
            sage: sage_input((2, 3.5, 'Hi'), verify=True)
            # Verified
            (2, RR(3.5000000000000000), 'Hi')
            sage: sage_input(lambda x: x)
            Traceback (most recent call last):
            ...
            ValueError: Can't convert <function <lambda> at 0x...> to sage_input form
            sage: sage_input(lambda x: x, allow_locals=True, verify=True)
            LOCALS:
              _sil1: <function <lambda> at 0x...>
            # Verified
            _sil1
        """
        # We want to look up x in our cache, to see if we've seen it before.
        # However, we don't want to assume that hashing x is always
        # efficient, so we only try the lookup if some value of the same
        # type as x has been cached.
        if type(x) in self._cached_types:
            v = self._cache.get((parent(x), x))
            if v is not None: return v

        if isinstance(x, SageInputExpression):
            return x

        if hasattr(x, '_sage_input_'):
            return x._sage_input_(self, coerced)

        if isinstance(x, bool):
            return SIE_literal_stringrep(self, str(x))

        if isinstance(x, int) or \
                (isinstance(x, long) and isinstance(int(x), long)):
            # For longs that don't fit in an int, we just use the int
            # code; it will get extended to long automatically.
            if self._preparse == True:
                if x < 0:
                    return -SIE_literal_stringrep(self, str(-x) + 'r')
                else:
                    return SIE_literal_stringrep(self, str(x) + 'r')
            elif self._preparse == False:
                return self.int(x)
            else:
                tyname = 'int' if isinstance(x, int) else 'long'
                if x < 0:
                    return -self.name(tyname)(self.int(-x))
                else:
                    return self.name(tyname)(self.int(x))

        if isinstance(x, long):
            # This must be a long that does fit in an int, so we need either
            # long(x) or an 'L' suffix.
            # With the current preparser, 1Lr does not work.
            # 1rL does work; but that's just ugly, so I don't use it.
            if self._preparse == False:
                if x < 0:
                    return -SIE_literal_stringrep(self, str(-x) + 'L')
                else:
                    return SIE_literal_stringrep(self, str(x) + 'L')
            else:
                if x < 0:
                    return -self.name('long')(self.int(-x))
                else:
                    return self.name('long')(self.int(x))

        if isinstance(x, str):
            return SIE_literal_stringrep(self, repr(x))

        if isinstance(x, tuple):
            return SIE_tuple(self, map(self, x))

        if self._allow_locals:
            loc = self._next_local
            self._next_local += 1
            loc_name = '_sil%d' % loc
            self._locals[loc_name] = x
            return SIE_literal_stringrep(self, loc_name)
        else:
            raise ValueError, "Can't convert %r to sage_input form"%x

    def preparse(self):
        r"""
        Checks the preparse status of this \class{SageInputBuilder}.
        (\code{True} if the preparser will be enabled, \code{False}
        if it will be disabled, and \code{None} if the result must
        work whether or not the preparser is enabled.)

        For example, this is useful in the \method{_sage_input_}
        methods of \class{Integer} and \class{RealNumber}; but most
        \method{_sage_input_} methods will not need to examine this.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder
            sage: SageInputBuilder().preparse()
            True
            sage: SageInputBuilder(preparse=False).preparse()
            False
        """
        return self._preparse

    def int(self, n):
        r"""
        Given an integer (an \class{Integer}, an \class{int}, or a
        \class{long}), produce a \class{SageInputExpression} that displays
        the integer with no marking for what kind of integer it is
        (so it may read back as an \class{Integer}, an \class{int}, or
        a \class{long}, depending on its size and whether the preparser
        is enabled).

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.int(-3^50))
            -717897987691852588770249

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.int(long(2^65)))
            36893488147419103232

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.int(-42r))
            -42
        """
        if n < 0:
            return -SIE_literal_stringrep(self, -n)
        else:
            return SIE_literal_stringrep(self, n)

    def float_str(self, n):
        r"""
        Given a string representing a floating-point number,
        produces a \class{SageInputExpression} that formats as that
        string.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.float_str(repr(RR(e))))
            2.71828182845905
        """
        return SIE_literal_stringrep(self, n)

    def name(self, n):
        r"""
        Given a string representing a Python name,
        produces a \class{SageInputExpression} for that name.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.name('pi') + sib.name('e'))
            pi + e
        """
        return SIE_literal_stringrep(self, n)

    def cache(self, x, sie, name):
        r"""
        INPUTS:
            x -- an arbitrary value
            sie -- a \class{SageInputExpression}
            name -- a requested variable name

        Enters \var{x} and \var{sie} in a cache, so that subsequent calls
        \code{self(x)} will directly return \var{sie}.  Also, marks the
        requested name of this \var{sie} to be \var{name}.

        This should almost always be called as part of the
        \method{_sage_input_} method of a parent.  It may also be called
        on values of an arbitrary type, which may be useful if the values
        are both large and likely to be used multiple times in a single
        expression.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie42 = sib(GF(101)(42))
            sage: sib.cache(GF(101)(42), sie42, 'the_ultimate_answer')
            sage: sib.result(sib(GF(101)(42)) + sib(GF(101)(42)))
            the_ultimate_answer = GF(101)(42)
            the_ultimate_answer + the_ultimate_answer

        Note that we don't assign the result to a variable if the value
        is only used once.
            sage: sib = SageInputBuilder()
            sage: sie42 = sib(GF(101)(42))
            sage: sib.cache(GF(101)(42), sie42, 'the_ultimate_answer')
            sage: sib.result(sib(GF(101)(42)) + sib(GF(101)(43)))
            GF_101 = GF(101)
            GF_101(42) + GF_101(43)
        """
        self._cached_types.add(type(x))
        self._cache[(parent(x), x)] = sie
        sie._sie_preferred_varname = name

    def empty_subscript(self, parent):
        r"""
        Given a \class{SageInputExpression} representing \code{foo},
        produces a \class{SageInputExpression} representing \code{foo[]}.
        Since this is not legal Python syntax, it is useful only for
        producing the \sage generator syntax for a polynomial ring.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.empty_subscript(sib(2) + sib(3)))
            (2 + 3)[]

        The following calls this method indirectly.
            sage: sage_input(polygen(ZZ['y']))
            R.<x> = ZZ['y'][]
            x
        """
        return SIE_subscript(self, parent, None)

    def parent_with_gens(self, parent, sie, gen_names, name, gens_syntax=None):
        r"""
        This method is used for parents with generators, to manage the
        \sage preparser generator syntax (like \code{K.<x> = QQ[]}).

        The \method{_sage_input_} method of a parent class with
        generators should construct a \class{SageInputExpression} for
        the parent, and then call this method with the parent itself,
        the constructed SIE, a sequence containing the names of the
        generators, and (optionally) another SIE to use if the \sage
        generator syntax is used; typically this will be the same as
        the first SIE except omitting a \var{names} parameter.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder


            sage: def test_setup(use_gens=True, preparse=True):
            ...       sib = SageInputBuilder(preparse=preparse)
            ...       gen_names=('foo', 'bar')
            ...       parent = "some parent"
            ...       normal_sie = sib.name('make_a_parent')(names=gen_names)
            ...       if use_gens:
            ...           gens_sie = sib.name('make_a_parent')()
            ...       else:
            ...           gens_sie = None
            ...       name = 'the_thing'
            ...       result = sib.parent_with_gens(parent, normal_sie,
            ...                                     gen_names, name,
            ...                                     gens_syntax=gens_sie)
            ...       return sib, result

            sage: sib, par_sie = test_setup()
            sage: sib.result(par_sie)
            make_a_parent(names=('foo', 'bar'))

            sage: sib, par_sie = test_setup()
            sage: sib.result(sib(3) * sib.gen("some parent", 0))
            the_thing.<foo,bar> = make_a_parent()
            3*foo

            sage: sib, par_sie = test_setup(preparse=False)
            sage: sib.result(par_sie)
            make_a_parent(names=('foo', 'bar'))

            sage: sib, par_sie = test_setup(preparse=False)
            sage: sib.result(sib(3) * sib.gen("some parent", 0))
            the_thing = make_a_parent(names=('foo', 'bar'))
            foo,bar = the_thing.gens()
            ZZ(3)*foo

            sage: sib, par_sie = test_setup(use_gens=False)
            sage: sib.result(par_sie)
            make_a_parent(names=('foo', 'bar'))

            sage: sib, par_sie = test_setup(use_gens=False)
            sage: sib.result(sib(3) * sib.gen("some parent", 0))
            the_thing = make_a_parent(names=('foo', 'bar'))
            foo,bar = the_thing.gens()
            3*foo

            sage: sib, par_sie = test_setup()
            sage: sib.result(par_sie - sib.gen("some parent", 1))
            the_thing.<foo,bar> = make_a_parent()
            the_thing - bar
        """
        v = SIE_gens_constructor(self, sie, gen_names, gens_syntax=gens_syntax)
        self.cache(parent, v, name)
        gens = [SIE_gen(self, v, n) for n in gen_names]
        self._parent_gens[parent] = gens
        v._sie_gens = gens
        return v

    def gen(self, parent, n=0):
        r"""
        Given a parent, returns a \class{SageInputExpression} for
        the $n$th (default 0) generator of the parent.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.gen(ZZ['y']))
            R.<y> = ZZ[]
            y
        """
        if not parent in self._parent_gens:
            self(parent)
            if not parent in self._parent_gens:
                raise ValueError, "%s did not register generators for sage_input" % parent

        gens = self._parent_gens[parent]

        if n > len(gens):
            raise ValueError, "%s registered only %d generators for sage_input" % (parent, len(gens))

        return gens[n]

    def prod(self, factors, simplify=False):
        r"""
        Given a sequence, returns a \class{SageInputExpression}
        for the product of the elements.

        With \code{simplify=True}, performs some simplifications
        first.  If any element is formatted as a string \code{'0'},
        then that element is returned directly.  If any element is
        formatted as a string \code{'1'}, then it is removed
        from the sequence (unless it is the only element in the sequence).
        And any negations are removed from the elements and moved to the
        outside of the product.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 0, 1, -2]))
            -1*0*1*-2

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 0, 1, 2], simplify=True))
            0

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 2, -3, -4], simplify=True))
            -2*3*4

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([-1, 1, -1, -1], simplify=True))
            -1

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.prod([1, 1, 1], simplify=True))
            1
        """
        neg = False
        factors = [self(factor) for factor in factors]
        if simplify:
            i = 0
            while i < len(factors):
                factor = factors[i]
                while isinstance(factor, SIE_unary) and factor._sie_op == '-':
                    neg = not neg
                    factor = factor._sie_operand
                    factors[i] = factor
                if isinstance(factor, SIE_literal_stringrep) and factor._sie_value == '0':
                    factors = [factor]
                    neg = False
                    break
                if isinstance(factor, SIE_literal_stringrep) and factor._sie_value == '1':
                    factors[i:i+1] = []
                else:
                    i += 1
            if len(factors) == 0:
                factors.append(SIE_literal_stringrep(self, '1'))

        prod = factors[0]
        for factor in factors[1:]:
            prod = prod * factor
        if neg:
            prod = -prod
        return prod

    def sum(self, terms, simplify=False):
        r"""
        Given a sequence, returns a \class{SageInputExpression}
        for the product of the elements.

        With \code{simplify=True}, performs some simplifications
        first.  If any element is formatted as a string \code{'0'},
        then it is removed from the sequence (unless it is the only
        element in the sequence); and any instances of \code{a + -b}
        are changed to \code{a - b}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.sum([-1, 0, 1, 0, -1]))
            -1 + 0 + 1 + 0 + -1

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.sum([-1, 0, 1, 0, -1], simplify=True))
            -1 + 1 - 1

            sage: sib = SageInputBuilder()
            sage: sib.result(sib.sum([0, 0, 0], simplify=True))
            0
        """
        terms = [self(term) for term in terms]
        if simplify:
            i = 0
            while i < len(terms):
                term = terms[i]
                if isinstance(term, SIE_literal_stringrep) and term._sie_value == '0':
                    terms[i:i+1] = []
                else:
                    i += 1
            if len(terms) == 0:
                terms.append(SIE_literal_stringrep(self, '0'))

        sum = terms[0]
        for term in terms[1:]:
            if simplify and isinstance(term, SIE_unary) and term._sie_op == '-':
                sum = sum - term._sie_operand
            else:
                sum = sum + term
        return sum

    def result(self, e):
        r"""
        Given a \class{SageInputExpression} constructed using \code{self},
        returns a tuple of a list of commands and an expression
        (and possibly a dictionary of local variables) suitable for
        \function{sage_eval}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: r = sib.result(sib(6) * sib(7)); r
            6*7
            sage: tuple(r)
            ('', '6*7')
        """
        sif = SageInputFormatter()

        e._sie_prepare(sif)

        s = sif.format(e, 0)

        locals = self._locals
        if len(locals):
            return SageInputAnswer(sif._commands, sif.format(e, 0), locals)
        else:
            return SageInputAnswer(sif._commands, sif.format(e, 0))

# Python's precedence levels.  Hand-transcribed from section 5.14 of
# the Python reference manual.
_prec_lambda = 2
_prec_or = 4
_prec_and = 6
_prec_not = 8
_prec_membership = 10
_prec_identity = 12
_prec_comparison = 14
_prec_bitor = 16
_prec_bitxor = 18
_prec_bitand = 20
_prec_shift = 22
_prec_addsub = 24
_prec_muldiv = 26
_prec_negate = 28
_prec_not = 30
_prec_exponent = 32
_prec_attribute = 34
_prec_subscript = 36
_prec_slicing = 38
_prec_funcall = 40
_prec_atomic = 42

class SageInputExpression(object):
    r"""
    Subclasses of this class represent expressions for \function{sage_input}.
    \sage classes should define a \method{_sage_input_} method, which
    will return an instance of \class{SageInputExpression}, created using
    methods of \class{SageInputBuilder}.

    To the extent possible, operations on \class{SageInputExpression} objects
    construct a new \class{SageInputExpression} representing that operation.
    That is, if \var{a} is a \class{SageInputExpression}, then \code{a + b}
    constructs a \class{SageInputExpression} representing this sum.
    This also works for attribute access, function calls, subscripts, etc.
    Since arbitrary attribute accesses might be used to construct a new
    attribte-access expression, all internal attributes and methods
    have names that begin with \code{_sie_} to reduce the chance of
    collisions.

    It is expected that instances of this class will not be directly
    created outside this module; instead, instances will be created
    using methods of \class{SageInputBuilder} and \class{SageInputExpression}.

    Values of type \class{SageInputExpression} print in a fairly ugly
    way, that reveals the internal structure of the expression tree.
    """

    def __init__(self, sib):
        r"""
        Initialize a \class{SageInputExpression}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3) # indirect doctest
            sage: sie
            {atomic:3}
            sage: sie._sie_builder is sib
            True
        """
        self._sie_refcount = 0
        self._sie_builder = sib
        self._sie_context = None
        self._sie_preferred_varname = None
        self._sie_varname = None
        self._sie_use_var = False
        self._sie_requested_varname = False

    def _sie_is_simple(self):
        r"""
        Returns \code{True} if this \class{SageInputExpression} is simple
        enough that duplicate uses are not worth caching.  Normally
        this will be true if the expression represents a single token.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.name('QQ')._sie_is_simple()
            True
            sage: sib(GF(2))._sie_is_simple()
            False
        """
        return False

    def _sie_referenced(self):
        r"""
        Returns a list of the immediate subexpressions of this
        \class{SageInputExpression}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: len(sib(GF(2))._sie_referenced())
            2
            sage: sib(5)._sie_referenced()
            []
        """
        return []

    def _sie_prepare(self, sif):
        r"""
        We traverse the entire expression DAG to prepare for printing.
        Here, we notice nodes with more than one parent, and mark them
        to replace with a variable (rather than generating the value
        multiple times).

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: pair = sib((GF(2), GF(2)))
            sage: single = sib(GF(2))
            sage: single._sie_refcount
            0
            sage: single._sie_use_var
            False
            sage: sib((GF(2), GF(2)))._sie_prepare(sif)
            sage: single._sie_refcount
            2
            sage: single._sie_use_var
            True
        """
        if self._sie_context is not sif:
            self._sie_context = sif
            self._sie_refcount = 0
        self._sie_refcount += 1
        if not self._sie_is_simple():
            if self._sie_refcount == 2:
                self._sie_require_varname(sif)
                self._sie_use_var = True
        if self._sie_refcount == 1:
            for r in self._sie_referenced():
                r._sie_prepare(sif)

    def _sie_require_varname(self, sif):
        r"""
        Mark this \class{SageInputExpression} as requiring a variable name,
        and register it with a \class{SageInputFormatter} (which will
        allocate a variable name at the end of the preparatory phase).

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(3)
            sage: sie._sie_require_varname(sif)
            sage: sie._sie_requested_varname
            True
        """
        if not self._sie_requested_varname:
            sif.register_name(self._sie_preferred_varname)
            self._sie_requested_varname = True
            self._sie_generated = False

    def _sie_get_varname(self, sif):
        r"""
        Get the variable name that the \class{SageInputFormatter} allocated
        for this \class{SageInputExpression}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(3)
            sage: sie._sie_require_varname(sif)
            sage: sie._sie_get_varname(sif)
            'si'
        """
        if self._sie_varname is None:
            self._sie_varname = sif.get_name(self._sie_preferred_varname)

        return self._sie_varname

    def __call__(self, *args, **kwargs):
        r"""
        Given a \class{SageInputExpression}, build a new
        \class{SageInputExpression} representing a function call node
        (with \code{self} as the function).

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie(4)
            {call: {atomic:3}({atomic:4})}
        """
        args = map(self._sie_builder, args)
        for k in kwargs:
            kwargs[k] = self._sie_builder(kwargs[k])
        return SIE_call(self._sie_builder, self, args, kwargs)

    def __getitem__(self, key):
        r"""
        Given a \class{SageInputExpression}, build a new
        \class{SageInputExpression} representing a subscript expression
        (with \code{self} as the value being subscripted).

        Currently, slices are not supported.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie[4]
            {subscr: {atomic:3}[{atomic:4}]}
            sage: sie[sib.name('x'), sib.name('y')]
            {subscr: {atomic:3}[{tuple: ({atomic:x}, {atomic:y})}]}
        """
        skey = self._sie_builder(key)
        return SIE_subscript(self._sie_builder, self, skey)

    def __pow__(self, other):
        r"""
        Compute an expression tree for \code{self ** other}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie ^ 4
            {binop:** {atomic:3} {atomic:4}}
        """
        return self._sie_binop('**', other)

    def __mul__(self, other):
        r"""
        Compute an expression tree for \code{self * other}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie * 4
            {binop:* {atomic:3} {atomic:4}}
        """
        return self._sie_binop('*', other)

    def __div__(self, other):
        r"""
        Compute an expression tree for \code{self / other}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie / 4
            {binop:/ {atomic:3} {atomic:4}}
        """
        return self._sie_binop('/', other)

    def __add__(self, other):
        r"""
        Compute an expression tree for \code{self + other}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie + 4
            {binop:+ {atomic:3} {atomic:4}}
        """
        return self._sie_binop('+', other)

    def __sub__(self, other):
        r"""
        Compute an expression tree for \code{self - other}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie - 4
            {binop:- {atomic:3} {atomic:4}}
        """
        return self._sie_binop('-', other)

    def _sie_binop(self, op, other):
        r"""
        Compute an expression tree for \code{self OP other},
        where OP is a string representing a binary operator (such as
        '+' or '**').

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: v = sib.name('x')._sie_binop('%', sib.name('y'))
            sage: type(v)
            <class 'sage.misc.sage_input.SIE_binary'>
            sage: (v)._sie_op
            '%'
            sage: v
            {binop:% {atomic:x} {atomic:y}}
        """
        return SIE_binary(self._sie_builder, op, self, self._sie_builder(other))

    def __neg__(self):
        r"""
        Compute an expression tree for \code{-self}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: -sie
            {unop:- {atomic:3}}
        """
        return self._sie_unop('-')

    def _sie_unop(self, op):
        r"""
        Compute an expression tree for \code{OP self},
        where OP is a string representing a unary operator (such as
        '-' or '~').

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: v = sie._sie_unop('~')
            sage: type(v)
            <class 'sage.misc.sage_input.SIE_unary'>
            sage: (v)._sie_op
            '~'
            sage: v
            {unop:~ {atomic:3}}
        """
        return SIE_unary(self._sie_builder, op, self)

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and the
        precedence of the top-level operator in the expression.

        EXAMPLES:

        Actually, all of these are examples of the \method{_sie_format}
        method on subclasses of \class{SageInputExpression};
        \class{SageInputExpression} itself is an abstract base class
        (that cannot be instantiated).

            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(3)

            sage: for v in (sie, sie+7, sie/5):
            ...       v._sie_prepare(sif)
            ...       v._sie_format(sif)
            ('3', 42)
            ('3 + 7', 24)
            ('3/5', 26)
        """
        raise NotImplementedError

class SIE_literal(SageInputExpression):
    r"""
    An abstract base class for ``literals'' (basically, values which
    consist of a single token).

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder, SIE_literal

        sage: sib = SageInputBuilder()
        sage: sie = sib(3)
        sage: sie
        {atomic:3}
        sage: isinstance(sie, SIE_literal)
        True
    """

    def _sie_is_simple(self):
        r"""
        Report that \class{SIE_literal} values are not worth replacing by
        variables (for ``common subexpression elimination'') even if they
        occur multiple times in an expression.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib(3)
            sage: sie._sie_is_simple()
            True
        """
        # Perhaps this should actually look at the formatted length of self,
        # and sometimes return false?  If some 50-digit integer occurs multiple
        # times in an expression, it might be better to do the replacement.
        return True

class SIE_literal_stringrep(SIE_literal):
    r"""
    Values in this class are leaves in a \function{sage_input} expression
    tree.  Typically they represent a single token, and consist of the
    string representation of that token.  They are used for integer,
    floating-point, and string literals, and for name expressions.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder, SIE_literal_stringrep

        sage: sib = SageInputBuilder()
        sage: isinstance(sib(3), SIE_literal_stringrep)
        True
        sage: isinstance(sib(3.14159, True), SIE_literal_stringrep)
        True
        sage: isinstance(sib.name('pi'), SIE_literal_stringrep)
        True
        sage: isinstance(sib(False), SIE_literal_stringrep)
        True
        sage: sib(False)
        {atomic:False}
    """

    def __init__(self, sib, n):
        r"""
        Initialize a \class{SIE_literal_stringrep} value.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            n -- a string; the value to be printed for this expression

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(3)
            {atomic:3}
            sage: sib(3)._sie_value
            '3'
        """
        super(SIE_literal_stringrep, self).__init__(sib)
        self._sie_value = str(n)

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_literal_stringrep}
        value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(3)
            {atomic:3}
            sage: sib("\n")
            {atomic:'\n'}
        """
        return "{atomic:%s}" % self._sie_value

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and an indication
        that it is ``atomic'' (never needs to be parenthesized).

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(True)
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ('True', 42)
        """
        return self._sie_value, _prec_atomic

class SIE_call(SageInputExpression):
    r"""
    This class represents a function-call node in a \function{sage_input}
    expression tree.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sie = sib.name('GF')
        sage: sie(49)
        {call: {atomic:GF}({atomic:49})}
    """

    def __init__(self, sib, func, args, kwargs):
        r"""
        Initialize an instance of \class{SIE_call}.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            func -- a \class{SageInputExpression} representing a function
            args -- a list of \class{SageInputExpression}s representing the
                    positional arguments
            kwargs -- a dictionary mapping strings to
                    \class{SageInputExpression}s representing the keyword
                    arguments

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib('RealField')(53, rnd='RNDZ')
        """
        super(SIE_call, self).__init__(sib)
        self._sie_func = func
        self._sie_args = args
        self._sie_kwargs = kwargs

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_call} value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib('RealField')(53, rnd='RNDZ')
        """
        func = repr(self._sie_func)
        args = [repr(arg) for arg in self._sie_args]
        kwargs = sorted(k + '=' + repr(v) for k, v in self._sie_kwargs.iteritems())
        all_args = ', '.join(args + kwargs)
        return "{call: %s(%s)}" % (func, all_args)

    def _sie_referenced(self):
        r"""
        Returns a list of the immediate subexpressions of this \class{SIE_call}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib('RealField')(53, rnd='RNDZ')
            sage: sie._sie_referenced()
            [{atomic:53}, {atomic:'RealField'}, {atomic:'RNDZ'}]
        """
        refs = self._sie_args[:]
        refs.append(self._sie_func)
        refs.extend(self._sie_kwargs.itervalues())
        return refs

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and an indication
        that it is a function call.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.name('RealField')(53, rnd='RNDZ')
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ("RealField(53, rnd='RNDZ')", 40)
        """
        func = sif.format(self._sie_func, _prec_attribute)
        args = [sif.format(arg, 0) for arg in self._sie_args]
        kwargs = sorted(k + '=' + sif.format(v, 0) for k, v in self._sie_kwargs.iteritems())
        all_args = ', '.join(args + kwargs)
        return ('%s(%s)' % (func, all_args), _prec_funcall)

class SIE_subscript(SageInputExpression):
    r"""
    This class represents a subscript node in a \function{sage_input}
    expression tree.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sie = sib.name('QQ')['x,y']
        sage: sie
        {subscr: {atomic:QQ}[{atomic:'x,y'}]}
    """

    def __init__(self, sib, coll, key):
        r"""
        Initialize an instance of \class{SIE_subscript}.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            func -- a \class{SageInputExpression} representing a collection
            key -- a \class{SageInputExpression} representing the subscript/key

        As a special case, \var{key} may be \code{None}; this represents an
        empty subscript.  This is not legal Python syntax, but it is legal
        in the \sage preparser in examples like \code{K.<x> = QQ[]}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.name('QQ')['x']
            {subscr: {atomic:QQ}[{atomic:'x'}]}
            sage: sib.name('x')[1,2,3]
            {subscr: {atomic:x}[{tuple: ({atomic:1}, {atomic:2}, {atomic:3})}]}
            sage: sib.empty_subscript(sib.name('QQ'))
            {subscr: {atomic:QQ}[]}
        """
        super(SIE_subscript, self).__init__(sib)
        self._sie_coll = coll
        self._sie_key = key

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_subscript} value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.name('ZZ')['x,y']
            {subscr: {atomic:ZZ}[{atomic:'x,y'}]}
        """
        coll = repr(self._sie_coll)
        if self._sie_key is None:
            key = ''
        else:
            key = repr(self._sie_key)
        return "{subscr: %s[%s]}" % (coll, key)

    def _sie_referenced(self):
        r"""
        Returns a list of the immediate subexpressions of this
        \class{SIE_subscript}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib.name('GF')(5)['x,y']
            sage: sie._sie_referenced()
            [{call: {atomic:GF}({atomic:5})}, {atomic:'x,y'}]
        """
        refs = [self._sie_coll]
        if self._sie_key is not None:
            refs.append(self._sie_key)
        return refs

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression, and an
        indication that it is a subscript.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.name('QQ')['x']
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ("QQ['x']", 36)
        """
        coll = sif.format(self._sie_coll, _prec_attribute)
        if self._sie_key is None:
            key = ''
        else:
            key = sif.format(self._sie_key, 0)
        return '%s[%s]' % (coll, key), _prec_subscript

class SIE_tuple(SageInputExpression):
    r"""
    This class represents a tuple node in a \function{sage_input}
    expression tree.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib((1, 'howdy'))
        {tuple: ({atomic:1}, {atomic:'howdy'})}
    """

    def __init__(self, sib, values):
        r"""
        Initialize an instance of \class{SIE_tuple}.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            values -- a list of \class{SageInputExpression}s representing
                 the elements of this tuple

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib((3.5, -2))
            {tuple: ({call: {atomic:RR}({atomic:3.5000000000000000})}, {unop:- {atomic:2}})}
        """
        super(SIE_tuple, self).__init__(sib)
        self._sie_values = values

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_tuple} value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib((2,3,5))
            {tuple: ({atomic:2}, {atomic:3}, {atomic:5})}
        """
        return "{tuple: (%s)}" % ', '.join([repr(v) for v in self._sie_values])

    def _sie_referenced(self):
        r"""
        Returns a list of the immediate subexpressions of this
        \class{SIE_tuple}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib((ZZ, GF(5)))
            sage: sie._sie_referenced()
            [{atomic:ZZ}, {call: {atomic:GF}({atomic:5})}]
        """
        return self._sie_values

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this tuple, and an
        indication that it is atomic (never needs to be parenthesized).

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: for v in ((), (1,), (1,2)):
            ...        sie = sib(v)
            ...        sie._sie_prepare(sif)
            ...        sie._sie_format(sif)
            ('()', 42)
            ('(1,)', 42)
            ('(1, 2)', 42)
        """
        values = [sif.format(val, 0) for val in self._sie_values]
        if len(values) == 1:
            return '(%s,)' % values[0], _prec_atomic
        else:
            return '(%s)' % ', '.join(values), _prec_atomic

class SIE_binary(SageInputExpression):
    r"""
    This class represents an arithmetic expression with a binary operator
    and its two arguments, in a \function{sage_input} expression tree.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib(3)+5
        {binop:+ {atomic:3} {atomic:5}}
    """

    def __init__(self, sib, op, lhs, rhs):
        r"""
        Initialize an instance of \class{SIE_binary}.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            op -- a string representing a binary operator, such as '*' or '%'
            lhs -- a \class{SageInputExpression}
            rhs -- a \class{SageInputExpression}

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(3)*5
            {binop:* {atomic:3} {atomic:5}}

        """
        super(SIE_binary, self).__init__(sib)
        self._sie_op = op
        self._sie_operands = (lhs, rhs)

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_binary} value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib(7)/9
            {binop:/ {atomic:7} {atomic:9}}
        """
        return "{binop:%s %s %s}" % (self._sie_op, repr(self._sie_operands[0]), repr(self._sie_operands[1]))

    def _sie_referenced(self):
        r"""
        Returns a tuple of the immediate subexpressions of this
        \class{SIE_binary}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = sib.name('x') + 5
            sage: sie._sie_referenced()
            ({atomic:x}, {atomic:5})
        """
        return self._sie_operands

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression,
        and the precedence of the top-level operator in the expression.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: x = sib.name('x')
            sage: y = sib.name('y')
            sage: for v in (x+y, x*y, x**y):
            ...       v._sie_prepare(sif)
            ...       v._sie_format(sif)
            ('x + y', 24)
            ('x*y', 26)
            ('x^y', 32)

        Note that the printing for $x^y$ varies depending on whether the
        preparser is enabled.
            sage: sibnp = SageInputBuilder(preparse=False)
            sage: sif = SageInputFormatter()
            sage: v = x**y
            sage: v._sie_prepare(sif)
            sage: v._sie_format(sif)
            ('x^y', 32)

        TESTS:
            sage: x = sib.name('x')
            sage: y = sib.name('y')
            sage: z = sib.name('z')
            sage: sib.result((x+y)+z)
            x + y + z
            sage: sib.result(x+(y+z))
            x + (y + z)
            sage: sib.result((x*y)*z)
            x*y*z
            sage: sib.result(x*(y*z))
            x*(y*z)
            sage: sib.result(x+(y*z))
            x + y*z
            sage: sib.result((x+y)*z)
            (x + y)*z
            sage: sib.result((x^y)^z)
            (x^y)^z
            sage: sib.result(x^(y^z))
            x^y^z
        """
        op = self._sie_op
        fop = op
        if op == '**':
            lhs = sif.format(self._sie_operands[0], _prec_exponent+1)
            rhs = sif.format(self._sie_operands[1], _prec_exponent)
            if self._sie_builder.preparse():
                return '%s^%s' % (lhs, rhs), _prec_exponent
            else:
                return '%s**%s' % (lhs, rhs), _prec_exponent

        if op == '*':
            prec = _prec_muldiv
        elif op == '/':
            prec = _prec_muldiv
        elif op == '+':
            fop = ' + '
            prec = _prec_addsub
        elif op == '-':
            fop = ' - '
            prec = _prec_addsub
        else:
            raise ValueError, 'Unhandled op %s in SIE_binary' % op

        lhs = sif.format(self._sie_operands[0], prec)
        rhs = sif.format(self._sie_operands[1], prec+1)
        return '%s%s%s' % (lhs, fop, rhs), prec

class SIE_unary(SageInputExpression):
    r"""
    This class represents an arithmetic expression with a unary operator
    and its argument, in a \function{sage_input} expression tree.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: -sib(256)
        {unop:- {atomic:256}}
    """

    def __init__(self, sib, op, operand):
        r"""
        Initialize an instance of \class{SIE_unary}.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            op -- a string representing a unary operator, such as '-'
            operand -- a \class{SageInputExpression}

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: -sib(3)
            {unop:- {atomic:3}}
        """
        super(SIE_unary, self).__init__(sib)
        self._sie_op = op
        self._sie_operand = operand

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_unary} value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: -sib(15)
            {unop:- {atomic:15}}
        """
        return "{unop:%s %s}" % (self._sie_op, repr(self._sie_operand))

    def _sie_referenced(self):
        r"""
        Returns a list of the immediate subexpressions of this
        \class{SIE_unary}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sie = -sib.name('x')
            sage: sie._sie_referenced()
            [{atomic:x}]
        """

        return [self._sie_operand]

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this expression,
        and the precedence of the top-level operator in the expression.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: x = sib.name('x')
            sage: v = -x
            sage: v._sie_prepare(sif)
            sage: v._sie_format(sif)
            ('-x', 28)

        TESTS:
            sage: x = sib.name('x')
            sage: y = sib.name('y')
            sage: sib.result((-x)+y)
            -x + y
            sage: sib.result(x+(-y))
            x + -y
            sage: sib.result(-(x+y))
            -(x + y)
            sage: sib.result(-(-x))
            --x
            sage: sib.result(x-(-y))
            x - -y

        We assume that -(x*y) is always equal to (-x)*y.  Using this
        assumption, we print -(x*y) as -x*y, which parses as (-x)*y.
            sage: sib.result(-(x*y))
            -x*y
            sage: sib.result((-x)*y)
            -x*y
            sage: sib.result(x*(-y))
            x*-y
        """
        op = self._sie_op
        fop = op
        rprec = None
        if op == '-':
            # We print -(a*b) as -a*b, even though that will parse as
            # (-a)*b.
            prec = _prec_muldiv
            rprec = _prec_negate
        else:
            raise ValueError, 'Unhandled op %s in SIE_unary' % op

        if rprec is None: rprec = prec

        return '%s%s' % (fop, sif.format(self._sie_operand, prec)), rprec

class SIE_gens_constructor(SageInputExpression):
    r"""
    This class represents an expression that can create a \sage parent
    with named generators, optionally using the \sage preparser
    generators syntax (like \code{K.<x> = QQ[]}).

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: qq = sib.name('QQ')
        sage: sib.parent_with_gens("some parent", qq['x'],
        ...                        ('x',), 'QQx',
        ...                        gens_syntax=sib.empty_subscript(qq))
        {constr_parent: {subscr: {atomic:QQ}[{atomic:'x'}]} with gens: ('x',)}
    """

    def __init__(self, sib, constr, gen_names, gens_syntax=None):
        r"""
        Initialize an instance of \class{SIE_gens_constructor}.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            constr -- a \class{SageInputExpression} for constructing
                this parent ``normally''
            gen_names -- a tuple of generator names
            gens_syntax -- an optional \class{SageInputExpression}
                for constructing this parent using the \sage preparser
                generators syntax

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: qq = sib.name('QQ')
            sage: sib.parent_with_gens("some parent", qq['x'],
            ...                        ('x',), 'QQx',
            ...                        gens_syntax=sib.empty_subscript(qq))
            {constr_parent: {subscr: {atomic:QQ}[{atomic:'x'}]} with gens: ('x',)}
        """
        super(SIE_gens_constructor, self).__init__(sib)
        self._sie_constr = constr
        self._sie_gen_names = gen_names
        self._sie_gens = None # will be overwritten from .parent_with_gens()
        self._sie_gens_constr = gens_syntax
        self._sie_assign_gens = False
        self._sie_generated = False

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_gens_constructor} value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: qq = sib.name('QQ')
            sage: sib.parent_with_gens("some parent", qq['x'],
            ...                        ('x',), 'QQx',
            ...                        gens_syntax=sib.empty_subscript(qq))
            {constr_parent: {subscr: {atomic:QQ}[{atomic:'x'}]} with gens: ('x',)}
        """
        return "{constr_parent: %s with gens: %s}" % (repr(self._sie_constr), self._sie_gen_names)

    def _sie_referenced(self):
        r"""
        Returns a list of the immediate subexpressions of this
        \class{SIE_gens_constructor}.


        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ...                             ('x',), 'QQx',
            ...                             gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_referenced()
            [{subscr: {atomic:QQ}[{atomic:'x'}]}]
        """
        # This is used to determine if some expressions should be replaced
        # by variables (if the expression has more than one parent in
        # the expression DAG).  We assume that all expressions in
        # self._sie_gens_constr also occur in self._sie_constr.
        return [self._sie_constr]

    def _sie_gens_referenced(self, sif):
        r"""
        Mark that at least one of the generators in this
        \class{SIE_gens_constructor} is used.  (This means we will actually
        construct all of the generators.)

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ...                             ('x',), 'QQx',
            ...                             gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_assign_gens
            False
            sage: gc._sie_gens_referenced(sif)
            sage: gc._sie_assign_gens
            True
        """
        self._sie_assign_gens = True
        self._sie_require_varname(sif)
        for gen in self._sie_gens:
            gen._sie_require_varname(sif)

    def _sie_add_command(self, sif):
        r"""
        Build commands to construct this parent and (if necessary)
        its associated generators.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ...                             ('x',), 'QQx',
            ...                             gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_gens_referenced(sif)
            sage: gc._sie_prepare(sif)
            sage: gc._sie_add_command(sif)
            sage: sif._commands
            'QQx.<x> = QQ[]\n'

        TESTS:
        There are several tricky cases here.

        We prefer the \sage preparser generators syntax:
            sage: sage_input(polygen(ZZ))
            R.<x> = ZZ[]
            x

        But of course we can't use that without the preparser:
            sage: sage_input(polygen(ZZ), preparse=False)
            R = ZZ['x']
            x = R.gen()
            x

        We also can't use the preparser syntax if there is a conflict
        between generator names.  For example, this works:
            sage: sage_input((polygen(ZZ), polygen(GF(17), 'y')))
            R1.<x> = ZZ[]
            R2.<y> = GF(17)[]
            (x, y)

        but this can't use the preparser syntax.
            sage: sage_input((polygen(ZZ), polygen(GF(17))))
            R1 = ZZ['x']
            x1 = R1.gen()
            R2 = GF(17)['x']
            x2 = R2.gen()
            (x1, x2)

        If we never use the generators, then we don't bother with the
        preparser syntax.
            sage: sage_input((ZZ['x'], ZZ['x'], GF(17)['y']))
            R = ZZ['x']
            (R, R, GF(17)['y'])
        """
        if not self._sie_generated:
            if self._sie_builder.preparse() and \
                    self._sie_gens_constr is not None and \
                    all(g._sie_got_preferred(sif) for g in self._sie_gens):
                s, _ = self._sie_gens_constr._sie_format(sif)
                sif._commands += '%s.<%s> = %s\n' % (self._sie_get_varname(sif), ','.join(self._sie_gen_names), s)
            else:
                s, _ = self._sie_constr._sie_format(sif)
                sif._commands += '%s = %s\n' % (self._sie_get_varname(sif), s)
                if self._sie_assign_gens:
                    if len(self._sie_gens) == 1:
                        sif._commands += '%s = %s.gen()\n' % (self._sie_gens[0]._sie_get_varname(sif), self._sie_get_varname(sif))
                    else:
                        sif._commands += '%s = %s.gens()\n' % (','.join([g._sie_get_varname(sif) for g in self._sie_gens]), self._sie_get_varname(sif))
            self._sie_generated = True

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this parent-construction
        expression, and its precedence.

        As a side effect, if the generators of this parent are used,
        this adds commands to assign the generators to names.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: qq = sib.name('QQ')
            sage: gc = sib.parent_with_gens("some parent", qq['x'],
            ...                             ('x',), 'QQx',
            ...                             gens_syntax=sib.empty_subscript(qq))
            sage: gc._sie_gens_referenced(sif)
            sage: gc._sie_prepare(sif)
            sage: gc._sie_format(sif)
            ('QQx', 42)
            sage: sif._commands
            'QQx.<x> = QQ[]\n'
        """
        if self._sie_assign_gens:
            self._sie_add_command(sif)
            return self._sie_get_varname(sif), _prec_atomic

        return self._sie_constr._sie_format(sif)

class SIE_gen(SageInputExpression):
    r"""
    This class represents a named generator of a parent with named
    generators.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputBuilder

        sage: sib = SageInputBuilder()
        sage: sib.gen(ZZ['x'])
        {gen:x {constr_parent: {subscr: {atomic:ZZ}[{atomic:'x'}]} with gens: ('x',)}}
    """

    def __init__(self, sib, parent, name):
        r"""
        Initializes an instance of \class{SIE_gen}.

        INPUTS:
            sib -- a \class{SageInputBuilder}
            parent -- a \class{SIE_gens_constructor}
            name -- a string with the name of this generator

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.gen(ZZ['x']) # indirect doctest
            {gen:x {constr_parent: {subscr: {atomic:ZZ}[{atomic:'x'}]} with gens: ('x',)}}
        """
        super(SIE_gen, self).__init__(sib)
        self._sie_parent = parent
        self._sie_preferred_varname = name

    def __repr__(self):
        r"""
        Returns a string representing this \class{SIE_gen} value.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.gen(ZZ['x']) # indirect doctest
            {gen:x {constr_parent: {subscr: {atomic:ZZ}[{atomic:'x'}]} with gens: ('x',)}}
        """
        return "{gen:%s %s}" % (self._sie_preferred_varname, repr(self._sie_parent))

    def _sie_is_simple(self):
        r"""
        Report that \class{SIE_gen} values are single tokens.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder

            sage: sib = SageInputBuilder()
            sage: sib.gen(ZZ['x'])._sie_is_simple()
            True
        """
        return True

    def _sie_prepare(self, sif):
        r"""
        We override the \method{_sie_prepare} method from
        \class{SageInputExpression} to additionally mark the parent of this
        generator that the generator names must be assigned.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.gen(GF(13)['z'])
            sage: sie._sie_parent._sie_assign_gens
            False
            sage: sie._sie_prepare(sif)
            sage: sie._sie_parent._sie_assign_gens
            True
        """
        super(SIE_gen, self)._sie_prepare(sif)
        self._sie_parent._sie_gens_referenced(sif)

    def _sie_format(self, sif):
        r"""
        Return the formatted string value of this named generator,
        and an indication that it is atomic.

        As a side effect, this generates commands to assign the generators
        of the parent to variables.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib.gen(GF(41)['x'])
            sage: sie._sie_prepare(sif)
            sage: sie._sie_format(sif)
            ('x', 42)
            sage: sif._commands
            'R.<x> = GF(41)[]\n'
        """
        self._sie_parent._sie_add_command(sif)
        return self._sie_get_varname(sif), _prec_atomic

    def _sie_got_preferred(self, sif):
        r"""
        Check whether the \class{SageInputFormatter} assigned us a
        variable name which is the same as the name of the generator
        name.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

        First we verify that if we use two generators with different
        names, then they get their preferred names.
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: v = sib.gen(GF(2)['x']); w = sib.gen(GF(3)['y'])
            sage: v._sie_prepare(sif); w._sie_prepare(sif)
            sage: v._sie_got_preferred(sif)
            True
            sage: w._sie_got_preferred(sif)
            True

        Now, we repeat the experiment, except that the generators now
        have the same names.  In this case, the \class{SageInputFormatter}
        will not use the generator name as the variable name, because
        of this conflict.
            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: v = sib.gen(GF(2)['x']); w = sib.gen(GF(3)['x'])
            sage: v._sie_prepare(sif); w._sie_prepare(sif)
            sage: v._sie_got_preferred(sif)
            False
            sage: w._sie_got_preferred(sif)
            False
        """
        return self._sie_get_varname(sif) == self._sie_preferred_varname

class SageInputFormatter:
    r"""
    An instance of this class is used to keep track of variable names
    and a sequence of generated commands during the \function{sage_input}
    formatting process.
    """

    def __init__(self):
        r"""
        Initialize an instance of \class{SageInputFormatter}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputFormatter

            sage: sif = SageInputFormatter()
        """
        self._commands = ''
        self._names = set()
        self._dup_names = {}

    def format(self, e, prec):
        r"""
        Format a \class{SageInputExpression} into a string.

        INPUTS:
            e -- a \class{SageInputExpression}
            prec -- an integer representing a precedence level

        First, we check to see if \var{e} should be replaced by a variable.
        If so, we generate the command to assign the variable, and return
        the name of the variable.

        Otherwise, we format the expression by calling its \method{_sie_format}
        method, and add parentheses if necessary.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputBuilder, SageInputFormatter

            sage: sib = SageInputBuilder()
            sage: sif = SageInputFormatter()
            sage: sie = sib(GF(5))

        Here we ``cheat'' by calling \method{_sie_prepare} twice, to make it
        use a variable.
            sage: sie._sie_prepare(sif)
            sage: sie._sie_prepare(sif)
            sage: sif._commands
            ''
            sage: sif.format(sie, 0)
            'GF_5'
            sage: sif._commands
            'GF_5 = GF(5)\n'
        """
        if e._sie_use_var:
            if not e._sie_generated:
                s, _ = e._sie_format(self)
                self._commands += '%s = %s\n' % (e._sie_get_varname(self), s)
                e._sie_generated = True

            return e._sie_get_varname(self)

        s, iprec = e._sie_format(self)
        if iprec < prec:
            s = '(' + s + ')'
        return s

    def register_name(self, name):
        r"""
        Register that some value would like to use a given name.
        If only one request for a name is received, then we will use the
        requested name; otherwise, we will add numbers to the end of the
        name to make it unique.

        If the input name is \code{None}, then it is treated as a name of
        \code{'si'}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputFormatter

            sage: sif = SageInputFormatter()
            sage: sif._names, sif._dup_names
            (set([]), {})
            sage: sif.register_name('x')
            sage: sif.register_name('y')
            sage: sif._names, sif._dup_names
            (set(['y', 'x']), {})
            sage: sif.register_name('x')
            sage: sif._names, sif._dup_names
            (set(['y', 'x']), {'x': 0})
        """
        if name is None: name = 'si'

        if name in self._names:
            self._dup_names[name] = 0
        else:
            self._names.add(name)

    def get_name(self, name):
        r"""
        Return a name corresponding to a given requested name.
        If only one request for a name is received, then we will use the
        requested name; otherwise, we will add numbers to the end of the
        name to make it unique.

        If the input name is \code{None}, then it is treated as a name of
        \code{'si'}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputFormatter

            sage: sif = SageInputFormatter()
            sage: names = ('x', 'x', 'y', 'z')
            sage: for n in names: sif.register_name(n)
            sage: for n in names: sif.get_name(n)
            'x1'
            'x2'
            'y'
            'z'
        """
        if name is None: name = 'si'

        if name in self._dup_names:
            next = self._dup_names[name] + 1
            self._dup_names[name] = next
            return name + str(next)
        else:
            return name

def verify_same(a, b):
    r"""
    Verify that two \sage values are the same.  This is an extended equality
    test; it checks that the values are equal and that their parents are equal.
    (For values which are not Elements, the types are checked instead.)

    If the values are the same, we return \code{None}; otherwise,
    we raise an exception.

    EXAMPLES:
        sage: from sage.misc.sage_input import verify_same
        sage: verify_same(1, 1)
        sage: verify_same(1, 2)
        Traceback (most recent call last):
        ...
            assert(a == b)
        AssertionError
        sage: verify_same(1, 1r)
        Traceback (most recent call last):
        ...
        AttributeError: 'int' object has no attribute 'parent'
        sage: verify_same(1r, 1)
        Traceback (most recent call last):
        ...
            assert(type(a) == type(b))
        AssertionError
        sage: verify_same(5, GF(7)(5))
        Traceback (most recent call last):
        ...
            assert(a.parent() == b.parent())
        AssertionError
    """
    from sage.structure.element import is_Element
    if is_Element(a):
        assert(a.parent() == b.parent())
    else:
        assert(type(a) == type(b))
    assert(a == b)

def verify_si_answer(x, answer, preparse):
    r"""
    INPUTS:
        x -- an arbitrary \sage value
        answer -- a string, or a \class{SageInputAnswer}
        preparse -- \code{True}, \code{False}, or \code{None}

    Verify that evaluating \var{answer} gives a value equal to \var{x}
    (with the same parent/type).  If \var{preparse} is \code{True} or
    \code{False}, then we evaluate \var{answer} with the preparser
    enabled or disabled, respectively; if \var{preparse} is \code{None},
    then we evaluate \var{answer} both with the preparser enabled and
    disabled and check both results.

    On success, we return \code{None}; on failure, we raise an exception.

    EXAMPLES:
        sage: from sage.misc.sage_input import verify_si_answer
        sage: verify_si_answer(1, '1', True)
        sage: verify_si_answer(1, '1', False)
        Traceback (most recent call last):
        ...
        AttributeError: 'int' object has no attribute 'parent'
        sage: verify_si_answer(1, 'ZZ(1)', None)
    """
    from sage.misc.sage_eval import sage_eval
    if preparse is None:
        verify_same(x, sage_eval(answer, preparse=True))
        verify_same(x, sage_eval(answer, preparse=False))
    else:
        verify_same(x, sage_eval(answer, preparse=preparse))

class SageInputAnswer(tuple):
    r"""
    This class inherits from tuple, so it acts like a tuple when passed
    to \function{sage_eval}; but it prints as a sequence of commands.

    EXAMPLES:
        sage: from sage.misc.sage_input import SageInputAnswer
        sage: v = SageInputAnswer('x = 22\n', 'x/7'); v
        x = 22
        x/7
        sage: isinstance(v, tuple)
        True
        sage: v[0]
        'x = 22\n'
        sage: v[1]
        'x/7'
        sage: len(v)
        2
        sage: v = SageInputAnswer('', 'sin(3.14)', {'sin': math.sin}); v
        LOCALS:
          sin: <built-in function sin>
        sin(3.14)
        sage: v[0]
        ''
        sage: v[1]
        'sin(3.14)'
        sage: v[2]
        {'sin': <built-in function sin>}
    """

    def __new__(cls, cmds, expr, locals=None):
        r"""
        Construct an instance of \class{SageInputAnswer}.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputAnswer
            sage: v = SageInputAnswer('', 'sin(3.14)', {'sin': math.sin}); v
            LOCALS:
              sin: <built-in function sin>
            sin(3.14)
            sage: v[0]
            ''
            sage: v[1]
            'sin(3.14)'
            sage: v[2]
            {'sin': <built-in function sin>}
        """
        if locals:
            return tuple.__new__(cls, (cmds, expr, locals))
        else:
            return tuple.__new__(cls, (cmds, expr))

    def __repr__(self):
        r"""
        Return a string representation for a \class{SageInputAnswer},
        such that if you evaluate this \class{SageInputAnswer} at the
        \sage command line, you get a result in a nice form ready to
        copy-and-paste.

        EXAMPLES:
            sage: from sage.misc.sage_input import SageInputAnswer
            sage: v = SageInputAnswer('', 'sin(3.14)', {'sin': math.sin}); v
            LOCALS:
              sin: <built-in function sin>
            sin(3.14)
        """
        if len(self) == 2:
            return self[0] + self[1]

        locals = self[2]
        locals_text = ''.join('  %s: %r\n' % (k, v) for k, v in locals.iteritems())
        return 'LOCALS:\n' + locals_text + self[0] + self[1]

