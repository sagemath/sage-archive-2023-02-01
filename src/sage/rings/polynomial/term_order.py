r"""
Term Orderings

\SAGE supports the following term orderings:

\begin{description}
\item[Lexicographic (\emph{lex})]

$x^a < x^b \Leftrightarrow \exists\; 0 \le i < n : a_0 = b_0, \ldots, a_{i-1} = b_{i-1}, a_i < b_i$

EXAMPLES:
    sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='lex')
    sage: x > y
    True
    sage: x > y^2
    True
    sage: x > 1
    True
    sage: x^1*y^2 > y^3*z^4
    True
    sage: x^3*y^2*z^4 < x^3*y^2*z^1
    False

This term ordering is called 'lp' in \Singular.

\item[Degree reverse lexicographic (\emph{degrevlex})]

Let $deg(x^a) = a_0 + \cdots + a_{n-1},$ then
$x^a < x^b \Leftrightarrow deg(x^a) < deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 \le i < n: a_{n-1} = b_{n-1}, \ldots, a_{i+1} = b_{i+1}, a_i > b_i.$

EXAMPLES:
    sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='degrevlex')
    sage: x > y
    True
    sage: x > y^2*z
    False
    sage: x > 1
    True
    sage: x^1*y^5*z^2 > x^4*y^1*z^3
    True
    sage: x^2*y*z^2 > x*y^3*z
    False

This term ordering is called 'dp' in \Singular.

\item[Degree lexicographic (\emph{deglex})]

Let $deg(x^a) = a_0 + \cdots + a_{n-1},$ then
$x^a < x^b \Leftrightarrow deg(x^a) < deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 \le i < n:a_0 = b_0, \ldots, a_{i-1} = b_{i-1}, a_i < b_i.$

EXAMPLES:
    sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='deglex')
    sage: x > y
    True
    sage: x > y^2*z
    False
    sage: x > 1
    True
    sage: x^1*y^2*z^3 > x^3*y^2*z^0
    True
    sage: x^2*y*z^2 > x*y^3*z
    True

This term order is called 'Dp' in \Singular.

\item[Inverse lexicographic (\emph{invlex})]

$x^a < x^b \Leftrightarrow \exists\; 0 \le i < n : a_{n-1} = b_{n-1}, \ldots, a_{i+1} = b_{i+1}, a_i < b_i.$

EXAMPLES:
    sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='invlex')
    sage: x > y
    False
    sage: y > x^2
    True
    sage: x > 1
    True
    sage: x*y > z
    False

This term ordering only makes sense in a non-commutative setting
because if P is the ring $k[x_1, \dots, x_n]$ and term ordering
'invlex' then it is equivalent to the ring $k[x_n, \dots, x_1]$ with
term ordering 'lex'.

This ordering is called 'rp' in \Singular.

\item[Negative lexicographic (\emph{neglex})]

$x^a < x^b \Leftrightarrow \exists\; 0 \le i < n : a_0 = b_0, \ldots, a_{i-1} = b_{i-1}, a_i > b_i$

EXAMPLES:
    sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='neglex')
    sage: x > y
    False
    sage: x > 1
    False
    sage: x^1*y^2 > y^3*z^4
    False
    sage: x^3*y^2*z^4 < x^3*y^2*z^1
    True
    sage: x*y > z
    False

This term ordering is called 'ls' in \Singular.

\item[Negative degree reverse lexicographic (\emph{negdegrevlex})]

Let $deg(x^a) = a_0 + \cdots + a_{n-1},$ then
$x^a < x^b \Leftrightarrow deg(x^a) > deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 \le i < n: a_{n-1} = b_{n-1}, \ldots, a_{i+1} = b_{i+1}, a_i > b_i.$

EXAMPLES:
    sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='negdegrevlex')
    sage: x > y
    True
    sage: x > x^2
    True
    sage: x > 1
    False
    sage: x^1*y^2 > y^3*z^4
    True
    sage: x^2*y*z^2 > x*y^3*z
    False

This term ordering is called 'ds' in \Singular.

\item[Negative degree lexicographic (\emph{negdeglex})]

Let $deg(x^a) = a_0 + \cdots + a_{n-1},$ then
$x^a < x^b \Leftrightarrow deg(x^a) > deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 \le i < n: a_0 = b_0, \ldots, a_{i-1} = b_{i-1}, a_i < b_i.$

EXAMPLES:
    sage: P.<x,y,z> = PolynomialRing(QQ, 3, order='negdeglex')
    sage: x > y
    True
    sage: x > x^2
    True
    sage: x > 1
    False
    sage: x^1*y^2 > y^3*z^4
    True
    sage: x^2*y*z^2 > x*y^3*z
    True

This term ordering is called 'Ds' in \Singular.

\end{description}

Of these, only 'degrevlex', 'deglex', 'invlex' and 'lex' are global orderings.

Additionally all these monomial orderings may be combined to product
or block orderings, defined as:

Let $x = (x_0, \ldots, x_{n-1})$ and $y = (y_0, \ldots, y_{m-1})$ be
two ordered sets of variables, $<_1$ a monomial ordering on $k[x]$ and
$<_2$ a monomial ordering on $k[y]$.

The product ordering (or block ordering) $<\ := (<_1,<_2)$ on $k[x,y]$
is defined as: $x^a y^b < x^A y^B \Leftrightarrow x^a <_1 x^A$ or
$(x^a =x^A \textrm{ and } y^b <_2 y^B)$.

These block orderings are constructed in \SAGE by giving a comma
separated list of monomial orderings with the length of each block
attached to them.

EXAMPLE:

   As an example, consider constructing a block ordering where the
   first four variables are compared using the degree reverse
   lexicographical ordering while the last two variables in the second
   block are compared using negative lexicographical ordering.

   sage: P.<a,b,c,d,e,f> = PolynomialRing(QQ, 6,order='degrevlex(4),neglex(2)')
   sage: a > c^4
   False
   sage: a > e^4
   True
   sage: e > f^2
   False

   The same result can be achieved by:

   sage: T1 = TermOrder('degrevlex',4)
   sage: T2 = TermOrder('neglex',2)
   sage: T = T1 + T2
   sage: P.<a,b,c,d,e,f> = PolynomialRing(QQ, 6, order=T)
   sage: a > c^4
   False
   sage: a > e^4
   True

If any other unsupported term ordering is given the provided string is
passed through as is to \Singular, \textsc{Macaulay2}, and
\Magma. This ensures that it is for example possible to calculated a
Groebner basis with respect to some term ordering \Singular supports
but \SAGE doesn't. However a warning is issued to make the user aware
of the situation and potential typos:

   sage: P.<a,b,c,d,e,f> = PolynomialRing(QQ, 6, order=T)
   sage: a > c^4
   False
   sage: a > e^4
   True

If any other unsupported term ordering is given the provided string is
passed through as is to \Singular, \textsc{Macaulay2}, and
\Magma. This ensures that it is for example possible to calculated a
Groebner basis with respect to some term ordering \Singular supports
but \SAGE doesn't. However a warning is issued to make the user aware
of the situation and potential typos:

   sage: T = TermOrder("royalorder")
   verbose 0 (...: term_order.py, __init__) Term ordering 'royalorder' unknown.

AUTHORS:
    -- David Joyner and William Stein: initial version multi_polynomial_ring
    -- Kiran S. Kedlaya: added macaulay2 interface
    -- Martin Albrecht: implemented native term orderings, refactoring
"""

import re
from sage.structure.sage_object import SageObject

print_name_mapping =     {'lex'          :'Lexicographic',
                          'invlex'       :'Inverse Lexicographic',
                          'degrevlex'    :'Degree reverse lexicographic',
                          'deglex'       :'Degree lexicographic',
                          'neglex'       :'Negative lexicographic',
                          'negdegrevlex' :'Negative degree reverse lexicographic',
                          'negdeglex'    :'Negative degree lexicographic'}

singular_name_mapping =  {'lex'          :'lp',
                          'invlex'       :'rp',
                          'degrevlex'    :'dp',
                          'deglex'       :'Dp',
                          'neglex'       :'ls',
                          'negdegrevlex' :'ds',
                          'negdeglex'    :'Ds'}

macaulay2_name_mapping = {'lex'          :'Lex',
                          'revlex'       :'RevLex, Global=>false',
                          'degrevlex'    :'GRevLex',
                          'deglex'       :'GLex'}

inv_macaulay2_name_mapping = { 'Lex':'lex',
                               'RevLex, Global=>false':'revlex',
                               'GRevLex':"degrevlex",
                               'GLex':'deglex'}

magma_name_mapping =     {'lex'          :'"lex"',
                          'degrevlex'    :'"grevlex"',
                          'deglex'       :'"glex"'}


inv_singular_name_mapping ={'lp':'lex'          ,
                            'rp':'invlex'       ,
                            'dp':'degrevlex'    ,
                            'Dp':'deglex'       ,
                            'ls':'neglex'       ,
                            'ds':'negdegrevlex' ,
                            'Ds':'negdeglex'    }


lp_description = """
Lexicographic (lex) term ordering.

$x^a < x^b <=> \exists\; 0 <= i < n : a_0 = b_0, ..., a_{i-1} = b_{i-1}, a_i < b_i$
"""

dp_description = """
Degree reverse lexicographic (degrevlex) term ordering.

Let $deg(x^a) = a_0 + ... + a_{n-1},$ then $x^a < x^b <=> deg(x^a) < deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 <= i < n: a_{n-1} = b_{n-1}, ..., a_{i+1} = b_{i+1}, a_i > b_i.$
"""

Dp_description = """
Degree lexicographic (deglex) term ordering.

Let $deg(x^a) = a_0 + \cdots + a_{n-1},$ then $x^a < x^b \Leftrightarrow deg(x^a) < deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 \le i < n:a_0 = b_0, \ldots, a_{i-1} = b_{i-1}, a_i < b_i.$
"""

rp_description = """
Inverse lexicographic (invlex) term ordering.

$x^a < x^b \Leftrightarrow \exists\; 0 \le i < n : a_{n-1} = b_{n-1}, \ldots, a_{i+1} = b_{i+1}, a_i < b_i.$
"""

ls_description = """
Negative lexicographic (neglex) term ordering.

$x^a < x^b \Leftrightarrow \exists\; 0 \le i < n : a_0 = b_0, \ldots, a_{i-1} = b_{i-1}, a_i > b_i$
"""

ds_description = """
Negative degree reverse lexicographic (negdegrevlex) term ordering.

Let $deg(x^a) = a_0 + \cdots + a_{n-1},$ then $x^a < x^b \Leftrightarrow deg(x^a) > deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 \le i < n: a_{n-1} = b_{n-1}, \ldots, a_{i+1} = b_{i+1}, a_i > b_i.$
"""

Ds_description = """
Negative degree lexicographic (negdeglex) term ordering.

Let $deg(x^a) = a_0 + \cdots + a_{n-1},$ then $x^a < x^b \Leftrightarrow deg(x^a) > deg(x^b)$ or
$deg(x^a) = deg(x^b)$ and $\exists\ 0 \le i < n: a_0 = b_0, \ldots, a_{i-1} = b_{i-1}, a_i < b_i.$
"""

description_mapping = {
    "lp":lp_description,
    "dp":dp_description,
    "Dp":Dp_description,
    "rp":rp_description,
    "ls":ls_description,
    "Ds":Ds_description,
    "ds":ds_description}

class TermOrder(SageObject):
    r"""
    A term order.

    See \code{sage.rings.polynomial.term_order} for details on
    supported term orderings.
    """
    def __init__(self, name='lex', n = 0, blocks=True):
        """
        Construct a new term ordering object.

        INPUT:
            name -- name of the term ordering (default: lex)
            n -- number of variables in the polynomial ring (default: 0)
            blocks -- controls whether a list of blocks is maintained
                      (internal use only, default:True)

        See the \code{sage.rings.polynomial.term_order} module for
        help which names and orderings are available.

        EXAMPLES:

            sage: t = TermOrder('lex')
            sage: t
            Lexicographic term order
            sage: loads(dumps(t)) == t
            True

        We can construct block orderings directly as

            sage: TermOrder('degrevlex(3),neglex(2)')
            degrevlex(3),neglex(2) term order

        or by adding together the blocks:

            sage: t1 = TermOrder('degrevlex',3)
            sage: t2 = TermOrder('neglex',2)
            sage: t1 + t2
            degrevlex(3),neglex(2) term order
            sage: t2 + t1
            neglex(2),degrevlex(3) term order

        NOTE: The optional $n$ parameter is not necessary if only
        non-block orderings like $deglex$ are constructed. However, it
        is useful if block orderings are to be constructed from this
        \code{TermOrder} object later.

        """
        if isinstance(name, TermOrder):
            if n == 0 and len(name) > 0:
                n = len(name)
            name = name.name()
        name = name.lower()

        #Block Orderings
        if "," in name:
            pattern  = re.compile("\(([0-9]+)\)$")

            self.__length = 0
            self.__name = ""
            self.__singular_str = "("
            self.__macaulay2_str = "{"
            self.__magma_str = "" # I (malb) believe MAGMA doesn't support this

            length = 0
            blocks = []
            name_str = []
            singular_str = []
            macaulay2_str = []

            for block in name.split(","):
                try:
                    block_name, block_length, _ = re.split(pattern,block)
                except ValueError:
                    raise TypeError, "%s is not a valid term ordering"%(name,)

                block_length = int(block_length)

                blocks.append( TermOrder(block_name, block_length) )
                name_str.append("%s(%d)"%(block_name, block_length))
                singular_str.append("%s(%d)"%(singular_name_mapping.get(block_name, block_name), block_length))
                macaulay2_str.append("%s => %d"%(macaulay2_name_mapping.get(block_name, block_name), block_length))

                length += block_length

            self.blocks = tuple(blocks)
            self.__length = length
            self.__name = ",".join(name_str)
            self.__singular_str = "(" + ",".join(singular_str) + ")"
            self.__macaulay2_str = "(" + ",".join(macaulay2_str) + ")"
            self.__magma_str = ""

            if n!=0 and self.__length != n:
                raise TypeError, "Term order length does not match number of generators"

        else:
            from sage.misc.misc import verbose
            if blocks:
                if name not in singular_name_mapping.keys() and \
                        name not in singular_name_mapping.values():
                    verbose("Term ordering '%s' unknown."%name,level=0)
                self.blocks = (TermOrder(name,n,blocks=False),)
            else:
                self.blocks = tuple()
            self.__length = n
            self.__name = name
            self.__singular_str = singular_name_mapping.get(name,name)
            self.__macaulay2_str = macaulay2_name_mapping.get(name,name)
            self.__magma_str = magma_name_mapping.get(name,name)

            if self.__singular_str in description_mapping:
                self.__doc__ = description_mapping[self.__singular_str]

    def __getattr__(self,name):
        """
        Return the correct compare_tuples/greater_tuple function.

        EXAMPLE:
            sage: TermOrder('lex').compare_tuples
            <bound method TermOrder.compare_tuples_lp of Lexicographic term order>

            sage: TermOrder('deglex').compare_tuples
            <bound method TermOrder.compare_tuples_Dp of Degree lexicographic term order>
        """
        if name=='compare_tuples':
            if len(self.blocks) <= 1:
                return getattr(self,'compare_tuples_'+self.__singular_str)
            else:
                return self.compare_tuples_block
        elif name=='greater_tuple':
            if len(self.blocks) <= 1:
                return getattr(self,'greater_tuple_'+self.__singular_str)
            else:
                return self.greater_tuple_block
        else:
            raise AttributeError,name

    def compare_tuples_lp(self,f,g):
        """
        Compares two exponent tuples with respect to the
        lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='lex')
            sage: x > y^2 # indirect doctest
            True
            sage: x > 1
            True
        """

        if f>g:
            return 1
        elif f<g:
            return -1
        else:
            return 0

    def compare_tuples_rp(self,f,g):
        """
        Compares two exponent tuples with respect to the inversed
        lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='invlex')
            sage: x > y^2 # indirect doctest
            False
            sage: x > 1
            True
        """
        return self.compare_tuples_lp(f.reversed(),g.reversed())

    def compare_tuples_Dp(self,f,g):
        """
        Compares two exponent tuples with respect to the
        degree lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='deglex')
            sage: x > y^2 # indirect doctest
            False
            sage: x > 1
            True
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        if sf > sg:
            return 1
        elif sf<sg:
            return -1
        elif sf == sg:
            return self.compare_tuples_lp(f,g)

    def compare_tuples_dp(self,f,g):
        """
        Compares two exponent tuples with respect to the degree
        reversed lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='degrevlex')
            sage: x > y^2 # indirect doctest
            False
            sage: x > 1
            True
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        if sf > sg:
            return 1
        elif sf<sg:
            return -1
        elif sf == sg:
            return (-1)*self.compare_tuples_lp(f.reversed(),g.reversed())


    def compare_tuples_ls(self,f,g):
        """
        Compares two exponent tuples with respect to the
        negative lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='neglex')
            sage: x > y^2 # indirect doctest
            False
            sage: x > 1
            False
        """
        return -self.compare_tuples_lp(f,g)

    def compare_tuples_ds(self,f,g):
        """
        Compares two exponent tuples with respect to the
        negative degree reverse lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='negdegrevlex')
            sage: x > y^2 # indirect doctest
            True
            sage: x > 1
            False
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        if sf > sg:
            return -1
        elif sf<sg:
            return 1
        elif sf == sg:
            return (-1)*self.compare_tuples_lp(f.reversed(),g.reversed())

    def compare_tuples_Ds(self,f,g):
        """
        Compares two exponent tuples with respect to the
        negative degree lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y> = PolynomialRing(QQbar, 2, order='negdeglex')
            sage: x > y^2 # indirect doctest
            True
            sage: x > 1
            False
        """
        sf = sum(f.nonzero_values(sort=False))
        sg = sum(g.nonzero_values(sort=False))
        if sf > sg:
            return -1
        elif sf<sg:
            return 1
        elif sf == sg:
            return self.compare_tuples_lp(f,g)

    def compare_tuples_block(self, f,g):
        """
        Compares two exponent tuple with respec to the block ordering
        as specified when constructing this element.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<a,b,c,d,e,f>=PolynomialRing(QQbar, 6, order='degrevlex(3),degrevlex(3)')
            sage: a > c^4 # indirect doctest
            False
            sage: a > e^4
            True
        """
        n = 0
        for order in self:
            r = getattr(self,"compare_tuples_"+order.singular_str())(f[n:n+len(order)],g[n:n+len(order)])
            if r != 0:
                return r
            n += len(order)
        return 0

    def greater_tuple_lp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the
        lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLES:
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='lex')
            sage: f = x + y^2; f.lm() # indirect doctest
            x

        This method is called by the lm/lc/lt methods of
        \code{MPolynomial_polydict}.
        """
        return f > g and f or g

    def greater_tuple_rp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the
        inversed lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='invlex')
            sage: f = x + y; f.lm() # indirect doctest
            y
            sage: f = y + x^2; f.lm()
            y

        This method is called by the lm/lc/lt methods of
        \code{MPolynomial_polydict}.
        """
        return f.reversed() > g.reversed()   and f or g

    def greater_tuple_Dp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the total
        degree lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='deglex')
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + y^2*z; f.lm()
            y^2*z

        This method is called by the lm/lc/lt methods of
        \code{MPolynomial_polydict}.
        """
        return (sum(f.nonzero_values(sort=False))>sum(g.nonzero_values(sort=False))
                or (sum(f.nonzero_values(sort=False))==sum(g.nonzero_values(sort=False)) and f  > g )) and f or g

    def greater_tuple_dp(self,f,g):
        """
        Returns the greater exponent tuple with respect to the total
        degree reversed lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLES:
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='degrevlex')
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + y^2*z; f.lm()
            y^2*z

        This method is called by the lm/lc/lt methods of
        \code{MPolynomial_polydict}.

        """
        return (sum(f.nonzero_values(sort=False))>sum(g.nonzero_values(sort=False))
                or (sum(f.nonzero_values(sort=False))==sum(g.nonzero_values(sort=False)) and f.reversed() < g.reversed())) and f or g

    def greater_tuple_ds(self,f,g):
        """
        Returns the greater exponent tuple with respect to the
        negative degree reverse lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='negdegrevlex')
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + x^2; f.lm()
            x
            sage: f = x^2*y*z^2 + x*y^3*z; f.lm()
            x*y^3*z

        This method is called by the lm/lc/lt methods of
        \code{MPolynomial_polydict}.
        """
        if self.compare_tuples_ds(f,g) >= 0:
            return f
        else:
            return g

    def greater_tuple_Ds(self,f,g):
        """
        Returns the greater exponent tuple with respect to the
        negative degree lexicographical term order.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<x,y,z> = PolynomialRing(QQbar, 3, order='negdeglex')
            sage: f = x + y; f.lm() # indirect doctest
            x
            sage: f = x + x^2; f.lm()
            x
            sage: f = x^2*y*z^2 + x*y^3*z; f.lm()
            x^2*y*z^2

        This method is called by the lm/lc/lt methods of
        \code{MPolynomial_polydict}.
        """
        if self.compare_tuples_Ds(f,g) >= 0:
            return f
        else:
            return g

    def greater_tuple_ls(self,f,g):
        """
        Returns the greater exponent tuple with respect to the
        negative lexicographical term order.

        This method is called by the lm/lc/lt methods of
        \code{MPolynomial_polydict}.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple


        EXAMPLE:
            sage: P.<a,b,c,d,e,f>=PolynomialRing(QQbar, 6, order='degrevlex(3),degrevlex(3)')
            sage: f = a + c^4; f.lm() # indirect doctest
            c^4
            sage: g = a + e^4; g.lm()
            a
        """
        if self.greater_tuple_lp(f,g) is f:
            return g
        else:
            return f

    def greater_tuple_block(self, f,g):
        """
        Compares two exponent tuple with respec to the block ordering
        as specified when constructing this element.

        This method is called by the lm/lc/lt methods of \code{MPolynomial_polydict}.

        INPUT:
            f -- exponent tuple
            g -- exponent tuple

        EXAMPLE:
            sage: P.<a,b,c,d,e,f>=PolynomialRing(QQbar, 6, order='degrevlex(3),degrevlex(3)')
            sage: f = a + c^4; f.lm() # indirect doctest
            c^4
            sage: g = a + e^4; g.lm()
            a
        """
        n = 0
        for order in self:
            r = getattr(self,"compare_tuples_"+order.singular_str())(f[n:n+len(order)],g[n:n+len(order)])
            if r != 0:
                if r < 0:
                    return g
                else:
                    return f
            n += len(order)
        return f

    def name(self):
        """
        EXAMPLE:
            sage: TermOrder('lex').name()
            'lex'
        """
        return self.__name

    def _repr_(self):
        """
        EXAMPLE:
            sage: TermOrder('lex') # indirect doctest
            Lexicographic term order
        """
        s = print_name_mapping.get(self.__name,self.__name)
        return '%s term order'%s

    def singular_str(self):
        """
        Return a SINGULAR representation of self.

        Used to convert polynomial rings to their SINGULAR
        representation.

        EXAMPLE:
            sage: P = PolynomialRing(GF(127),10,names='x',order='lex(3),deglex(5),lex(2)')
            sage: T = P.term_order()
            sage: T.singular_str()
            '(lp(3),Dp(5),lp(2))'
            sage: P._singular_()
            //   characteristic : 127
            //   number of vars : 10
            //        block   1 : ordering lp
            //                  : names    x0 x1 x2
            //        block   2 : ordering Dp
            //                  : names    x3 x4 x5 x6 x7
            //        block   3 : ordering lp
            //                  : names    x8 x9
            //        block   4 : ordering C

        """
        return self.__singular_str

    def macaulay2_str(self):
        """
        Return a Macaulay2 representation of self.

        Used to convert polynomial rings to their Macaulay2
        representation.

        EXAMPLE:
            sage: P = PolynomialRing(GF(127), 8,names='x',order='degrevlex(3),lex(5)')
            sage: T = P.term_order()
            sage: T.macaulay2_str()
            '(GRevLex => 3,Lex => 5)'
            sage: P._macaulay2_() # optional - macaulay2
            ZZ/127 [x0, x1, x2, x3, x4, x5, x6, x7, MonomialOrder => {GRevLex => 3, Lex => 5}, MonomialSize => 16]
        """
        return self.__macaulay2_str

    def magma_str(self):
        """
        Return a MAGMA representation of self.

        Used to convert polynomial rings to their MAGMA
        representation.

        EXAMPLE:
            sage: P = PolynomialRing(GF(127), 10,names='x',order='degrevlex')
            sage: magma(P)                                                        # optional - magma
            Polynomial ring of rank 10 over GF(127)
            Graded Reverse Lexicographical Order
            Variables: x0, x1, x2, x3, x4, x5, x6, x7, x8, x9

            sage: T = P.term_order()
            sage: T.magma_str()
            '"grevlex"'
        """
        return self.__magma_str

    def __cmp__(self, other):
        """
        Only equality testing makes sense here.

        EXAMPLE:
            sage: TermOrder('lex') == TermOrder('lex',3)
            True

            sage: TermOrder('degrevlex') == TermOrder('lex')
            False

            sage: T1 = TermOrder('lex',2)+TermOrder('lex',3)
            sage: T2 = TermOrder('lex',3)+TermOrder('lex',2)
            sage: T1 == T2
            False

            sage: T1 = TermOrder('lex',2)+TermOrder('neglex',3)
            sage: T2 = TermOrder('lex',2)+TermOrder('neglex',3)
            sage: T1 == T2
            True
        """
        if not isinstance(other, TermOrder):
            if isinstance(other, str):
                other = TermOrder(other)
            else:
                return cmp(type(self), type(other))
        return cmp(self.singular_str(), other.singular_str())


    def __add__(self, other):
        """
        Block ordering constructor.

        INPUT:
            other -- a term order

        OUTPUT:
            a block ordering

        EXAMPLE:
            sage: from sage.rings.polynomial.term_order import TermOrder
            sage: TermOrder('deglex',2) + TermOrder('degrevlex(3),neglex(3)')
            deglex(2),degrevlex(3),neglex(3) term order
        """
        if other is None:
            return self

        if not isinstance(other,TermOrder):
            other = TermOrder(other)
        if len(self) == 0 or len(other) == 0:
            raise ArithmeticError, "Can only concatenate term orders with length attribute."

        name = []
        for block_order in self.blocks + other.blocks:
            nom = block_order.singular_str()
            name.append("%s(%d)"%(inv_singular_name_mapping.get(nom,nom), len(block_order)))

        name = ",".join(name)
        return TermOrder(name, len(self)+len(other))

    def __len__(self):
        """
        Return the length of this term ordering, i.e. the number of
        variables it covers. This may be zero for undefinitely many
        variables.

        EXAMPLE:
            sage: T = TermOrder('lex')
            sage: len(T)
            0
            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: len(T)
            5
        """
        return self.__length

    def __getitem__(self, i):
        r"""
        Return the i-th block of this term ordering.

        INPUT:
            i -- index

        EXAMPLE:
            sage: T = TermOrder('lex')
            sage: T[0]
            Lexicographic term order

            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: T[1]
            Degree reverse lexicographic term order

        Note that \code{len(self)} does not count blocks but variables.

            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: T[len(T)-1]
            Traceback (most recent call last):
            ...
            IndexError: tuple index out of range
        """
        return self.blocks[i]

    def __iter__(self):
        r"""
        Iterate over the blocks of this term ordering.

        EXAMPLE:
            sage: T = TermOrder('lex')
            sage: list(T) # indirect doctest
            [Lexicographic term order]

            sage: T = TermOrder('lex', 2) + TermOrder('degrevlex', 3)
            sage: list(T)
            [Lexicographic term order, Degree reverse lexicographic term order]

        Note that \code{len(self)} and \code{len(list(self))} are not
        the same. The former counts the number of variables in
        \code{self} while the latter counts the number of blocks.
        """
        return iter(self.blocks)

    def is_global(self):
        r"""
        Return \code{True} if this term ordering is definitely
        global. Return \code{False} otherwise, which includes unknown
        term orderings.

        EXAMPLE:
            sage: T = TermOrder('lex')
            sage: T.is_global()
            True
            sage: T = TermOrder('degrevlex', 3) + TermOrder('degrevlex', 3)
            sage: T.is_global()
            True
            sage: T = TermOrder('degrevlex', 3) + TermOrder('negdegrevlex', 3)
            sage: T.is_global()
            False
        """
        if len(self.blocks) == 1:
            if self.singular_str() in ('lp','dp','Dp'):
                return True
            else:
                return False
        else:
            return all([t.is_global() for t in self.blocks])

    def is_local(self):
        r"""
        Return \code{True} if this term ordering is definitely
        local. Return \code{False} otherwise, which includes unknown
        term orderings.

        EXAMPLE:
            sage: T = TermOrder('lex')
            sage: T.is_local()
            False
            sage: T = TermOrder('negdeglex', 3) + TermOrder('negdegrevlex', 3)
            sage: T.is_local()
            True
            sage: T = TermOrder('degrevlex', 3) + TermOrder('negdegrevlex', 3)
            sage: T.is_local()
            False
        """
        if len(self.blocks) == 1:
            if self.singular_str() in ('rp', 'ls','ds','Ds'):
                return True
            else:
                return False
        else:
            return all([t.is_local() for t in self.blocks])
