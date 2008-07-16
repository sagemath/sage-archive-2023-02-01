r"""
A Sample Session using Sympy

In this first part, we do all of the examples in the Sympy
tutorial (\url{http://code.google.com/p/sympy/wiki/Tutorial}), but using
\sage instead of Sympy.

sage: a = Rational((1,2))
sage: a
1/2
sage: a*2
1
sage: Rational(2)^50 / Rational(10)^50
1/88817841970012523233890533447265625
sage: 1.0/2
0.500000000000000
sage: 1/2
1/2
sage: pi^2
pi^2
sage: float(pi)
3.1415926535897931
sage: RealField(200)(pi)
3.1415926535897932384626433832795028841971693993751058209749
sage: float(pi + exp(1))
5.85987448204883...
sage: oo != 2
True


sage: var('x y')
(x, y)
sage: x + y + x - y
2*x
sage: (x+y)^2
(y + x)^2
sage: ((x+y)^2).expand()
y^2 + 2*x*y + x^2
sage: ((x+y)^2).subs(x=1)
(y + 1)^2
sage: ((x+y)^2).subs(x=y)
4*y^2


sage: limit(sin(x)/x, x=0)
1
sage: limit(x, x=oo)
+Infinity
sage: limit((5^x + 3^x)^(1/x), x=oo)
5



sage: diff(sin(x), x)
cos(x)
sage: diff(sin(2*x), x)
2*cos(2*x)
sage: diff(tan(x), x)
sec(x)^2
sage: limit((tan(x+y) - tan(x))/y, y=0)
1/cos(x)^2
sage: diff(sin(2*x), x, 1)
2*cos(2*x)
sage: diff(sin(2*x), x, 2)
-4*sin(2*x)
sage: diff(sin(2*x), x, 3)
-8*cos(2*x)


sage: cos(x).taylor(x,0,10)
1 - x^2/2 + x^4/24 - x^6/720 + x^8/40320 - x^10/3628800
sage: (1/cos(x)).taylor(x,0,10)
1 + x^2/2 + 5*x^4/24 + 61*x^6/720 + 277*x^8/8064 + 50521*x^10/3628800


sage: matrix([[1,0], [0,1]])
[1 0]
[0 1]
sage: var('x y')
(x, y)
sage: A = matrix([[1,x], [y,1]])
sage: A
[1 x]
[y 1]
sage: A^2
[x*y + 1     2*x]
[    2*y x*y + 1]
sage: R.<x,y> = QQ[]
sage: A = matrix([[1,x], [y,1]])
sage: print A^10
[x^5*y^5 + 45*x^4*y^4 + 210*x^3*y^3 + 210*x^2*y^2 + 45*x*y + 1     10*x^5*y^4 + 120*x^4*y^3 + 252*x^3*y^2 + 120*x^2*y + 10*x]
[    10*x^4*y^5 + 120*x^3*y^4 + 252*x^2*y^3 + 120*x*y^2 + 10*y x^5*y^5 + 45*x^4*y^4 + 210*x^3*y^3 + 210*x^2*y^2 + 45*x*y + 1]
sage: var('x y')
(x, y)

And here are some actual tests of sympy:
sage: from sympy import Symbol, cos, sympify, pprint
sage: from sympy.abc import x

sage: e = sympify(1)/cos(x)**3; e
cos(x)**(-3)
sage: f = e.series(x, 0, 10); f
1 + (3/2)*x**2 + (11/8)*x**4 + (241/240)*x**6 + (8651/13440)*x**8 + O(x**10)

And the pretty-printer:
sage: pprint(e)
       1
    -------
       3
    cos (x)
sage: pprint(f)
           2       4        6         8
        3*x    11*x    241*x    8651*x
    1 + ---- + ----- + ------ + ------- + O(x**10)
         2       8      240      13440

And the functionality to convert from sympy format to \sage format:
sage: e._sage_()
1/cos(x)^3
sage: print e._sage_()
                                           1
                                        -------
                                           3
                                        cos (x)
sage: e._sage_().taylor(x._sage_(), 0, 8)
1 + 3*x^2/2 + 11*x^4/8 + 241*x^6/240 + 8651*x^8/13440
sage: f._sage_()
8651*x^8/13440 + 241*x^6/240 + 11*x^4/8 + 3*x^2/2 + 1



Mixing SymPy with \sage:
sage: import sympy
sage: sympy.sympify(var("y"))+sympy.Symbol("x")
x + y
sage: o = var("omega")
sage: s = sympy.Symbol("x")
sage: t1 = s + o
sage: t2 = o + s
sage: print type(t1)
<class 'sage.calculus.calculus.SymbolicArithmetic'>
sage: print type(t2)
<class 'sage.calculus.calculus.SymbolicArithmetic'>
sage: print t1, t2
                                       x + omega
                                       x + omega
sage: e=sympy.sin(var("y"))+sage.all.cos(Symbol("x"))
sage: print type(e)
<class 'sage.calculus.calculus.SymbolicArithmetic'>
sage: print e
                                    sin(y) + cos(x)
sage: e=e._sage_()
sage: print type(e)
<class 'sage.calculus.calculus.SymbolicArithmetic'>
sage: print e
                                sin(y) + cos(x)
sage: e = sage.all.cos(var("y")**3)**4+var("x")**2
sage: e = e._sympy_()
sage: print e
    x**2 + cos(y**3)**4

sage: a = sympy.Matrix([1, 2, 3])
sage: a[1]
2

sage: sympify(1.5)
1.5
sage: sympify(2)
2
sage: sympify(-2)
-2
"""
