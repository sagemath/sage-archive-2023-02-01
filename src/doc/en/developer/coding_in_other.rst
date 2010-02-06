.. _chapter-cython:

=========================
Coding in Other Languages
=========================

When writing code for Sage, use Python for the basic structure and
interface. For speed, efficiency, or convenience, you can implement
parts of the code using any of the following languages: Cython, C/C++,
Fortran 95, GAP, Common Lisp, Singular, and GP/PARI. You can also use
all C/C++ libraries included with Sage  [3]_. (And if you are okay
with your code depending on optional Sage packages, you can use
Octave, or even Magma, Mathematica, or Maple.)

The first section of this chapter discusses Cython, which is a
compiled language based on Python. Many components of Sage are written
in Cython. Later sections discuss the interfaces between Sage and
PARI, GAP, and Singular.


Cython
======

Cython is a compiled version of Python. It is based on Pyrex
(http://www.cosc.canterbury.ac.nz/greg.ewing/python/Pyrex/). To a
large degree, Cython has changed based on what Sage's developers
needed; Cython has been developed in concert with Sage. However, it is
an independent project now, which is used beyond the scope of Sage.

As such, it is a young, but developing language, with young, but
developing documentation. See its web page,
http://www.cython.org/, for the most up-to-date information.

Python is an interpreted language and has no declared data types for
variables. These features make it easy to write and debug, but Python
code can sometimes be slow. Cython code can look a lot like Python,
but it gets translated into C code (often very efficient C code) and
then compiled. Thus it offers a language which is familiar to Python
developers, but with the potential for much greater speed.

There are several ways to create and build Cython code in Sage.

#. In the Sage Notebook, begin any cell with ``%cython``. When you
   evaluate that cell,

   #. It is saved to a file.

   #. Cython is run on it with all the standard Sage libraries
      automatically linked if necessary.

   #. The resulting ``.so`` file is then loaded into your running
      instance of Sage.

   #. The functionality defined in that cell is now available for you
      to use in the notebook. Also, the output cell has a link to the C
      program that was compiled to create the ``.so`` file.

#. Create an ``.spyx`` file and attach or load it from the command
   line. This is similar to creating a ``%cython`` cell in the
   notebook but works completely from the command line (and not from
   the notebook).

#. Create a ``.pyx`` file and add it to the Sage library.

   #. First, add a listing for the Cython extension to the variable
      ``ext_modules`` in the file
      ``SAGE_ROOT/devel/sage/module_list.py``. See the
      ``distutils.extension.Extension`` class for more information on
      creating a new Cython extension.

   #. Then, if you created a new directory for your ``.pyx`` file, add
      the directory name to the ``packages`` list in the file
      ``SAGE_ROOT/devel/sage/setup.py``.  (See also the section on
      "Creating a new directory" in :ref:`chapter-python`.)

   #. Run ``sage -b`` to rebuild Sage.

   For example, the file
   ``SAGE_ROOT/devel/sage/sage/graphs/chrompoly.pyx`` has the lines

   ::

     Extension('sage.graphs.chrompoly',
               sources = ['sage/graphs/chrompoly.pyx']),

   in ``module_list.py``. In addition, ``sage.graphs`` is included in
   the ``packages`` list under the Distutils section of ``setup.py``
   since ``chrompoly.pyx`` is contained in the directory
   ``sage/graphs``.


Special pragmas
---------------

If Cython code is either attached or loaded as a ``.spyx`` file or
loaded from the notebook as a ``%cython`` block, the following
pragmas are available:

* clang --- may be either c or c++ indicating whether a C or C++
  compiler should be used.

* clib --- additional libraries to be linked in, the space separated
  list is split and passed to distutils.

* cinclude --- additional directories to search for header files. The
  space separated list is split and passed to distutils.

For example:

::

    #clang C++
    #clib givaro
    #cinclude /usr/local/include/


Attaching or loading .spyx
--------------------------

The easiest way to try out Cython without having to learn anything
about distutils, etc., is to create a file with the extension
``spyx``, which stands for "Sage Pyrex":

#. Create a file ``power2.spyx``.

#. Put the following in it:

   ::

       def is2pow(n):
           while n != 0 and n%2 == 0:
               n = n >> 1
           return n == 1

#. Start the Sage command line interpreter and load the ``spyx`` file
   (this will fail if you do not have a C compiler installed).

   .. skip

   ::

       sage: load "power2.spyx"
       Compiling power2.spyx...
       sage: is2pow(12)
       False

Note that you can change ``power2.spyx``, then load it again and it
will be recompiled on the fly. You can also attach ``power2.spyx`` so
it is reloaded whenever you make changes:

.. skip

::

    sage: attach "power2.spyx"

Cython is used for its speed. Here is a timed test on a 2.6 GHz
Opteron:

.. skip

::

    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 0.60 s, sys: 0.00 s, total: 0.60 s
    Wall time: 0.60 s

Now, the code in the file ``power2.spyx`` is valid Python, and if we
copy this to a file ``powerslow.py`` and load that, we get the
following:

.. skip

::

    sage: load "powerslow.py"
    sage: %time [n for n in range(10^5) if is2pow(n)]
    [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048, 4096, 8192, 16384, 32768, 65536]
    CPU times: user 1.01 s, sys: 0.04 s, total: 1.05 s
    Wall time: 1.05 s

By the way, we could gain even a little more speed with the Cython
version with a type declaration, by changing ``def is2pow(n):`` to
``def is2pow(unsigned int n):``.


Other languages
===============

Since Sage is based on Python, it interfaces with C and C++, as well
as other languages. See the Python documentation at
http://www.python.org/doc/ for more details. In particular, the
section "Extending and Embedding the Python Interpreter", available at
http://docs.python.org/ext/ext.html, describes how to write C or
C++ modules for use in Python.


The PARI C library interface
============================

(This chapter was written by Martin Albrecht.)

Here is a step-by-step guide to adding new PARI functions to Sage. We
use the Frobenius form of a matrix as an example.

Some heavy lifting for matrices over integers is implemented using
the PARI library. To compute the Frobenius form in PARI, the
``matfrobenius`` function is used.

There are two ways to interact with the PARI library from Sage. The
gp interface uses the gp interpreter. The PARI interface uses
direct calls to the PARI C functions---this is the preferred way
as it is much faster. Thus this section focuses on using PARI.

We will add a new method to the gen class. This is the abstract
representation of all PARI library objects. That means that once we
add a method to this class, every PARI object, whether it is a number,
polynomial or matrix, will have our new method. So you can do
``pari(1).matfrobenius()``, but since PARI wants to apply
``matfrobenius`` to matrices, not numbers, you will receive a
PariError in this case.

The gen class is defined in
``SAGE_ROOT/devel/sage/sage/libs/pari/gen.pyx``, and this is where we
add the method ``matfrobenius``:

::

        def matfrobenius(self, flag=0):
            """
            matfrobenius(M,{flag}): Return the Frobenius form of the
            square matrix M. If flag is 1, return only the elementary
            divisors. If flag is 2, return a two-components vector [F,B]
            where F is the Frobenius form and B is the basis change
            so that M=B^-1*F*B.
            """
            _sig_on
            return self.new_gen(matfrobenius(self.g, flag))

The ``_sig_on`` statement is some magic for catching segfault signals.
In this way, it prevents SIGSEGVs from the PARI C library crashing the
Sage interpreter. Note that ``self.new_gen()`` calls a closing
``_sig_off`` macro. These two *must always* come in pairs, i.e. every
``_sig_on`` must be matched by a closing ``_sig_off``. The
``self.new_gen()`` call constructs a new Sage-python-gen object from a
given pari-C-gen where the pari-C-gen is stored as the
Sage-python-gen.g attribute. The ``matfrobenius`` call is just a call
to the PARI C library function ``matfrobenius`` with the appropriate
parameters.

The information about which function to call and how to call it can be
retrieved from the PARI user's manual (note: Sage includes the
development version of PARI, so check that version of the user's
manual). Looking for ``matfrobenius`` you can find:
``"The library syntax is matfrobenius(M,flag)"``.

In case you are familiar with gp, please note that the PARI C function
may have a name that is different from the corresponding gp function
(for example, see ``mathnf``), so always check the manual.

We can also add a ``frobenius(flag)`` method to the ``matrix_integer``
class where we call the ``matfrobenius()`` method on the PARI object
associated to the matrix after doing some sanity checking. Then we
convert output from PARI to Sage objects:

::

        def frobenius(self,flag=0):
            """
            If flag is 0 (the default value), return the Frobenius
                form of this matrix.
            If flag is 1, return only the elementary divisors.
            If flag is 2, return a two-component vector [F,B]
                where F is the Frobenius form and B is the basis change
                so that M=B^-1*F*B.

            INPUT:
               flag -- 0,1 or 2 as described above

            ALGORITHM: uses pari's matfrobenius()

            EXAMPLE:
               sage: A = MatrixSpace(IntegerRing(), 3)(range(9))
               sage: A.frobenius(0)
               [ 0  0  0]
               [ 1  0 18]
               [ 0  1 12]
               sage: A.frobenius(1)
               [x3 - 12*x2 - 18*x]
               sage: A.frobenius(2)
               ([ 0  0  0]
               [ 1  0 18]
               [ 0  1 12],
               [    -1      2     -1]
               [     0  23/15 -14/15]
               [     0  -2/15   1/15])
            """
            if self.nrows()!=self.ncols():
                raise ArithmeticError, \
                "frobenius matrix of non-square matrix not defined."
            v = self._pari_().matfrobenius(flag)
            if flag==0:
                return self.matrix_space()(v.python())
            elif flag==1:
                r = polynomial_ring.PolynomialRing(self.base_ring())
                #BUG: this should be handled in PolynomialRing not here
                return [eval(str(x).replace("^","**"),{},r.gens_dict())
                        for x in v.python_list()]
            elif flag==2:
                F = matrix_space.MatrixSpace(rational_field.RationalField(),
                                             self.nrows())(v[0].python())
                B = matrix_space.MatrixSpace(rational_field.RationalField(),
                                             self.nrows())(v[1].python())
                return F,B


GAP
===

(The first version of this chapter was written by David Joyner.)

Wrapping a GAP function in Sage is a matter of writing a program in
Python that uses the pexpect interface to pipe various commands to GAP
and read back the input into Sage. This is sometimes easy, sometimes
hard.

For example, suppose we want to make a wrapper for the computation of
the Cartan matrix of a simple Lie algebra. The Cartan matrix of `G_2`
is available in GAP using the commands

::

    gap> L:= SimpleLieAlgebra( "G", 2, Rationals );
    <Lie algebra of dimension 14 over Rationals>
    gap> R:= RootSystem( L );
    <root system of rank 2>
    gap> CartanMatrix( R );

(Incidentally, most of the GAP Lie algebra implementation was written
by Thomas Breuer, Willem de Graaf and Craig Struble.)

In Sage, one can access these commands by typing

::

    sage: L = gap.SimpleLieAlgebra('"G"', 2, 'Rationals'); L
    Algebra( Rationals, [ v.1, v.2, v.3, v.4, v.5, v.6, v.7, v.8, v.9, v.10,
      v.11, v.12, v.13, v.14 ] )
    sage: R = L.RootSystem(); R
    <root system of rank 2>
    sage: R.CartanMatrix()
    [ [ 2, -1 ], [ -3, 2 ] ]

Note the ``'"G"'`` which is evaluated in GAP as the string ``"G"``.

The purpose of this section is to use this example to show how one
might write a Python/Sage program whose input is, say, ``('G',2)`` and
whose output is the matrix above (but as a Sage Matrix---see the code
in the directory ``SAGE_ROOT/devel/sage/sage/matrix/`` and the
corresponding parts of the Sage reference manual).

First, the input must be converted into strings consisting of legal
GAP commands. Then the GAP output, which is also a string, must be
parsed and converted if possible to a corresponding Sage/Python
object.

::

    def cartan_matrix(type, rank):
        """
        Return the Cartan matrix of given Chevalley type and rank.

        INPUT:
            type -- a Chevalley letter name, as a string, for
                    a family type of simple Lie algebras
            rank -- an integer (legal for that type).

        EXAMPLES:
            sage: cartan_matrix("A",5)
            [ 2 -1  0  0  0]
            [-1  2 -1  0  0]
            [ 0 -1  2 -1  0]
            [ 0  0 -1  2 -1]
            [ 0  0  0 -1  2]
            sage: cartan_matrix("G",2)
            [ 2 -1]
            [-3  2]
        """

        L = gap.SimpleLieAlgebra('"%s"'%type, rank, 'Rationals')
        R = L.RootSystem()
        sM = R.CartanMatrix()
        ans = eval(str(sM))
        MS  = MatrixSpace(QQ, rank)
        return MS(ans)

The output ``ans`` is a Python list. The last two lines convert that
list to an instance of the Sage class ``Matrix``.

Alternatively, one could replace the first line of the above function
with this:

::

        L = gap.new('SimpleLieAlgebra("%s", %s, Rationals);'%(type, rank))

Defining "easy" and "hard" is subjective, but here is one definition.
Wrapping a GAP function is "easy" if there is already a corresponding
class in Python or Sage for the output data type of the GAP function
you are trying to wrap. For example, wrapping any GUAVA (GAP's
error-correcting codes package) function is "easy" since
error-correcting codes are vector spaces over finite fields and GUAVA
functions return one of the following data types:

- vectors over finite fields,

- polynomials over finite fields,

- matrices over finite fields,

- permutation groups or their elements,

- integers.


Sage already has classes for each of these.

A "hard" example is left as an exercise! Here are a few ideas.

- Write a wrapper for GAP's ``FreeLieAlgebra`` function (or, more
  generally, all the finitely presented Lie algebra functions in
  GAP). This would require creating new Python objects.

- Write a wrapper for GAP's ``FreeGroup`` function (or, more
  generally, all the finitely presented groups functions in GAP). This
  would require writing some new Python objects.

- Write a wrapper for GAP's character tables. Though this could be
  done without creating new Python objects, to make the most use of
  these tables, it probably would be best to have new Python objects
  for this.


Singular
========

(The first version of this chapter was written by David Joyner.)

Using Singular functions from Sage is not much different conceptually
from using GAP functions from Sage. As with GAP, this can range from
easy to hard, depending on how much of the data structure of the
output of the Singular function is already present in Sage.

First, some terminology. For us, a *curve* `X` over a finite field `F`
is an equation of the form `f(x,y) = 0`, where `f \in F[x,y]` is a
polynomial. It may or may not be singular. A *place of degree* `d` is
a Galois orbit of `d` points in `X(E)`, where `E/F` is of degree
`d`. For example, a place of degree `1` is also a place of degree `3`,
but a place of degree `2` is not since no degree `3` extension of `F`
contains a degree `2` extension. Places of degree `1` are also called
`F`-rational points.

As an example of the Sage/Singular interface, we will explain how to
wrap Singular's ``NSplaces``, which computes places on a curve over a
finite field. (The command ``closed_points`` also does this in some
cases.) This is "easy" since no new Python classes are needed in Sage
to carry this out.

Here is an example on how to use this command in Singular:

::

     A Computer Algebra System for Polynomial Computations   /   version 3-0-0
                                                           0<
         by: G.-M. Greuel, G. Pfister, H. Schoenemann        \   May 2005
    FB Mathematik der Universitaet, D-67653 Kaiserslautern    \
    > LIB "brnoeth.lib";
    [...]
    > ring s=5,(x,y),lp;
    > poly f=y^2-x^9-x;
    > list X1=Adj_div(f);
    Computing affine singular points ...
    Computing all points at infinity ...
    Computing affine singular places ...
    Computing singular places at infinity ...
    Computing non-singular places at infinity ...
    Adjunction divisor computed successfully

    The genus of the curve is 4
    > list X2=NSplaces(1,X1);
    Computing non-singular affine places of degree 1 ...
    > list X3=extcurve(1,X2);

    Total number of rational places : 6

    > def R=X3[1][5];
    > setring R;
    > POINTS;
    [1]:
       [1]:
          0
       [2]:
          1
       [3]:
          0
    [2]:
       [1]:
          -2
       [2]:
          1
       [3]:
          1
    [3]:
       [1]:
          -2
       [2]:
          1
       [3]:
          1
    [4]:
       [1]:
          -2
       [2]:
          -1
       [3]:
          1
    [5]:
       [1]:
          2
       [2]:
          -2
       [3]:
          1
    [6]:
       [1]:
          0
       [2]:
          0
       [3]:
          1

Here is another way of doing this same calculation in the Sage
interface to Singular:

::

    sage: singular.LIB("brnoeth.lib")
    sage: singular.ring(5,'(x,y)','lp')
        //   characteristic : 5
        //   number of vars : 2
        //        block   1 : ordering lp
        //                  : names    x y
        //        block   2 : ordering C
    sage: f = singular('y^2-x^9-x')
    sage: print singular.eval("list X1=Adj_div(%s);"%f.name())
    Computing affine singular points ...
    Computing all points at infinity ...
    Computing affine singular places ...
    Computing singular places at infinity ...
    Computing non-singular places at infinity ...
    Adjunction divisor computed successfully
    <BLANKLINE>
    The genus of the curve is 4
    sage: print singular.eval("list X2=NSplaces(1,X1);")
    Computing non-singular affine places of degree 1 ...
    sage: print singular.eval("list X3=extcurve(1,X2);")
    <BLANKLINE>
    Total number of rational places : 6
    <BLANKLINE>
    sage: singular.eval("def R=X3[1][5];")
    'def R=X3[1][5];'
    sage: singular.eval("setring R;")
    'setring R;'
    sage: L = singular.eval("POINTS;")

.. link

::

    sage: print L
    [1]:
       [1]:
          0
       [2]:
          1
       [3]:
          0
    [2]:
       [1]:
          0    # 32-bit
          -2   # 64-bit
       [2]:
          0    # 32-bit
          -1   # 64-bit
       [3]:
          1
    ...

From looking at the output, notice that our wrapper function will need
to parse the string represented by `L` above, so let us write a
separate function to do just that. This requires figuring out how to
determine where the coordinates of the points are placed in the string
`L`. Python has some very useful string manipulation commands to do
just that.

::

    def points_parser(string_points,F):
        """
        This function will parse a string of points
        of X over a finite field F returned by Singular's NSplaces
        command into a Python list of points with entries from F.

        EXAMPLES:
            sage: F = GF(5)
            sage: points_parser(L,F)
            ((0, 1, 0), (3, 4, 1), (0, 0, 1), (2, 3, 1), (3, 1, 1), (2, 2, 1))
        """
        Pts=[]
        n=len(L)
        #print n
        #start block to compute a pt
        L1=L
        while len(L1)>32:
            idx=L1.index("     ")
            pt=[]
            ## start block1 for compute pt
            idx=L1.index("     ")
            idx2=L1[idx:].index("\n")
            L2=L1[idx:idx+idx2]
            #print L2
            pt.append(F(eval(L2)))
            # end block1 to compute pt
            L1=L1[idx+8:] # repeat block 2 more times
            #print len(L1)
            ## start block2 for compute pt
            idx=L1.index("     ")
            idx2=L1[idx:].index("\n")
            L2=L1[idx:idx+idx2]
            pt.append(F(eval(L2)))
            # end block2 to compute pt
            L1=L1[idx+8:] # repeat block 1 more time
            ## start block3 for compute pt
            idx=L1.index("     ")
            if "\n" in L1[idx:]:
                idx2=L1[idx:].index("\n")
            else:
                idx2=len(L1[idx:])
            L2=L1[idx:idx+idx2]
            pt.append(F(eval(L2)))
            #print pt
            # end block3 to compute pt
            #end block to compute a pt
            Pts.append(tuple(pt))  # repeat until no more pts
            L1=L1[idx+8:] # repeat block 2 more times
        return tuple(Pts)

Now it is an easy matter to put these ingredients together into a Sage
function which takes as input a triple `(f,F,d)`: a polynomial `f` in
`F[x,y]` defining `X:\  f(x,y)=0` (note that the variables `x,y` must
be used), a finite field `F` *of prime order*, and the degree `d`. The
output is the number of places in `X` of degree `d=1` over `F`. At the
moment, there is no "translation" between elements of `GF(p^d)` in
Singular and Sage unless `d=1`. So, for this reason, we restrict
ourselves to points of degree one.

::

    def places_on_curve(f,F):
        """
        INPUT:
            f -- element of F[x,y], defining X: f(x,y)=0
            F -- a finite field of *prime order*

        OUTPUT:
            integer -- the number of places in X of degree d=1 over F

        EXAMPLES:
            sage: F=GF(5)
            sage: R=MPolynomialRing(F,2,names=["x","y"])
            sage: x,y=R.gens()
            sage: f=y^2-x^9-x
            sage: places_on_curve(f,F)
            ((0, 1, 0), (3, 4, 1), (0, 0, 1), (2, 3, 1), (3, 1, 1), (2, 2, 1))
        """
        d = 1
        p = F.characteristic()
        singular.eval('LIB "brnoeth.lib";')
        singular.eval("ring s="+str(p)+",(x,y),lp;")
        singular.eval("poly f="+str(f))
        singular.eval("list X1=Adj_div(f);")
        singular.eval("list X2=NSplaces("+str(d)+",X1);")
        singular.eval("list X3=extcurve("+str(d)+",X2);")
        singular.eval("def R=X3[1][5];")
        singular.eval("setring R;")
        L = singular.eval("POINTS;")
        return points_parser(L,F)

Note that the ordering returned by this Sage function is exactly the
same as the ordering in the Singular variable ``POINTS``.

One more example (in addition to the one in the docstring):

.. skip

::

    sage: F = GF(2)
    sage: R = MPolynomialRing(F,2,names = ["x","y"])
    sage: x,y = R.gens()
    sage: f = x^3*y+y^3+x
    sage: places_on_curve(f,F)
    ((0, 1, 0), (1, 0, 0), (0, 0, 1))


Singular: Another approach
==========================

There is also a more Python-like interface to Singular. Using this,
the code is much simpler, as illustrated below. First, we demonstrate
computing the places on a curve in a particular case.

::

    sage: singular.lib('brnoeth.lib')
    sage: R = singular.ring(5, '(x,y)', 'lp')
    sage: f = singular.new('y^2 - x^9 - x')
    sage: X1 = f.Adj_div()
    sage: X2 = singular.NSplaces(1, X1)
    sage: X3 = singular.extcurve(1, X2)
    sage: R = X3[1][5]
    sage: singular.set_ring(R)
    sage: L = singular.new('POINTS')

.. link

::

    sage: [(L[i][1], L[i][2], L[i][3]) for i in range(1,7)]
    [(0, 1, 0), (2, 2, 1), (0, 0, 1), (-2, -1, 1), (-2, 1, 1), (2, -2, 1)]  # 32-bit
    [(0, 1, 0), (-2, 1, 1), (-2, -1, 1), (2, 2, 1), (0, 0, 1), (2, -2, 1)]  # 64-bit

Next, we implement the general function (for brevity we omit the
docstring, which is the same as above). Note that the ``point_parser``
function is not required.

::

    def places_on_curve(f,F):
        p = F.characteristic()
        if F.degree() > 1:
            raise NotImplementedError
        singular.lib('brnoeth.lib')
        R = singular.ring(5, '(x,y)', 'lp')
        f = singular.new('y^2 - x^9 - x')
        X1 = f.Adj_div()
        X2 = singular.NSplaces(1, X1)
        X3 = singular.extcurve(1, X2)
        R = X3[1][5]
        singular.setring(R)
        L = singular.new('POINTS')
        return [(int(L[i][1]), int(L[i][2]), int(L[i][3])) \
                 for i in range(1,int(L.size())+1)]

This code is much shorter, nice, and more readable. However, it
depends on certain functions, e.g. ``singular.setring`` having been
implemented in the Sage/Singular interface, whereas the code in the
previous section used only the barest minimum of that interface.


Creating a new pseudo-tty interface
===================================

You can create Sage pseudo-tty interfaces that allow Sage to work with
almost any command line program, and which do not require any
modification or extensions to that program. They are also surprisingly
fast and flexible (given how they work!), because all I/O is buffered,
and because interaction between Sage and the command line program can
be non-blocking (asynchronous). A pseudo-tty Sage interface is
asynchronous because it derives from the Sage class ``Expect``, which
handles the communication between Sage and the external process.

For example, here is part of the file
``SAGE_ROOT/devel/sage/sage/interfaces/octave.py``, which
defines an interface between Sage and Octave, an open source program
for doing numerical computations, among other things.

::

    import os
    from expect import Expect, ExpectElement

    class Octave(Expect):
        ...

The first two lines import the library ``os``, which contains
operating system routines, and also the class ``Expect``, which is the
basic class for interfaces. The third line defines the class
``Octave``; it derives from ``Expect`` as well. After this comes a
docstring, which we omit here (see the file for details). Next comes:

::

        def __init__(self, maxread=100, script_subdirectory="", logfile=None,
                     server=None, server_tmpdir=None):
            Expect.__init__(self,
                            name = 'octave',
                            prompt = '>',
                            command = "octave --no-line-editing --silent",
                            maxread = maxread,
                            server = server,
                            server_tmpdir = server_tmpdir,
                            script_subdirectory = script_subdirectory,
                            restart_on_ctrlc = False,
                            verbose_start = False,
                            logfile = logfile,
                            eval_using_file_cutoff=100)

This uses the class ``Expect`` to set up the Octave interface.

::

        def set(self, var, value):
            """
            Set the variable var to the given value.
            """
            cmd = '%s=%s;'%(var,value)
            out = self.eval(cmd)
            if out.find("error") != -1:
                raise TypeError, "Error executing code in Octave\nCODE:\n\t%s\nOctave ERROR:\n\t%s"%(cmd, out)

        def get(self, var):
            """
            Get the value of the variable var.
            """
            s = self.eval('%s'%var)
            i = s.find('=')
            return s[i+1:]

        def console(self):
            octave_console()

These let users type ``octave.set('x', 3)``, after which
``octave.get('x')`` returns ``' 3'``. Running ``octave.console()``
dumps the user into an Octave interactive shell.

::

        def solve_linear_system(self, A, b):
            """
            Use octave to compute a solution x to A*x = b, as a list.

            INPUT:
                A -- mxn matrix A with entries in QQ or RR
                b -- m-vector b entries in QQ or RR (resp)

            OUTPUT:
                An list x (if it exists) which solves M*x = b

            EXAMPLES:
                sage: M33 = MatrixSpace(QQ,3,3)
                sage: A   = M33([1,2,3,4,5,6,7,8,0])
                sage: V3  = VectorSpace(QQ,3)
                sage: b   = V3([1,2,3])
                sage: octave.solve_linear_system(A,b)    # requires optional octave
                [-0.33333299999999999, 0.66666700000000001, -3.5236600000000002e-18]

            AUTHOR: David Joyner and William Stein
            """
            m = A.nrows()
            n = A.ncols()
            if m != len(b):
                raise ValueError, "dimensions of A and b must be compatible"
            from sage.matrix.all import MatrixSpace
            from sage.rings.all import QQ
            MS = MatrixSpace(QQ,m,1)
            b  = MS(list(b)) # converted b to a "column vector"
            sA = self.sage2octave_matrix_string(A)
            sb = self.sage2octave_matrix_string(b)
            self.eval("a = " + sA )
            self.eval("b = " + sb )
            soln = octave.eval("c = a \\ b")
            soln = soln.replace("\n\n ","[")
            soln = soln.replace("\n\n","]")
            soln = soln.replace("\n",",")
            sol  = soln[3:]
            return eval(sol)

This code defines the method ``solve_linear_system``, which works as
documented.

These are only excerpts from ``octave.py``; check that file for more
definitions and examples. Look at other files in the directory
``SAGE_ROOT/devel/sage/sage/interfaces/`` for examples of interfaces
to other software packages.


.. [3] See http://www.sagemath.org/links-components.html for a list
