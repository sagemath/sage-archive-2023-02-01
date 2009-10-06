.. _chapter-groups:

******
Groups
******

.. index::
   pair: group; permutation

.. _section-permutation:

Permutation groups
==================

A permutation group is a subgroup of some symmetric group
:math:`S_n`. Sage has a Python class ``PermutationGroup``, so you
can work with such groups directly:

::

    sage: G = PermutationGroup(['(1,2,3)(4,5)'])
    sage: G
    Permutation Group with generators [(1,2,3)(4,5)]
    sage: g = G.gens()[0]; g
    (1,2,3)(4,5)
    sage: g*g
    (1,3,2)
    sage: G = PermutationGroup(['(1,2,3)'])
    sage: g = G.gens()[0]; g
    (1,2,3)
    sage: g.order()
    3

For the example of the Rubik's cube group (a permutation subgroup
of :math:`S_{48}`, where the non-center facets of the Rubik's
cube are labeled :math:`1,2,...,48` in some fixed way), you can
use the GAP-Sage interface as follows.

.. index::
   pair: group; Rubik's cube

.. skip

::

    sage: cube = "cubegp := Group(
    ( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
    ( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
    (17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
    (25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
    (33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
    (41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) )"
    sage: gap(cube)
    'permutation group with 6 generators'
    sage: gap("Size(cubegp)")
    43252003274489856000'

Another way you can choose to do this:


-  Create a file ``cubegroup.py`` containing the lines

   ::

       cube = "cubegp := Group(
       ( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)(11,35,27,19),
       ( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)( 6,22,46,35),
       (17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)( 8,30,41,11),
       (25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)( 8,33,48,24),
       (33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)( 1,14,48,27),
       (41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)(16,24,32,40) )"

-  Place the file in the subdirectory
   ``$SAGE_ROOT/local/lib/python2.4/site-packages/sage`` of your Sage home
   directory.

-  Read (i.e.,``{import``) it into Sage:

   .. skip

::

       sage: import sage.cubegroup
       sage: sage.cubegroup.cube
       'cubegp := Group(( 1, 3, 8, 6)( 2, 5, 7, 4)( 9,33,25,17)(10,34,26,18)
       (11,35,27,19),( 9,11,16,14)(10,13,15,12)( 1,17,41,40)( 4,20,44,37)
       ( 6,22,46,35),(17,19,24,22)(18,21,23,20)( 6,25,43,16)( 7,28,42,13)
       ( 8,30,41,11),(25,27,32,30)(26,29,31,28)( 3,38,43,19)( 5,36,45,21)
       ( 8,33,48,24),(33,35,40,38)(34,37,39,36)( 3, 9,46,32)( 2,12,47,29)
       ( 1,14,48,27),(41,43,48,46)(42,45,47,44)(14,22,30,38)(15,23,31,39)
       (16,24,32,40) )'
       sage: gap(sage.cubegroup.cube)
       'permutation group with 6 generators'
       sage: gap("Size(cubegp)")
       '43252003274489856000'

   (You will have line wrap instead of the above carriage returns in
   your Sage output.)

-  Use the ``CubeGroup`` class

   ::

       sage: rubik = CubeGroup()
       sage: rubik
       The PermutationGroup of all legal moves of the Rubik's cube.
       sage: print rubik
       The Rubik's cube group with generators R,L,F,B,U,D in SymmetricGroup(48).
       sage: G = rubik.group()
       sage: G.order()
       43252003274489856000


    (1) has implemented classical groups (such as
    :math:`GU(3,GF(5))`) and matrix groups over a finite field with
    user-defined generators.

    (2) also has implemented finite and infinite (but finitely
    generated) abelian groups.

.. index::
   pair: group; conjugacy classes

.. _section-conjugacy:

Conjugacy classes
=================

You can compute conjugacy classes of a finite group using
"natively":

::

    sage: G = PermutationGroup(['(1,2,3)', '(1,2)(3,4)', '(1,7)'])
    sage: CG = G.conjugacy_classes_representatives()
    sage: gamma = CG[2]
    sage: CG; gamma
    [(), (1,2), (1,2)(3,4), (1,2,3), (1,2,3)(4,7), (1,2,3,4), (1,2,3,4,7)]
    (1,2)(3,4)

You can use the Sage-GAP interface.

::

    sage: gap.eval("G := Group((1,2)(3,4),(1,2,3))")
    'Group([ (1,2)(3,4), (1,2,3) ])'
    sage: gap.eval("CG := ConjugacyClasses(G)")
    '[ ()^G, (1,2)(3,4)^G, (1,2,3)^G, (1,2,4)^G ]'
    sage: gap.eval("gamma := CG[3]")
    '(1,2,3)^G'
    sage: gap.eval("g := Representative(gamma)")
    '(1,2,3)'

Or, here's another (more "pythonic") way to do this type of
computation:

::

    sage: G = gap.Group('[(1,2,3), (1,2)(3,4), (1,7)]')
    sage: CG = G.ConjugacyClasses()
    sage: gamma = CG[2]
    sage: g = gamma.Representative()
    sage: CG; gamma; g
    [ ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), () ),
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2) ),
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2)(3,4) ),
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2,3) ),
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2,3)(4,7) ),
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2,3,4) ),
      ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2,3,4,7) ) ]
    ConjugacyClass( SymmetricGroup( [ 1, 2, 3, 4, 7 ] ), (1,2) )
    (1,2)

.. index::
   pair: group; normal subgroups

.. _section-normal:

Normal subgroups
================

If you want to find all the normal subgroups of a permutation group
:math:`G` (up to conjugacy), you can use 's interface to GAP:

::

    sage: G = AlternatingGroup( 5 )
    sage: gap(G).NormalSubgroups()
    [ Group( () ), AlternatingGroup( [ 1 .. 5 ] ) ]

or

::

    sage: G = gap("AlternatingGroup( 5 )")
    sage: G.NormalSubgroups()
    [ Group( () ), AlternatingGroup( [ 1 .. 5 ] ) ]

Here's another way, working more directly with GAP:

::

    sage: print gap.eval("G := AlternatingGroup( 5 )")
    Alt( [ 1 .. 5 ] )
    sage: print gap.eval("normal := NormalSubgroups( G )")
    [ Group(()), Alt( [ 1 .. 5 ] ) ]
    sage: G = gap.new("DihedralGroup( 10 )")
    sage: G.NormalSubgroups()
    [ Group( <identity> of ... ), Group( [ f2 ] ), Group( [ f1, f2 ] ) ]
    sage: print gap.eval("G := SymmetricGroup( 4 )")
    Sym( [ 1 .. 4 ] )
    sage: print gap.eval("normal := NormalSubgroups( G );")
    [ Group(()), Group([ (1,4)(2,3), (1,3)(2,4) ]),
      Group([ (2,4,3), (1,4)(2,3), (1,3)(2,4) ]), Sym( [ 1 .. 4 ] ) ]

.. index::
   pair: groups; center

.. _section-center:

Centers
=======

How do you compute the center of a group in Sage?

Although Sage calls GAP to do the computation of the group center,
``center`` is "wrapped" (i.e., Sage has a class PermutationGroup with
associated class method "center"), so the user does not need to use
the ``gap`` command. Here's an example:

::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.center()
    Permutation Group with generators [()]

A similar syntax for matrix groups also works:

::

    sage: G = SL(2, GF(5) )
    sage: G.center()
    Matrix group over Finite Field of size 5 with 1 generators:
     [[[4, 0], [0, 4]]]
    sage: G = PSL(2, 5 )
    sage: G.center()
    Permutation Group with generators [()]

Note: ``center`` can be spelled either way in GAP, not so in Sage.

The group id database
=====================

The function ``group_id`` requires that the Small Groups Library of
E. A. O'Brien, B. Eick, and H. U. Besche be installed (you can do
this by typing ``./sage -i database_gap-4.4.9`` in the Sage home
directory).

::

    sage: G = PermutationGroup(['(1,2,3)(4,5)', '(3,4)'])
    sage: G.order()
    120
    sage: G.group_id()      # requires optional GAP database package
    [120, 34]

Another example of using the small groups database: ``group_id``

.. skip

::

    sage: gap_console()
    GAP4, Version: 4.4.6 of 02-Sep-2005, x86_64-unknown-linux-gnu-gcc
    gap> G:=Group((4,6,5)(7,8,9),(1,7,2,4,6,9,5,3));
    Group([ (4,6,5)(7,8,9), (1,7,2,4,6,9,5,3) ])
    gap> StructureDescription(G);
    "(((C3 x C3) : Q8) : C3) : C2"
