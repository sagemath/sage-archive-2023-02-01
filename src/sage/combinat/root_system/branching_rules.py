"""
Branching Rules
"""
#*****************************************************************************
#  Copyright (C) 2014 Daniel Bump <bump at match.stanford.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.root_system.weyl_characters import branching_rule
from sage.combinat.root_system.cartan_type import CartanType

def maximal_subgroups(ct, mode="print_rules"):
    """
    Given a classical Cartan type (of rank less than or equal to 8)
    this prints the Cartan types of maximal subgroups, with a method
    of obtaining the branching rule. The string to the right of the
    colon in the output is a command to create a branching rule.

    INPUT:

    - ``ct`` -- a classical irreducible Cartan type

    EXAMPLES::

        sage: maximal_subgroups("D4")
        B3:branching_rule("D4","B3","symmetric")
        A2:branching_rule("D4","A2(1,1)","plethysm")
        A1xC2:branching_rule("D4","C1xC2","tensor")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])
        A1xA1xA1xA1:branching_rule("D4","D2xD2","orthogonal_sum")*branching_rule("D2xD2","A1xA1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])

    You may use the optional argument ``mode="get_rule"`` to extract one branching rule
    from this list. This gives the same result as typing the given string
    at the command line::

        sage: maximal_subgroups("D4",mode="get_rule")["A2"]
        plethysm (along A2(1,1)) branching rule D4 => A2

    It is believed that the list of maximal subgroups is complete, except that some
    subgroups may be not be invariant under outer automorphisms. It is reasonable
    to want a list of maximal subgroups that is complete up to conjugation,
    but to obtain such a list you may have to apply outer automorphisms.
    The group of outer automorphisms modulo inner automorphisms is isomorphic
    to the group of symmetries of the Dynkin diagram, and these are available

    For example, we verify that the maximal subgroup `A_1\times C_2` from
    the above list is not invariant under inner automorphisms, by checking
    that applying the triality automorphism changes it to another. We
    check this by observing that the first spin representation has different
    branching. There are (up to conjugacy) three maximal subgroups of
    `D_4` of type `A_1\times C_2`::

        sage: b = maximal_subgroups("D4",mode="get_rule")["A1xC2"]; b
        composite branching rule D4 => (tensor) C1xC2 => ([isomorphic branching rule C1 => A1, 'identity']) A1xC2
        sage: b1 = branching_rule("D4","D4","triality")*b
        sage: [D4,A1xC2]=[WeylCharacterRing(x,style="coroots") for x in ["D4","A1xC2"]]
        sage: D4(0,0,1,0).branch(A1xC2,rule=b)
        A1xC2(1,1,0)
        sage: D4(0,0,1,0).branch(A1xC2,rule=b1)
        A1xC2(2,0,0) + A1xC2(0,0,1)
    """

    if CartanType(ct) == CartanType("A2"):
        rul = ["""A1:branching_rule("A2","A1","levi")"""]
    elif CartanType(ct) == CartanType("A3"):
        rul = ["""A2:branching_rule("A3","A2","levi")""", \
               """A1xA1:branching_rule("A3","A1xA1","tensor")""", \
               """C2:branching_rule("A3","C2","symmetric")""", \
               """A1xA1:branching_rule("A3","A1xA1","levi")"""]
    elif CartanType(ct) == CartanType("A4"):
        rul = ["""A3:branching_rule("A4","A3","levi")""", \
               """B2:branching_rule("A4","B2","symmetric")""", \
               """A1xA2:branching_rule("A4","A1xA2","levi")"""]
    elif CartanType(ct) == CartanType("A5"):
        rul = ["""A4:branching_rule("A5","A4","levi")""", \
               """A3:branching_rule("A5","D3","symmetric")*branching_rule("D3","A3","isomorphic")""", \
               """A3:branching_rule("A5","A3(0,1,0)","plethysm") # alternative""", \
               """C3:branching_rule("A5","C3","symmetric")""", \
               """A2:branching_rule("A5","A2(2,0)","plethysm")""", \
               """A1xA2:branching_rule("A5","A1xA2","tensor")""", \
               """A1xA3:branching_rule("A5","A1xA3","levi")""", \
               """A2xA2:branching_rule("A5","A2xA2","levi")"""]
    elif CartanType(ct) == CartanType("A6"):
        rul = ["""A5:branching_rule("A6","A5","levi")""", \
               """B3:branching_rule("A6","B3","symmetric")""", \
               """A1xA4:branching_rule("A6","A1xA4","levi")""", \
               """A2xA3:branching_rule("A6","A2xA3","levi")"""]
    elif CartanType(ct) == CartanType("A7"):
        rul = ["""A6:branching_rule("A7","A6","levi")""", \
               """C4:branching_rule("A7","C4","symmetric")""", \
               """D4:branching_rule("A7","D4","symmetric")""", \
               """A1xA3:branching_rule("A7","A1xA3","tensor")""", \
               """A1xA5:branching_rule("A7","A1xA5","levi")""", \
               """A2xA4:branching_rule("A7","A2xA4","levi")""", \
               """A3xA3:branching_rule("A7","A3xA3","levi")"""]
    elif CartanType(ct) == CartanType("A8"):
        rul = ["""A7:branching_rule("A8","A7","levi")""", \
               """B4:branching_rule("A8","B4","symmetric")""", \
               """A2xA2:branching_rule("A8","A2xA2","tensor")""", \
               """A1xA6:branching_rule("A8","A1xA6","levi")""", \
               """A2xA5:branching_rule("A8","A2xA5","levi")""", \
               """A3xA4:branching_rule("A8","A3xA4","levi")"""]
    elif CartanType(ct) == CartanType("B3"):
        rul = ["""G2:branching_rule("B3","G2","miscellaneous")""", \
               """A3:branching_rule("B3","D3","extended")*branching_rule("D3","A3","isomorphic")""", \
               """A1xA1xA1:branching_rule("B3","D2xB1","orthogonal_sum")*branching_rule("D2xB1","A1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("B1","A1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("B4"):
        rul = ["""D4:branching_rule("B4","D4","extended")""", \
               """A1:branching_rule("B4","A1","symmetric_power")""", \
               """A1xA1:branching_rule("B4","B1xB1","tensor")*branching_rule("B1xB1","A1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("B1","A1","isomorphic")])""", \
               """A1xA1xB2:branching_rule("B4","D2xB2","extended")*branching_rule("D2xB2","A1xA1xB2",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """A1xA3:branching_rule("B4","B1xD3","extended")*branching_rule("B1xD3","A1xA3",[branching_rule("B1","A1","isomorphic"),branching_rule("D3","A3","isomorphic")])"""]
    elif CartanType(ct) == CartanType("B5"):
        rul = ["""D5:branching_rule("B5","D5","extended")""", \
               """A1:branching_rule("B5","A1","symmetric_power")""", \
               """A1xA2xB3:branching_rule("B5","D2xB3","extended")*branching_rule("D2xB3","A1xA2xB3",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """A1xD4:branching_rule("B5","B1xD4","orthogonal_sum")*branching_rule("B1xD4","A1xD4",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """A3xB2:branching_rule("B5","D3xB2","orthogonal_sum")*branching_rule("D3xB2","A3xB2",[branching_rule("D3","A3","isomorphic"),"identity"])"""]
    elif CartanType(ct) == CartanType("B6"):
        rul = ["""D6:branching_rule("B6","D6","extended")""", \
               """A1:branching_rule("B6","A1","symmetric_power")""", \
               """A1xA1xB4:branching_rule("B6","D2xB4","orthogonal_sum")*branching_rule("D2xB4","A1xA1xB4",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """A1xD5:branching_rule("B6","B1xD5","orthogonal_sum")*branching_rule("B1xD5","A1xD5",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """A3xB3:branching_rule("B6","D3xB3","orthogonal_sum")*branching_rule("D3xB3","A3xB3",[branching_rule("D3","A3","isomorphic"),"identity"])""", \
               """B2xD4:branching_rule("B6","B2xD4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("B7"):
        rul = ["""D7:branching_rule("B7","D7","extended")""", \
               """A3:branching_rule("B7","A3(1,0,1)","plethysm")""", \
               """A1:branching_rule("B7","A1","symmetric_power")""", \
               """A1xB2:branching_rule("B7","B1xB2","tensor")*branching_rule("B1xB2","A1xB2",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """A1xD6:branching_rule("B7","B1xD6","extended")*branching_rule("B1xD6","A1xD6",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """A1xA1xB5:branching_rule("B7","D2xB5","extended")*branching_rule("D2xB5","A1xA1xB5",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """B2xD5:branching_rule("B7","B2xD5","orthogonal_sum")""", \
               """A3xB4:branching_rule("B7","D3xB4","orthogonal_sum")*branching_rule("D3xB4","A3xB4",[branching_rule("D3","A3","isomorphic"),"identity"])""", \
               """B3xD4:branching_rule("B7","B3xD4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("B8"):
        rul = ["""D8:branching_rule("B8","D8","extended")""", \
               """A1:branching_rule("B8","A1","symmetric_power")""", \
               """A1xD7:branching_rule("B8","B1xD7","orthogonal_sum")*branching_rule("B1xD7","A1xD7",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """A1xA1xB6:branching_rule("B8","D2xB6","orthogonal_sum")*branching_rule("D2xB6","A1xA1xB6",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """B2xD6:branching_rule("B8","B2xD6","orthogonal_sum")""", \
               """A3xB5:branching_rule("B8","D3xB5","orthogonal_sum")*branching_rule("D3xB5","A3xB5",[branching_rule("D3","A3","isomorphic"),"identity"])""", \
               """B3xD5:branching_rule("B8","B3xD5","orthogonal_sum")""", \
               """B4xD4:branching_rule("B8","B4xD4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("C2"):
        rul = ["""A1:branching_rule("C2","A1","symmetric_power")""", \
               """A1xA1:branching_rule("C2","C1xC1","orthogonal_sum")*branching_rule("C1xC1","A1xA1",[branching_rule("C1","A1","isomorphic"),branching_rule("C1","A1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("C3"):
        rul = ["""A2:branching_rule("C3","A2","levi")""", \
               """A1:branching_rule("C3","A1","symmetric_power")""", \
               """A1xA1:branching_rule("C3","B1xC1","tensor")*branching_rule("B1xC1","A1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("C1","A1","isomorphic")])""", \
               """A1xC2:branching_rule("C3","C1xC2","orthogonal_sum")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])"""]
    elif CartanType(ct) == CartanType("C4"):
        rul = ["""A3:branching_rule("C4","A3","levi")""", \
               """A1:branching_rule("C4","A1","symmetric_power")""", \
               """A1xA3:branching_rule("C4","C1xC3","orthogonal_sum")*branching_rule("C1xC3","A1xA3",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """C2xC2:branching_rule("C4","C2xC2","orthogonal_sum")""", \
               """A1xA1xA1:branching_rule("C4","C1xD2","tensor")*branching_rule("C1xD2","A1xA1xA1",[branching_rule("C1","A1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("C5"):
        rul = ["""A4:branching_rule("C5","A4","levi")""", \
               """A1:branching_rule("C5","A1","symmetric_power")""", \
               """A1xC4:branching_rule("C5","C1xC4","orthogonal_sum")*branching_rule("C1xC4","A1xC4",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """C2xC3:branching_rule("C5","C2xC3","orthogonal_sum")""", \
               """A1xB2:branching_rule("C5","C1xB2","tensor")*branching_rule("C1xB2","A1xB2",[branching_rule("C1","A1","isomorphic"),"identity"])"""]
    elif CartanType(ct) == CartanType("C6"):
        rul = ["""A5:branching_rule("C6","A5","levi")""", \
               """A1:branching_rule("C6","A1","symmetric_power")""", \
               """A1xA3:branching_rule("C6","C1xD3","tensor")*branching_rule("C1xD3","A1xA3",[branching_rule("C1","A1","isomorphic"),branching_rule("D3","A3","isomorphic")])""", \
               """A1xC2:branching_rule("C6","B1xC2","tensor")*branching_rule("B1xC2","A1xC2",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """A1xC5:branching_rule("C6","C1xC5","orthogonal_sum")*branching_rule("C1xC5","A1xC5",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """C2xC4:branching_rule("C6","C2xC4","orthogonal_sum")""", \
               """C3xC3:branching_rule("C6","C3xC3","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("C7"):
        rul = ["""A6:branching_rule("C7","A6","levi")""", \
               """A1:branching_rule("C7","A1","symmetric_power")""", \
               """A1xB3:branching_rule("C7","C1xB3","tensor")*branching_rule("C1xB3","A1xB3",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """A1xC6:branching_rule("C7","C1xC6","orthogonal_sum")*branching_rule("C1xC6","A1xC6",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """C2xC5:branching_rule("C7","C2xC5","orthogonal_sum")""", \
               """C3xC4:branching_rule("C7","C3xC4","orthogonal_sum")""", \
               """C3:branching_rule("C7","C3(0,0,1)","plethysm") # overlooked by Patera and McKay"""]
    elif CartanType(ct) == CartanType("C8"):
        rul = ["""A7:branching_rule("C8","A7","levi")""", \
               """A1:branching_rule("C8","A1","symmetric_power")""", \
               """C2:branching_rule("C8","C2(1,1)","plethysm")""", \
               """A1xD4:branching_rule("C8","C1xD4","tensor")*branching_rule("C1xD4","A1xD4",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """A1xC7:branching_rule("C8","C1xC7","orthogonal_sum")*branching_rule("C1xC7","A1xC7",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """C2xC6:branching_rule("C8","C2xC6","orthogonal_sum")""", \
               """C3xC5:branching_rule("C8","C3xC5","orthogonal_sum")""", \
               """C4xC4:branching_rule("C8","C4xC4","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("D4"):
        rul = ["""B3:branching_rule("D4","B3","symmetric")""", \
               """A2:branching_rule("D4","A2(1,1)","plethysm")""", \
               """A1xC2:branching_rule("D4","C1xC2","tensor")*branching_rule("C1xC2","A1xC2",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """A1xA1xA1xA1:branching_rule("D4","D2xD2","orthogonal_sum")*branching_rule("D2xD2","A1xA1xA1xA1",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("D5"):
        rul = ["""A4:branching_rule("D5","A4","levi")""", \
               """B4:branching_rule("D5","B4","symmetric")""", \
               """C2:branching_rule("D5","C2(2,0)","plethysm")""", \
               """A1xA1xA3:branching_rule("D5","D2xD3","orthogonal_sum")*branching_rule("D2xD3","A1xA1xA3",[branching_rule("D2","A1xA1","isomorphic"),branching_rule("D3","A3","isomorphic")])""", \
               """A1xA3:branching_rule("D5","B1xB3","orthogonal_sum")*branching_rule("B1xB3","A1xA3",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """B2xB2:branching_rule("D5","B2xB2","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("D6"):
        rul = ["""A5:branching_rule("D6","A5","levi")""", \
               """B5:branching_rule("D6","B5","symmetric")""", \
               """A1xA3:branching_rule("D6","C1xC3","tensor")*branching_rule("C1xC3","A1xA3",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """A1xA1xD4:branching_rule("D6","D2xD4","orthogonal_sum")*branching_rule("D2xD4","A1xA1xD4",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """A3xA3:branching_rule("D6","D3xD3","orthogonal_sum")*branching_rule("D3xD3","A3xA3",[branching_rule("D3","A3","isomorphic"),branching_rule("D3","A3","isomorphic")])""", \
               """A1xB4:branching_rule("D6","B1xB4","orthogonal_sum")*branching_rule("B1xB4","A1xB4",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """B2xB3:branching_rule("D6","B2xB3","orthogonal_sum")""", \
               """A1xA1xA1:branching_rule("D6","B1xD2","tensor")*branching_rule("B1xD2","A1xA1xA1",[branching_rule("B1","A1","isomorphic"),branching_rule("D2","A1xA1","isomorphic")])"""]
    elif CartanType(ct) == CartanType("D7"):
        rul = ["""A6:branching_rule("D7","A6","levi")""", \
               """B6:branching_rule("D7","B6","symmetric")""", \
               """C3:branching_rule("D7","C3(0,1,0)","plethysm")""", \
               """C2:branching_rule("D7","C2(0,2)","plethysm")""", \
               """G2:branching_rule("D7","G2(0,1)","plethysm")""", \
               """A1xA1xD5:branching_rule("D7","D2xD5","orthogonal_sum")*branching_rule("D2xD5","A1xA1xD5",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """A3xD4:branching_rule("D7","D3xD4","orthogonal_sum")*branching_rule("D3xD4","A3xD4",[branching_rule("D3","A3","isomorphic"),"identity"])""", \
               """A1xB5:branching_rule("D7","B1xB5","orthogonal_sum")*branching_rule("B1xB5","A1xB5",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """B2xB4:branching_rule("D7","B2xB4","orthogonal_sum")""", \
               """B3xB3:branching_rule("D7","B3xB3","orthogonal_sum")"""]
    elif CartanType(ct) == CartanType("D8"):
        rul = ["""A7:branching_rule("D8","A7","levi")""", \
               """B7:branching_rule("D8","B7","symmetric")""", \
               """B4:branching_rule("D8","B4(0,0,0,1)","plethysm")""", \
               """A1xC4:branching_rule("D8","C1xC4","tensor")*branching_rule("C1xC4","A1xC4",[branching_rule("C1","A1","isomorphic"),"identity"])""", \
               """A1xA1xD6:branching_rule("D8","D2xD6","orthogonal_sum")*branching_rule("D2xD6","A1xA1xD6",[branching_rule("D2","A1xA1","isomorphic"),"identity"])""", \
               """A3xD5:branching_rule("D8","D3xD5","orthogonal_sum")*branching_rule("D3xD5","A3xD5",[branching_rule("D3","A3","isomorphic"),"identity"])""", \
               """D4xD4:branching_rule("D8","D4xD4","orthogonal_sum")""", \
               """A1xB6:branching_rule("D8","B1xB6","orthogonal_sum")*branching_rule("B1xB6","A1xB6",[branching_rule("B1","A1","isomorphic"),"identity"])""", \
               """B2xB5:branching_rule("D8","B2xB5","orthogonal_sum")""", \
               """B3xB4:branching_rule("D8","B3xB4","orthogonal_sum")""", \
               """C2xC2:branching_rule("D8","C2xC2","tensor")"""]
    elif CartanType(ct) == CartanType("G2"):
        rul = ["""A2:branching_rule("G2","A2","extended")""", \
               """A1:branching_rule("G2","A1","i")""", \
               """A1xA1:branching_rule("G2","A1xA1","extended")"""]
    elif CartanType(ct) == CartanType("F4"):
        rul = ["""B4:branching_rule("F4","B4","extended")""", \
               """A1:branching_rule("F4","A1","ii")""", \
               """A1xG2:branching_rule("F4","A1xG2","miscellaneous")""", \
               """A1xC3:branching_rule("F4","A1xC3","extended")""", \
               """A2xA2:branching_rule("F4","A2xA2","extended")"""]
    elif CartanType(ct) == CartanType("E6"):
        rul = ["""D5:branching_rule("E6","D5","levi")""", \
               """C4:branching_rule("E6","C4","symmetric")""", \
               """F4:branching_rule("E6","F4","symmetric")""", \
               """A2:branching_rule("E6","A2","miscellaneous")""", \
               """G2:branching_rule("E6","G2","miscellaneous")""", \
               """A2xG2:branching_rule("E6","A2xG2","miscellaneous")""", \
               """A1xA5:branching_rule("E6","A1xA5","extended")""", \
               """A2xA2xA2:branching_rule("E6","A2xA2xA2","extended")"""]
    elif CartanType(ct) == CartanType("E7"):
        rul = ["""A7:branching_rule("E7","A7","extended")""", \
               """E6:branching_rule("E7","E6","levi")""", \
               """A2:branching_rule("E7","A2","miscellaneous")""", \
               """A1:branching_rule("E7","A1","iii")""", \
               """A1:branching_rule("E7","A1","iv")""", \
               """A1xF4:branching_rule("E7","A1xF4","miscellaneous")""", \
               """G2xC3:branching_rule("E7","G2xC3","miscellaneous")""", \
               """A1xG2:branching_rule("E7","A1xG2","miscellaneous")""", \
               """A1xA1:branching_rule("E7","A1xA1","miscellaneous")""", \
               """A1xD6:branching_rule("E7","A1xD6","extended")""", \
               """A5xA2:branching_rule("E7","A5xA2","extended")"""]
    elif CartanType(ct) == CartanType("E8"):
        rul = ["""A4xA4:branching_rule("E8","A4xA4","extended")""", \
               """G2xF4:branching_rule("E8","G2xF4","miscellaneous")""", \
               """E6xA2:branching_rule("E8","E6xA2","extended")""", \
               """E7xA1:branching_rule("E8","E7xA1","extended")""", \
               """D8:branching_rule("E8","D8","extended")""", \
               """A8:branching_rule("E8","A8","extended")""", \
               """B2:branching_rule("E8","B2","miscellaneous")""", \
               """A1xA2:branching_rule("E8","A1xA2","miscellaneous")""", \
               """A1:branching_rule("E8","A1","v")""", \
               """A1:branching_rule("E8","A1","vi")""", \
               """A1:branching_rule("E8","A1","vii")"""]
    else:
        raise ValueError("Argument must be an irreducible classical Cartan Type with rank less than or equal to 8")
    if mode == "print_rules":
        for line in rul:
            print line
    elif mode == "get_rule":
        d = {}
        for line in rul:
            [k, br] = line.split(":")
            br = eval(br)
            if d.has_key(k):
                if type(d[k]) is not list:
                    d[k] = [d[k]]
                d[k].append(br)
            else:
                d[k] = br
        return d
            
