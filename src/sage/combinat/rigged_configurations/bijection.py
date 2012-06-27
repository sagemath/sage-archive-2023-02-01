r"""
Bijection between rigged configurations and KR tableaux

Functions which are big switch statements to create the bijection class of the
correct type.

AUTHORS:

- Travis Scrimshaw (2011-04-15): Initial version
"""

#*****************************************************************************
#       Copyright (C) 2011, 2012 Travis Scrimshaw <tscrim@ucdavis.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.combinat.rigged_configurations.bij_type_A import KRTToRCBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_A import RCToKRTBijectionTypeA
from sage.combinat.rigged_configurations.bij_type_D import KRTToRCBijectionTypeD
from sage.combinat.rigged_configurations.bij_type_D import RCToKRTBijectionTypeD

def KRTToRCBijection(krt_elt):
    r"""
    Return the correct KR tableaux to rigged configuration bijection helper class.

    TESTS::

        sage: KRT = TensorProductOfKirillovReshetikhinTableaux(['A', 4, 1], [[2,1]])
        sage: from sage.combinat.rigged_configurations.bijection import KRTToRCBijection
        sage: bijection = KRTToRCBijection(KRT(pathlist=[[4,3]]))
        sage: bijection.next_state(3, 0)
        sage: bijection.ret_rig_con
        <BLANKLINE>
        -1[ ]-1
        <BLANKLINE>
        -1[ ]-1
        <BLANKLINE>
        (/)
        <BLANKLINE>
        (/)
        <BLANKLINE>
    """
    ct = krt_elt.cartan_type()
    if ct.letter == 'A':
        return KRTToRCBijectionTypeA(krt_elt)
    elif ct.letter == 'D':
        return KRTToRCBijectionTypeD(krt_elt)
    else:
        raise NotImplementedError

def RCToKRTBijection(rigged_configuration_elt):
    r"""
    Return the correct rigged configuration to KR tableaux bijection helper class.

    TESTS::

        sage: RC = RiggedConfigurations(['A', 4, 1], [[2, 1]])
        sage: from sage.combinat.rigged_configurations.bijection import RCToKRTBijection
        sage: bijection = RCToKRTBijection(RC(partition_list=[[1],[1],[1],[1]]))
        sage: bijection.tj(1)
        sage: bijection.next_state()
        5
        sage: bijection.cur_partitions
        [(/)
        , (/)
        , (/)
        , (/)
        ]
    """
    ct = rigged_configuration_elt.cartan_type()
    if ct.letter == 'A':
        return RCToKRTBijectionTypeA(rigged_configuration_elt)
    elif ct.letter == 'D':
        return RCToKRTBijectionTypeD(rigged_configuration_elt)
    else:
        raise NotImplementedError
