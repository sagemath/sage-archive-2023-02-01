
#*****************************************************************************
#      Copyright (C) 2006 - 2008 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../../../ext/cdefs.pxi'
include '../../../ext/stdsage.pxi'
include '../../../misc/bitset_pxd.pxi'

cdef struct OrbitPartition:
    # Disjoint set data structure for representing the orbits of the generators
    # so far found. Also keeps track of the minimum elements of the cells and
    # their sizes.
    int degree
    int *parent
    int *rank
    int *mcr # minimum cell representatives - only valid at the root of a cell
    int *size # also only valid at the root of a cell

cdef struct PartitionStack:
    # Representation of a node of the search tree. A sequence of partitions of
    # length depth + 1, each of which is finer than the last. Partition k is
    # represented as PS.entries in order, broken up immediately after each
    # entry of levels which is at most k.
    int *entries
    int *levels
    int depth
    int degree


