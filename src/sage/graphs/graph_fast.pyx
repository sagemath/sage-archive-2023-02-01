"""
Graph Theory SageX functions

AUTHOR:
    -- Robert L. Miller (2007-02-13): initial version
"""

#*****************************************************************************
#           Copyright (C) 2007 Robert L. Miller <rlmillster@gmail.com>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

include '../ext/cdefs.pxi'

def binary(n, length=None):
    """
    A quick python int to binary string conversion.

    sage.: timeit sage.graphs.graph_fast.binary(389)
    100000 loops, best of 3: 11.4 [micro]s per loop
    sage.: timeit Integer(389).binary()
    10000 loops, best of 3: 16.8 [micro]s per loop

    EXAMPLE:
    sage: sage.graphs.graph_fast.binary(2007)
    '11111010111'
    """
    cdef mpz_t i
    mpz_init(i)
    mpz_set_ui(i,n)
    return mpz_get_str(NULL, 2, i)

def R(x):
    """
    A helper function for the graph6 format. Described in [McK]

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    # pad on the right to make a multiple of 6
    x = x + ( '0' * ((6 - len(x))%6) )

    # split into groups of 6, and convert numbers to decimal, adding 63
    six_bits = ''
    cdef int i
    for i from 0 <= i < len(x)/6:
        six_bits += chr( int( x[6*i:6*(i+1)], 2) + 63 )
    return six_bits

def N(n):
    """
    A helper function for the graph6 format. Described in [McK]

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    if n < 63:
        return chr(n + 63)
    else:
        # get 18-bit rep of n
        n = binary(n)
        n = '0'*(18-len(n)) + n
        return chr(126) + R(n)

def N_inverse(s):
    """
    A helper function for the graph6 format. Described in [McK]

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    if s[0] == chr(126): # first four bytes are N
        a = binary(ord(s[1]) - 63)
        b = binary(ord(s[2]) - 63)
        c = binary(ord(s[3]) - 63)
        n = int(a + b + c,2)
        s = s[4:]
    else: # only first byte is N
        n = ord(s[0]) - 63
        s = s[1:]
    return n, s

def R_inverse(s, n):
    """
    A helper function for the graph6 format. Described in [McK]

    REFERENCES:
    McKay, Brendan. 'Description of graph6 and sparse6 encodings.'
    http://cs.anu.edu.au/~bdm/data/formats.txt (2007-02-13)
    """
    l = []
    cdef int i
    for i from 0 <= i < len(s):
        a = binary(ord(s[i])-63)
        l.append( '0'*(6-len(a)) + a )
    m = "".join(l)
    return m[:(n*(n-1)/2)]










