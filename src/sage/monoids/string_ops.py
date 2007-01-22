#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.rings.integer import Integer
from sage.rings.real_mpfr import RealField
from sage.rings.rational_field import RationalField
from sage.probability.random_variable import DiscreteProbabilitySpace

def strip_encoding(S):
    """
    The upper case string of S stripped of all non-alphabetic characters.
    EXAMPLES:
        sage: S = "The cat in the hat."
        sage: strip_encoding(S)
	'THECATINTHEHAT'
    """
    X = ''
    for i in range(len(S)):
	C = S[i]
	if C.isalpha():
	    X += S[i].upper()
    return X

def frequency_distribution(S, n=1, field=None):
    """
    The probability space of frequencies of n-character substrings of S.
    """
    N = len(S)-n+1
    if field is None:
        field = RealField()
    P = {}
    alph = []
    eps = field(1)/N
    for i in range(N):
        c = S[i:i+n]
        if P.has_key(c):
	    P[c] += eps
        else:
	    P[c] = eps
	    alph.append(c)
    return DiscreteProbabilitySpace(alph,P,field)

def coincidence_index(S):
    """
    The coincidence index of the string S.
    EXAMPLES:
	sage: S = strip_encoding("The cat in the hat.")
	sage: coincidence_index(S)
	0.120879120879120
    """
    S = strip_encoding(S)
    n = len(S)
    X = [ 0 for i in range(26) ]
    AZ = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
    for i in range(n):
        X[AZ.index(S[i])] += 1
    RR = RealField()
    return sum([ RR(m*(m-1)) for m in X ])/RR(n*(n-1))

def coincidence_discriminant(S):
    """
    Input
        A sequence of 2-character strings, e.g. produced as decimation
        of transposition ciphertext, or of adjacent characters in some
        sample plaintext.
    Output
        A measure of the difference of probability of association of
        two characters, relative to their independent probabilities.
    EXAMPLES:
	sage: S = strip_encoding("The cat in the hat.")
	sage: T = [ S[i:i+2] for i in range(len(S)-1) ]
	sage: coincidence_discriminant(T)
        0.0300925925925925
    """
    AZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    AA = [ AZ[i] + AZ[j] for i in range(26) for j in range(26) ]
    X1 = frequency_distribution(''.join([ s[0] for s in S ]))
    X2 = frequency_distribution(''.join([ s[1] for s in S ]))
    F2 = {}
    RR = RealField()
    for XY in AA:
        F2[XY] = RR(0)
    eps = RR(1/len(S))
    for AB in S:
        F2[AB] += eps
    return sum([ (F2[AZ[i]+AZ[j]]-X1[AZ[i]]*X2[AZ[j]])**2 for i in range(26) for j in range(26) ])
