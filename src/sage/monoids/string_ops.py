#*****************************************************************************
#       Copyright (C) 2007 David Kohel <kohel@maths.usyd.edu.au>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

from sage.misc.misc import prod
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

def coincidence_index(S,n=1):
    """
    The coincidence index of the string S.
    EXAMPLES:
	sage: S = strip_encoding("The cat in the hat.")
	sage: coincidence_index(S)
	0.120879120879120
    """
    S = strip_encoding(S)
    N = len(S)-n+1
    X = {}
    for i in range(N):
        c = S[i:i+n]
        if X.has_key(c):
            X[c] += 1
        else:
            X[c] = 1
    RR = RealField()
    return RR(sum([ m*(m-1) for m in X.values() ]))/RR(N*(N-1))

def coincidence_discriminant(S,n=2):
    """
    Input: A string, e.g. produced as decimation of transposition ciphertext,
    or a sample plaintext.
    Output: A measure of the difference of probability of association of
    character pairs, relative to their independent one-character probabilities.

    EXAMPLES:
	sage: S = strip_encoding("The cat in the hat.")
	sage: coincidence_discriminant(S)
        0.0827001855677322
    """
    if n != 2:
        raise ValueError, "Argument n (= %s) is only implemented for n = 2" % n
    S = strip_encoding(S)
    X1 = [ frequency_distribution(S[i:len(S)-n+i+1]) for i in range(n) ]
    XX = frequency_distribution(S,n)
    AZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    return sum([ (XX(AZ[i]+AZ[j])-X1[0](AZ[i])*X1[1](AZ[j]))**2 for i in range(26) for j in range(26) ])
