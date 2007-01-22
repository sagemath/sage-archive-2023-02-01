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

## def coincidence_discriminant(S):
##     """

##     INPUT:
##         S -- A sequence of 2-character strings, e.g. produced as
##              decimation of transposition ciphertext, or of adjacent
##              characters in some sample plaintext.
##     OUTPUT:
##         A measure of the difference of probability of association of
##         two characters, relative to their independent probabilities.

##     EXAMPLES:
## 	sage: S = strip_encoding("The cat in the hat.")
## 	sage: T = [ S[i:i+2] for i in range(len(S)-1) ]
## 	sage: coincidence_discriminant(T)
##         0.0300925925925925
##     """
##     AZ = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
##     AA = [ AZ[i] + AZ[j] for i in range(26) for j in range(26) ]
##     X1 = frequency_distribution(''.join([ s[0] for s in S ]))
##     X2 = frequency_distribution(''.join([ s[1] for s in S ]))
##     F2 = {}
##     RR = RealField()
##     for XY in AA:
##         F2[XY] = RR(0)
##     eps = RR(1/len(S))
##     for AB in S:
##         F2[AB] += eps

## def translation_correlations(X1,X2):
##      """
##      The sequence of correlations of the sequence X1 with the cyclic
##      translations of the sequence X2.
##      """
##      n = len(X1)
##      if n != len(X2):
##          raise ValueError, "Arguments must be of the same length."

##      # Compute the mean value of each sequence:
##      mu1 = sum(X1)/n
##      mu2 = sum(X2)/n

##      # Compute the standard deviations of each sequence:
##      sig1 = sqrt(sum([ (S1[k]-mu1)^2 for k in range(n) ]))
##      sig2 = sqrt(sum([ (S2[k]-mu2)^2 for k in range(n) ]))

##      sig = sig1*sig2
##      corr_dict = { }
##      for j in range(n):
##          corr_dict[j] = sum([
##              (S1[i] - mu1) * (S2[(i+j)%n] - mu2) / sig for i in range(n) ])
##      return corr_dict

## def translation_matches(S,X,r):
##     """
##     Input
##        S : Test string.
##        X : Sequence of standard frequencies for the language.
##        r : A real number between 0 and 1.
##     Output
##        A list of integers k such that affine translation of S by k
##        has correlation at least r with the standard frequencies
##        given by the real sequence X.
##     """
##     Y = S.frequency_distribution()
##     corr_dict = translation_correlations(X,Y)
##     I = []
##     for i in keys(X):
##         if x[2] > r:
##             I.append(x[1])
##     return I
