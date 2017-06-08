from __future__ import absolute_import
from sage.misc.lazy_import import lazy_import

lazy_import('sage.crypto.classical', 'AffineCryptosystem')
lazy_import('sage.crypto.classical', 'HillCryptosystem')
lazy_import('sage.crypto.classical', 'SubstitutionCryptosystem')
lazy_import('sage.crypto.classical', 'ShiftCryptosystem')
lazy_import('sage.crypto.classical', 'TranspositionCryptosystem')
lazy_import('sage.crypto.classical', 'VigenereCryptosystem')

lazy_import('sage.crypto.stream', 'LFSRCryptosystem')
lazy_import('sage.crypto.stream', 'ShrinkingGeneratorCryptosystem')

lazy_import('sage.crypto.lfsr', 'lfsr_sequence')
lazy_import('sage.crypto.lfsr', 'lfsr_autocorrelation')
lazy_import('sage.crypto.lfsr', 'lfsr_connection_polynomial')
