import gen as _gen

pari = _gen.pari

pari_gen = _gen.gen

allocatemem = pari.allocatemem

PariError = _gen.PariError

from gen import (prec_dec_to_words, prec_dec_to_bits,
                 prec_bits_to_words, prec_bits_to_dec,
                 prec_words_to_bits, prec_words_to_dec)
