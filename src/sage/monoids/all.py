from monoid import is_Monoid

from free_monoid import FreeMonoid, is_FreeMonoid
from string_monoid import (BinaryStrings, OctalStrings, HexadecimalStrings,
                           Radix64Strings, AlphabeticStrings)

from free_abelian_monoid import FreeAbelianMonoid, is_FreeAbelianMonoid

from free_monoid_element import is_FreeMonoidElement
from string_ops import (
    strip_encoding,
    frequency_distribution,
    coincidence_index,
    coincidence_discriminant)

from string_monoid_element import (
    is_BinaryStringMonoidElement,
    is_OctalStringMonoidElement,
    is_HexadecimalStringMonoidElement,
    is_Radix64StringMonoidElement)

from free_abelian_monoid_element import is_FreeAbelianMonoidElement

