r"""
Iterable of the keys of a Mapping associated with nonzero values
"""

from collections.abc import MappingView, Sequence, Set


class SupportView(MappingView, Sequence, Set):

    def __init__(self, mapping, *, zero=None):
        self._mapping = mapping
        self._zero = zero

    def __len__(self):
        length = 0
        for key in self:
            length += 1
        return length

    def __getitem__(self, index):
        for i, key in enumerate(self):
            if i == index:
                return key
        raise IndexError

    def __iter__(self):
        zero = self._zero
        if zero is None:
            for key, value in self._mapping.items():
                if value:
                    yield key
        else:
            for key, value in self._mapping.items():
                if value != zero:
                    yield key

    def __contains__(self, key):
        try:
            value = self._mapping[key]
        except KeyError:
            return False
        if zero is None:
            return bool(value)
        return value != zero
