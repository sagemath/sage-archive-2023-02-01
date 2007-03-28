from __builtin__ import min as python_min
from __builtin__ import max as python_max
from sage.rings.infinity import infinity

def min(*L):
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_min(L)
    except ValueError:
        return infinity

def max(*L):
    if len(L) == 1 and isinstance(L[0], (list, tuple)):
        L = L[0]
    try:
        return python_max(L)
    except ValueError:
        return -infinity
