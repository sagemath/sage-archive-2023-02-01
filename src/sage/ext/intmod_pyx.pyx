"""
Integers modulo n.
"""

def Intmod(x, m):
    if abs(m) < 50000:
        return Intmod_cint(x,m)
    raise NotImplementedError

cdef class Intmod_cint:
    cdef int x
    cdef int m

    def __init__(self, int value, int modulus):
        self.x = value % modulus
        self.m = modulus
        if self.x < 0:
            self.x = self.x + modulus

    def __add__(Intmod_cint self, Intmod_cint other):
        return Intmod_cint((self.x + other.x)%self.m, self.m)

    def __mul__(Intmod_cint self, Intmod_cint other):
        return Intmod_cint((self.x * other.x)%self.m, self.m)

    def __repr__(self):
        return str(self.x)

    def modulus(self):
        return self.m
