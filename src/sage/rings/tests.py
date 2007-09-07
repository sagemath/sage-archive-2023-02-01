"""
TESTS:
    sage: K.<x>=FractionField(QQ['x'])
    sage: V.<z> = K[]
    sage: x+z
    z + x

    sage: (1/2)^(2^100)
    Traceback (most recent call last):
    ...
    RuntimeError: exponent must be at most 9223372036854775807     # 64-bit
    RuntimeError: exponent must be at most 2147483647              # 32-bit
"""
