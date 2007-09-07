"""
TESTS:
    sage: K.<x>=FractionField(QQ['x'])
    sage: V.<z> = K[]
    sage: x+z
    z + x

    sage: (1/2)^(2^100)
    Traceback (most recent call last):
    ...
    RuntimeError: exponent must be at most 18446744073709551614    # 64-bit
    RuntimeError: exponent must be at most 4294967294              # 32-bit
"""
