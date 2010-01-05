#!/usr/bin/env python
# as suggested by 'javier' on sage-devel.

import sys
b = 1
x = sys.maxint

while x:
    x = x >> 1
    b = b+1

print(b)

