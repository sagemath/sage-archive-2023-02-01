#!/usr/bin/env python
import os
for x in zip(open(os.sys.argv[1]).readlines(),open(os.sys.argv[2]).readlines()):
    if x[0]!=x[1]:
        exit(1)
