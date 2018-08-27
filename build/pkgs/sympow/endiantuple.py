import sys
import platform

if sys.byteorder == "little":
    endian = "le"
else:
    endian = "be"

bits = platform.architecture()[0][:-3]

# for example "le64" or "be32"
print(endian + bits)
