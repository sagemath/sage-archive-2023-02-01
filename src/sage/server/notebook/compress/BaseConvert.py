# BaseConvert class,
#    converts an unsigned base 10 integer to a different base and vice versa.
# ______________________________________________________________
# BaseConvert
#    |
#    |________ constructor(newBase:string)
#    |             uses newBase string var for convertion
#    |                 [i.e. "0123456789abcdef" for an hex convertion]
#    |
#    |________ toBase(unsignedInteger:uint):string
#    |             return base value of input
#    |
#    |________ fromBase(baseString:string):uint
#                  return base 10 integer value of base input
# --------------------------------------------------------------
# @Compatibility	>= Python 2.3
# @Author		Andrea Giammarchi
# @Site			http://www.devpro.it/
# @Date			2006/06/05
# @Version		1.0
# @License              GNU General Public License (GPL)

class BaseConvert:

	__base = ""
	__baseLength = 0

	def __init__(self, __base):
		self.__base = __base
		self.__baseLength = len(__base)

	def toBase(self, num):
		module = 0
		result = ""
		while num != 0:
			module = num % self.__baseLength
			result = self.__base[module] + result
			num = int((num - module) / self.__baseLength)
		if result == "":
			result = self.__base[0]
		return result


	def fromBase(self, str):
		pos = 0
		strLen = len(str) - 1
		result = 0
		while pos < strLen:
			result = result + (pow(self.__baseLength, strLen - pos) * self.__base.find(str[pos]))
			pos = pos + 1
		if strLen >= 0:
			result = result + self.__base.find(str[pos])
		return result;
