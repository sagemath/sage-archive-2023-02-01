# SourceMap class,
#	reads a generic language source code and returns its map.
# ______________________________________________________________
# The SourceMap goals is to create a map of a generic script/program language.
# The getMap method returns an array/list of arrays/dictionary/objects
# of source map using delimeters variable to map correctly:
#  - multi line comments
#  - single line comments
#  - double quoted strings
#  - single quoted strings
#  - pure code
#  - everything else (for example regexp [/re/] with javascript), just adding a correct delimeter
# --------------------------------------------------------------
# What about the delimeter
# 	It's an array/list of arrays/dictionary/obects with some properties to find what you're looking for.
#
# parameters are:
#  - name, the name of the delimeter (i.e. "doublequote")
#  - start, one or mode chars to find as start delimeter (i.e. " for double quoted string)
#  - end, one or mode chars to find as end delimeter (i.e. " for double quoted string) [end should be an array/list too]
#
# optional parameters are:
#  - noslash, if true find the end of the delimeter only if last char is not slashed (i.e. "string\"test" find " after test)
#  - match, if choosed language has regexp, verify if string from start to end matches used regexp (i.e. /^\/[^\n\r]+\/$/ for JavaScript regexp)
#
# If end parameter is an array, match and noslash are not supported (i.e. ["\n", "\r"] for end delimeter of a single line comment)
# --------------------------------------------------------------
# What about SourceMap usage
# 	It should be a good solution to create sintax highlighter, parser,
# 	verifier or some other source code parsing procedure
# --------------------------------------------------------------
# What about SourceMap performance script/languages
# 	I've created different version of this class to test each script/program language performance too.
# Python with or without Psyco is actually the faster parser.
# --------------------------------------------------------------
# @Compatibility	>= PHP 4
# @Author		Andrea Giammarchi
# @Site		http://www.devpro.it/
# @Date		2006/08/01
# @LastMOd		2006/08/01
# @Version		0.1
# @Application		Last version of JavaScriptCompressor class use this one to map source code.
# @License              GNU General Public License (GPL)

import re
class SourceMap:

	# public method
        # 	getMap(&$source:string, &$delimeters:array):array
	# Maps the source code using $delimeters rules and returns map as an array
        # NOTE: read comments to know more about map and delimeter
        #
        # @param	string		generic source code
        # @param	array		array with nested array with code rules

	def getMap(self, source, delimeters):

		# "unsigned" integer variables
		sourcePosition = 0
		delimetersPosition = 0
		findLength = 0
		templen = 0
		tempIndex = 0
		sourceLength = len(source)
		delimetersLength = len(delimeters)

		# integer variables
		tempPosition = -1
		endPosition = -1

		# list variables
		codeMap = []
		tempMap = []

		# dictionary variable
		tempDelimeter = {}

		while sourcePosition < sourceLength:
			endPosition = -1
			for delimetersPosition in range(0, delimetersLength):
				tempPosition = source.find(delimeters[delimetersPosition]["start"], sourcePosition)
				if tempPosition != -1 and (tempPosition < endPosition or endPosition == -1):
					endPosition = tempPosition
					tempIndex = delimetersPosition
			if endPosition != -1:
				sourcePosition = endPosition
				tempDelimeter = delimeters[tempIndex]
				findLength = len(tempDelimeter["start"])
				if type(tempDelimeter["end"]) == type([]):
					endPosition = -1
					for delimetersPosition in range(0, len(tempDelimeter["end"])):
						tempPosition = source.find(tempDelimeter["end"][delimetersPosition], sourcePosition + findLength)
						if tempPosition != -1 and (tempPosition < endPosition or endPosition == -1):
							endPosition = tempPosition
							tempIndex = delimetersPosition
					if endPosition != -1:
						endPosition = endPosition + len(tempDelimeter["end"][tempIndex])
					else:
						endPosition = sourceLength
					codeMap.append({"name":tempDelimeter["name"], "start":sourcePosition, "end":endPosition})
					sourcePosition = endPosition - 1
				elif self.__has(tempDelimeter, "match"):
					tempPosition = source.find(tempDelimeter["end"], sourcePosition + findLength)
					templen = len(tempDelimeter["end"])
					if tempPosition != -1 and re.match(tempDelimeter["match"], source[sourcePosition:tempPosition+templen]) != None:
						if self.__has(tempDelimeter, "noslash") and tempDelimeter["noslash"] == True:
							endPosition = self.__endCharNoSlash(source, sourcePosition, tempDelimeter["end"], sourceLength)
						else:
							endPosition = tempPosition + len
						codeMap.append({"name":tempDelimeter["name"], "start":sourcePosition, "end":endPosition})
						sourcePosition = endPosition - 1
				else:
					if self.__has(tempDelimeter, "noslash") and tempDelimeter["noslash"] == True:
						endPosition = self.__endCharNoSlash(source, sourcePosition, tempDelimeter["end"], sourceLength)
					else:
						tempPosition = source.find(tempDelimeter["end"], sourcePosition + findLength)
						if tempPosition != -1:
							endPosition = tempPosition + len(tempDelimeter["end"])
						else:
							endPosition = sourceLength
					codeMap.append({"name":tempDelimeter["name"], "start":sourcePosition, "end":endPosition})
					sourcePosition = endPosition - 1
			else:
				sourcePosition = sourceLength - 1
			sourcePosition = sourcePosition + 1
		templen = len(codeMap)
		if templen == 0:
			tempMap.append({"name":"code", "start":0, "end":sourceLength})
		else:
			for tempIndex in range(0, templen):
				if tempIndex == 0 and codeMap[tempIndex]["start"] > 0:
					tempMap.append({"name":"code", "start":0, "end":codeMap[tempIndex]["start"]});
				elif tempIndex > 0 and codeMap[tempIndex]["start"] > codeMap[tempIndex-1]["end"]:
					tempMap.append({"name":"code", "start":codeMap[tempIndex-1]["end"], "end":codeMap[tempIndex]["start"]});
				tempMap.append({"name":codeMap[tempIndex]["name"], "start":codeMap[tempIndex]["start"], "end":codeMap[tempIndex]["end"]});
				if tempIndex + 1 == templen and codeMap[tempIndex]["end"] < sourceLength:
					tempMap.append({"name":"code", "start":codeMap[tempIndex]["end"], "end":sourceLength});
		return tempMap

	def __has(self, dict, name):
		return dict.get(name, None) != None

	def __endCharNoSlash(self, source, position, find, sourceLen):
		loop = True
		temp = len(find)
		while loop:
			position = source.find(find, position + 1)
			if not (position != -1 and not self.__charNoSlash(source, position)):
				loop = False
		if position == -1:
			position = sourceLen - temp
		return position + temp

	def __charNoSlash(self, source, position):
		next = 1
		sourceLen = position - next
		while sourceLen > 0 and source[sourceLen] == '\\':
			next = next + 1
			sourceLen = position - next
		return ((next - 1) % 2 == 0)
