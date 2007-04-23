# JavaScriptCompressor class,
#	removes comments or pack JavaScript source[s] code.
# ______________________________________________________________
# JavaScriptCompressor (just 2 public methods)
#    |
#    |________ self.getClean(jsSource:mixed):string
#    |         	returns one or more JavaScript code without comments,
#    |         	by default removes some spaces too
#    |
#    |________ self.getPacked(jsSource:mixed):string
#              	returns one or more JavaScript code packed,
#	        	using getClean and obfuscating output
# --------------------------------------------------------------
# Note about jsSource input varible:
# 	this var should be a string (i.e. jsSource = "source file string")
#      should be a list of strings (i.e. ["source file string 1", "source file string 2"])
#      should be a dictionary with 1 or 2 keys:
#      	(i.e. {"code":"source file string 1"})
#              (i.e. {"code":"source file string 1", "name":"mySource"})
#      ... and should be a list of dictionaries created with theese rules
#      [
#		"source javascript 1",
#              {"code":"source file string 2"},
#              {"code":"source file string 3", "name":"JSApplication 1.0"}
#      ]
#
#      The name used on dedicated key, will be write on parsed source header
# --------------------------------------------------------------
# Note about returned strings:
# 	Your browser should wrap very long strings, then don't use
#      cut and paste from your browser, save output into your database or directly
#      in a file or print them only inside <script> and </script> tags
# --------------------------------------------------------------
# Note about JavaScript packed compatibility:
# 	To be sure about compatibility include before every script JSL Library:
#      http://www.devpro.it/JSL/
# JSL library add some features for old or buggy browsers, one of
# those functions is String.replace with function as second argument,
# used by JavaScript generated packed code to rebuild original code.
#
# Remember that KDE 3.5, Safari and IE5 will not work correctly with packed version
# if you'll not include JSL.
# --------------------------------------------------------------
# @Compatibility	>= Python 2.3
# @Author		Andrea Giammarchi
# @Site			http://www.devpro.it/
# @Date			2006/08/02
# @LastMOd		2006/08/03 [add __getSize method to have stats like PHP version]
# @Version		0.1b
# @Dependencies		Python: BaseConvert.py
#			Python: SourceMap.py	[ ... I'm waiting that this site will approve SourceMap class ... ]
#			Client: JSL.js (http://www.devpro.it/JSL/)
# @Browsers		Convertion is supported by every browser with JSL Library (FF 1+ Opera 8+ and IE5.5+ are supported without JSL too)
# @Credits		Dean Edwards for his originally idea [dean.edwards.name] and his JavaScript packer
# @License		GNU General Public License (GPL)

import re, time, string, SourceMap, BaseConvert
class JavaScriptCompressor:

	# public variables
        # 	stats:string		after every compression has some informations
        #      version:string		version of this class
	stats = "",
	version = "0.1";

	# private variables, any comment sorry
	__startTime = 0
	__sourceLength = 0
	__sourceNewLength = 0
	__totalSources = 0
	__sources = []
	__delimeter = [{"name":"doublequote", "start":'"', "end":'"', "noslash":True}, {"name":"singlequote", "start":"'", "end":"'", "noslash":True}, {"name":"singlelinecomment", "start":"//", "end":["\n", "\r"]}, {"name":"multilinecomment", "start":"/*", "end":"*/"}, {"name":"regexp", "start":"/", "end":"/", "match":"^/[^\n\r]+/$", "noslash":True}]
	__cleanFinder = ["(\n|\r)+", "( |\t)+", "(\n )|( \n)|( \n )", "\s+(\)|})", "(\(|{)\s+", "\s*(;|,|:|<|>|\&|\||\=|\?|\+|\-|\%)\s*", "\)\s+{", "}\s+\("]
	__cleanReplacer = ["\n", " ", "\n", "\\1", "\\1", "\\1", "){", "}("]
	__container = []
	__BC = BaseConvert.BaseConvert("0123456789abcdefghijklmnopqrstuvwxyz")
	__SourceMap = SourceMap.SourceMap()

	# public method
        # 	string self.getClean(jsSource:mixed)
        #      compress JavaScript removing comments and somespaces
        # @param	mixed		view example and notes on class comments
	def getClean(self, jsSource):
		return self.__commonInitMethods(jsSource, False)

	# public method
        # 	string self.getClean(jsSource:mixed)
        #      compress JavaScript replaceing words and removing comments and some spaces
        # @param	mixed		view example and notes on class comments
	def getPacked(self, jsSource):
		return self.__commonInitMethods(jsSource, True)

	# private methods, any comment sorry
	def __addCleanCode(self, str):
		for a in range(0, len(self.__cleanFinder)):
			str = re.sub(self.__cleanFinder[a], self.__cleanReplacer[a], str)
		return str
	def __addslashes(self, str):
		return str.replace("\\", "\\\\").replace("'", "\'").replace("\"", "\\\"")
	def __clean(self, str):
		type = ""
		clean = []
		map = self.__SourceMap.getMap(str, self.__delimeter)
		for a in range(0, len(map)):
			type = map[a]["name"]
			if type == "code":
				clean.append(self.__addCleanCode(str[map[a]["start"]:map[a]["end"]]))
			elif type == "regexp" or type == "doublequote" or type == "singlequote":
				clean.append(str[map[a]["start"]:map[a]["end"]])
			if type != "regexp":
				clean.append("\n")
		return re.sub("/(\n)+/", "\n", re.sub("/^\s*|\s*$/", "", string.join(clean, "")))
	def __commonInitMethods(self, jsSource, packed):
		header = ""
		sources = ""
		tempSources = []
		self.__startTime = time.clock()
		self.__sourceLength = 0
		self.__sourceManager(jsSource)
		for a in range(0, self.__totalSources):
			self.__sources[a]["code"] = self.__clean(self.__sources[a]["code"])
		header = self.__getHeader()
		for a in range(0, self.__totalSources):
			tempSources.append(self.__sources[a]["code"])
		sources = string.join(tempSources, ";")
		if packed:
			sources = self.__pack(sources)
		self.__sourceNewLength = len(sources)
		self.__setStats()
		return header + sources
	def __doReplacement(self, matchobj):
		return self.__BC.toBase(self.__wordsParser(matchobj.group(0)))
	def __getHeader(self):
		return string.join([
			"/* ", self.__getScriptNames(), "JavaScriptCompressor ", self.version, " [www.devpro.it], ",
			"thanks to Dean Edwards for idea [dean.edwards.name]",
			" */\n"
		], "")
	def __getScriptNames(self):
		a = 0
		result = []
		for a in range(0, self.__totalSources):
			if self.__sources[a]["name"] != "":
				result.append(self.__sources[a]["name"])
		a = len(result)
		if a > 0:
			a = a - 1
			result[a] = result[a] + " with ";
		return string.join(result, ", ")
	def __getSize(self, size):
		times = 0
		fsize = float(size)
		sizeType = ["bytes", "Kb", "Mb", "Gb", "Tb", "Zb"]
		sizeTypeLen = len(sizeType)
		resultSize = ""
		while fsize > 1024 and times < sizeTypeLen:
			fsize = fsize / 1024
			times = times + 1
		if times > 0:
			resultSize = "%.2f" % (fsize)
		else:
			resultSize = str(size)
		return resultSize + " " + sizeType[times]
	def __pack(self, str):
		self.__container = []
		str = self.__addslashes(re.sub("\w+", self.__doReplacement, self.__clean(str))).replace("\n", "\\n")
		return 'eval(function(A,G){return A.replace(/(\\w+)/g,function(a,b){return G[parseInt(b,36)]})}("' + str + '","' + string.join(self.__container, ",") + '".split(",")));'
	def __setStats(self):
		endTime = "%.3f" % ((time.clock() - self.__startTime) / 1000)
		self.stats = string.join([
			self.__getSize(self.__sourceLength),
			"to",
			self.__getSize(self.__sourceNewLength),
			"in",
			endTime,
			"seconds"
		], " ")
	def __sourceManager(self, jsSource):
		b = len(jsSource)
		dictType = type({})
		striType = type("")
		self.__sources = []
		if type(jsSource) == striType:
			self.__sourcePusher(jsSource, "")
		elif type(jsSource) == dictType:
			self.__sourcePusher(jsSource.get("code", ""), jsSource.get("name", ""))
		elif type(jsSource) == type([]) and b > 0:
			for a in range(0, b):
				if type(jsSource[a]) == dictType and jsSource.get("code", None) != None and jsSource.get("name", None) != None:
					self.__sourcePusher(jsSource.get("code", None), jsSource.get("name", None))
				elif type(jsSource[a]) == striType:
					self.__sourcePusher(jsSource[a], "")
		self.__totalSources = len(self.__sources)
	def __sourcePusher(self, code, name):
		self.__sourceLength = self.__sourceLength + len(code)
		self.__sources.append({"code":code, "name":name})
	def __wordsParser(self, str):
		key = 0
		try:
			key = self.__container.index(str)
		except:
			key = len(self.__container)
			self.__container.append(str)
		return key
