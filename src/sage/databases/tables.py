"""
Access to external tables of modular forms and elliptic curves
"""

#*****************************************************************************
#
#      Sage: Copyright (C) 2004 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import os

import sage.misc.db as db
import sage.misc.misc as misc

PATH = "%s/src/tables/"%misc.SAGE_ROOT

class __Cremona:
    def __init__(self):
        self.__path = "%s/cremona/"%PATH

    def __repr__(self):
        return "Cremona's tables of elliptic curves."

    def allcurves(self):
        """
        Returns the entire allcurves tables from levels 1 through 25000
        as a single Python dictionary , where allcurves[N] is a list
        of the data for level N as lists:
             [isogeny_class (a letter), curve number,
              [a1,a2,a3,a4,a6],  rank, torsion order]
        """
        try:
            return self.__allcurves
        except AttributeError:
            # not defined yet, so load it and return
            self.__allcurves = {}
            # first try to load from previously saved data
            file = "%s/allcurves.pickle"%self.__path
            if os.path.exists(file):
                misc.verbose("Loading Cremona database.")
                self.__allcurves = db.load(file)
                return self.__allcurves
            # If previously saved version is missing, load and parse
            # from whatever cremona data files are in the directory.
            misc.verbose("Parsing Cremona allcurves files and saving as database.")
            for F in os.listdir(self.__path):
                if F[:9] != "allcurves":
                    continue
                name = "%s/%s"%(self.__path, F)
                misc.verbose("Parsing %s..."%name)
                for line in open(name).readlines():
                    level, iso, num, ainvs, rank, tors = line.split()
                    level = int(level)
                    num = int(num)
                    ainvs = eval(ainvs)
                    rank = int(rank)
                    tors = int(tors)
                    if not self.__allcurves.has_key(level):
                        self.__allcurves[level] = []
                    self.__allcurves[level].append([level, iso, num, ainvs, rank, tors])
            db.save(self.__allcurves, file)
        return self.__allcurves

    def curves(self, level):
        """
        Return the data for level N as a list of lists
             [isogeny_class (a letter), curve number,
              [a1,a2,a3,a4,a6],  rank, torsion order]
        """
        level = int(level)
        #assert isinstance(level, int)
        if level <= 10:
            return []
        try:
            # if we happen to have the full database in memory, use it.
            return self.__allcurves[level]
        except AttributeError:
            # load directly from whatever cremona data files are in the directory.
            misc.verbose("Parsing Cremona allcurves file.")
            ans = []
            for F in os.listdir(self.__path):
                if F[:9] != "allcurves":
                    continue
                name = "%s/%s"%(self.__path, F)
                file = open(name).read()
                i = file.find("\n%s "%level)
                if i == -1:
                    continue
                for line in file[i:].split("\n"):
                    X = line.split()
                    if len(X) < 6:
                        continue
                    N, iso, num, ainvs, rank, tors = line.split()
                    N = int(N)
                    if N != level:
                        return ans
                    num = int(num)
                    ainvs = eval(ainvs)
                    rank = int(rank)
                    tors = int(tors)
                    ans.append([N, iso, num, ainvs, rank, tors])
                return ans
            return []

    def curve(self, level, iso="A", num=1):
        if isinstance(level,str):
            label = level
            i=0
            while i<len(label) and label[i]>='0' and label[i]<='9':
                i += 1
            j=i+1
            if j>len(label):
                label += "A"
            while j<len(label) and not (label[j]>='0' and label[j]<='9'):
                j += 1
            if j>=len(label):
                label += "1"
            level, iso, num = int(label[:i]), label[i:j], int(label[j:])
        C = self.curves(level)
        for x in C:
            if x[1]==iso.upper() and x[2]==num:
                return x
        raise RuntimeError, "There is no elliptic curve %s%s%s"%(level,iso,num)

Cremona = __Cremona()

def zeta_zeros():
    """
    List of the first 10000 imaginary parts of nontrivial zeros of the Riemann
    zeta function.

    NOTE: Andrew Odlyzko computed these to precision within 3*10^(-9).  See
              http://www.dtc.umn.edu/~odlyzko/zeta_tables/
    """
    path = "%s/odlyzko"%PATH
    file = "%s/zeros1"%path
    if os.path.exists(file+".pickle"):
        misc.verbose("Loading Odlyzko database from " + file + ".pickle")
        return db.load(file+".pickle")
    misc.verbose("Creating Odlyzko Database.")
    F = [eval(x) for x in open(file).read().split()]
    db.save(F, file+".pickle")
    return F

