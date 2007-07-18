#*****************************************************************************
#      Copyright (C) 2007 Michael W. Hansen <mwhansen@odin.ac.hmc.edu>
#
# Distributed  under  the  terms  of  the  GNU  General  Public  License (GPL)
#                         http://www.gnu.org/licenses/
#*****************************************************************************

import networkx as NX
import sys

class GenerateLinearExtensionsStruct:
   def __init__(self,le,a,b,dag):
       self.le = le[:]
       self.a  = a[:]
       self.b  = b[:]
       self.dag = dag
       self.mrb = 0
       self.mra = 0
       self.IsPlus = True
       self.linearExtensions = []

def Switch(i, gles):
   if i == -1:
       gles.IsPlus = not gles.IsPlus
   if i >= 0:
       aIndex = gles.le.index(gles.a[i])
       bIndex = gles.le.index(gles.b[i])
       gles.le[aIndex] = gles.b[i]
       gles.le[bIndex] = gles.a[i]

       temp = gles.b[i]
       gles.b[i] = gles.a[i]
       gles.a[i] = temp

   if gles.IsPlus:
       gles.linearExtensions.append(gles.le[:])


def Move(element, gles, direction):
   index = gles.le.index(element)
   if direction == "right":
       gles.le[index] = gles.le[index+1]
       gles.le[index+1] = element
   elif direction == "left":
       gles.le[index] = gles.le[index-1]
       gles.le[index-1] = element
   else:
       print "Bad direction!"
       sys.exit()
   if gles.IsPlus:
       gles.linearExtensions.append(gles.le[:])


def incomparable(x, y, dag):
   if (not NX.path.shortest_path(dag, x, y)) and (not
NX.path.shortest_path(dag, y, x)):
       return True
   return False


def Right(i, gles, letter):
   if letter == "a":
       x = gles.a[i]
       yindex = gles.le.index(x) + 1
       if yindex >= len(gles.le):
           return False
       y = gles.le[ yindex ]
       if incomparable(x,y,gles.dag) and y != gles.b[i]:
           return True
       return False
   elif letter == "b":
       x = gles.b[i]
       yindex = gles.le.index(x) + 1
       if yindex >= len(gles.le):
           return False
       y = gles.le[ yindex ]
       if incomparable(x,y,gles.dag):
           return True
       return False
   else:
       print "Bad letter!"
       sys.exit()


def GenerateLinearExtensions(i, gles):
   if i >= 0:
       #print "GenerateLinearExtensions(%s):" % str(i)
       GenerateLinearExtensions(i-1, gles)
       mrb = 0
       typical = False
       while Right(i, gles, "b"):
           mrb += 1
           Move(gles.b[i], gles, "right")
           GenerateLinearExtensions(i-1, gles)
           mra = 0
           if Right(i, gles, "a"):
               typical = True
               cont = True
               while cont:
                   mra += 1
                   Move(gles.a[i], gles, "right")
                   GenerateLinearExtensions(i-1, gles)
                   cont = Right(i, gles, "a")
           if typical:
               Switch(i-1, gles)
               GenerateLinearExtensions(i-1, gles)
               if mrb % 2 == 1:
                   mla = mra -1
               else:
                   mla = mra + 1
               for x in range(mla):
                   Move(gles.a[i], gles, "left")
                   GenerateLinearExtensions(i-1, gles)

       if typical and (mrb % 2 == 1):
           Move(gles.a[i], gles, "left")
       else:
           Switch(i-1, gles)
       GenerateLinearExtensions(i-1, gles)
       for x in range(mrb):
           Move(gles.b[i], gles, "left")
           GenerateLinearExtensions(i-1, gles)


#Compute all linear extensions of a directed acyclic DiGraph
def linearExtensions(dag):
   ################
   #Precomputation#
   ################
   i = 0
   j = 0
   dagCopy = dag.copy()
   le = []
   a  = []
   b  = []
   while dagCopy.number_of_nodes() != 0:
       #Find all the minimal elements of dagCopy
       minimalElements = []
       for node in dagCopy.nodes():
           if len(dagCopy.in_edges(node)) == 0:
               minimalElements.append(node)
       if len(minimalElements) == 1:
           le.append(minimalElements[0])
           dagCopy.delete_node(minimalElements[0])
       else:
           ap = minimalElements[0]
           bp = minimalElements[1]
           a.append(ap)
           b.append(bp)
           le.append(ap)
           le.append(bp)
           dagCopy.delete_node(ap)
           dagCopy.delete_node(bp)
   maxPair = len(a) - 1

   gles = GenerateLinearExtensionsStruct(le, a, b, dag)
   gles.linearExtensions.append(le[:])
   GenerateLinearExtensions(maxPair , gles)
   Switch(maxPair, gles)
   GenerateLinearExtensions(maxPair , gles)
   gles.linearExtensions.sort()
   return gles.linearExtensions[:]
