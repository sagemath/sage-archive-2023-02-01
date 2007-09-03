"""nodoctest
Here we list bugs.  Enable doctesting of this file to find all kinds of problems.

Modular forms with character over a finite field just breaks:

   sage: m = ModularForms(DirichletGroup(8).1,2,GF(7)); m
   Modular Forms space of dimension 2, character [1, -1] and weight 2 over Finite Field of size 7
   sage: m.basis()   # this just goes into infinite loop (???)
   boom

"""
