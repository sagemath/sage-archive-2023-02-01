def flatten(in_list, ltypes=(list, tuple)):
   """
   Flattens a nested list.

   INPUT:
       in_list -- a list or tuple
       ltypes -- optional list of particular types to flatten

   OUTPUT:
       a flat list of the entries of in_list

   EXAMPLES:
       sage: flatten([[1,1],[1],2])
       [1, 1, 1, 2]
       sage: flatten([[1,2,3], (4,5), [[[1],[2]]]])
       [1, 2, 3, 4, 5, 1, 2]

   In the following example, the vector isn't flattened because
   it is not given in the ltypes input.
       sage: flatten((['Hi',2,vector(QQ,[1,2,3])],(4,5,6)))
       ['Hi', 2, (1, 2, 3), 4, 5, 6]

   We give the vector type and then even the vector gets flattened:
       sage: flatten((['Hi',2,vector(QQ,[1,2,3])], (4,5,6)), ltypes=(list, tuple,sage.modules.vector_rational_dense.Vector_rational_dense))
       ['Hi', 2, 1, 2, 3, 4, 5, 6]

   We flatten a finite field.
       sage: flatten(GF(5))
       [0, 1, 2, 3, 4]
       sage: flatten([GF(5)])
       [Finite Field of size 5]
       sage: flatten([GF(5)], ltypes = (list, tuple, sage.rings.finite_field.FiniteField_prime_modn))
       [0, 1, 2, 3, 4]

   Degenerate cases:
      sage: flatten([[]])
      []
      sage: flatten([[[]]])
      []
   """
   index = 0
   new_list = [x for x in in_list]
   while index < len(new_list):
      while isinstance(new_list[index], ltypes):
         v = list(new_list[index])
         if len(v) != 0:
            new_list[index : index + 1] = v
         else:
            new_list.pop(index)
            break
      index += 1
   return new_list
