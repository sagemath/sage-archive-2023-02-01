.. _cython_interface:
.. nodoctest
================================
How to call a C code from Sage ?
================================

If you have some C/C++ code that you would like to call from Sage for your own
use, this document is for you.

- Do you want to **contibute** to Sage by adding your interface to its code? The
  (more complex) instructions are `available here
  <http://www.sagemath.org/doc/developer/index.html#packaging-third-party-code>`_

Calling "hello_world()" from hello.c
------------------------------------

Let us suppose that you have in your current directory a file named hello.c::

  ~/a$ cat hello.c
  #include <stdio.h>

  void hello_world(){
  printf("Hello World\n");
  }

  void main(){
  hello_world();
  }
  ~/a$ gcc hello.c -o hello; ./hello
  Hello World

In order to call this function from Sage, you must create a Cython file (i.e. a
file whose extension is .pyx). This file contains a header containing the
signature of the function that you want to call::

  ~/a$ cat hello_sage.pyx
  cdef extern from "hello.c":
  void hello_world()

  def my_interface_function():
      hello_world() # This is the C function from hello.c

You can now load this file in Sage, and call the C code though the Python
function ``my_interface_function``::

  sage: %runfile hello_sage.pyx
  Compiling ./hello_sage.pyx...
  sage: my_bridge_function()
  Hello World

Arguments and return value
--------------------------

Calling function with more complex arguments and return values works the same
way. To learn more about the Cython language, `click here
<http://docs.cython.org/src/reference/language_basics.html>`_

The following example defines a function taking and returning ``int *``
pointers, and involves some memory allocation. The C code defines a function
whose purpose is to return the sum of two vectors as a third vector.

**The C file**::

  #include <string.h>

  int * sum_of_two_vectors(int n, int * vec1, int * vec2){
    /*
       INPUT : two arrays vec1,vec2 of n integers
       OUTPUT: an array of size n equal to vec1+vec2
    */
    int * sum = (int *) malloc(n*sizeof(int));
    int i;

    for(i=0;i<n;i++)
      sum[i] = vec1[i] + vec2[i];
    return sum;
  }

**The Cython file**::

  cdef extern from "double_vector.c":
      int * sum_of_two_vectors(int n, int * vec1, int * vec2)

  from libc.stdlib cimport malloc, free

  def sage_sum_of_vectors(n, list1, list2):
      cdef int * vec1 = <int *> malloc(n*sizeof(int))
      cdef int * vec2 = <int *> malloc(n*sizeof(int))

      # Fill the vectors
      for i in range(n):
          vec1[i] = list1[i]
          vec2[i] = list2[i]

      # Call the C function
      cdef int * vec3 = sum_of_two_vectors(n,vec1,vec2)

      # Save the answer in a Python object
      answer = [vec3[i] for i in range(n)]

      free(vec1)
      free(vec2)
      free(vec3)

      return answer


**Call from Sage**::

  sage: %runfile double_vector_sage.pyx
  Compiling ./double_vector_sage.pyx...
  sage: sage_sum_of_vectors(3,[1,1,1],[2,3,4])
  [3, 4, 5]
