Numerical Optimization
======================

.. toctree::
   :maxdepth: 1

   sage/numerical/knapsack
   sage/numerical/mip
   sage/numerical/sdp
   sage/numerical/linear_functions
   sage/numerical/linear_tensor
   sage/numerical/linear_tensor_element
   sage/numerical/linear_tensor_constraints
   sage/numerical/optimize
   sage/numerical/interactive_simplex_method

Linear Optimization (LP) and Mixed Integer Linear Optimization (MIP) Solver backends
------------------------------------------------------------------------------------

.. toctree::
   :maxdepth: 1

   sage/numerical/backends/generic_backend
   sage/numerical/backends/interactivelp_backend
   sage/numerical/backends/glpk_backend
   sage/numerical/backends/glpk_exact_backend
   sage/numerical/backends/glpk_graph_backend
   sage/numerical/backends/ppl_backend
   sage/numerical/backends/cvxopt_backend

.. backends depending on optional packages
..   sage/numerical/backends/coin_backend
..   sage/numerical/backends/cplex_backend
..   sage/numerical/backends/gurobi_backend

Sage also supports, via optional packages, CBC (COIN-OR), CPLEX (ILOG), and Gurobi. In order to find out
how to use them in Sage, please refer to the `Thematic Tutorial on Linear
Programming
<http://www.sagemath.org/doc/thematic_tutorials/linear_programming.html>`_.

Semidefinite Optimization (SDP) Solver backends
-----------------------------------------------

.. toctree::
   :maxdepth: 1

   sage/numerical/backends/cvxopt_sdp_backend

For more details on CVXOPT, see `CVXOPT documentation <http://cvxopt.org/documentation/index.html>`_.

.. include:: ../footer.txt
