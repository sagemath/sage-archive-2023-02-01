cdef extern from *:
    ctypedef double* const_double_ptr "const double*"

cdef extern from "../../local/include/coin/CoinPackedVector.hpp":
     ctypedef struct c_CoinPackedVector "CoinPackedVector":
         void insert(float, float)
     c_CoinPackedVector *new_c_CoinPackedVector "new CoinPackedVector" ()
     void del_CoinPackedVector "delete" (c_CoinPackedVector *)
cdef extern from "../../local/include/coin/CoinPackedMatrix.hpp":
     ctypedef struct c_CoinPackedMatrix "CoinPackedMatrix":
         void setDimensions(int, int)
         void appendRow(c_CoinPackedVector)
     c_CoinPackedMatrix *new_c_CoinPackedMatrix "new CoinPackedMatrix" (bool, double, double)
     void del_CoinPackedMatrix "delete" (c_CoinPackedMatrix *)

cdef extern from "../../local/include/coin/OsiCpxSolverInterface.hpp":
     ctypedef struct c_OsiCpxSolverInterface "OsiCpxSolverInterface":
         double getInfinity()
         void loadProblem(c_CoinPackedMatrix, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr)
         void assignProblem(c_CoinPackedMatrix *, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr)
         void writeMps(char *, char *, double)
         void initialSolve()
         void branchAndBound()
         void readMps(string)
         float getObjValue()
         double * getColSolution()
         void setObjSense (double )
         void setLogLevel(int)
         void setInteger(int)
         void setContinuous(int)
         int isAbandoned ()
         int isProvenOptimal ()
         int isProvenPrimalInfeasible ()
         int isProvenDualInfeasible ()
         int isPrimalObjectiveLimitReached ()
         int isDualObjectiveLimitReached ()
         int isIterationLimitReached ()
         void setMaximumSolutions(int)
         int getMaximumSolutions()
     c_OsiCpxSolverInterface *new_c_OsiCpxSolverInterface "new OsiCpxSolverInterface" ()
     void del_OsiCpxSolverInterface "delete" (c_OsiCpxSolverInterface *)



from sage.numerical.osi_interface cimport Osi_interface
from sage.numerical.osi_interface cimport c_OsiSolverInterface
