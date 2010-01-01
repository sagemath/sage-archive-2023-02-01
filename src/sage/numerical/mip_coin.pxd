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

cdef extern from "../../local/include/coin/CoinMessageHandler.hpp":
     ctypedef struct c_CoinMessageHandler "CoinMessageHandler":
         void setLogLevel (int)
     c_CoinMessageHandler *new_c_CoinMessageHandler "new CoinMessageHandler" ()
     void del_CoinMessageHandler "delete" (c_CoinMessageHandler *)


cdef extern from "../../local/include/coin/CbcModel.hpp":
     ctypedef struct c_CbcModel "CbcModel":
         c_CoinMessageHandler * messageHandler ()
         void setNumberThreads (int)
         int getSolutionCount()

     c_CbcModel *new_c_CbcModel "new CbcModel" ()
     void del_CbcModel "delete" (c_CbcModel *)

cdef extern from "../../local/include/coin/OsiCbcSolverInterface.hpp":
     ctypedef struct c_OsiCbcSolverInterface "OsiCbcSolverInterface":
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
         c_CoinMessageHandler * messageHandler ()
         c_CbcModel * getModelPtr  ()
         int isAbandoned ()
         int isProvenOptimal ()
         int isProvenPrimalInfeasible ()
         int isProvenDualInfeasible ()
         int isPrimalObjectiveLimitReached ()
         int isDualObjectiveLimitReached ()
         int isIterationLimitReached ()
         void setMaximumSolutions(int)
         int getMaximumSolutions()
     c_OsiCbcSolverInterface *new_c_OsiCbcSolverInterface "new OsiCbcSolverInterface" ()
     void del_OsiCbcSolverInterface "delete" (c_OsiCbcSolverInterface *)
