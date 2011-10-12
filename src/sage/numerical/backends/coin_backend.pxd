##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from sage.numerical.backends.generic_backend cimport GenericBackend

include '../../../../../devel/sage/sage/ext/stdsage.pxi'
include '../../ext/cdefs.pxi'

cdef extern from *:
    ctypedef double* const_double_ptr "const double*"

cdef extern from "../../local/include/coin/CoinPackedVector.hpp":
     ctypedef struct c_CoinPackedVector "CoinPackedVector":
         void insert(float, float)
     c_CoinPackedVector *new_c_CoinPackedVector "new CoinPackedVector" ()
     void del_CoinPackedVector "delete" (c_CoinPackedVector *)

cdef extern from "../../local/include/coin/CoinShallowPackedVector.hpp":
     ctypedef struct c_CoinShallowPackedVector "CoinShallowPackedVector":
         void insert(float, float)
         int * getIndices ()
         double * getElements ()
         int getNumElements ()
     c_CoinShallowPackedVector *new_c_CoinShallowPackedVector "new CoinShallowPackedVector" ()
     void del_CoinShallowPackedVector "delete" (c_CoinShallowPackedVector *)

cdef extern from "../../local/include/coin/CoinPackedMatrix.hpp":
     ctypedef struct c_CoinPackedMatrix "CoinPackedMatrix":
         void setDimensions(int, int)
         void appendRow(c_CoinPackedVector)
         c_CoinShallowPackedVector getVector(int)
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

cdef extern from "../../local/include/coin/OsiSolverInterface.hpp":
     cdef cppclass OsiSolverInterface:
         pass

cdef extern from "../../local/include/coin/OsiCbcSolverInterface.hpp":
     cdef cppclass c_OsiCbcSolverInterface "OsiCbcSolverInterface":
         double getInfinity()
         void loadProblem(c_CoinPackedMatrix, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr)
         void assignProblem(c_CoinPackedMatrix *, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr, const_double_ptr)
         void writeMps(char *, char *, double)
         void initialSolve()
         void branchAndBound()
         void readMps(string)
         double getObjValue()
         int getObjSense()
         double * getColSolution()
         void setObjSense (double )
         void setObjCoeff(int, double)
         double * getObjCoefficients ()
         void setLogLevel(int)
         void setInteger(int)
         void setContinuous(int)
         c_CoinMessageHandler * messageHandler ()
         c_CbcModel * getModelPtr  ()
         int isContinuous (int)
         int isAbandoned ()
         int isProvenOptimal ()
         int isProvenPrimalInfeasible ()
         int isProvenDualInfeasible ()
         int isPrimalObjectiveLimitReached ()
         int isDualObjectiveLimitReached ()
         int isIterationLimitReached ()
         void setMaximumSolutions(int)
         int getMaximumSolutions()
         int getNumCols()
         int getNumRows()
         void setColLower (int elementIndex, double elementValue)
         void setColUpper (int elementIndex, double elementValue)
         double * getRowLower()
         double * getRowUpper()
         double * getColLower()
         double * getColUpper()
         void addCol (int numberElements, int *rows, double *elements, double collb, double colub, double obj)
         void addRow (c_CoinPackedVector vec, double rowlb, double rowub)
         c_CoinPackedMatrix * getMatrixByRow()
         c_OsiCbcSolverInterface * c_OsiCbcSolverInterface(OsiSolverInterface * solver)
     #c_OsiCbcSolverInterface *new_c_OsiCbcSolverInterface "new OsiCbcSolverInterface" ()
     #c_OsiCbcSolverInterface * new2_c_OsiCbcSolverInterface "new OsiCbcSolverInterface" (OsiSolverInterface * solver)
     #void del_OsiCbcSolverInterface "delete" (c_OsiCbcSolverInterface *)

cdef class CoinBackend(GenericBackend):
    cdef c_OsiCbcSolverInterface * si
    cpdef CoinBackend copy(self)
