##############################################################################
#       Copyright (C) 2010 Nathann Cohen <nathann.cohen@gmail.com>
#  Distributed under the terms of the GNU General Public License (GPL)
#  The full text of the GPL is available at:
#                  http://www.gnu.org/licenses/
##############################################################################

from generic_backend cimport GenericBackend
include '../../../../../devel/sage/sage/ext/stdsage.pxi'


cdef extern from *:
    ctypedef double* const_double_ptr "const double*"
    ctypedef char * const_char_ptr "const char*"

cdef extern from "float.h":
    cdef double DBL_MAX

cdef extern from "../../../local/include/glpk.h":
     ctypedef struct c_glp_prob "glp_prob":
         pass
     ctypedef struct c_glp_iocp "glp_iocp":
         int presolve
         int msg_lev
         int gmi_cuts
         int fp_heur
         int mir_cuts
     c_glp_iocp * new_c_glp_iocp "new glp_iocp" ()
     void glp_init_iocp(c_glp_iocp *)
     c_glp_prob * glp_create_prob()
     void glp_set_prob_name(c_glp_prob *, char *)
     const_char_ptr glp_get_prob_name(c_glp_prob *)
     void glp_set_obj_dir(c_glp_prob *, int)
     void glp_add_rows(c_glp_prob *, int)
     void glp_add_cols(c_glp_prob *, int)
     void glp_set_row_name(c_glp_prob *, int, char *)
     void glp_set_col_name(c_glp_prob *, int, char *)
     void glp_set_row_bnds(c_glp_prob *, int, int, double, double)
     void glp_set_col_bnds(c_glp_prob *, int, int, double, double)
     void glp_set_obj_coef(c_glp_prob *, int, double)
     void glp_load_matrix(c_glp_prob *, int, int *, int *, double *)
     void glp_simplex(c_glp_prob *, int)
     int glp_intopt(c_glp_prob *, c_glp_iocp *)
     int lpx_intopt(c_glp_prob *)
     void glp_delete_prob(c_glp_prob *)
     double glp_get_col_prim(c_glp_prob *, int)
     double glp_get_obj_val(c_glp_prob *)
     int glp_get_num_rows(c_glp_prob *)
     int glp_get_num_cols(c_glp_prob *)
     double glp_mip_col_val(c_glp_prob *, int)
     double glp_mip_obj_val(c_glp_prob *)
     void glp_set_col_kind(c_glp_prob *, int, int)
     int glp_write_mps(c_glp_prob *lp, int fmt, void *parm, char *fname)
     int glp_write_lp(c_glp_prob *lp, void *parm, char *fname)
     void glp_set_prob_name(c_glp_prob *lp, char *name)
     void glp_set_obj_name(c_glp_prob *lp, char *name)
     void glp_set_row_name(c_glp_prob *lp, int i, char *name)
     void glp_set_col_name(c_glp_prob *lp, int i, char *name)

     double glp_get_row_ub(c_glp_prob *lp, int i)
     double glp_get_row_lb(c_glp_prob *lp, int i)

     double glp_get_col_ub(c_glp_prob *lp, int i)
     double glp_get_col_lb(c_glp_prob *lp, int i)
     void glp_set_col_ub(c_glp_prob *lp, int i, double value)
     void glp_set_col_lb(c_glp_prob *lp, int i, double value)


     const_char_ptr glp_get_row_name(c_glp_prob *lp, int i)
     const_char_ptr glp_get_col_name(c_glp_prob *lp, int i)

     void glp_create_index(c_glp_prob *lp)

     double glp_get_col_lb(c_glp_prob *lp, int i)
     double glp_get_col_ub(c_glp_prob *lp, int i)

     int glp_mip_status(c_glp_prob *lp)
     int glp_set_mat_row(c_glp_prob *lp, int, int, int *, double * )
     int glp_set_mat_col(c_glp_prob *lp, int, int, int *, double * )
     int glp_get_mat_row(c_glp_prob *lp, int, int *, double * )
     double glp_get_row_ub(c_glp_prob *lp, int)
     double glp_get_row_lb(c_glp_prob *lp, int)
     int glp_get_col_kind(c_glp_prob *lp, int)
     double glp_get_obj_coef(c_glp_prob *lp, int)
     int glp_get_obj_dir(c_glp_prob *lp)
     void glp_copy_prob(c_glp_prob *dst, c_glp_prob *src, int names)



     int GLP_ON
     int GLP_MAX
     int GLP_MIN
     int GLP_UP
     int GLP_FR
     int GLP_DB
     int GLP_FX
     int GLP_LO
     int GLP_CV
     int GLP_IV
     int GLP_BV
     int GLP_MSG_OFF
     int GLP_MSG_ERR
     int GLP_MSG_ON
     int GLP_MSG_ALL
     int GLP_MPS_DECK
     int GLP_MPS_FILE

     int GLP_UNDEF
     int GLP_OPT
     int GLP_FEAS
     int GLP_NOFEAS

cdef class GLPKBackend(GenericBackend):
    cdef c_glp_prob * lp
    cdef c_glp_iocp * iocp
    cpdef GLPKBackend copy(self)
