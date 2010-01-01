cdef extern from *:
    ctypedef double* const_double_ptr "const double*"

cdef extern from "../../local/include/glpk.h":
     ctypedef struct c_glp_prob "glp_prob":
         pass
     ctypedef struct c_glp_iocp "glp_iocp":
         int presolve
         int msg_lev
     c_glp_iocp *new_c_glp_iocp "new glp_iocp" ()
     void glp_init_iocp(c_glp_iocp *)
     c_glp_prob * glp_create_prob()
     void glp_set_prob_name(c_glp_prob *, char *)
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
     void glp_intopt(c_glp_prob *, c_glp_iocp *)
     int lpx_intopt(c_glp_prob *)
     void glp_delete_prob(c_glp_prob *)
     double glp_get_col_prim(c_glp_prob *, int)
     double glp_get_obj_val(c_glp_prob *)
     double glp_mip_col_val(c_glp_prob *, int)
     double glp_mip_obj_val(c_glp_prob *)
     void glp_set_col_kind(c_glp_prob *, int, int)
     int glp_write_mps(c_glp_prob *lp, int fmt, void *parm, char *fname)
     int glp_write_lp(c_glp_prob *lp, void *parm, char *fname)
     void glp_set_prob_name(c_glp_prob *lp, char *name)
     void glp_set_obj_name(c_glp_prob *lp, char *name)
     void glp_set_row_name(c_glp_prob *lp, int i, char *name)
     void glp_set_col_name(c_glp_prob *lp, int i, char *name)
     int glp_mip_status(c_glp_prob *lp)

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
