include 'mipGlpk.pxd'
include '../../../../devel/sage/sage/ext/stdsage.pxi'

def solveGlpk(self,log=False,objective_only=False):
  r"""
  Solves the MIP using GLPK. Use solve() instead.
  """

  if self.objective==None:
    raise Exception("No objective function has been defined !")
  cdef c_glp_prob * lp
  k=1
  for c in self.constraints:
    k=k+c["card"]
  lp = glp_create_prob()
  glp_add_rows(lp, len(self.constraints))
  glp_add_cols(lp, len(self.variables))


  # Gives an ID to each one of the variable. It will be their index in the problem Matrix
  #self.setVariablesId();


  # Definition of the objective function
  for (id,coeff) in [(id,coeff) for (id,coeff) in self._NormalForm(self.objective).items() if id!=-1]:
    #glp_set_obj_coef(lp, self.variables[i],self.objective[i] )
    glp_set_obj_coef(lp, id,coeff )

  # Upper and lower bounds on the variables

  for i in self.variables.keys():
      if self.getmin(i)!=None and self.getmax(i)!=None:
          glp_set_col_bnds(lp, self.variables[i], GLP_DB, self.getmin(i), self.getmax(i))
      elif self.getmin(i)==None and self.getmax(i)!=None:
          glp_set_col_bnds(lp, self.variables[i], GLP_UP, 0, self.getmax(i))
      elif self.getmin(i)!=None and self.getmax(i)==None:
          glp_set_col_bnds(lp, self.variables[i], GLP_LO, self.getmin(i),0)
      else:
          glp_set_col_bnds(lp, self.variables[i], GLP_FR, 0,0)



  # Definition of the problem matrix

  cdef int * ia = <int*> sage_malloc(k*sizeof(int))
  cdef int * ja = <int*> sage_malloc(k*sizeof(int))
  cdef double * ar = <double*> sage_malloc(k*sizeof(double))
  a=1
  b=1
  for c in self.constraints:
      minimum=c["min"]
      maximum=c["max"]

      function=c["function"]
      for (id,coeff) in [(id,coeff) for (id,coeff) in self._NormalForm(function).items() if id!=-1]:
          ar[b]=coeff
          ia[b]=a
          ja[b]=id
          b=b+1

      if minimum!=None and maximum!=None:
          if minimum==maximum:
               glp_set_row_bnds(lp, a, GLP_FX, minimum, maximum)
          else:
               glp_set_row_bnds(lp, a, GLP_DB, minimum, maximum)
      elif minimum==None and maximum!=None:
          glp_set_row_bnds(lp, a, GLP_UP, 0, maximum)
      elif minimum!=None and maximum==None:
          glp_set_row_bnds(lp, a, GLP_LO, minimum,0)
      else:
          glp_set_row_bnds(lp, a, GLP_FR, 0,0)
      a=a+1

  # Direction of the problem : maximization, minimization...
  if self.sense==-1:
      glp_set_obj_dir(lp, GLP_MIN)
  else:
      glp_set_obj_dir(lp, GLP_MAX)

  # Sets variables as integers when needed

  for i in self.variables.keys():
     if self.isinteger(i):
          glp_set_col_kind(lp,self.variables[i],GLP_IV)
     elif self.isbinary(i):
          glp_set_col_kind(lp,self.variables[i],GLP_BV)
     else:
          glp_set_col_kind(lp,self.variables[i],GLP_CV)

  glp_load_matrix(lp, k-1, ia, ja, ar)
  cdef c_glp_iocp * iocp
  iocp=new_c_glp_iocp()
  glp_init_iocp(iocp)

  iocp.presolve=GLP_ON
  # stuff related to the loglevel
  if log==0:
      iocp.msg_lev=GLP_MSG_OFF
  elif log==1:
      iocp.msg_lev=GLP_MSG_ERR
  elif log==2:
      iocp.msg_lev=GLP_MSG_ON
  else:
      iocp.msg_lev=GLP_MSG_ALL
  glp_intopt(lp,iocp)
  obj = glp_mip_obj_val(lp)

  if objective_only==False:
       for (v,id) in self.variables.items():
            self.values[v]=glp_mip_col_val(lp, id) if self.isreal(v) else round(glp_mip_col_val(lp, id))
  glp_delete_prob(lp)
  return obj
