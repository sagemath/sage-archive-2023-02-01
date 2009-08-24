include 'mipGlpk.pxd'
include '../../../devel/sage/sage/ext/stdsage.pxi'

def solveGlpk(self,log=False,objective_only=False):
  r"""
  Solves the MIP using GLPK. Use solve() instead.
  """

  if self.objective==None:
    raise Exception("No objective function has been defined !")
  cdef c_glp_prob * lp
  k=1
  for c in self.constraints:
    k=k+c.card()
  lp = glp_create_prob()
  glp_add_rows(lp, len(self.constraints))
  glp_add_cols(lp, len(self.variables))


  # Gives an ID to each one of the variable. It will be their index in the problem Matrix
  self.setVariablesId();


  # Definition of the objective function
  for i in self.objective.keys():
    glp_set_obj_coef(lp, self.variables[i].getid()+1,self.objective[i] )

  # Upper and lower bounds on the variables

  for i in self.variables.keys():
      if self.variables[i].getmin()!=None and self.variables[i].getmax()!=None:
          glp_set_col_bnds(lp, 1+self.variables[i].getid(), GLP_DB, self.variables[i].getmin(), self.variables[i].getmax())
      elif self.variables[i].getmin()==None and self.variables[i].getmax()!=None:
          glp_set_col_bnds(lp, 1+self.variables[i].getid(), GLP_UP, 0, self.variables[i].getmax())
      elif self.variables[i].getmin()!=None and self.variables[i].getmax()==None:
          glp_set_col_bnds(lp, 1+self.variables[i].getid(), GLP_LO, self.variables[i].getmin(),0)
      else:
          glp_set_col_bnds(lp, 1+self.variables[i].getid(), GLP_FR, 0,0)


  # Definition of the problem matrix

  cdef int * ia = <int*> sage_malloc(k*sizeof(int))
  cdef int * ja = <int*> sage_malloc(k*sizeof(int))
  cdef double * ar = <double*> sage_malloc(k*sizeof(double))
  a=1
  b=1
  for c in self.constraints:
      for v in c.dict.keys():
          ar[b]=c.dict[v]
          ia[b]=a
          ja[b]=self.variables[v].getid()+1
          b=b+1
      if c.getmin()!=None and c.getmax()!=None:
          if c.getmin()==c.getmax():
               glp_set_row_bnds(lp, a, GLP_FX, c.getmin(), c.getmax())
          else:
               glp_set_row_bnds(lp, a, GLP_DB, c.getmin(), c.getmax())
      elif c.getmin()==None and c.getmax()!=None:
          glp_set_row_bnds(lp, a, GLP_UP, 0, c.getmax())
      elif c.getmin()!=None and c.getmax()==None:
          glp_set_row_bnds(lp, a, GLP_LO, c.getmin(),0)
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
     if self.variables[i].isinteger():
          glp_set_col_kind(lp,self.variables[i].getid()+1,GLP_IV)
     elif self.variables[i].isbinary():
          glp_set_col_kind(lp,self.variables[i].getid()+1,GLP_BV)
     else:
          glp_set_col_kind(lp,self.variables[i].getid()+1,GLP_CV)



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

  if objective_only==True:
       glp_delete_prob(lp)
       return obj
  else:
       list={}
       for v in self.variables.keys():
            list[v]=glp_mip_col_val(lp, self.variables[v].getid()+1)
       glp_delete_prob(lp)
       return [obj,list]
