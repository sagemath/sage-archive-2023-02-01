#include "dgs.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>

static inline void _dgs_disc_gauss_dp_init_bexp(dgs_disc_gauss_dp_t *self, double sigma, long upper_bound) {
  self->f = (2*sigma*sigma);
  size_t l = 2*ceil(log2(upper_bound));
  self->Bexp = dgs_bern_exp_dp_init(self->f, l);
}

dgs_disc_gauss_dp_t *dgs_disc_gauss_dp_init(double sigma, long c, size_t tau, dgs_disc_gauss_alg_t algorithm) {
  assert(tau > 0); 
  size_t upper_bound;
  
  dgs_disc_gauss_dp_t *self = (dgs_disc_gauss_dp_t*)calloc(sizeof(dgs_disc_gauss_dp_t),1);
  if (!self) abort();

  self->sigma = sigma;
  self->c = c;  
  self->tau = tau;
  
  switch(algorithm) {

  case DGS_DISC_GAUSS_UNIFORM_ONLINE:
    self->call = dgs_disc_gauss_dp_call_uniform_online;
    upper_bound = ceil(self->sigma*tau);
    self->upper_bound = upper_bound;
    self->two_upper_bound_plus_one = 2*upper_bound + 1;
    self->f = -1.0/(2.0*(self->sigma*self->sigma));
    break;

  case DGS_DISC_GAUSS_UNIFORM_TABLE:
    self->call = dgs_disc_gauss_dp_call_uniform_table;
    upper_bound = ceil(self->sigma*tau);
    self->upper_bound = upper_bound;
    self->two_upper_bound_plus_one = 2*upper_bound + 1;
    self->B = dgs_bern_uniform_init(0);
    self->f = -1.0/(2.0*(sigma*sigma));
    
    self->rho = (double*)malloc(sizeof(double)*self->upper_bound);
    if (!self->rho) abort();
    for(unsigned long x=0; x<self->upper_bound; x++) {
      self->rho[x] = exp(((double)(x*x))*self->f);
    }
    self->rho[0]/= 2.0;
    break;

  case DGS_DISC_GAUSS_UNIFORM_LOGTABLE:
    self->call = dgs_disc_gauss_dp_call_uniform_logtable;
    upper_bound = ceil(self->sigma*tau);
    self->upper_bound = upper_bound;
    self->two_upper_bound_plus_one = 2*upper_bound + 1;
    _dgs_disc_gauss_dp_init_bexp(self, self->sigma, self->upper_bound);
   break;
 
  case DGS_DISC_GAUSS_SIGMA2_LOGTABLE: {
    self->call = dgs_disc_gauss_dp_call_sigma2_logtable;
    double sigma2 = sqrt(1.0/(2*log(2.0)));
    double k = sigma/sigma2;
    self->k = round(k);
    self->sigma = self->k * sigma2;
    upper_bound = ceil(self->sigma*tau);
    self->upper_bound = upper_bound;
    self->two_upper_bound_plus_one = 2*upper_bound + 1;
    _dgs_disc_gauss_dp_init_bexp(self, self->sigma, self->upper_bound);
    self->B = dgs_bern_uniform_init(0);
    self->D2 = dgs_disc_gauss_sigma2p_init();
    break;
  }
    
  default:
    abort();
  }  
  return self;
}

long dgs_disc_gauss_dp_call_uniform_online(dgs_disc_gauss_dp_t *self) {
  long x;
  double y, z;
  long c = self->c;
  do {
    x = _dgs_randomm_libc(self->two_upper_bound_plus_one) - self->upper_bound + c;
    z = exp(((double)(x-c)*(x-c))*self->f);
    y = drand48();
  } while (y >= z);

  return x;
}

long dgs_disc_gauss_dp_call_uniform_table(dgs_disc_gauss_dp_t *self) {
  long x;
  double y;
  do {
    x = _dgs_randomm_libc(self->upper_bound);
    y = drand48();
  } while (y >= self->rho[x]);

  if(dgs_bern_uniform_call_libc(self->B))
    x = -x;
  return x + self->c;
}

long dgs_disc_gauss_dp_call_uniform_logtable(dgs_disc_gauss_dp_t *self) {
  long x;
  do {
    x = _dgs_randomm_libc(self->two_upper_bound_plus_one) - self->upper_bound;
  } while (dgs_bern_exp_dp_call(self->Bexp, x*x) == 0);
  return x + self->c;
}

// TODO
long dgs_disc_gauss_dp_call_sigma2_logtable(dgs_disc_gauss_dp_t *self) {
  long x, y, z;
  long k = self->k;

  do {
    do {
      x = dgs_disc_gauss_sigma2p_dp_call(self->D2);
      y = _dgs_randomm_libc(self->k);
    } while (dgs_bern_exp_dp_call(self->Bexp, y*(y + 2*k*x)) == 0);
    z = k*x + y;
    if (!z) {
      if (dgs_bern_uniform_call_libc(self->B))
        break;
    } else {
      break;
    }
  } while (1);
  if(dgs_bern_uniform_call_libc(self->B))
    z = -z;
  return z + self->c;
}


void dgs_disc_gauss_dp_clear(dgs_disc_gauss_dp_t *self) {
  if (self->B) dgs_bern_uniform_clear(self->B);
  if (self->Bexp) dgs_bern_exp_dp_clear(self->Bexp);
  if (self->rho) free(self->rho);
  free(self);
}



