/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: longrat.cc,v 1.23 2007/07/27 07:37:53 Singular Exp $ */
/*
* ABSTRACT: computation with long rational numbers (Hubert Grassmann)
*/

#ifndef LONGRAT_CC
#define LONGRAT_CC

#include <string.h>
#include <float.h>
#include "mod2.h"
#include "structs.h"
#include "omalloc.h"
#include "febase.h"
#include "numbers.h"
#include "modulop.h"
#include "ring.h"
#include "shortfl.h"
#include "mpr_complex.h"
#include "longrat.h"


#ifndef BYTES_PER_MP_LIMB
#define BYTES_PER_MP_LIMB sizeof(mp_limb_t)
#endif

/*-----------------------------------------------------------------*/
/*3
* parameter s in number:
* 0 (or FALSE): not normalised rational
* 1 (or TRUE):  normalised rational
* 3          :  integer with n==NULL
*/
/*3
**  'SR_INT' is the type of those integers small enough to fit into  29  bits.
**  Therefor the value range of this small integers is: $-2^{28}...2^{28}-1$.
**
**  Small integers are represented by an immediate integer handle, containing
**  the value instead of pointing  to  it,  which  has  the  following  form:
**
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**      | guard | sign  | bit   | bit   |       | bit   | tag   | tag   |
**      | bit   | bit   | 27    | 26    |       | 0     | 0     | 1     |
**      +-------+-------+-------+-------+- - - -+-------+-------+-------+
**
**  Immediate integers handles carry the tag 'SR_INT', i.e. the last bit is 1.
**  This distuingishes immediate integers from other handles which  point  to
**  structures aligned on 4 byte boundaries and therefor have last bit  zero.
**  (The second bit is reserved as tag to allow extensions of  this  scheme.)
**  Using immediates as pointers and dereferencing them gives address errors.
**
**  To aid overflow check the most significant two bits must always be equal,
**  that is to say that the sign bit of immediate integers has a  guard  bit.
**
**  The macros 'INT_TO_SR' and 'SR_TO_INT' should be used to convert  between
**  a small integer value and its representation as immediate integer handle.
**
**  Large integers and rationals are represented by z and n
**  where n may be undefined (if s==3)
**  NULL represents only deleted values
*/
#define SR_HDL(A) ((long)(A))
/*#define SR_INT    1*/
/*#define INT_TO_SR(INT)  ((number) (((long)INT << 2) + SR_INT))*/
// #define SR_TO_INT(SR)   (((long)SR) >> 2)

#define MP_SMALL 1
//#define mpz_isNeg(A) (mpz_cmp_si(A,(long)0)<0)
#define mpz_isNeg(A) ((A)->_mp_size<0)
#define mpz_limb_size(A) ((A)->_mp_size)
#define mpz_limb_d(A) ((A)->_mp_d)
#define MPZ_DIV(A,B,C) mpz_tdiv_q((A),(B),(C))
#define MPZ_EXACTDIV(A,B,C) mpz_divexact((A),(B),(C))


/***************************************************************
 *
 * Routines which are never inlined by p_Numbers.h
 *
 *******************************************************************/
#ifndef P_NUMBERS_H

omBin rnumber_bin = omGetSpecBin(sizeof(snumber));

number nlOne=INT_TO_SR(1);

#if (__GNU_MP_VERSION*10+__GNU_MP_VERSION_MINOR < 31)
void mpz_mul_si (mpz_ptr r, mpz_srcptr s, long int si)
{
  if (si>=0)
    mpz_mul_ui(r,s,si);
  else
  {
    mpz_mul_ui(r,s,-si);
    mpz_neg(r,r);
  }
}
#endif

static int nlPrimeM;
static number nlMapP(number from)
{
  number to;
  int save=npPrimeM;
  npPrimeM=nlPrimeM;
  to = nlInit(npInt(from));
  npPrimeM=save;
  return to;
}

static number nlMapLongR(number from);
static number nlMapR(number from);

nMapFunc nlSetMap(ring src, ring dst)
{
  if (rField_is_Q(src))
  {
    return nlCopy;
  }
  if (rField_is_Zp(src))
  {
    nlPrimeM=rChar(src);
    return nlMapP;
  }
  if (rField_is_R(src))
  {
    return nlMapR;
  }
  if (rField_is_long_R(src))
  {
    return nlMapLongR; /* long R -> Q */
  }
  return NULL;
}

#ifdef LDEBUG
BOOLEAN nlDBTest(number a, char *f,int l)
{
  if (a==NULL)
  {
    Print("!!longrat: NULL in %s:%d\n",f,l);
    return FALSE;
  }
  //if ((int)a==1) Print("!! 0x1 as number ? %s %d\n",f,l);
  if ((((long)a)&3L)==3L)
  {
    Print(" !!longrat:ptr(3) in %s:%d\n",f,l);
    return FALSE;
  }
  if ((((long)a)&3L)==1L)
  {
    if (((((long)a)<<1)>>1)!=((long)a))
    {
      Print(" !!longrat:arith:%x in %s:%d\n",(long)a, f,l);
      return FALSE;
    }
    return TRUE;
  }
  omCheckIf(omCheckAddrSize(a,sizeof(*a)), return FALSE);
  if (a->debug!=123456)
  {
    Print("!!longrat:debug:%d in %s:%d\n",a->debug,f,l);
    a->debug=123456;
    return FALSE;
  }
  if ((a->s<0)||(a->s>4))
  {
    Print("!!longrat:s=%d in %s:%d\n",a->s,f,l);
    return FALSE;
  }
  omCheckAddrSize(a->z._mp_d,a->z._mp_alloc*BYTES_PER_MP_LIMB);
  if (a->z._mp_alloc==0)
    Print("!!longrat:z->alloc=0 in %s:%l\n",f,l);

  if (a->s<2)
  {
    omCheckIf(omCheckAddrSize(a->n._mp_d,a->n._mp_alloc*BYTES_PER_MP_LIMB), return FALSE);
    if (a->z._mp_alloc==0)
      Print("!!longrat:n->alloc=0 in %s:%l\n",f,l);
    if ((mpz_size1(&a->n) ==1) && (mpz_cmp_si(&a->n,(long)1)==0))
    {
      Print("!!longrat:integer as rational in %s:%d\n",f,l);
      return FALSE;
    }
    else if (mpz_isNeg(&a->n))
    {
      Print("!!longrat:div. by negative in %s:%d\n",f,l);
      mpz_neg(&a->z,&a->z);
      mpz_neg(&a->n,&a->n);
      return FALSE;
    }
    return TRUE;
  }
  if (a->s==2)
  {
    Print("!!longrat:s=2 in %s:%d\n",f,l);
    return FALSE;
  }
  if (mpz_size1(&a->z)>MP_SMALL) return TRUE;
  int ui=(int)mpz_get_si(&a->z);
  if ((((ui<<3)>>3)==ui)
  && (mpz_cmp_si(&a->z,(long)ui)==0))
  {
    Print("!!longrat:im int %d in %s:%d\n",ui,f,l);
    f=NULL;
    return FALSE;
  }
  return TRUE;
}
#endif

number nlRInit (int i);

static number nlMapR(number from)
{
  double f=nrFloat(from);
  if (f==0.0) return INT_TO_SR(0);
  int f_sign=1;
  if (f<0.0)
  {
    f_sign=-1;
    f=-f;
  }
  int i=0;
  mpz_t h1;
  mpz_init_set_ui(h1,1);
  while((FLT_RADIX*f) < DBL_MAX && i<DBL_MANT_DIG)
  {
    f*=FLT_RADIX;
    mpz_mul_ui(h1,h1,FLT_RADIX);
    i++;
  }
  number r=nlRInit(1);
  mpz_set_d(&(r->z),f);
  memcpy4(&(r->n),&h1,sizeof(h1));
  r->s=0; /* not normalized */
  nlNormalize(r);
  return r;
}

static number nlMapLongR(number from)
{
  gmp_float *ff=(gmp_float*)from;
  mpf_t *f=ff->_mpfp();
  number res;
  lint *dest,*ndest;
  int size, i,negative;
  int e,al,bl;
  mp_ptr qp,dd,nn;

  size = (*f)[0]._mp_size;
  if (size == 0)
    return INT_TO_SR(0);
  if(size<0)
  {
    negative = 1;
    size = -size;
  }
  else
    negative = 0;

  qp = (*f)[0]._mp_d;
  while(qp[0]==0)
  {
    qp++;
    size--;
  }

  e=(*f)[0]._mp_exp-size;
  res = (number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  res->debug=123456;
#endif
  dest = &(res->z);

  if (e<0)
  {
    al = dest->_mp_size = size;
    if (al<2) al = 2;
    dd = (mp_ptr)omAlloc(sizeof(mp_limb_t)*al);
    for (i=0;i<size;i++) dd[i] = qp[i];
    bl = 1-e;
    nn = (mp_ptr)omAlloc(sizeof(mp_limb_t)*bl);
    nn[bl-1] = 1;
    for (i=bl-2;i>=0;i--) nn[i] = 0;
    ndest = &(res->n);
    ndest->_mp_d = nn;
    ndest->_mp_alloc = ndest->_mp_size = bl;
    res->s = 0;
  }
  else
  {
    al = dest->_mp_size = size+e;
    if (al<2) al = 2;
    dd = (mp_ptr)omAlloc(sizeof(mp_limb_t)*al);
    for (i=0;i<size;i++) dd[i+e] = qp[i];
    for (i=0;i<e;i++) dd[i] = 0;
    res->s = 3;
  }

  dest->_mp_d = dd;
  dest->_mp_alloc = al;
  if (negative) dest->_mp_size = -dest->_mp_size;

  if (res->s==0)
    nlNormalize(res);
  else if (mpz_size1(&res->z)<=MP_SMALL)
  {
    // res is new, res->ref is 1
    int ui=(int)mpz_get_si(&res->z);
    if ((((ui<<3)>>3)==ui)
       && (mpz_cmp_si(&res->z,(long)ui)==0))
    {
       mpz_clear(&res->z);
#if defined(LDEBUG)
       res->debug=654324;
#endif
       omFreeBin((ADDRESS)res, rnumber_bin);
       return INT_TO_SR(ui);
    }
  }
#if defined(LDEBUG)
  nlTest(res);
#endif
  return res;
}

//static number nlMapLongR(number from)
//{
//  gmp_float *ff=(gmp_float*)from;
//  const mpf_t *f=ff->mpfp();
//  int f_size=ABS((*f)[0]._mp_size);
//  if (f_size==0)
//    return nlInit(0);
//  int f_sign=1;
//  number work=ngcCopy(from);
//  if (!ngcGreaterZero(work))
//  {
//    f_sign=-1;
//    work=ngcNeg(work);
//  }
//  int i=0;
//  mpz_t h1;
//  mpz_init_set_ui(h1,1);
//  while((FLT_RADIX*f) < DBL_MAX && i<DBL_MANT_DIG)
//  {
//    f*=FLT_RADIX;
//    mpz_mul_ui(h1,h1,FLT_RADIX);
//    i++;
//  }
//  number r=nlRInit(1);
//  mpz_set_d(&(r->z),f);
//  memcpy4(&(r->n),&h1,sizeof(h1));
//  r->s=0; /* not normalized */
//  nlNormalize(r);
//  return r;
//
//
//  number r=nlRInit(1);
//  int f_shift=f_size+(*f)[0]._mp_exp;
//  if ( f_shift > 0)
//  {
//    r->s=0;
//    mpz_init(&r->n);
//    mpz_setbit(&r->n,f_shift*BYTES_PER_MP_LIMB*8);
//    mpz_setbit(&r->z,f_size*BYTES_PER_MP_LIMB*8-1);
//    // now r->z has enough space
//    memcpy4(mpz_limb_d(&r->z),((*f)[0]._mp_d),f_size*BYTES_PER_MP_LIMB);
//    nlNormalize(r);
//  }
//  else
//  {
//    r->s=3;
//    if (f_shift==0)
//    {
//      mpz_setbit(&r->z,f_size*BYTES_PER_MP_LIMB*8-1);
//      // now r->z has enough space
//      memcpy4(mpz_limb_d(&r->z),((*f)[0]._mp_d),f_size*BYTES_PER_MP_LIMB);
//    }
//    else /* f_shift < 0 */
//    {
//      mpz_setbit(&r->z,(f_size-f_shift)*BYTES_PER_MP_LIMB*8-1);
//      // now r->z has enough space
//      memcpy4(mpz_limb_d(&r->z)-f_shift,((*f)[0]._mp_d),
//        f_size*BYTES_PER_MP_LIMB);
//    }
//  }
//  if ((*f)[0]._mp_size<0);
//    r=nlNeg(r);
//  return r;
//}

int nlSize(number a)
{
  if (a==INT_TO_SR(0))
     return 0; /* rational 0*/
  if (SR_HDL(a) & SR_INT)
     return 1; /* immidiate int */
  int s=a->z._mp_alloc;
//  while ((s>0) &&(a->z._mp_d[s]==0L)) s--;
//#if SIZEOF_LONG == 8
//  if (a->z._mp_d[s] < (unsigned long)0x100000000L) s=s*2-1;
//  else s *=2;
//#endif
//  s++;
  if (a->s<2)
  {
    int d=a->n._mp_alloc;
//    while ((d>0) && (a->n._mp_d[d]==0L)) d--;
//#if SIZEOF_LONG == 8
//    if (a->n._mp_d[d] < (unsigned long)0x100000000L) d=d*2-1;
//    else d *=2;
//#endif
    s+=d;
  }
  return s;
}

/*2
* convert number to int
*/
int nlInt(number &i)
{
#ifdef LDEBUG
  nlTest(i);
#endif
  nlNormalize(i);
  if (SR_HDL(i) &SR_INT) return SR_TO_INT(i);
  if (i->s==3)
  {
    if(mpz_size1(&i->z)>MP_SMALL) return 0;
    int ul=(int)mpz_get_si(&i->z);
    if (mpz_cmp_si(&i->z,(long)ul)!=0) return 0;
    return ul;
  }
  lint tmp;
  int ul;
  mpz_init(&tmp);
  MPZ_DIV(&tmp,&i->z,&i->n);
  if(mpz_size1(&tmp)>MP_SMALL) ul=0;
  else
  {
    ul=(int)mpz_get_si(&tmp);
    if (mpz_cmp_si(&tmp,(long)ul)!=0) ul=0;
  }
  mpz_clear(&tmp);
  return ul;
}

/*
* 1/a
*/
number nlInvers(number a)
{
#ifdef LDEBUG
  nlTest(a);
#endif
  number n;
  if (SR_HDL(a) & SR_INT)
  {
    if ((a==INT_TO_SR(1)) || (a==INT_TO_SR(-1)))
    {
      return a;
    }
    if (nlIsZero(a))
    {
      WerrorS(nDivBy0);
      return INT_TO_SR(0);
    }
    n=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
    n->debug=123456;
#endif
    n->s=1;
    if ((long)a>0L)
    {
      mpz_init_set_si(&n->z,(long)1);
      mpz_init_set_si(&n->n,(long)SR_TO_INT(a));
    }
    else
    {
      mpz_init_set_si(&n->z,(long)-1);
      mpz_init_set_si(&n->n,(long)-SR_TO_INT(a));
    }
#ifdef LDEBUG
    nlTest(n);
#endif
    return n;
  }
  n=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  n->debug=123456;
#endif
  {
    n->s=a->s;
    mpz_init_set(&n->n,&a->z);
    switch (a->s)
    {
      case 0:
      case 1:
              mpz_init_set(&n->z,&a->n);
              if (mpz_isNeg(&n->n)) /* && n->s<2*/
              {
                mpz_neg(&n->z,&n->z);
                mpz_neg(&n->n,&n->n);
              }
              if (mpz_cmp_si(&n->n,(long)1)==0)
              {
                mpz_clear(&n->n);
                n->s=3;
                if (mpz_size1(&n->z)<=MP_SMALL)
                {
                  int ui=(int)mpz_get_si(&n->z);
                  if ((((ui<<3)>>3)==ui)
                  && (mpz_cmp_si(&n->z,(long)ui)==0))
                  {
                    mpz_clear(&n->z);
#if defined(LDEBUG)
                    n->debug=654324;
#endif
                    omFreeBin((ADDRESS)n, rnumber_bin);
                    n=INT_TO_SR(ui);
                  }
                }
              }
              break;
      case 3:
              n->s=1;
              if (mpz_isNeg(&n->n)) /* && n->s<2*/
              {
                mpz_neg(&n->n,&n->n);
                mpz_init_set_si(&n->z,(long)-1);
              }
              else
              {
                mpz_init_set_si(&n->z,(long)1);
              }
              break;
    }
  }
#ifdef LDEBUG
  nlTest(n);
#endif
  return n;
}


/*2
* u := a / b in Z, if b | a (else undefined)
*/
number   nlExactDiv(number a, number b)
{
  if (b==INT_TO_SR(0))
  {
    WerrorS(nDivBy0);
    return INT_TO_SR(0);
  }
  if (a==INT_TO_SR(0))
    return INT_TO_SR(0);
  number u;
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    /* the small int -(1<<28) divided by -1 is the large int (1<<28)   */
    if ((a==INT_TO_SR(-(1<<28)))&&(b==INT_TO_SR(-1)))
    {
      return nlRInit(1<<28);
    }
    int aa=SR_TO_INT(a);
    int bb=SR_TO_INT(b);
    return INT_TO_SR(aa/bb);
  }
  number bb=NULL;
  if (SR_HDL(b) & SR_INT)
  {
    bb=nlRInit(SR_TO_INT(b));
    b=bb;
  }
  u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  u->debug=123456;
#endif
  mpz_init(&u->z);
  /* u=a/b */
  u->s = 3;
  MPZ_EXACTDIV(&u->z,&a->z,&b->z);
  if (bb!=NULL)
  {
    mpz_clear(&bb->z);
#if defined(LDEBUG)
    bb->debug=654324;
#endif
    omFreeBin((ADDRESS)bb, rnumber_bin);
  }
  if (mpz_size1(&u->z)<=MP_SMALL)
  {
    int ui=(int)mpz_get_si(&u->z);
    if ((((ui<<3)>>3)==ui)
    && (mpz_cmp_si(&u->z,(long)ui)==0))
    {
      mpz_clear(&u->z);
#if defined(LDEBUG)
      u->debug=654324;
#endif
      omFreeBin((ADDRESS)u, rnumber_bin);
      u=INT_TO_SR(ui);
    }
  }
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

/*2
* u := a / b in Z
*/
number nlIntDiv (number a, number b)
{
  if (b==INT_TO_SR(0))
  {
    WerrorS(nDivBy0);
    return INT_TO_SR(0);
  }
  if (a==INT_TO_SR(0))
    return INT_TO_SR(0);
  number u;
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    /* the small int -(1<<28) divided by -1 is the large int (1<<28)   */
    if ((a==INT_TO_SR(-(1<<28)))&&(b==INT_TO_SR(-1)))
    {
      return nlRInit(1<<28);
    }
    int aa=SR_TO_INT(a);
    int bb=SR_TO_INT(b);
    return INT_TO_SR(aa/bb);
  }
  if (SR_HDL(a) & SR_INT)
  {
    /* the small int -(1<<28) divided by 2^28 is 1   */
    if (a==INT_TO_SR(-(1<<28)))
    {
      if(mpz_cmp_si(&b->z,(long)(1<<28))==0)
      {
        return INT_TO_SR(-1);
      }
    }
    /* a is a small and b is a large int: -> 0 */
    return INT_TO_SR(0);
  }
  number bb=NULL;
  if (SR_HDL(b) & SR_INT)
  {
    bb=nlRInit(SR_TO_INT(b));
    b=bb;
  }
  u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  u->debug=123456;
#endif
  assume(a->s==3);
  assume(b->s==3);
  mpz_init_set(&u->z,&a->z);
  /* u=u/b */
  u->s = 3;
  MPZ_DIV(&u->z,&u->z,&b->z);
  if (bb!=NULL)
  {
    mpz_clear(&bb->z);
#if defined(LDEBUG)
    bb->debug=654324;
#endif
    omFreeBin((ADDRESS)bb, rnumber_bin);
  }
  if (mpz_size1(&u->z)<=MP_SMALL)
  {
    int ui=(int)mpz_get_si(&u->z);
    if ((((ui<<3)>>3)==ui)
    && (mpz_cmp_si(&u->z,(long)ui)==0))
    {
      mpz_clear(&u->z);
#if defined(LDEBUG)
      u->debug=654324;
#endif
      omFreeBin((ADDRESS)u, rnumber_bin);
      u=INT_TO_SR(ui);
    }
  }
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

/*2
* u := a mod b in Z, u>=0
*/
number nlIntMod (number a, number b)
{
  if (b==INT_TO_SR(0))
  {
    WerrorS(nDivBy0);
    return INT_TO_SR(0);
  }
  if (a==INT_TO_SR(0))
    return INT_TO_SR(0);
  number u;
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    if ((long)a>0L)
    {
      if ((long)b>0L)
        return INT_TO_SR(SR_TO_INT(a)%SR_TO_INT(b));
      else
        return INT_TO_SR(SR_TO_INT(a)%(-SR_TO_INT(b)));
    }
    else
    {
      if ((long)b>0L)
      {
        int i=(-SR_TO_INT(a))%SR_TO_INT(b);
        if ( i != 0 ) i = (SR_TO_INT(b))-i;
        return INT_TO_SR(i);
      }
      else
      {
        int i=(-SR_TO_INT(a))%(-SR_TO_INT(b));
        if ( i != 0 ) i = (-SR_TO_INT(b))-i;
        return INT_TO_SR(i);
      }
    }
  }
  if (SR_HDL(a) & SR_INT)
  {
    /* a is a small and b is a large int: -> a or (a+b) or (a-b) */
    if ((long)a<0L)
    {
      if (mpz_isNeg(&b->z))
        return nlSub(a,b);
      /*else*/
        return nlAdd(a,b);
    }
    /*else*/
      return a;
  }
  number bb=NULL;
  if (SR_HDL(b) & SR_INT)
  {
    bb=nlRInit(SR_TO_INT(b));
    b=bb;
  }
  u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  u->debug=123456;
#endif
  mpz_init(&u->z);
  u->s = 3;
  mpz_mod(&u->z,&a->z,&b->z);
  if (bb!=NULL)
  {
    mpz_clear(&bb->z);
#if defined(LDEBUG)
    bb->debug=654324;
#endif
    omFreeBin((ADDRESS)bb, rnumber_bin);
  }
  if (mpz_isNeg(&u->z))
  {
    if (mpz_isNeg(&b->z))
      mpz_sub(&u->z,&u->z,&b->z);
    else
      mpz_add(&u->z,&u->z,&b->z);
  }
  if (mpz_size1(&u->z)<=MP_SMALL)
  {
    int ui=(int)mpz_get_si(&u->z);
    if ((((ui<<3)>>3)==ui)
    && (mpz_cmp_si(&u->z,(long)ui)==0))
    {
      mpz_clear(&u->z);
#if defined(LDEBUG)
      u->debug=654324;
#endif
      omFreeBin((ADDRESS)u, rnumber_bin);
      u=INT_TO_SR(ui);
    }
  }
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

/*2
* u := a / b
*/
number nlDiv (number a, number b)
{
  number u;
  if (nlIsZero(b))
  {
    WerrorS(nDivBy0);
    return INT_TO_SR(0);
  }
  u=(number)omAllocBin(rnumber_bin);
  u->s=0;
#if defined(LDEBUG)
  u->debug=123456;
#endif
// ---------- short / short ------------------------------------
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    int i=SR_TO_INT(a);
    int j=SR_TO_INT(b);
    int r=i%j;
    if (r==0)
    {
      omFreeBin((ADDRESS)u, rnumber_bin);
      return INT_TO_SR(i/j);
    }
    mpz_init_set_si(&u->z,(long)i);
    mpz_init_set_si(&u->n,(long)j);
  }
  else
  {
    mpz_init(&u->z);
// ---------- short / long ------------------------------------
    if (SR_HDL(a) & SR_INT)
    {
      // short a / (z/n) -> (a*n)/z
      if (b->s<2)
      {
        mpz_mul_si(&u->z,&b->n,SR_TO_INT(a));
      }
      else
      // short a / long z -> a/z
      {
        mpz_set_si(&u->z,SR_TO_INT(a));
      }
      if (mpz_cmp(&u->z,&b->z)==0)
      {
        mpz_clear(&u->z);
        omFreeBin((ADDRESS)u, rnumber_bin);
        return INT_TO_SR(1);
      }
      mpz_init_set(&u->n,&b->z);
    }
// ---------- long / short ------------------------------------
    else if (SR_HDL(b) & SR_INT)
    {
      mpz_set(&u->z,&a->z);
      // (z/n) / b -> z/(n*b)
      if (a->s<2)
      {
        mpz_init_set(&u->n,&a->n);
        if ((long)b>0L)
          mpz_mul_ui(&u->n,&u->n,SR_TO_INT(b));
        else
        {
          mpz_mul_ui(&u->n,&u->n,-SR_TO_INT(b));
          mpz_neg(&u->z,&u->z);
        }
      }
      else
      // long z / short b -> z/b
      {
        //mpz_set(&u->z,&a->z);
        mpz_init_set_si(&u->n,SR_TO_INT(b));
      }
    }
// ---------- long / long ------------------------------------
    else
    {
      mpz_set(&u->z,&a->z);
      mpz_init_set(&u->n,&b->z);
      if (a->s<2) mpz_mul(&u->n,&u->n,&a->n);
      if (b->s<2) mpz_mul(&u->z,&u->z,&b->n);
    }
  }
  if (mpz_isNeg(&u->n))
  {
    mpz_neg(&u->z,&u->z);
    mpz_neg(&u->n,&u->n);
  }
  if (mpz_cmp_si(&u->n,(long)1)==0)
  {
    mpz_clear(&u->n);
    if (mpz_size1(&u->z)<=MP_SMALL)
    {
      int ui=(int)mpz_get_si(&u->z);
      if ((((ui<<3)>>3)==ui)
      && (mpz_cmp_si(&u->z,(long)ui)==0))
      {
        mpz_clear(&u->z);
        omFreeBin((ADDRESS)u, rnumber_bin);
        return INT_TO_SR(ui);
      }
    }
    u->s=3;
  }
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

/*2
* u:= x ^ exp
*/
void nlPower (number x,int exp,number * u)
{
  if (!nlIsZero(x))
  {
#ifdef LDEBUG
    nlTest(x);
#endif
    number aa=NULL;
    if (SR_HDL(x) & SR_INT)
    {
      aa=nlRInit(SR_TO_INT(x));
      x=aa;
    }
    else if (x->s==0)
      nlNormalize(x);
    *u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
    (*u)->debug=123456;
#endif
    mpz_init(&(*u)->z);
    mpz_pow_ui(&(*u)->z,&x->z,(unsigned long)exp);
    if (x->s<2)
    {
      if (mpz_cmp_si(&x->n,(long)1)==0)
      {
        x->s=3;
        mpz_clear(&x->n);
      }
      else
      {
        mpz_init(&(*u)->n);
        mpz_pow_ui(&(*u)->n,&x->n,(unsigned long)exp);
      }
    }
    (*u)->s = x->s;
    if (((*u)->s==3) && (mpz_size1(&(*u)->z)<=MP_SMALL))
    {
      int ui=(int)mpz_get_si(&(*u)->z);
      if ((((ui<<3)>>3)==ui)
      && (mpz_cmp_si(&(*u)->z,(long)ui)==0))
      {
        mpz_clear(&(*u)->z);
        omFreeBin((ADDRESS)*u, rnumber_bin);
        *u=INT_TO_SR(ui);
      }
    }
    if (aa!=NULL)
    {
      mpz_clear(&aa->z);
      omFreeBin((ADDRESS)aa, rnumber_bin);
    }
  }
  else
    *u = INT_TO_SR(0);
#ifdef LDEBUG
  nlTest(*u);
#endif
}


/*2
* za >= 0 ?
*/
BOOLEAN nlGreaterZero (number a)
{
#ifdef LDEBUG
  nlTest(a);
#endif
  if (SR_HDL(a) & SR_INT) return SR_HDL(a)>1 /* represents number(0) */;
  return (!mpz_isNeg(&a->z));
}

/*2
* a > b ?
*/
BOOLEAN nlGreater (number a, number b)
{
#ifdef LDEBUG
  nlTest(a);
  nlTest(b);
#endif
  number r;
  BOOLEAN rr;
  r=nlSub(a,b);
  rr=(!nlIsZero(r)) && (nlGreaterZero(r));
  nlDelete(&r,currRing);
  return rr;
}

/*2
* a == -1 ?
*/
BOOLEAN nlIsMOne (number a)
{
#ifdef LDEBUG
  if (a==NULL) return FALSE;
  nlTest(a);
#endif
  if (SR_HDL(a) & SR_INT) return (a==INT_TO_SR(-1));
  return FALSE;
}

/*2
* result =gcd(a,b)
*/
number nlGcd(number a, number b, const ring r)
{
  number result;
#ifdef LDEBUG
  nlTest(a);
  nlTest(b);
#endif
  //nlNormalize(a);
  //nlNormalize(b);
  if ((a==INT_TO_SR(1))||(a==INT_TO_SR(-1))
  ||  (b==INT_TO_SR(1))||(b==INT_TO_SR(-1)))
    return INT_TO_SR(1);
  if (a==INT_TO_SR(0)) /* gcd(0,b) ->b */
    return nlCopy(b);
  if (b==INT_TO_SR(0)) /* gcd(a,0) -> a */
    return nlCopy(a);
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    int i=SR_TO_INT(a);
    int j=SR_TO_INT(b);
    if((i==0)||(j==0))
      return INT_TO_SR(1);
    int l;
    i=ABS(i);
    j=ABS(j);
    do
    {
      l=i%j;
      i=j;
      j=l;
    } while (l!=0);
    return INT_TO_SR(i);
  }
  if (((!(SR_HDL(a) & SR_INT))&&(a->s<2))
  ||  ((!(SR_HDL(b) & SR_INT))&&(b->s<2))) return INT_TO_SR(1);
  result=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  result->debug=123456;
#endif
  mpz_init(&result->z);
  if (SR_HDL(a) & SR_INT)
  {
    unsigned long t=mpz_gcd_ui(NULL,&b->z,ABS(SR_TO_INT(a)));
    return INT_TO_SR((int)t);
  }
  else
  if (SR_HDL(b) & SR_INT)
  {
    unsigned long t=mpz_gcd_ui(NULL,&a->z,ABS(SR_TO_INT(b)));
    return INT_TO_SR((int)t);
  }
  else
  {
    mpz_gcd(&result->z,&a->z,&b->z);
  }
  result->s = 3;
  if (mpz_size1(&result->z)<=MP_SMALL)
  {
    int ui=(int)mpz_get_si(&result->z);
    if ((((ui<<3)>>3)==ui)
    && (mpz_cmp_si(&result->z,(long)ui)==0))
    {
      mpz_clear(&result->z);
      omFreeBin((ADDRESS)result, rnumber_bin);
      result=INT_TO_SR(ui);
    }
  }
#ifdef LDEBUG
  nlTest(result);
#endif
  return result;
}

/*2
* simplify x
*/
void nlNormalize (number &x)
{
  if ((SR_HDL(x) & SR_INT) ||(x==NULL))
    return;
#ifdef LDEBUG
  if (!nlTest(x)) { x->s=1; return; }
#endif
  if (x->s==3)
  {
    if (mpz_cmp_ui(&x->z,(long)0)==0)
    {
      nlDelete(&x,currRing);
      x=INT_TO_SR(0);
      return;
    }
    if (mpz_size1(&x->z)<=MP_SMALL)
    {
      int ui=(int)mpz_get_si(&x->z);
      if ((((ui<<3)>>3)==ui)
      && (mpz_cmp_si(&x->z,(long)ui)==0))
      {
        mpz_clear(&x->z);
        omFreeBin((ADDRESS)x, rnumber_bin);
        x=INT_TO_SR(ui);
        return;
      }
    }
  }
  else if (x->s==0)
  {
    if (mpz_cmp_si(&x->n,(long)1)==0)
    {
      mpz_clear(&x->n);
      if (mpz_size1(&x->z)<=MP_SMALL)
      {
        int ui=(int)mpz_get_si(&x->z);
        if ((((ui<<3)>>3)==ui)
        && (mpz_cmp_si(&x->z,(long)ui)==0))
        {
          mpz_clear(&x->z);
#if defined(LDEBUG)
          x->debug=654324;
#endif
          omFreeBin((ADDRESS)x, rnumber_bin);
          x=INT_TO_SR(ui);
          return;
        }
      }
      x->s=3;
    }
    else
    {
      MP_INT gcd;
      mpz_init(&gcd);
      mpz_gcd(&gcd,&x->z,&x->n);
      x->s=1;
      if (mpz_cmp_si(&gcd,(long)1)!=0)
      {
        MP_INT r;
        mpz_init(&r);
        MPZ_EXACTDIV(&r,&x->z,&gcd);
        mpz_set(&x->z,&r);
        MPZ_EXACTDIV(&r,&x->n,&gcd);
        mpz_set(&x->n,&r);
        mpz_clear(&r);
        if (mpz_cmp_si(&x->n,(long)1)==0)
        {
          mpz_clear(&x->n);
          if (mpz_size1(&x->z)<=MP_SMALL)
          {
            int ui=(int)mpz_get_si(&x->z);
            if ((((ui<<3)>>3)==ui)
            && (mpz_cmp_si(&x->z,(long)ui)==0))
            {
              mpz_clear(&x->z);
#if defined(LDEBUG)
              x->debug=654324;
#endif
              omFreeBin((ADDRESS)x, rnumber_bin);
              x=INT_TO_SR(ui);
              return;
            }
          }
          x->s=3;
        }
      }
      mpz_clear(&gcd);
    }
  }
#ifdef LDEBUG
  nlTest(x);
#endif
}

/*2
* returns in result->z the lcm(a->z,b->n)
*/
number nlLcm(number a, number b, const ring r)
{
  number result;
#ifdef LDEBUG
  nlTest(a);
  nlTest(b);
#endif
  if ((SR_HDL(b) & SR_INT)
  || (b->s==3))
  {
    // b is 1/(b->n) => b->n is 1 => result is a
    return nlCopy(a);
  }
  result=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  result->debug=123456;
#endif
  result->s=3;
  MP_INT gcd;
  mpz_init(&gcd);
  mpz_init(&result->z);
  if (SR_HDL(a) & SR_INT)
    mpz_gcd_ui(&gcd,&b->n,ABS(SR_TO_INT(a)));
  else
    mpz_gcd(&gcd,&a->z,&b->n);
  if (mpz_cmp_si(&gcd,(long)1)!=0)
  {
    MP_INT bt;
    mpz_init_set(&bt,&b->n);
    MPZ_EXACTDIV(&bt,&bt,&gcd);
    if (SR_HDL(a) & SR_INT)
      mpz_mul_si(&result->z,&bt,SR_TO_INT(a));
    else
      mpz_mul(&result->z,&bt,&a->z);
    mpz_clear(&bt);
  }
  else
    if (SR_HDL(a) & SR_INT)
      mpz_mul_si(&result->z,&b->n,SR_TO_INT(a));
    else
      mpz_mul(&result->z,&b->n,&a->z);
  mpz_clear(&gcd);
  if (mpz_size1(&result->z)<=MP_SMALL)
  {
    int ui=(int)mpz_get_si(&result->z);
    if ((((ui<<3)>>3)==ui)
    && (mpz_cmp_si(&result->z,(long)ui)==0))
    {
      mpz_clear(&result->z);
      omFreeBin((ADDRESS)result, rnumber_bin);
      return INT_TO_SR(ui);
    }
  }
#ifdef LDEBUG
  nlTest(result);
#endif
  return result;
}

int nlModP(number n, int p)
{
  if (SR_HDL(n) & SR_INT)
  {
    int i=SR_TO_INT(n);
    if (i<0) return (p-((-i)%p));
    return i%p;
  }
  int iz=(int)mpz_mmod_ui(NULL,&n->z,(unsigned long)p);
  if (n->s!=3)
  {
    int in=mpz_mmod_ui(NULL,&n->n,(unsigned long)p);
    #ifdef NV_OPS
    if (npPrimeM>NV_MAX_PRIME)
    return (int)((long)nvDiv((number)iz,(number)in));
    #endif
    return (int)((long)npDiv((number)iz,(number)in));
  }
  return iz;
}

/*2
* acces to denominator, other 1 for integers
*/
number   nlGetDenom(number &n, const ring r)
{
  if (!(SR_HDL(n) & SR_INT))
  {
    if (n->s==0)
    {
      nlNormalize(n);
    }
    if (!(SR_HDL(n) & SR_INT))
    {
      if (n->s!=3)
      {
        number u=(number)omAllocBin(rnumber_bin);
        u->s=3;
#if defined(LDEBUG)
        u->debug=123456;
#endif
        mpz_init_set(&u->z,&n->n);
        {
          int ui=(int)mpz_get_si(&u->z);
          if ((((ui<<3)>>3)==ui)
          && (mpz_cmp_si(&u->z,(long)ui)==0))
          {
            mpz_clear(&u->z);
            omFreeBin((ADDRESS)u, rnumber_bin);
            return INT_TO_SR(ui);
          }
        }
        return u;
      }
    }
  }
  return INT_TO_SR(1);
}

/*2
* acces to Nominator, nlCopy(n) for integers
*/
number   nlGetNom(number &n, const ring r)
{
  if (!(SR_HDL(n) & SR_INT))
  {
    if (n->s==0)
    {
      nlNormalize(n);
    }
    if (!(SR_HDL(n) & SR_INT))
    {
      number u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
      u->debug=123456;
#endif
      u->s=3;
      mpz_init_set(&u->z,&n->z);
      if (n->s!=3)
      {
        int ui=(int)mpz_get_si(&u->z);
        if ((((ui<<3)>>3)==ui)
        && (mpz_cmp_si(&u->z,(long)ui)==0))
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(ui);
        }
      }
      return u;
    }
  }
  return n; // imm. int
}

/***************************************************************
 *
 * routines which are needed by Inline(d) routines
 *
 *******************************************************************/
BOOLEAN _nlEqual_aNoImm_OR_bNoImm(number a, number b)
{
  assume(! (SR_HDL(a) & SR_HDL(b) & SR_INT));
//  long - short
  BOOLEAN bo;
  if (SR_HDL(b) & SR_INT)
  {
    if (a->s!=0) return FALSE;
    number n=b; b=a; a=n;
  }
//  short - long
  if (SR_HDL(a) & SR_INT)
  {
    if (b->s!=0)
      return FALSE;
    if (((long)a > 0L) && (mpz_isNeg(&b->z)))
      return FALSE;
    if (((long)a < 0L) && (!mpz_isNeg(&b->z)))
      return FALSE;
    MP_INT  bb;
    mpz_init_set(&bb,&b->n);
    mpz_mul_si(&bb,&bb,(long)SR_TO_INT(a));
    bo=(mpz_cmp(&bb,&b->z)==0);
    mpz_clear(&bb);
    return bo;
  }
// long - long
  if (((a->s==1) && (b->s==3))
  ||  ((b->s==1) && (a->s==3)))
    return FALSE;
  if (mpz_isNeg(&a->z)&&(!mpz_isNeg(&b->z)))
    return FALSE;
  if (mpz_isNeg(&b->z)&&(!mpz_isNeg(&a->z)))
    return FALSE;
  MP_INT  aa;
  MP_INT  bb;
  mpz_init_set(&aa,&a->z);
  mpz_init_set(&bb,&b->z);
  if (a->s<2) mpz_mul(&bb,&bb,&a->n);
  if (b->s<2) mpz_mul(&aa,&aa,&b->n);
  bo=(mpz_cmp(&aa,&bb)==0);
  mpz_clear(&aa);
  mpz_clear(&bb);
  return bo;
}

// copy not immediate number a
number _nlCopy_NoImm(number a)
{
  assume(!((SR_HDL(a) & SR_INT)||(a==NULL)));
#ifdef LDEBUG
  nlTest(a);
#endif
  number b=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  b->debug=123456;
#endif
  switch (a->s)
  {
    case 0:
    case 1:
            mpz_init_set(&b->n,&a->n);
    case 3:
            mpz_init_set(&b->z,&a->z);
            break;
  }
  b->s = a->s;
#ifdef LDEBUG
  nlTest(b);
#endif
  return b;
}

void _nlDelete_NoImm(number *a, const ring r)
{
  {
    switch ((*a)->s)
    {
      case 0:
      case 1:
        mpz_clear(&(*a)->n);
      case 3:
        mpz_clear(&(*a)->z);
#ifdef LDEBUG
        (*a)->s=2;
#endif
    }
    omFreeBin((ADDRESS) *a, rnumber_bin);
  }
}

number _nlNeg_NoImm(number a)
{
  {
    mpz_neg(&a->z,&a->z);
    if ((a->s==3) && (mpz_size1(&a->z)<=MP_SMALL))
    {
      int ui=(int)mpz_get_si(&a->z);
      if ((((ui<<3)>>3)==ui)
        && (mpz_cmp_si(&a->z,(long)ui)==0))
      {
        mpz_clear(&a->z);
        omFreeBin((ADDRESS)a, rnumber_bin);
        a=INT_TO_SR(ui);
      }
    }
  }
#ifdef LDEBUG
  nlTest(a);
#endif
  return a;
}

number _nlAdd_aNoImm_OR_bNoImm(number a, number b)
{
  number u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  u->debug=123456;
#endif
  mpz_init(&u->z);
  if (SR_HDL(b) & SR_INT)
  {
    number x=a;
    a=b;
    b=x;
  }
  if (SR_HDL(a) & SR_INT)
  {
    switch (b->s)
    {
      case 0:
      case 1:/* a:short, b:1 */
      {
        MP_INT x;
        mpz_init(&x);
        mpz_mul_si(&x,&b->n,SR_TO_INT(a));
        mpz_add(&u->z,&b->z,&x);
        mpz_clear(&x);
        if (mpz_cmp_ui(&u->z,(long)0)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(0);
        }
        if (mpz_cmp(&u->z,&b->n)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(1);
        }
        mpz_init_set(&u->n,&b->n);
        u->s = 0;
        break;
      }
      case 3:
      {
        if ((long)a>0L)
          mpz_add_ui(&u->z,&b->z,SR_TO_INT(a));
        else
          mpz_sub_ui(&u->z,&b->z,-SR_TO_INT(a));
        if (mpz_cmp_ui(&u->z,(long)0)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(0);
        }
        //u->s = 3;
        if (mpz_size1(&u->z)<=MP_SMALL)
        {
          int ui=(int)mpz_get_si(&u->z);
          if ((((ui<<3)>>3)==ui)
          && (mpz_cmp_si(&u->z,(long)ui)==0))
          {
            mpz_clear(&u->z);
            omFreeBin((ADDRESS)u, rnumber_bin);
            return INT_TO_SR(ui);
          }
        }
        u->s = 3;
        break;
      }
    }
  }
  else
  {
    switch (a->s)
    {
      case 0:
      case 1:
      {
        switch(b->s)
        {
          case 0:
          case 1:
          {
            MP_INT x;
            MP_INT y;
            mpz_init(&x);
            mpz_init(&y);
            mpz_mul(&x,&b->z,&a->n);
            mpz_mul(&y,&a->z,&b->n);
            mpz_add(&u->z,&x,&y);
            mpz_clear(&x);
            mpz_clear(&y);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            mpz_init(&u->n);
            mpz_mul(&u->n,&a->n,&b->n);
            if (mpz_cmp(&u->z,&u->n)==0)
            {
               mpz_clear(&u->z);
               mpz_clear(&u->n);
               omFreeBin((ADDRESS)u, rnumber_bin);
               return INT_TO_SR(1);
            }
            u->s = 0;
            break;
          }
          case 3: /* a:1 b:3 */
          {
            MP_INT x;
            mpz_init(&x);
            mpz_mul(&x,&b->z,&a->n);
            mpz_add(&u->z,&a->z,&x);
            mpz_clear(&x);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            if (mpz_cmp(&u->z,&a->n)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(1);
            }
            mpz_init_set(&u->n,&a->n);
            u->s = 0;
            break;
          }
        } /*switch (b->s) */
        break;
      }
      case 3:
      {
        switch(b->s)
        {
          case 0:
          case 1:/* a:3, b:1 */
          {
            MP_INT x;
            mpz_init(&x);
            mpz_mul(&x,&a->z,&b->n);
            mpz_add(&u->z,&b->z,&x);
            mpz_clear(&x);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            if (mpz_cmp(&u->z,&b->n)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(1);
            }
            mpz_init_set(&u->n,&b->n);
            u->s = 0;
            break;
          }
          case 3:
          {
            mpz_add(&u->z,&a->z,&b->z);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            if (mpz_size1(&u->z)<=MP_SMALL)
            {
              int ui=(int)mpz_get_si(&u->z);
              if ((((ui<<3)>>3)==ui)
              && (mpz_cmp_si(&u->z,(long)ui)==0))
              {
                mpz_clear(&u->z);
                omFreeBin((ADDRESS)u, rnumber_bin);
                return INT_TO_SR(ui);
              }
            }
            u->s = 3;
            break;
          }
        }
        break;
      }
    }
  }
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

number _nlSub_aNoImm_OR_bNoImm(number a, number b)
{
  number u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  u->debug=123456;
#endif
  mpz_init(&u->z);
  if (SR_HDL(a) & SR_INT)
  {
    switch (b->s)
    {
      case 0:
      case 1:/* a:short, b:1 */
      {
        MP_INT x;
        mpz_init(&x);
        mpz_mul_si(&x,&b->n,SR_TO_INT(a));
        mpz_sub(&u->z,&x,&b->z);
        mpz_clear(&x);
        if (mpz_cmp_ui(&u->z,(long)0)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(0);
        }
        if (mpz_cmp(&u->z,&b->n)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(1);
        }
        mpz_init_set(&u->n,&b->n);
        u->s = 0;
        break;
      }
      case 3:
      {
        if ((long)a>0L)
        {
          mpz_sub_ui(&u->z,&b->z,SR_TO_INT(a));
          mpz_neg(&u->z,&u->z);
        }
        else
        {
          mpz_add_ui(&u->z,&b->z,-SR_TO_INT(a));
          mpz_neg(&u->z,&u->z);
        }
        if (mpz_cmp_ui(&u->z,(long)0)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(0);
        }
        if (mpz_size1(&u->z)<=MP_SMALL)
        {
          int ui=(int)mpz_get_si(&u->z);
          if ((((ui<<3)>>3)==ui)
          && (mpz_cmp_si(&u->z,(long)ui)==0))
          {
            mpz_clear(&u->z);
            omFreeBin((ADDRESS)u, rnumber_bin);
            return INT_TO_SR(ui);
          }
        }
        u->s = 3;
        break;
      }
    }
  }
  else if (SR_HDL(b) & SR_INT)
  {
    switch (a->s)
    {
      case 0:
      case 1:/* b:short, a:1 */
      {
        MP_INT x;
        mpz_init(&x);
        mpz_mul_si(&x,&a->n,SR_TO_INT(b));
        mpz_sub(&u->z,&a->z,&x);
        mpz_clear(&x);
        if (mpz_cmp_ui(&u->z,(long)0)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(0);
        }
        if (mpz_cmp(&u->z,&a->n)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(1);
        }
        mpz_init_set(&u->n,&a->n);
        u->s = 0;
        break;
      }
      case 3:
      {
        if ((long)b>0L)
        {
          mpz_sub_ui(&u->z,&a->z,SR_TO_INT(b));
        }
        else
        {
          mpz_add_ui(&u->z,&a->z,-SR_TO_INT(b));
        }
        if (mpz_cmp_ui(&u->z,(long)0)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(0);
        }
        //u->s = 3;
        if (mpz_size1(&u->z)<=MP_SMALL)
        {
          int ui=(int)mpz_get_si(&u->z);
          if ((((ui<<3)>>3)==ui)
          && (mpz_cmp_si(&u->z,(long)ui)==0))
          {
            mpz_clear(&u->z);
            omFreeBin((ADDRESS)u, rnumber_bin);
            return INT_TO_SR(ui);
          }
        }
        u->s = 3;
        break;
      }
    }
  }
  else
  {
    switch (a->s)
    {
      case 0:
      case 1:
      {
        switch(b->s)
        {
          case 0:
          case 1:
          {
            MP_INT x;
            MP_INT y;
            mpz_init(&x);
            mpz_init(&y);
            mpz_mul(&x,&b->z,&a->n);
            mpz_mul(&y,&a->z,&b->n);
            mpz_sub(&u->z,&y,&x);
            mpz_clear(&x);
            mpz_clear(&y);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            mpz_init(&u->n);
            mpz_mul(&u->n,&a->n,&b->n);
            if (mpz_cmp(&u->z,&u->n)==0)
            {
              mpz_clear(&u->z);
              mpz_clear(&u->n);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(1);
            }
            u->s = 0;
            break;
          }
          case 3: /* a:1, b:3 */
          {
            MP_INT x;
            mpz_init(&x);
            mpz_mul(&x,&b->z,&a->n);
            mpz_sub(&u->z,&a->z,&x);
            mpz_clear(&x);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            if (mpz_cmp(&u->z,&a->n)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(1);
            }
            mpz_init_set(&u->n,&a->n);
            u->s = 0;
            break;
          }
        }
        break;
      }
      case 3:
      {
        switch(b->s)
        {
          case 0:
          case 1: /* a:3, b:1 */
          {
            MP_INT x;
            mpz_init(&x);
            mpz_mul(&x,&a->z,&b->n);
            mpz_sub(&u->z,&x,&b->z);
            mpz_clear(&x);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            if (mpz_cmp(&u->z,&b->n)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(1);
            }
            mpz_init_set(&u->n,&b->n);
            u->s = 0;
            break;
          }
          case 3: /* a:3 , b:3 */
          {
            mpz_sub(&u->z,&a->z,&b->z);
            if (mpz_cmp_ui(&u->z,(long)0)==0)
            {
              mpz_clear(&u->z);
              omFreeBin((ADDRESS)u, rnumber_bin);
              return INT_TO_SR(0);
            }
            //u->s = 3;
            if (mpz_size1(&u->z)<=MP_SMALL)
            {
              int ui=(int)mpz_get_si(&u->z);
              if ((((ui<<3)>>3)==ui)
              && (mpz_cmp_si(&u->z,(long)ui)==0))
              {
                mpz_clear(&u->z);
                omFreeBin((ADDRESS)u, rnumber_bin);
                return INT_TO_SR(ui);
              }
            }
            u->s = 3;
            break;
          }
        }
        break;
      }
    }
  }
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

// a and b are intermediate, but a*b not
number _nlMult_aImm_bImm_rNoImm(number a, number b)
{
  number u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  u->debug=123456;
#endif
  u->s=3;
  mpz_init_set_si(&u->z,(long)SR_TO_INT(a));
  mpz_mul_si(&u->z,&u->z,(long)SR_TO_INT(b));
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

// a or b are not immediate
number _nlMult_aNoImm_OR_bNoImm(number a, number b)
{
  assume(! (SR_HDL(a) & SR_HDL(b) & SR_INT));
  number u=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  u->debug=123456;
#endif
  mpz_init(&u->z);
  if (SR_HDL(b) & SR_INT)
  {
    number x=a;
    a=b;
    b=x;
  }
  if (SR_HDL(a) & SR_INT)
  {
    u->s=b->s;
    if (u->s==1) u->s=0;
    if ((long)a>0L)
    {
      mpz_mul_ui(&u->z,&b->z,(unsigned long)SR_TO_INT(a));
    }
    else
    {
      if (a==INT_TO_SR(-1))
      {
        mpz_set(&u->z,&b->z);
        mpz_neg(&u->z,&u->z);
        u->s=b->s;
      }
      else
      {
        mpz_mul_ui(&u->z,&b->z,(unsigned long)-SR_TO_INT(a));
        mpz_neg(&u->z,&u->z);
      }
    }
    if (u->s<2)
    {
      if (mpz_cmp(&u->z,&b->n)==0)
      {
        mpz_clear(&u->z);
        omFreeBin((ADDRESS)u, rnumber_bin);
        return INT_TO_SR(1);
      }
      mpz_init_set(&u->n,&b->n);
    }
    else //u->s==3
    {
      if (mpz_size1(&u->z)<=MP_SMALL)
      {
        int ui=(int)mpz_get_si(&u->z);
        if ((((ui<<3)>>3)==ui)
            && (mpz_cmp_si(&u->z,(long)ui)==0))
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(ui);
        }
      }
    }
  }
  else
  {
    mpz_mul(&u->z,&a->z,&b->z);
    u->s = 0;
    if(a->s==3)
    {
      if(b->s==3)
      {
        u->s = 3;
      }
      else
      {
        if (mpz_cmp(&u->z,&b->n)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(1);
        }
        mpz_init_set(&u->n,&b->n);
      }
    }
    else
    {
      if(b->s==3)
      {
        if (mpz_cmp(&u->z,&a->n)==0)
        {
          mpz_clear(&u->z);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(1);
        }
        mpz_init_set(&u->n,&a->n);
      }
      else
      {
        mpz_init(&u->n);
        mpz_mul(&u->n,&a->n,&b->n);
        if (mpz_cmp(&u->z,&u->n)==0)
        {
          mpz_clear(&u->z);
          mpz_clear(&u->n);
          omFreeBin((ADDRESS)u, rnumber_bin);
          return INT_TO_SR(1);
        }
      }
    }
  }
#ifdef LDEBUG
  nlTest(u);
#endif
  return u;
}

/*2
* z := i
*/
number nlRInit (int i)
{
  number z=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  z->debug=123456;
#endif
  mpz_init_set_si(&z->z,(long)i);
  z->s = 3;
  return z;
}

/*2
* z := i/j
*/
number nlInit2 (int i, int j)
{
  number z=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  z->debug=123456;
#endif
  mpz_init_set_si(&z->z,(long)i);
  mpz_init_set_si(&z->n,(long)j);
  z->s = 0;
  nlNormalize(z);
  return z;
}

number nlInit2gmp (mpz_t i, mpz_t j)
{
  number z=(number)omAllocBin(rnumber_bin);
#if defined(LDEBUG)
  z->debug=123456;
#endif
  mpz_init_set(&z->z,i);
  mpz_init_set(&z->n,j);
  z->s = 0;
  nlNormalize(z);
  return z;
}


#else // DO_LINLINE

// declare immedate routines
number nlRInit (int i);
BOOLEAN _nlEqual_aNoImm_OR_bNoImm(number a, number b);
number  _nlCopy_NoImm(number a);
void    _nlDelete_NoImm(number *a, const ring r);
number  _nlNeg_NoImm(number a);
number  _nlAdd_aNoImm_OR_bNoImm(number a, number b);
number  _nlSub_aNoImm_OR_bNoImm(number a, number b);
number  _nlMult_aNoImm_OR_bNoImm(number a, number b);
number  _nlMult_aImm_bImm_rNoImm(number a, number b);

#endif


/***************************************************************
 *
 * Routines which might be inlined by p_Numbers.h
 *
 *******************************************************************/
#if defined(DO_LINLINE) || !defined(P_NUMBERS_H)

// routines which are always inlined/static

/*2
* a = b ?
*/
LINLINE BOOLEAN nlEqual (number a, number b)
{
#ifdef LDEBUG
  nlTest(a);
  nlTest(b);
#endif
// short - short
  if (SR_HDL(a) & SR_HDL(b) & SR_INT) return a==b;
  return _nlEqual_aNoImm_OR_bNoImm(a, b);
}


LINLINE number nlInit (int i)
{
  number n;
  if ( ((i << 3) >> 3) == i ) n=INT_TO_SR(i);
  else                        n=nlRInit(i);
#ifdef LDEBUG
  nlTest(n);
#endif
  return n;
}


/*2
* a == 1 ?
*/
LINLINE BOOLEAN nlIsOne (number a)
{
#ifdef LDEBUG
  if (a==NULL) return FALSE;
  nlTest(a);
#endif
  if (SR_HDL(a) & SR_INT) return (a==INT_TO_SR(1));
  return FALSE;
}

LINLINE BOOLEAN nlIsZero (number a)
{
  return (a==INT_TO_SR(0));
}

/*2
* copy a to b
*/
LINLINE number nlCopy(number a)
{
  if ((SR_HDL(a) & SR_INT)||(a==NULL))
  {
    return a;
  }
  return _nlCopy_NoImm(a);
}


LINLINE void nlNew (number * r)
{
  *r=NULL;
}

/*2
* delete a
*/
LINLINE void nlDelete (number * a, const ring r)
{
  if (*a!=NULL)
  {
#ifdef LDEBUG
    nlTest(*a);
#endif
    if ((SR_HDL(*a) & SR_INT)==0)
    {
      _nlDelete_NoImm(a,r);
    }
    *a=NULL;
  }
}

/*2
* za:= - za
*/
LINLINE number nlNeg (number a)
{
#ifdef LDEBUG
  nlTest(a);
#endif
  if(SR_HDL(a) &SR_INT)
  {
    int r=SR_TO_INT(a);
    if (r==(-(1<<28))) a=nlRInit(1<<28);
    else               a=INT_TO_SR(-r);
    return a;
  }
  return _nlNeg_NoImm(a);
}

/*2
* u:= a + b
*/
LINLINE number nlAdd (number a, number b)
{
  number u;
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    int r=SR_HDL(a)+SR_HDL(b)-1;
    if ( ((r << 1) >> 1) == r )
      return (number)r;
    else
      return nlRInit(SR_TO_INT(r));
  }
  return _nlAdd_aNoImm_OR_bNoImm(a, b);
}

LINLINE number nlMult (number a, number b)
{
#ifdef LDEBUG
  nlTest(a);
  nlTest(b);
#endif
  if (a==INT_TO_SR(0)) return INT_TO_SR(0);
  if (b==INT_TO_SR(0)) return INT_TO_SR(0);
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    int r=(SR_HDL(a)-1)*(SR_HDL(b)>>1);
    if ((r/(SR_HDL(b)>>1))==(SR_HDL(a)-1))
    {
      number u=((number) ((r>>1)+SR_INT));
      if (((SR_HDL(u)<<1)>>1)==SR_HDL(u)) return (u);
      return nlRInit(SR_HDL(u)>>2);
    }
    return _nlMult_aImm_bImm_rNoImm(a, b);
  }
  return _nlMult_aNoImm_OR_bNoImm(a, b);
}


/*2
* u:= a - b
*/
LINLINE number nlSub (number a, number b)
{
  if (SR_HDL(a) & SR_HDL(b) & SR_INT)
  {
    int r=SR_HDL(a)-SR_HDL(b)+1;
    if ( ((r << 1) >> 1) == r )
    {
      return (number)r;
    }
    else
      return nlRInit(SR_TO_INT(r));
  }
  return _nlSub_aNoImm_OR_bNoImm(a, b);
}

#endif // DO_LINLINE

#ifndef P_NUMBERS_H

void nlInpGcd(number &a, number b, ring r)
{
  if ((SR_HDL(b)|SR_HDL(a))&SR_INT)
  {
    number n=nlGcd(a,b,r);
    nlDelete(&a,r);
    a=n;
  }
  else
  {
    mpz_gcd(&a->z,&a->z,&b->z);
    if (mpz_size1(&a->z)<=MP_SMALL)
    {
      int ui=(int)mpz_get_si(&a->z);
      if ((((ui<<3)>>3)==ui)
      && (mpz_cmp_si(&a->z,(long)ui)==0))
      {
        mpz_clear(&a->z);
        omFreeBin((ADDRESS)a, rnumber_bin);
        a=INT_TO_SR(ui);
      }
    }
  }
}
void nlInpIntDiv(number &a, number b, ring r)
{
  if ((SR_HDL(b)|SR_HDL(a))&SR_INT)
  {
    number n=nlIntDiv(a,b);
    nlDelete(&a,r);
    a=n;
  }
  else
  {
    if (mpz_isNeg(&a->z))
    {
      if (mpz_isNeg(&b->z))
      {
        mpz_add(&a->z,&a->z,&b->z);
      }
      else
      {
        mpz_sub(&a->z,&a->z,&b->z);
      }
      mpz_add_ui(&a->z,&a->z,1);
    }
    else
    {
      if (mpz_isNeg(&b->z))
      {
        mpz_sub(&a->z,&a->z,&b->z);
      }
      else
      {
        mpz_add(&a->z,&a->z,&b->z);
      }
      mpz_sub_ui(&a->z,&a->z,1);
    }
    MPZ_DIV(&a->z,&a->z,&b->z);
    if (mpz_size1(&a->z)<=MP_SMALL)
    {
      int ui=(int)mpz_get_si(&a->z);
      if ((((ui<<3)>>3)==ui)
      && (mpz_cmp_si(&a->z,(long)ui)==0))
      {
        mpz_clear(&a->z);
        omFreeBin((ADDRESS)a, rnumber_bin);
        a=INT_TO_SR(ui);
      }
    }
  }
}
void nlInpAdd(number &a, number b, ring r)
{
  // TODO
  if ((SR_HDL(b)|SR_HDL(a))&SR_INT)
  {
    number n=nlGcd(a,b,r);
    nlDelete(&a,r);
    a=n;
  }
  else
  {
    mpz_gcd(&a->z,&a->z,&b->z);
    if (mpz_size1(&a->z)<=MP_SMALL)
    {
      int ui=(int)mpz_get_si(&a->z);
      if ((((ui<<3)>>3)==ui)
      && (mpz_cmp_si(&a->z,(long)ui)==0))
      {
        mpz_clear(&a->z);
        omFreeBin((ADDRESS)a, rnumber_bin);
        a=INT_TO_SR(ui);
      }
    }
  }
}
void nlInpMult(number &a, number b, ring r)
{
  if (((SR_HDL(b)|SR_HDL(a))&SR_INT)
  )
  {
    number n=nlMult(a,b);
    nlDelete(&a,r);
    a=n;
  }
  else
  {
    mpz_mul(&a->z,&a->z,&b->z);
    if (a->s==3)
    {
      if(b->s!=3)
      {
        mpz_init_set(&a->n,&b->n);
        a->s=0;
      }
    }
    else
    {
      if(b->s!=3)
      {
        mpz_mul(&a->n,&a->n,&b->n);
      }
      a->s=0;
    }
  }
}

#if 0
number nlMod(number a, number b)
{
  if (((SR_HDL(b)&SR_HDL(a))&SR_INT)
  {
    int bi=SR_TO_INT(b);
    int ai=SR_TO_INT(a);
    int bb=ABS(bi);
    int c=ai%bb;
    if (c<0)  c+=bb;
    return (INT_TO_SR(c));
  }
  number al;
  number bl;
  if (SR_HDL(a))&SR_INT)
    al=nlRInit(SR_TO_INT(a));
  else
    al=nlCopy(a);
  if (SR_HDL(b))&SR_INT)
    bl=nlRInit(SR_TO_INT(b));
  else
    bl=nlCopy(b);
  number r=nlRInit(0);
  mpz_mod(r->z,al->z,bl->z);
  nlDelete(&al);
  nlDelete(&bl);
  if (mpz_size1(&r->z)<=MP_SMALL)
  {
    int ui=(int)mpz_get_si(&r->z);
    if ((((ui<<3)>>3)==ui)
    && (mpz_cmp_si(&x->z,(long)ui)==0))
    {
      mpz_clear(&r->z);
      omFreeBin((ADDRESS)r, rnumber_bin);
      r=INT_TO_SR(ui);
    }
  }
  return r;
}
#endif
#endif // not P_NUMBERS_H
#endif // LONGRAT_CC
