#line 2 "../src/kernel/gmp/mp.c"
/* $Id: mp.c 12385 2010-06-03 14:42:34Z kb $

Copyright (C) 2002-2003  The PARI group.

This file is part of the PARI/GP package.

PARI/GP is free software; you can redistribute it and/or modify it under the
terms of the GNU General Public License as published by the Free Software
Foundation. It is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY WHATSOEVER.

Check the License for details. You should have received a copy of it, along
with the package; see the file 'COPYING'. If not, write to the Free Software
Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA. */

/***********************************************************************/
/**                                                                   **/
/**                               GMP KERNEL                          **/
/** BA2002Sep24                                                       **/
/***********************************************************************/
/* GMP t_INT as just like normal t_INT, just the mantissa is the other way
 * round
 *
 *   `How would you like to live in Looking-glass House, Kitty?  I
 *   wonder if they'd give you milk in there?  Perhaps Looking-glass
 *   milk isn't good to drink--But oh, Kitty! now we come to the
 *   passage.  You can just see a little PEEP of the passage in
 *   Looking-glass House, if you leave the door of our drawing-room
 *   wide open:  and it's very like our passage as far as you can see,
 *   only you know it may be quite different on beyond.  Oh, Kitty!
 *   how nice it would be if we could only get through into Looking-
 *   glass House!  I'm sure it's got, oh! such beautiful things in it!
 *
 *  Through the Looking Glass,  Lewis Carrol
 *
 *  (pityful attempt to beat GN code/comments rate)
 *  */

#include <gmp.h>
#include "pari.h"
#include "paripriv.h"
#include "../src/kernel/none/tune-gen.h"

/*We need PARI invmod renamed to invmod_pari*/
#define INVMOD_PARI

static void *gmp_realloc(void *ptr, size_t old_size, size_t new_size) {
  (void)old_size; return (void *) pari_realloc(ptr,new_size);
}

static void gmp_free(void *ptr, size_t old_size){
  (void)old_size; pari_free(ptr);
}

int pari_kernel_init(void)
{
  /* Use pari_malloc instead of malloc */
  /* patch for Sage
  mp_set_memory_functions((void *(*)(size_t)) pari_malloc, gmp_realloc, gmp_free);
  */
  return 0;
}

#define LIMBS(x)  ((mp_limb_t *)((x)+2))
#define NLIMBS(x) (lgefint(x)-2)
/*This one is for t_REALs to emphasize they are not t_INTs*/
#define RLIMBS(x)  ((mp_limb_t *)((x)+2))
#define RNLIMBS(x) (lg(x)-2)

INLINE void
xmpn_copy(mp_limb_t *x, mp_limb_t *y, long n)
{
  while (--n >= 0) x[n]=y[n];
}

INLINE void
xmpn_mirror(mp_limb_t *x, long n)
{
  long i;
  for(i=0;i<(n>>1);i++)
  {
    ulong m=x[i];
    x[i]=x[n-1-i];
    x[n-1-i]=m;
  }
}

INLINE void
xmpn_mirrorcopy(mp_limb_t *z, mp_limb_t *x, long n)
{
  long i;
  for(i=0;i<n;i++)
    z[i]=x[n-1-i];
}

INLINE void
xmpn_zero(mp_ptr x, mp_size_t n)
{
  while (--n >= 0) x[n]=0;
}

INLINE GEN
icopy_ef(GEN x, long l)
{
  register long lx = lgefint(x);
  const GEN y = cgeti(l);

  while (--lx > 0) y[lx]=x[lx];
  return y;
}

/* NOTE: arguments of "spec" routines (muliispec, addiispec, etc.) aren't
 * GENs but pairs (long *a, long na) representing a list of digits (in basis
 * BITS_IN_LONG) : a[0], ..., a[na-1]. [ In ordre to facilitate splitting: no
 * need to reintroduce codewords ]
 * Use speci(a,na) to visualize the corresponding GEN.
 */

/***********************************************************************/
/**                                                                   **/
/**                     ADDITION / SUBTRACTION                        **/
/**                                                                   **/
/***********************************************************************/

GEN
setloop(GEN a)
{
  pari_sp av = avma - 2 * sizeof(long);
  (void)cgetg(lgefint(a) + 3, t_VECSMALL);
  return icopy_avma(a, av); /* two cells of extra space after a */
}

/* we had a = setloop(?), then some incloops. Reset a to b */
GEN
resetloop(GEN a, GEN b) {
  a[0] = evaltyp(t_INT) | evallg(lgefint(b));
  affii(b, a); return a;
}

/* assume a > 0, initialized by setloop. Do a++ */
static GEN
incpos(GEN a)
{
  long i, l = lgefint(a);
  for (i=2; i<l; i++)
    if (++a[i]) return a;
  a[l] = 1; l++;
  a[0]=evaltyp(t_INT) | _evallg(l);
  a[1]=evalsigne(1) | evallgefint(l);
  return a;
}

/* assume a < 0, initialized by setloop. Do a++ */
static GEN
incneg(GEN a)
{
  long i, l = lgefint(a);
  if (a[2]--)
  {
    if (l == 3 && !a[2])
    {
      a[0] = evaltyp(t_INT) | _evallg(2);
      a[1] = evalsigne(0) | evallgefint(2);
    }
    return a;
  }
  for (i=3; i<l; i++)
    if (a[i]--) break;
  l -= i - 2;
  a[0] = evaltyp(t_INT) | _evallg(l);
  a[1] = evalsigne(-1) | evallgefint(l);
  return a;
}

/* assume a initialized by setloop. Do a++ */
GEN
incloop(GEN a)
{
  switch(signe(a))
  {
    case 0:
      a[0]=evaltyp(t_INT) | _evallg(3);
      a[1]=evalsigne(1) | evallgefint(3);
      a[2]=1; return a;
    case -1: return incneg(a);
    default: return incpos(a);
  }
}

INLINE GEN
adduispec(ulong s, GEN x, long nx)
{
  GEN  zd;
  long lz;

  if (nx == 1) return adduu((ulong)x[0], s);
  lz = nx+3; zd = cgeti(lz);
  if (mpn_add_1(LIMBS(zd),(mp_limb_t *)x,nx,s))
    zd[lz-1]=1;
  else
    lz--;
  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

GEN
adduispec_offset(ulong s, GEN x, long offset, long nx)
{
  GEN xd=x+2+offset;
  while (nx && *(xd+nx-1)==0) nx--;
  if (!nx) return utoi(s);
  return adduispec(s,xd,nx);
}

INLINE GEN
addiispec(GEN x, GEN y, long nx, long ny)
{
  GEN zd;
  long lz;

  if (nx < ny) swapspec(x,y, nx,ny);
  if (ny == 1) return adduispec(*y,x,nx);
  lz = nx+3; zd = cgeti(lz);

  if (mpn_add(LIMBS(zd),(mp_limb_t *)x,nx,(mp_limb_t *)y,ny))
    zd[lz-1]=1;
  else
    lz--;

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

/* assume x >= y */
INLINE GEN
subiuspec(GEN x, ulong s, long nx)
{
  GEN zd;
  long lz;

  if (nx == 1) return utoi(x[0] - s);

  lz = nx + 2; zd = cgeti(lz);
  mpn_sub_1 (LIMBS(zd), (mp_limb_t *)x, nx, s);
  if (! zd[lz - 1]) { --lz; }

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

/* assume x > y */
INLINE GEN
subiispec(GEN x, GEN y, long nx, long ny)
{
  GEN zd;
  long lz;
  if (ny==1) return subiuspec(x,*y,nx);
  lz = nx+2; zd = cgeti(lz);

  mpn_sub (LIMBS(zd), (mp_limb_t *)x, nx, (mp_limb_t *)y, ny);
  while (lz >= 3 && zd[lz - 1] == 0) { lz--; }

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

static void
roundr_up_ip(GEN x, long l)
{
  long i = l;
  for(;;)
  {
    if (++x[--i]) break;
    if (i == 2) { x[2] = HIGHBIT; setexpo(x, expo(x)+1); break; }
  }
}

void
affir(GEN x, GEN y)
{
  const long s = signe(x), ly = lg(y);
  long lx, sh, i;

  if (!s)
  {
    y[1] = evalexpo(-bit_accuracy(ly));
    return;
  }
  lx = lgefint(x); sh = bfffo(*int_MSW(x));
  y[1] = evalsigne(s) | evalexpo(bit_accuracy(lx)-sh-1);
  if (sh) {
    if (lx <= ly)
    {
      for (i=lx; i<ly; i++) y[i]=0;
      mpn_lshift(LIMBS(y),LIMBS(x),lx-2,sh);
      xmpn_mirror(LIMBS(y),lx-2);
      return;
    }
    mpn_lshift(LIMBS(y),LIMBS(x)+lx-ly,ly-2,sh);
    y[2]|=((ulong) x[lx-ly+1])>>(BITS_IN_LONG-sh);
    xmpn_mirror(LIMBS(y),ly-2);
    /* lx > ly: round properly */
    if ((x[lx-ly+1]<<sh) & HIGHBIT) roundr_up_ip(y, ly);
  }
  else {
    GEN xd=int_MSW(x);
    if (lx <= ly)
    {
      for (i=2; i<lx; i++,xd=int_precW(xd)) y[i]=*xd;
      for (   ; i<ly; i++) y[i]=0;
      return;
    }
    for (i=2; i<ly; i++,xd=int_precW(xd)) y[i]=*xd;
    /* lx > ly: round properly */
    if (x[lx-ly+1] & HIGHBIT) roundr_up_ip(y, ly);
  }
}

GEN
shifti(GEN x, long n)
{
  const long s=signe(x);
  long lz,lx,m;
  GEN z;

  if (!s) return gen_0;
  if (!n) return icopy(x);
  lx = lgefint(x);
  if (n > 0)
  {
    long d = dvmdsBIL(n, &m);
    long i;

    lz = lx + d + (m!=0);
    z = cgeti(lz);
    for (i=0; i<d; i++) LIMBS(z)[i] = 0;

    if (!m) xmpn_copy(LIMBS(z)+d, LIMBS(x), NLIMBS(x));
    else
    {
      ulong carry = mpn_lshift(LIMBS(z)+d, LIMBS(x), NLIMBS(x), m);
      if (carry) z[lz - 1] = carry;
      else lz--;
    }
  }
  else
  {
    long d = dvmdsBIL(-n, &m);

    lz = lx - d;
    if (lz<3) return gen_0;
    z = cgeti(lz);

    if (!m) xmpn_copy(LIMBS(z), LIMBS(x) + d, NLIMBS(x) - d);
    else
    {
      mpn_rshift(LIMBS(z), LIMBS(x) + d, NLIMBS(x) - d, m);
      if (z[lz - 1] == 0)
      {
        if (lz == 3) { avma = (pari_sp)(z+3); return gen_0; }
        lz--;
      }
    }
  }
  z[1] = evalsigne(s)|evallgefint(lz);
  return z;
}

GEN
trunc2nr_lg(GEN x, long lx, long n)
{
  long ly, i, m, s = signe(x);
  GEN y;
  if (!s) return gen_0;
  if (!n)
  {
    y = cgeti(lx);
    y[1] = evalsigne(s) | evallgefint(lx);
    xmpn_mirrorcopy(LIMBS(y),RLIMBS(x),lx-2);
    return y;
  }
  if (n > 0)
  {
    GEN z = (GEN)avma;
    long d = dvmdsBIL(n, &m);

    ly = lx+d; y = new_chunk(ly);
    for ( ; d; d--) *--z = 0;
    if (!m) for (i=2; i<lx; i++) y[i]=x[i];
    else
    {
      register const ulong sh = BITS_IN_LONG - m;
      shift_left(y,x, 2,lx-1, 0,m);
      i = ((ulong)x[2]) >> sh;
      /* Extend y on the left? */
      if (i) { ly++; y = new_chunk(1); y[2] = i; }
    }
  }
  else
  {
    ly = lx - dvmdsBIL(-n, &m);
    if (ly<3) return gen_0;
    y = new_chunk(ly);
    if (m) {
      shift_right(y,x, 2,ly, 0,m);
      if (y[2] == 0)
      {
        if (ly==3) { avma = (pari_sp)(y+3); return gen_0; }
        ly--; avma = (pari_sp)(++y);
      }
    } else {
      for (i=2; i<ly; i++) y[i]=x[i];
    }
  }
  xmpn_mirror(LIMBS(y),ly-2);
  y[1] = evalsigne(s)|evallgefint(ly);
  y[0] = evaltyp(t_INT)|evallg(ly); return y;
}

GEN
truncr(GEN x)
{
  long s, e, d, m, i;
  GEN y;
  if ((s=signe(x)) == 0 || (e=expo(x)) < 0) return gen_0;
  d = nbits2prec(e+1); m = remsBIL(e);
  if (d > lg(x)) pari_err(precer, "truncr (precision loss in truncation)");

  y=cgeti(d); y[1] = evalsigne(s) | evallgefint(d);
  if (++m == BITS_IN_LONG)
    for (i=2; i<d; i++) y[d-i+1]=x[i];
  else
  {
    GEN z=cgeti(d);
    for (i=2; i<d; i++) z[d-i+1]=x[i];
    mpn_rshift(LIMBS(y),LIMBS(z),d-2,BITS_IN_LONG-m);
    avma=(pari_sp)y;
  }
  return y;
}

/* integral part */
GEN
floorr(GEN x)
{
  long e, d, m, i, lx;
  GEN y;
  if (signe(x) >= 0) return truncr(x);
  if ((e=expo(x)) < 0) return gen_m1;
  d = nbits2prec(e+1); m = remsBIL(e);
  lx=lg(x); if (d>lx) pari_err(precer, "floorr (precision loss in truncation)");
  y = cgeti(d+1);
  if (++m == BITS_IN_LONG)
  {
    for (i=2; i<d; i++) y[d-i+1]=x[i];
    i=d; while (i<lx && !x[i]) i++;
    if (i==lx) goto END;
  }
  else
  {
    GEN z=cgeti(d);
    for (i=2; i<d; i++) z[d-i+1]=x[i];
    mpn_rshift(LIMBS(y),LIMBS(z),d-2,BITS_IN_LONG-m);
    if (x[d-1]<<m == 0)
    {
      i=d; while (i<lx && !x[i]) i++;
      if (i==lx) goto END;
    }
  }
  if (mpn_add_1(LIMBS(y),LIMBS(y),d-2,1))
    y[d++]=1;
END:
  y[1] = evalsigne(-1) | evallgefint(d);
  return y;
}

INLINE int
absi_cmp_lg(GEN x, GEN y, long l)
{
  return mpn_cmp(LIMBS(x),LIMBS(y),l-2);
}

INLINE int
absi_equal_lg(GEN x, GEN y, long l)
{
  return !mpn_cmp(LIMBS(x),LIMBS(y),l-2);
}

/***********************************************************************/
/**                                                                   **/
/**                          MULTIPLICATION                           **/
/**                                                                   **/
/***********************************************************************/
/* assume ny > 0 */
INLINE GEN
muluispec(ulong x, GEN y, long ny)
{
  if (ny == 1)
    return muluu(x, *y);
  else
  {
    long lz = ny+3;
    GEN z = cgeti(lz);
    ulong hi = mpn_mul_1 (LIMBS(z), (mp_limb_t *)y, ny, x);
    if (hi) { z[lz - 1] = hi; } else lz--;
    z[1] = evalsigne(1) | evallgefint(lz);
    return z;
  }
}

/* a + b*|y| */
GEN
addumului(ulong a, ulong b, GEN y)
{
  GEN z;
  long i, lz;
  ulong hi;
  if (!signe(y)) return utoi(a);
  lz = lgefint(y)+1;
  z = cgeti(lz);
  z[2]=a;
  for(i=3;i<lz;i++) z[i]=0;
  hi=mpn_addmul_1(LIMBS(z), LIMBS(y), NLIMBS(y), b);
  if (hi) z[lz-1]=hi; else lz--;
  z[1] = evalsigne(1) | evallgefint(lz);
  avma=(pari_sp)z; return z;
}

GEN muliispec(GEN x, GEN y, long nx, long ny);

/* We must have nx>=ny. This lets garbage on the stack.
   This handle squares correctly (mpn_mul is optimized
   for squares).
*/

INLINE GEN
muliispec_mirror(GEN x, GEN y, long nx, long ny)
{
  GEN cx=new_chunk(nx),cy;
  GEN z;
  xmpn_mirrorcopy((mp_limb_t *)cx,(mp_limb_t *)x,nx);
  if (x==y) cy=cx; /*If nx<ny cy will be too short*/
  else
  {
    cy=new_chunk(ny);
    xmpn_mirrorcopy((mp_limb_t *)cy,(mp_limb_t *)y,ny);
  }
  z=muliispec(cx, cy, nx, ny);
  xmpn_mirror(LIMBS(z), NLIMBS(z));
  return z;
}

/***********************************************************************/
/**                                                                   **/
/**                          DIVISION                                 **/
/**                                                                   **/
/***********************************************************************/

ulong
umodiu(GEN y, ulong x)
{
  long sy=signe(y);
  ulong hi;
  if (!x) pari_err(gdiver);
  if (!sy) return 0;
  hi = mpn_mod_1(LIMBS(y),NLIMBS(y),x);
  if (!hi) return 0;
  return (sy > 0)? hi: x - hi;
}

/* return |y| \/ x */
GEN
diviu_rem(GEN y, ulong x, ulong *rem)
{
  long ly;
  GEN z;

  if (!x) pari_err(gdiver);
  if (!signe(y)) { *rem = 0; return gen_0; }

  ly = lgefint(y);
  if (ly == 3 && (ulong)x > (ulong)y[2]) { *rem = (ulong)y[2]; return gen_0; }

  z = cgeti(ly);
  *rem = mpn_divrem_1(LIMBS(z), 0, LIMBS(y), NLIMBS(y), x);
  if (z [ly - 1] == 0) ly--;
  z[1] = evallgefint(ly) | evalsigne(1);
  return z;
}

GEN
divis_rem(GEN y, long x, long *rem)
{
  long sy=signe(y),ly,s;
  GEN z;

  if (!x) pari_err(gdiver);
  if (!sy) { *rem = 0; return gen_0; }
  if (x<0) { s = -sy; x = -x; } else s = sy;

  ly = lgefint(y);
  if (ly == 3 && (ulong)x > (ulong)y[2]) { *rem = itos(y); return gen_0; }

  z = cgeti(ly);
  *rem = mpn_divrem_1(LIMBS(z), 0, LIMBS(y), NLIMBS(y), x);
  if (sy<0) *rem = -  *rem;
  if (z[ly - 1] == 0) ly--;
  z[1] = evallgefint(ly) | evalsigne(s);
  return z;
}

GEN
divis(GEN y, long x)
{
  long sy=signe(y),ly,s;
  GEN z;

  if (!x) pari_err(gdiver);
  if (!sy) return gen_0;
  if (x<0) { s = -sy; x = -x; } else s=sy;

  ly = lgefint(y);
  if (ly == 3 && (ulong)x > (ulong)y[2]) return gen_0;

  z = cgeti(ly);
  (void)mpn_divrem_1(LIMBS(z), 0, LIMBS(y), NLIMBS(y), x);
  if (z[ly - 1] == 0) ly--;
  z[1] = evallgefint(ly) | evalsigne(s);
  return z;
}

/* We keep llx bits of x and lly bits of y*/
static GEN
divrr_with_gmp(GEN x, GEN y)
{
  long lx=RNLIMBS(x),ly=RNLIMBS(y);
  long lw=minss(lx,ly);
  long lly=minss(lw+1,ly);
  GEN  w=cgetr(lw+2);
  long lu=lw+lly;
  long llx=minss(lu,lx);
  mp_limb_t *u=(mp_limb_t *)new_chunk(lu);
  mp_limb_t *z=(mp_limb_t *)new_chunk(lly);
  mp_limb_t *q,*r;
  long e=expo(x)-expo(y);
  long sx=signe(x),sy=signe(y);
  xmpn_mirrorcopy(z,RLIMBS(y),lly);
  xmpn_mirrorcopy(u+lu-llx,RLIMBS(x),llx);
  xmpn_zero(u,lu-llx);
  q = (mp_limb_t *)new_chunk(lw+1);
  r = (mp_limb_t *)new_chunk(lly);

  mpn_tdiv_qr(q,r,0,u,lu,z,lly);

  /*Round up: This is not exactly correct we should test 2*r>z*/
  if ((ulong)r[lly-1] > ((ulong)z[lly-1]>>1))
    mpn_add_1(q,q,lw+1,1);

  xmpn_mirrorcopy(RLIMBS(w),q,lw);

  if (q[lw] == 0) e--;
  else if (q[lw] == 1) { shift_right(w,w, 2,lw+2, 1,1); }
  else { w[2] = HIGHBIT; e++; }
  if (sy < 0) sx = -sx;
  w[1] = evalsigne(sx) | evalexpo(e);
  avma=(pari_sp) w;
  return w;
}

/* We keep llx bits of x and lly bits of y*/
static GEN
divri_with_gmp(GEN x, GEN y)
{
  long llx=RNLIMBS(x),ly=NLIMBS(y);
  long lly=minss(llx+1,ly);
  GEN  w=cgetr(llx+2);
  long lu=llx+lly, ld=ly-lly;
  mp_limb_t *u=(mp_limb_t *)new_chunk(lu);
  mp_limb_t *z=(mp_limb_t *)new_chunk(lly);
  mp_limb_t *q,*r;
  long sh=bfffo(y[ly+1]);
  long e=expo(x)-expi(y);
  long sx=signe(x),sy=signe(y);
  if (sh) mpn_lshift(z,LIMBS(y)+ld,lly,sh);
  else xmpn_copy(z,LIMBS(y)+ld,lly);
  xmpn_mirrorcopy(u+lu-llx,RLIMBS(x),llx);
  xmpn_zero(u,lu-llx);
  q = (mp_limb_t *)new_chunk(llx+1);
  r = (mp_limb_t *)new_chunk(lly);

  mpn_tdiv_qr(q,r,0,u,lu,z,lly);

  /*Round up: This is not exactly correct we should test 2*r>z*/
  if ((ulong)r[lly-1] > ((ulong)z[lly-1]>>1))
    mpn_add_1(q,q,llx+1,1);

  xmpn_mirrorcopy(RLIMBS(w),q,llx);

  if (q[llx] == 0) e--;
  else if (q[llx] == 1) { shift_right(w,w, 2,llx+2, 1,1); }
  else { w[2] = HIGHBIT; e++; }
  if (sy < 0) sx = -sx;
  w[1] = evalsigne(sx) | evalexpo(e);
  avma=(pari_sp) w;
  return w;
}

GEN
divri(GEN x, GEN y)
{
  long  s = signe(y);

  if (!s) pari_err(gdiver);
  if (!signe(x)) return real_0_bit(expo(x) - expi(y));
  if (!is_bigint(y)) {
    GEN z = divru(x, y[2]);
    if (s < 0) togglesign(z);
    return z;
  }
  return divri_with_gmp(x,y);
}

GEN
divrr(GEN x, GEN y)
{
  long sx=signe(x), sy=signe(y), lx,ly,lr,e,i,j;
  ulong y0,y1;
  GEN r, r1;

  if (!sy) pari_err(gdiver);
  e = expo(x) - expo(y);
  if (!sx) return real_0_bit(e);
  if (sy<0) sx = -sx;

  lx=lg(x); ly=lg(y);
  if (ly==3)
  {
    ulong k = x[2], l = (lx>3)? x[3]: 0;
    LOCAL_HIREMAINDER;
    if (k < (ulong)y[2]) e--;
    else
    {
      l >>= 1; if (k&1) l |= HIGHBIT;
      k >>= 1;
    }
    r = cgetr(3); r[1] = evalsigne(sx) | evalexpo(e);
    hiremainder=k; r[2]=divll(l,y[2]); return r;
  }

  if (ly>=DIVRR_GMP_LIMIT)
    return divrr_with_gmp(x,y);

  lr = minss(lx,ly); r = new_chunk(lr);
  r1 = r-1;
  r1[1] = 0; for (i=2; i<lr; i++) r1[i]=x[i];
  r1[lr] = (lx>ly)? x[lr]: 0;
  y0 = y[2]; y1 = y[3];
  for (i=0; i<lr-1; i++)
  { /* r1 = r + (i-1), OK up to r1[2] (accesses at most r[lr]) */
    ulong k, qp;
    LOCAL_HIREMAINDER;
    LOCAL_OVERFLOW;

    if ((ulong)r1[1] == y0)
    {
      qp = ULONG_MAX; k = addll(y0,r1[2]);
    }
    else
    {
      if ((ulong)r1[1] > y0) /* can't happen if i=0 */
      {
        GEN y1 = y+1;
        j = lr-i; r1[j] = subll(r1[j],y1[j]);
        for (j--; j>0; j--) r1[j] = subllx(r1[j],y1[j]);
        j=i; do r[--j]++; while (j && !r[j]);
      }
      hiremainder = r1[1]; overflow = 0;
      qp = divll(r1[2],y0); k = hiremainder;
    }
    j = lr-i+1;
    if (!overflow)
    {
      long k3, k4;
      k3 = mulll(qp,y1);
      if (j == 3) /* i = lr - 2 maximal, r1[3] undefined -> 0 */
        k4 = subll(hiremainder,k);
      else
      {
        k3 = subll(k3, r1[3]);
        k4 = subllx(hiremainder,k);
      }
      while (!overflow && k4) { qp--; k3=subll(k3,y1); k4=subllx(k4,y0); }
    }
    if (j<ly) (void)mulll(qp,y[j]); else { hiremainder = 0 ; j = ly; }
    for (j--; j>1; j--)
    {
      r1[j] = subll(r1[j], addmul(qp,y[j]));
      hiremainder += overflow;
    }
    if ((ulong)r1[1] != hiremainder)
    {
      if ((ulong)r1[1] < hiremainder)
      {
        qp--;
        j = lr-i-(lr-i>=ly); r1[j] = addll(r1[j], y[j]);
        for (j--; j>1; j--) r1[j] = addllx(r1[j], y[j]);
      }
      else
      {
        r1[1] -= hiremainder;
        while (r1[1])
        {
          qp++; if (!qp) { j=i; do r[--j]++; while (j && !r[j]); }
          j = lr-i-(lr-i>=ly); r1[j] = subll(r1[j],y[j]);
          for (j--; j>1; j--) r1[j] = subllx(r1[j],y[j]);
          r1[1] -= overflow;
        }
      }
    }
    *++r1 = qp;
  }
  /* i = lr-1 */
  /* round correctly */
  if ((ulong)r1[1] > (y0>>1))
  {
    j=i; do r[--j]++; while (j && !r[j]);
  }
  r1 = r-1; for (j=i; j>=2; j--) r[j]=r1[j];
  if (r[0] == 0) e--;
  else if (r[0] == 1) { shift_right(r,r, 2,lr, 1,1); }
  else { /* possible only when rounding up to 0x2 0x0 ... */
    r[2] = (long)HIGHBIT; e++;
  }
  r[0] = evaltyp(t_REAL)|evallg(lr);
  r[1] = evalsigne(sx) | evalexpo(e);
  return r;
}

/* Integer division x / y: such that sign(r) = sign(x)
 *   if z = ONLY_REM return remainder, otherwise return quotient
 *   if z != NULL set *z to remainder
 *   *z is the last object on stack (and thus can be disposed of with cgiv
 *   instead of gerepile)
 * If *z is zero, we put gen_0 here and no copy.
 * space needed: lx + ly
 */
GEN
dvmdii(GEN x, GEN y, GEN *z)
{
  long sx=signe(x),sy=signe(y);
  long lx, ly, lq;
  pari_sp av;
  GEN r,q;

  if (!sy) { if (z == ONLY_REM && !sx) return gen_0; pari_err(gdiver); }
  if (!sx)
  {
    if (!z || z == ONLY_REM) return gen_0;
    *z=gen_0; return gen_0;
  }
  lx=lgefint(x);
  ly=lgefint(y); lq=lx-ly;
  if (lq <= 0)
  {
    if (lq == 0)
    {
      long s=mpn_cmp(LIMBS(x),LIMBS(y),NLIMBS(x));
      if (s>0) goto DIVIDE;
      if (s==0)
      {
        if (z == ONLY_REM) return gen_0;
        if (z) *z = gen_0;
        if (sx < 0) sy = -sy;
        return stoi(sy);
      }
    }
    if (z == ONLY_REM) return icopy(x);
    if (z) *z = icopy(x);
    return gen_0;
  }
DIVIDE: /* quotient is non-zero */
  av=avma; if (sx<0) sy = -sy;
  if (ly==3)
  {
    ulong lq = lx;
    ulong si;
    q  = cgeti(lq);
    si = mpn_divrem_1(LIMBS(q), 0, LIMBS(x), NLIMBS(x), y[2]);
    if (q[lq - 1] == 0) lq--;
    if (z == ONLY_REM)
    {
      avma=av; if (!si) return gen_0;
      r=cgeti(3);
      r[1] = evalsigne(sx) | evallgefint(3);
      r[2] = si; return r;
    }
    q[1] = evalsigne(sy) | evallgefint(lq);
    if (!z) return q;
    if (!si) { *z=gen_0; return q; }
    r=cgeti(3);
    r[1] = evalsigne(sx) | evallgefint(3);
    r[2] = si; *z=r; return q;
  }
  if (z == ONLY_REM)
  {
    ulong lr = lgefint(y);
    ulong lq = lgefint(x)-lgefint(y)+3;
    GEN r = cgeti(lr);
    GEN q = cgeti(lq);
    mpn_tdiv_qr(LIMBS(q), LIMBS(r),0, LIMBS(x), NLIMBS(x), LIMBS(y), NLIMBS(y));
    if (!r[lr - 1])
    {
      while(lr>2 && !r[lr - 1]) lr--;
      if (lr == 2) {avma=av; return gen_0;} /* exact division */
    }
    r[1] = evalsigne(sx) | evallgefint(lr);
    avma = (pari_sp) r; return r;
  }
  else
  {
    ulong lq = lgefint(x)-lgefint(y)+3;
    ulong lr = lgefint(y);
    GEN q = cgeti(lq);
    GEN r = cgeti(lr);
    mpn_tdiv_qr(LIMBS(q), LIMBS(r),0, LIMBS(x), NLIMBS(x), LIMBS(y), NLIMBS(y));
    if (q[lq - 1] == 0) lq--;
    q[1] = evalsigne(sy) | evallgefint(lq);
    if (!z) { avma = (pari_sp)q; return q; }
    if (!r[lr - 1])
    {
      while(lr>2 && !r[lr - 1]) lr--;
      if (lr == 2) {avma=(pari_sp) q; *z=gen_0; return q;} /* exact division */
    }
    r[1] = evalsigne(sx) | evallgefint(lr);
    avma = (pari_sp) r; *z = r; return q;
  }
}

/* Montgomery reduction.
 * N has k words, assume T >= 0 has less than 2k.
 * Return res := T / B^k mod N, where B = 2^BIL
 * such that 0 <= res < T/B^k + N  and  res has less than k words */
GEN
red_montgomery(GEN T, GEN N, ulong inv)
{
  (void)T; (void)N; (void)inv;
  pari_err(impl, "Montgomery reduction in gmp kernel");
  return NULL; /* not reached */
}

/* EXACT INTEGER DIVISION */

#if 1 /* use undocumented GMP interface */
static void
GEN2mpz(mpz_t X, GEN x)
{
  long l = lgefint(x)-2;
  X->_mp_alloc = l;
  X->_mp_size = signe(x) > 0? l: -l;
  X->_mp_d = LIMBS(x);
}
static void
mpz2GEN(GEN z, mpz_t Z)
{
  long l = Z->_mp_size;
  z[1] = evalsigne(l > 0? 1: -1) | evallgefint(labs(l)+2);
}

/* assume y != 0 and the division is exact */
GEN
diviuexact(GEN x, ulong y)
{
  if (!signe(x)) return gen_0;
  {
    long l = lgefint(x);
    mpz_t X, Z;
    GEN z = cgeti(l);
    GEN2mpz(X, x);
    Z->_mp_alloc = l-2;
    Z->_mp_size  = l-2;
    Z->_mp_d = LIMBS(z);
    mpz_divexact_ui(Z, X, y);
    mpz2GEN(z, Z); return z;
  }
}

/* Find z such that x=y*z, knowing that y | x (unchecked) */
GEN
diviiexact(GEN x, GEN y)
{
  if (!signe(y)) pari_err(gdiver);
  if (lgefint(y) == 3)
  {
    GEN z = diviuexact(x, y[2]);
    if (signe(y) < 0) togglesign(z);
    return z;
  }
  if (!signe(x)) return gen_0;
  {
    long l = lgefint(x);
    mpz_t X, Y, Z;
    GEN z = cgeti(l);
    GEN2mpz(X, x);
    GEN2mpz(Y, y);
    Z->_mp_alloc = l-2;
    Z->_mp_size  = l-2;
    Z->_mp_d = LIMBS(z);
    mpz_divexact(Z, X, Y);
    mpz2GEN(z, Z); return z;
  }
}
#else
/* assume y != 0 and the division is exact */
GEN
diviuexact(GEN x, ulong y)
{
  /*TODO: implement true exact division.*/
  return divii(x,utoi(y));
}

/* Find z such that x=y*z, knowing that y | x (unchecked)
 * Method: y0 z0 = x0 mod B = 2^BITS_IN_LONG ==> z0 = 1/y0 mod B.
 *    Set x := (x - z0 y) / B, updating only relevant words, and repeat */
GEN
diviiexact(GEN x, GEN y)
{
  /*TODO: use mpn_bdivmod instead*/
  return divii(x,y);
}
#endif


/********************************************************************/
/**                                                                **/
/**               INTEGER MULTIPLICATION                           **/
/**                                                                **/
/********************************************************************/

/* nx >= ny = num. of digits of x, y (not GEN, see mulii) */
GEN
muliispec(GEN x, GEN y, long nx, long ny)
{
  GEN zd;
  long lz;
  ulong hi;

  if (nx < ny) swapspec(x,y, nx,ny);
  if (!ny) return gen_0;
  if (ny == 1) return muluispec((ulong)*y, x, nx);

  lz = nx+ny+2;
  zd = cgeti(lz);
  hi = mpn_mul(LIMBS(zd), (mp_limb_t *)x, nx, (mp_limb_t *)y, ny);
  if (!hi) lz--;
  /*else zd[lz-1]=hi; GH tell me it is not necessary.*/

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

GEN
sqrispec(GEN x, long nx)
{
  GEN zd;
  long lz;

  if (!nx) return gen_0;
  if (nx==1) return sqru(*x);

  lz = (nx<<1)+2;
  zd = cgeti(lz);
  mpn_mul_n(LIMBS(zd), (mp_limb_t *)x, (mp_limb_t *)x, nx);
  if (zd[lz-1]==0) lz--;

  zd[1] = evalsigne(1) | evallgefint(lz);
  return zd;
}

/* x % (2^n), assuming x, n >= 0 */
GEN
remi2n(GEN x, long n)
{
  ulong hi;
  long l, k, lx, ly, sx = signe(x);
  GEN z, xd, zd;

  if (!sx || !n) return gen_0;

  k = dvmdsBIL(n, &l);
  lx = lgefint(x);
  if (lx < k+3) return icopy(x);

  xd = x + (2 + k);
  /* x = |k|...|1|#|... : copy the last l bits of # and the first k words
   *              ^--- initial xd  */
  hi = ((ulong)*xd) & ((1UL<<l)-1); /* last l bits of # = top bits of result */
  if (!hi)
  { /* strip leading zeroes from result */
    xd--; while (k && !*xd) { k--; xd--; }
    if (!k) return gen_0;
    ly = k+2;
  }
  else
    ly = k+3;

  zd = z = cgeti(ly);
  *++zd = evalsigne(sx) | evallgefint(ly);
  xd = x+1; for ( ;k; k--) *++zd = *++xd;
  if (hi) *++zd = hi;
  return z;
}

/********************************************************************/
/**                                                                **/
/**                      INTEGER SQUARE ROOT                       **/
/**                                                                **/
/********************************************************************/

/* Return S (and set R) s.t S^2 + R = N, 0 <= R <= 2S.
 * As for dvmdii, R is last on stack and guaranteed to be gen_0 in case the
 * remainder is 0. R = NULL is allowed. */
GEN
sqrtremi(GEN a, GEN *r)
{
  long l, na = NLIMBS(a);
  mp_size_t nr;
  GEN S;
  if (!na) {
    if (r) *r = gen_0;
    return gen_0;
  }
  l = (na + 5) >> 1; /* 2 + ceil(na/2) */
  S = cgetipos(l);
  if (r) {
    GEN R = cgeti(2 + na);
    nr = mpn_sqrtrem(LIMBS(S), LIMBS(R), LIMBS(a), na);
    if (nr) R[1] = evalsigne(1) | evallgefint(nr+2);
    else    { avma = (pari_sp)S; R = gen_0; }
    *r = R;
  }
  else
    (void)mpn_sqrtrem(LIMBS(S), NULL, LIMBS(a), na);
  return S;
}

/* compute sqrt(|a|), assuming a != 0 */
GEN
sqrtr_abs(GEN a)
{
  GEN res;
  mp_limb_t *b, *c;
  long l = RNLIMBS(a), e = expo(a), er = e>>1;
  long n;
  res = cgetr(2 + l);
  res[1] = evalsigne(1) | evalexpo(er);
  if (e&1)
  {
    b = (mp_limb_t *) new_chunk(l<<1);
    xmpn_zero(b,l);
    xmpn_mirrorcopy(b+l, RLIMBS(a), l);
    c = (mp_limb_t *) new_chunk(l);
    n = mpn_sqrtrem(c,b,b,l<<1); /* c <- sqrt; b <- rem */
    if (n>l || (n==l && mpn_cmp(b,c,l) > 0)) mpn_add_1(c,c,l,1);
  }
  else
  {
    ulong u;
    b = (mp_limb_t *) trunc2nr_lg(a,l+2,-1);
    b[1] = ((ulong)a[l+1])<<(BITS_IN_LONG-1);
    b = (mp_limb_t *) new_chunk(l);
    xmpn_zero(b,l+1); /* overwrites the former b[0] */
    c = (mp_limb_t *) new_chunk(l + 1);
    n = mpn_sqrtrem(c,b,b,(l<<1)+2); /* c <- sqrt; b <- rem */
    u = (ulong)*c++;
    if ( u&HIGHBIT || (u == ~HIGHBIT &&
             (n>l || (n==l && mpn_cmp(b,c,l)>0))))
      mpn_add_1(c,c,l,1);
  }
  xmpn_mirrorcopy(RLIMBS(res),c,l);
  avma = (pari_sp)res; return res;
}

/* Normalize a non-negative integer */
GEN
int_normalize(GEN x, long known_zero_words)
{
  long i =  lgefint(x) - 1 - known_zero_words;
  for ( ; i > 1; i--)
    if (x[i]) { setlgefint(x, i+1); return x; }
  x[1] = evalsigne(0) | evallgefint(2); return x;
}

/********************************************************************
 **                                                                **
 **                           Base Conversion                      **
 **                                                                **
 ********************************************************************/

ulong *
convi(GEN x, long *l)
{
  long n = nchar2nlong(2 + (long)(NLIMBS(x) * (BITS_IN_LONG * LOG10_2)));
  GEN str = cgetg(n+1, t_VECSMALL);
  unsigned char *res = (unsigned char*) GSTR(str);
  long llz = mpn_get_str(res, 10, LIMBS(icopy(x)), NLIMBS(x));
  long lz;
  ulong *z;
  long i, j;
  unsigned char *t;
  while (!*res) {res++; llz--;} /*Strip leading zeros*/
  lz  = (8+llz)/9;
  z = (ulong*)new_chunk(1+lz);
  t=res+llz+9;
  for(i=0;i<llz-8;i+=9)
  {
    ulong s;
    t-=18;
    s=*t++;
    for (j=1; j<9;j++)
      s=10*s+*t++;
    *z++=s;
  }
  if (i<llz)
  {
    unsigned char *t = res;
    ulong s=*t++;
    for (j=i+1; j<llz;j++)
      s=10*s+*t++;
    *z++=s;
  }
  *l = lz;
  return z;
}
