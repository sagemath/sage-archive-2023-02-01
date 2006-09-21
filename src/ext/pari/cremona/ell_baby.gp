/*  Babystep-giantstep function and applications

Author: John Cremona (john.cremona@nottingham.ac.uk)
Last edited: 13/10/04

[Based on JEC's C++ code which in turn is based on LiDIA code by
Volker Muller]


function bg_algorithm takes a curve e and points P,Q and finds i with
lower<=i<=upper such that i*P=Q (in the group), returning i or -1 if
no such i exists.

Applications:  finding order of a point given (a) a multiple of the
order;  (b) same with the factorization of the multiple;  (c) same
with lower and uppoer bounds on a multiple of the order.  [So far
these can be used for an elliptic curve over any field.]  Finally, (d)
for curves over Z/pZ only, the order of a point given no information,
using the Hasse bounds.  The latter could trivially be extended to
more general finite fields:  just change e.disc.mod to some function
returning the field size.

*/


global(MAX_BG_STEPS = 3000000);

{
hasse_bounds(q)=
local(u,l);
  u=sqrt(q << 2);
  l = ceil(q + 1 - u);    \\ lower bound of Hasse interval
  if(l<0,l=1);
  u = floor(q + 1 + u);   \\ upper bound of Hasse interval
[l,u];
}

\\ Order of a point given a factorizaton of a multiple of the order
{
order_point_from_fact(e,P,fact)=
local(ans,fm,p,m,s);
ans=1; fm=fact~;
for(j=1,#fm,p=fm[1,j];
m=factorback(matextract(fm,concat("^",Str(j)))~);
s=ellpow(e,P,m);
while(s!=[0],s=ellpow(e,s,p);ans*=p));
ans;
}

\\ Order of a point given a multiple of the order
{
order_point_from_mult(e,P, m)=order_point_from_fact(e,P,factor(m));
}

\\ Order of a point given lower & upper bounds of a multple, via bg
{
order_point_from_bounds(e,P, lower, upper)=
order_point_from_mult(e,P,bg_algorithm(e,P,[0],lower,upper,0));
}

\\ Order of a point given nothing except (finite!) field size
{
order_point(e,P)=
lu=hasse_bounds(e.disc.mod);
order_point_from_bounds(e,P,lu[1],lu[2]);
}

{
bg_algorithm(e, PP, QQ, lower, upper, info)=
local(P,Q,H,H2,H3,i,l,nb,ng,j,h,H_table,step_size,isort,H_table_sorted);
  if ((PP==[0]) && (QQ!=[0]), return(-1));

  if ((lower==0) && (QQ==[0]), return(0));

  if ((lower<0) || (upper<lower),
      print("bg_algorithm: lower bound > upper bound");
      return(-1));

  if (info, print("\nBabystep Giantstep algorithm: "));

  P=PP; Q=QQ;

  if ((upper - lower < 30),    \\ for very small intervals
      if (info, print("\nTesting ",(upper - lower)," possibilities ... "));
      H=ellpow(e,P,lower);
      h=lower;
      if (H==Q, return(h));
      while (h<upper, H=elladd(e,H,P); h+=1; if(H==Q,  return(h)));
      return(-1));

  \\**** otherwise we use the Babystep Giantstep method **************

  h = 1 + floor(sqrt(upper - lower));  \\ compute number of babysteps
  if (h > MAX_BG_STEPS, h = MAX_BG_STEPS);
  nb=h;

  \\ H_table[i] stores Str(lift(i*P)) for later sorting
\\  H_table=listcreate(nb);
  H_table=vector(nb);
  H2 = ellsub(e,Q,ellpow(e,P,lower));
  H = [0];

  \\****** Babysteps, store [i*P,i]  *********************

  if (info, print(" #Babysteps = ", nb));

  for (i = 1, nb,
      H=elladd(e,H,P);  \\ H = i*P and H2 = Q-lower*P
      if (H==H2, return (lower + i));
        \\ i * P = Q - lower* P, solution = lower+i
      H_table[i]=Str(lift(H));
\\      listput(H_table,[H,i]);  \\ store [H,i] in table
   );

  \\ Now for all i up to nb we have a table of  i*P, which we sort:
  isort=vecsort(H_table,,1);
  H_table_sorted=vecextract(H_table,isort);
  \\ and H  =  nb*P
  \\ and H2 =  Q-lower*P

  \\ We will subtract H from H2 repeatedly, so H2=Q-lower*P-j*H in the loop

  \\****** Giantsteps ***************************************************

  ng = 1+floor((upper - lower)/nb);

  if (info, print(", #Giantsteps = " , ng));

  step_size = nb;
  for (j = 0, ng,
      \\ Here H2=Q-(lower+j*step_size)*P
      if (H2==[0], \\ on the nail, no need to check table
	  h = lower + j * step_size;
	  if (h <= upper, return(h), return(-1)));

      \\ look in table to see if H2= i*P for a suitable i
      i=setsearch(H_table_sorted,Str(lift(H2)));
      if(i!=0,
	  i = isort[i];
\\	  H3=ellpow(e,P,i);
\\	  if (H3==H2,  \\ H2 is in table
	      h = lower + i + j * step_size;
	      if (h <= upper, return(h),  return(-1))
\\)
);
      H2=ellsub(e,H2,H)); \\ loop on j
 -1
}
