/* ======================================================== */

/* ===                                                  === */

/* ===                M o r i C o n e . c               === */

/* ===                                                  === */

/* ===	Authors: Maximilian Kreuzer, Nils-Ole Walliser	=== */

/* ===	Last update: 22/06/11                           === */

/* ===                                                  === */

/* ======================================================== */





/* ======================================================== */

/* =========            H E A D E R s             ========= */



#include "Global.h"

#include "Rat.h"

#include "Mori.h"





/* ======================================================== */

/* =========            D E F I N I T I O N s     ========= */



/*** local for Moricone.c ***/

#define Inci64		unsigned long long

#define		naT		FACE_Nmax	/* allocate: triangulation */



/* change to options */

#define NewtonMonomCOORD	(1)	/* (1) t_i, (2) u,v,w,x,... */

#define PRINT_MONOMIALS		(0)

#define EXCEPT_DIV_CLASS_BASE	(1)	/* try 0: D_min, 1: D_max */



/* diagnostic stuff */

#define TRACE_TRIANGULATION	(0)	/* detailed triangulation info */





/* ====================================================== */

/* =========           T Y P E D E F s          ========= */



/*** from Polynf.c ***/

typedef struct {

	int v, d;

	Long **x;

} Matrix;



/* ====================================================== */

/* =========            P R O T O T Y P E s     ========= */



/* from Polynf.c */

void Init_Matrix(Matrix *M,int v, int d);



int Make_G_for_GxMT_UT(Matrix M,Matrix G);



void Free_Matrix(Matrix *M);



Long VxV(Long *X,Long *Y,int d);



void Aux_IPS_Print_WP(Long *W,int w,int cd);

void Print_LMatrix(Matrix M, char *s);

 void Print_QuotZ(int Z[][VERT_Nmax],int *M,int p,int n);

/*

Auxiliary print functions

*/



/* from LG.c */

int  Init_Multiloop(int *N, int *I, int *j, int *J);

int  Multiloop(int *N,int *I,int *j,int *J);





// ===========================================================================

//

// Examples of interest for the extension of the triangulation algorithm

//

// 	make poly && poly.x -gPV zfree.c

//	poly.x -e ~/h/pic4| cws.x -N -f| poly.x -f| grep " 5 H">pic4.v5

//	make poly && tail -n11 ~/h/pic4.v5| poly.x -fgPV|grep pri



// 6 2 1 1 1 1 0 0  3 1 1 0 0 0 1 0  3 0 0 1 0 0 1 1 M:74 7 N:8 7 H:3,63 [-120]

// -> non-convex-triangulation :: Singular crash because of norm=0



// 7 1 2 2 0 0 1 1  9 1 3 3 1 1 0 0 15 2 5 4 1 0 0 3 M:48 10 N:9 7 H:4,37 [-66]

// -> incompatible induced triangulations  11010111 & 11101011 = 11000011



// 6 1 2 1 1 0 1  12 3 4 1 0 2 2 M:59 6 N:9 6 H:4,52 [-96]

// -> "triangles+squares"



// 6 2 1 1 0 1 1  3 1 0 0 1 0 1 /Z3: 1 1 2 2 0 0 M:41 7 N:9 6 H:4,34 [-60]

// -> non-triangular chamber in 3d 2ndary fan (double-intersecting edges)



// 8 1 0 4 0 1 1 0 1  6 1 0 3 1 0 1 0 0  6 1 0 3 0 1 0 1 0  6 0 1 3 1 0 0 1 0

// M:106 8 N:10 8 H:4,82 [-156] -> coincident (triple) intersections of edges

//============================================================================





//	I N C I D E N C E S

Inci64 makeN(int N){return ((Inci64) 1)<<N;}

void putN(int N,Inci64 *I){*I |= 1<<N;}			/* make INCIDENCE */

void setN(int N,Inci64 *I){*I |= 1<<N;}			/* make INCIDENCE */

int  getN(int N,Inci64 I){return (I>>N)%2;}		/* read INCIDENCE */

void prnI(int N,Inci64 I){int i; for(i=0;i<N;i++) printf("%d",getN(i,I));}

void fprI(int N,Inci64 I){int i;

  for(i=0;i<N;i++) fprintf(outFILE,"%d",getN(i,I));}    /* print INCIDENCE */



int Inci64_LE(Inci64 A,Inci64 B){return (A&B)==A;}

int Inci64_LT(Inci64 A,Inci64 B){if((A&B)==A) return A!=B; else return 0;}

//int Inci64_LmR(Inci64 *A,Inci64 *B){return (*A==*B) ? 0 : ((*A>*B) ? 1:-1);}

//int Inci64_diff(const void *A,const void *B){return Inci64_LmR(A,B);}



void PRNtriang(triang *SR,const char *c){int i; printf("%d %s\n",SR->n,c);

  for(i=0;i<SR->n;i++){if(i)printf(" ");prnI(SR->v,SR->I[i]);}puts("");}



int MaxBit(Inci64 I,int p){assert(I!=0);while(0==getN(--p,I));return p;}



/*   do {;} while(Multiloop(N,I,j,J));  <-->  do forall 0<=I[j]<N[j], 0<=j<J */

int  Init_Multiloop(int *N, int *I, int *j, int *J);

int  Multiloop(int *N,int *I,int *j,int *J);



/*  do {;} while(Choose(n,k,C));       <--> forall 0 <= C[0]<...<C[k-1] < n */

int Init_Choose(int n,int k,int *C){int i;

  long X=1; assert(k>0); assert(k<=n); for(i=0;i<k;i++){

    C[i]=i; X*=(n-i); X/=(i+1);} return X;}	      /* return (n choose k) */

int Choose(int n, int k, int *C){int x;

  for(x=0;x<k-1;x++) if(C[x]<C[x+1]-1) {int i;	// hole => ++ and below -> min

    for(i=0;i<x;i++) C[i]=i; C[x]++; return 1;} assert(x==k-1);

  if(C[x]<n-1){int i;for(i=0;i<x;i++)C[i]=i;C[x]++;return 1;} else return 0;}



/*  ==============    I N T E R S E C T I O N   R I N G	  ================   */



/*  find basis for divisor classes, i.e. integral basis of intersection ring */

/*

 *   Simply trying to find T-divisors that span the lattice of DivClasses

 *   i.e. find subdeterminants=volumes=1 for elimination.

 */

/*---QUICK-FIX----------------------------------------------------------------*/



void DivClassBasis(FILE *SF,PolyPointList *P,int v,char *D,char *B){

  Long cdiv=0, sv, nok, *X[POLY_Dmax]; int d=P->n, C[VERT_Nmax];

  nok=Init_Choose(v,d,C);

  do { int c, CC[VERT_Nmax];



    if(EXCEPT_DIV_CLASS_BASE) for(c=0;c<d;c++) CC[c]=C[c];

    else for(c=0;c<d;c++) CC[c]=v-1-C[c];	/* eliminate D_max first */

    for(c=0;c<d;c++) X[c]=P->x[CC[c]];

    sv=SimplexVolume(X,P->n); cdiv=NNgcd(cdiv,sv);

    if(sv==1){int x=0; fprintf(SF,"\nideal DCbase=");

      Inci64 I=makeN(CC[0]); for(c=1;c<d;c++) I+=makeN(CC[c]);

      for(c=0;c<v;c++) if(!getN(c,I)) {if(x) fprintf(SF,","); /* complement */

	fprintf(SF,"%s%d-%s%d",B,++x,D,c+1);}

      assert(x==v-d); fprintf(SF,";\n"); return; }

    nok--;} while(Choose(v,d,C)); assert(nok==0);



  if(cdiv>1){

	  printf("Fundamental group = Z%d\n",(int)cdiv);

    nok=Init_Choose(v,d,C);

    do { int c, CC[VERT_Nmax];



      if(EXCEPT_DIV_CLASS_BASE) for(c=0;c<d;c++) CC[c]=C[c];

      else for(c=0;c<d;c++) CC[c]=v-1-C[c];	/* eliminate D_max first */

      for(c=0;c<d;c++) X[c]=P->x[CC[c]]; sv=SimplexVolume(X,P->n);

      if(sv==cdiv){int x=0; fprintf(SF,"\nideal DCbase=");

      	Inci64 I=makeN(CC[0]); for(c=1;c<d;c++) I+=makeN(CC[c]);

        for(c=0;c<v;c++) if(!getN(c,I)) {if(x)fprintf(SF,",");

	  fprintf(SF,"%s%d-%s%d",B,++x,D,c+1);}

        assert(x==v-d); fprintf(SF,";\n"); return; }

      nok--;} while(Choose(v,d,C));}



  puts("IMPROVE CODE: no Vol=1 simplex in DivClassBasis()");exit(0);}





void OLD_LinRelLatticeBasis(FILE *SF,PolyPointList *P,int v,char *D,char *B){

  int i,j,k,l; Long cdiv=0; if(v>8) {DivClassBasis(SF,P,v,D,B);return;}

  if(P->n != 4) {DivClassBasis(SF,P,v,D,B);return;}





  for(l=3;l<v;l++)for(k=2;k<l;k++)for(j=1;j<k;j++)for(i=0;i<j;i++){Long *X[4];

    Inci64 I=makeN(i)+makeN(j)+makeN(k)+makeN(l); int x=0,y; Long sv;

    for(y=0;y<v;y++) if(getN(y,I)) X[x++]=P->x[y]; assert(x==4);

/*---ENDE QUICK-FIX----------------------------------------------------------*/



//printf("SimpVol[%d%d%d%d]=%ld\n",i,j,k,l,SimplexVolume(X,P->n));

    sv=SimplexVolume(X,P->n); cdiv=NNgcd(cdiv,sv);

    if(sv==1) {x=0; fprintf(SF,"\nideal DCbase=");

      for(y=0;y<v;y++) if(!getN(y,I)) {		      /* integral basis B[] */

        if(x) fprintf(SF,","); fprintf(SF,"%s%d-%s%d",B,++x,D,y+1);}

      assert(x==v-4); fprintf(SF,";\n"); return;}}



  if(cdiv>1){printf("Fundamental group = Z%d\n",(int)cdiv);fflush(0);

   for(l=3;l<v;l++)for(k=2;k<l;k++)for(j=1;j<k;j++)for(i=0;i<j;i++){Long *X[4];

    Inci64 I=makeN(i)+makeN(j)+makeN(k)+makeN(l); int x=0,y; Long sv;

    for(y=0;y<v;y++) if(getN(y,I)) X[x++]=P->x[y]; assert(x==4);

    sv=SimplexVolume(X,P->n);

    if(sv==cdiv) {x=0; fprintf(SF,"\nideal DCbase=");

      for(y=0;y<v;y++) if(!getN(y,I)) {		      /* integral basis B[] */

        if(x)fprintf(SF,",");fprintf(SF,"%s%d-%s%d",B,++x,D,y+1);} /* OFFSET */

      assert(x==v-4); fprintf(SF,";\n"); return;}}

   }

  puts("IMPROVE CODE: no Vol=1 simplex in LinRelLatticeBasis!!!");exit(0);

  }





/*  ideal IRingNorm = (Q*U - Di.Dj.Dk); on the CY where

 *  Q = Di.Dj.Dk.Dcy / Di.Dj.Dk.Dl*Vol(\s_ijkl) and Di.Dj.Dk.Dl*Vol(\s) = U = 1

 *  on ambient space IP_\S, i.e. we introduce a formal variable U = Unit

 */

/* new version:	number Vol=1; poly norm=Vol*reduce(d1*d2*d3*d4*d5,chow);

 *		reduce(nef*d1*d1*d1,chow)/norm;  // D=d divisor, B=h basis

 */

void CatFile(char *fn){char CAT[30+L_tmpnam];strcpy(CAT,"cat ");

  strcat(CAT,fn);printf("======= FILE content of %s:\n",fn); fflush(0);

  assert(0==system(CAT));

  printf("====== End of FILE content of %s\n\n",fn); fflush(0);}









/*   ==============	    M O R I   C O N E		================   */



void MoriGen(Matrix T,Long *V){int i; Matrix G; Init_Matrix(&G,T.d,T.d);

  assert(T.d==T.v+1); assert(T.v==Make_G_for_GxMT_UT(T,G));

  for(i=0;i<T.d;i++) V[i]=G.x[T.v][i]; Free_Matrix(&G);}



void Print_CMatrix(Matrix M, char *s){int i,j;

  fprintf(outFILE,"%d %d CV %s\n",M.d,M.v,s);for(i=0;i<M.d;i++){for(j=0;j<M.v;

    j++)fprintf(outFILE,"%2d%s",(int) M.x[j][i],(j+1==M.v) ? "\n" : " ");}}



void Print_VNL(VertexNumList *V){int i; for(i=0;i<V->nv;i++)

  fprintf(outFILE,"%d ",V->v[i]); fprintf(outFILE,"#=%d\n",V->nv);}



int  Vdiff_LmR(Long *L,Long *R,int d){

  while(d--){Long D=L[d]-R[d]; if(D) return (D>0) ? 1 : -1;} return 0;}



void Inci64_2_VNL(Inci64 X, VertexNumList *V, int n){int i; V->nv=0;

  for(i=0;i<n;i++) if(getN(i,X)) V->v[V->nv++]=i;

//Print_VNL(V);Sort_VL(V);Print_VNL(V);exit(0);

  }



void Print_Inci64_list(int n,Inci64 *I,int p){

  printf("To be done: Print_Inci64_list n=%d p=%d I=%lld\n",n,p,*I);exit(0);}



#define	Inci64_AND(I,J)		((I)&(J))



int Inci64_abs(Inci64 X){int abs=X%2; while(X/=2) abs+=X%2; return abs;}



void IDerr(){puts("\n       ********       INPUT DATA ERROR    	  ********");}



int Make_triCD2F(triang *T,Inci64 *cd2I){int i,j,cd2n=0;

  for(i=1;i<T->n;i++)for(j=0;j<i;j++){int k; Inci64 S=T->I[i]&T->I[j];

    for(k=0;k<cd2n;k++) if(Inci64_LE(S,cd2I[k])) break;

      else if(Inci64_LE(cd2I[k],S)) cd2I[k--]=cd2I[--cd2n];

  	if(k==cd2n) cd2I[cd2n++]=S;

    }	for(i=0;i<cd2n;i++) if(Inci64_abs(cd2I[i])!=T->d-1) break;

  if((i<cd2n)||(cd2n*2!=T->n*T->d)){IDerr();PRNtriang(T,"Triangulation ERROR");

    T->n=cd2n; T->I=cd2I; PRNtriang(T,"Codim-2 faces:");exit(0);}

  return cd2n;}



int Check_Mori(PolyPointList *P,int p,triang *T){	// strongly convex(?)

  int nI=T->n; Inci64 *I=T->I, cd2F[CD2F_Nmax];

  int i,j,r=0,d=P->n,ngen=0, ng0,nv,np; Long Z[POLY_Dmax+1]; VertexNumList V;

			// int e0=0,nm=0,m[VERT_Nmax];; Inci64 IE[VERT_Nmax];

  PolyPointList *UT = (PolyPointList *) malloc(sizeof(PolyPointList));

  Matrix R,VT,G; EqList *E = (EqList*)malloc(sizeof(EqList));



  assert(E!=NULL); assert(UT!=NULL);

  for(i=1;i<nI;i++)for(j=0;j<i;j++)ngen+=(Inci64_abs((I[i])&(I[j]))==d-1);

  Init_Matrix(&VT,d,d+1);Init_Matrix(&R,ngen,p);Init_Matrix(&G,p,p); ng0=ngen;

assert(ngen==Make_triCD2F(T,cd2F));

  for(i=1;i<nI;i++)for(j=0;j<i;j++)if(Inci64_abs(Inci64_AND(I[i],I[j]))==d-1){

    int a,b; Inci64 Ia=I[i]|I[j]; Inci64_2_VNL(Ia,&V,p); assert(V.nv==d+1);

    for(a=0;a<=d;a++) for(b=0;b<d;b++) VT.x[b][a]=P->x[V.v[a]][b];

    MoriGen(VT,Z); for(a=0;a<p;a++)R.x[r][a]=0;

    Ia=I[i]^I[j]; for(b=0;b<d;b++)if(getN(V.v[b],Ia)) break;

//prnI(p,I[i]);printf("=Ii Ij=");prnI(p,I[j]); printf(" Ia=");prnI(p,Ia);

//printf(" -> ");	Print_VNL(&V);

assert(Inci64_abs(Ia)==2); assert(b<d);

    if(0==Z[b]){ printf("Error: Z[b]==0 for I[%d]&I[%d] !\n",j,i);exit(0);}

    if(Z[b]>0) for(a=0;a<=d;a++) R.x[r][V.v[a]]=Z[a];

    else for(a=0;a<=d;a++) R.x[r][V.v[a]]=-Z[a];

//for(a=0;a<p;a++)printf("%ld ",R.x[r][a]);printf("=Z::R[%d]\n",r);fflush(0);

    for(b=0;b<r;b++) if(!Vdiff_LmR(R.x[r],R.x[b],p)) break; if(b==r) r++; }

  R.v=UT->np=ngen=r; r=Make_G_for_GxMT_UT(R,G); 	// base change -> UT

  if(r!=p-d) { Matrix GR; R.v=ngen; Print_LMatrix(R,"R");

    Print_LMatrix(G,"GLZ"); Init_Matrix(&GR,R.v,p);

    for(i=0;i<R.v;i++) for(j=0;j<p;j++) GR.x[i][j]=VxV(G.x[j],R.x[i],p);

    Print_LMatrix(GR,"GR"); printf("rank=%d != p-d !!!\n",r); exit(0);}

//Print_LMatrix(R,"Matrix of all rays");

  if(UT->np>=POINT_Nmax){fprintf(outFILE,"need POINT_Nmax>=%d\n",UT->np+1);

    exit(0);} UT->n=r;

  if(r>POLY_Dmax){fprintf(outFILE,"need POLY_Dmax>=%d\n",UT->n);exit(0);}

  for(i=0;i<ngen;i++) for(j=0;j<r;j++) UT->x[i][j]=VxV(G.x[j],R.x[i],p);

  for(i=0;i<ngen;i++) for(j=r;j<p;j++) assert(0==VxV(G.x[j],R.x[i],p));

  for(i=0;i<UT->n;i++)UT->x[UT->np][i]=0; UT->np++; Find_Equations(UT,&V,E);

  np=UT->np-1; nv=V.nv-1; Sort_VL(&V);



  Free_Matrix(&VT); Free_Matrix(&R); Free_Matrix(&G); free(E); free(UT);



  if(np!=V.v[nv]) {PRNtriang(T,"Non-coherent Triangulation"); return 0;}

  else return 1;

  }







void Print_Mori(PolyPointList *P,int p,int nI, Inci64 *I){

  int i,j,r=0,d=P->n,ngen=0, ng0,e0=0,nm=0,nv,np;

  int m[VERT_Nmax]; Long Z[POLY_Dmax+1]; VertexNumList V; Inci64 IE[VERT_Nmax];

  PolyPointList *UT = (PolyPointList *) malloc(sizeof(PolyPointList));

  Matrix R,VT,G; EqList *E = (EqList*)malloc(sizeof(EqList));



  assert(E!=NULL); assert(UT!=NULL);

  for(i=1;i<nI;i++)for(j=0;j<i;j++)ngen+=(Inci64_abs((I[i])&(I[j]))==d-1);

  Init_Matrix(&VT,d,d+1);Init_Matrix(&R,ngen,p);Init_Matrix(&G,p,p); ng0=ngen;

  for(i=1;i<nI;i++)for(j=0;j<i;j++)if(Inci64_abs(Inci64_AND(I[i],I[j]))==d-1){

    int a,b; Inci64 Ia=I[i]|I[j]; Inci64_2_VNL(Ia,&V,p); assert(V.nv==d+1);

    for(a=0;a<=d;a++) for(b=0;b<d;b++) VT.x[b][a]=P->x[V.v[a]][b];

    MoriGen(VT,Z); for(a=0;a<p;a++)R.x[r][a]=0;

    Ia=I[i]^I[j]; for(b=0;b<d;b++)if(getN(V.v[b],Ia)) break;

//prnI(p,I[i]);printf("=Ii Ij=");prnI(p,I[j]); printf(" Ia=");prnI(p,Ia);

//printf(" -> ");	Print_VNL(&V);

assert(Inci64_abs(Ia)==2); assert(b<d);

    if(0==Z[b]){ printf("Error: Z[b]==0 for I[%d]&I[%d] !\n",j,i);exit(0);}

    if(Z[b]>0) for(a=0;a<=d;a++) R.x[r][V.v[a]]=Z[a];

    else for(a=0;a<=d;a++) R.x[r][V.v[a]]=-Z[a];

//for(a=0;a<p;a++)printf("%ld ",R.x[r][a]);printf("=Z::R[%d]\n",r);fflush(0);

    for(b=0;b<r;b++) if(!Vdiff_LmR(R.x[r],R.x[b],p)) break; if(b==r) r++; }

  R.v=UT->np=ngen=r; r=Make_G_for_GxMT_UT(R,G); 	// base change -> UT

  if(r!=p-d) { Matrix GR; R.v=ngen;Print_LMatrix(R,"R");

    Print_LMatrix(G,"GLZ"); Init_Matrix(&GR,R.v,p);

    for(i=0;i<R.v;i++) for(j=0;j<p;j++) GR.x[i][j]=VxV(G.x[j],R.x[i],p);

    Print_LMatrix(GR,"GR"); printf("rank=%d != p-d !!!\n",r); exit(0);}

//Print_LMatrix(R,"Matrix of all rays");

  if(UT->np>=POINT_Nmax){fprintf(outFILE,"need POINT_Nmax>=%d\n",UT->np+1);

    exit(0);} UT->n=r;

  if(r>POLY_Dmax){fprintf(outFILE,"need POLY_Dmax>=%d\n",UT->n);exit(0);}

  for(i=0;i<ngen;i++) for(j=0;j<r;j++) UT->x[i][j]=VxV(G.x[j],R.x[i],p);

  for(i=0;i<ngen;i++) for(j=r;j<p;j++) assert(0==VxV(G.x[j],R.x[i],p));

  for(i=0;i<UT->n;i++)UT->x[UT->np][i]=0; UT->np++; Find_Equations(UT,&V,E);

  np=UT->np-1; nv=V.nv-1; Sort_VL(&V);

  if(np!=V.v[nv]){IDerr();puts("MORI CONE not strictly convex:");

    Print_Inci64_list(nI,I,p);puts("... non-convex triangulation?\n");exit(0);}



/* The extremal rays of the Mori cone are those that have maximal incidences *

 * with faces of the cone, i.e. with equations containing the origin.	     */

  for(i=0;i<nv;i++) IE[i]=0;  	   /* igonore the origin = UT.x[np==V.v[nv]] */

  for(i=0;i<E->ne;i++) if(Eval_Eq_on_V(&E->e[i],UT->x[np],UT->n)==0){e0++;

    for(j=0;j<nv;j++) IE[j]=(2*IE[j]+

      !Eval_Eq_on_V(&E->e[i],UT->x[V.v[j]],r));}/* compute Eq(0)-INCIs for Vs */

  if(e0>VERT_Nmax){fprintf(outFILE,"need VERT_Nmax >= %d\n",e0);exit(0);}

//printf("p=%d nm=%d\n",p,nm);

  for(i=0;i<nv;i++){

	  for(j=0;j<nv;j++)

//    if(i!=j) if(INCI_LE(IE[i],IE[j]))			/* equivalent ??? */

	  if((((IE[i]) | (IE[j])) ==(IE[j]))&&!((IE[i])==(IE[j])))

		  break;

  if(j==nv) m[nm++]=i;

  }

  assert(nm>=r);

//  fprintf(outFILE,

//  "%d MORI GENERATORS / dim(cone)=%d   [#rays=%d<=%d #eq=%d<=%d #v=%d<=%d]\n",

//    nm,r,ngen,ng0,e0,E->ne,nm,V.nv);

  fprintf(outFILE,

  "%d MORI GENERATORS / dim(cone)=%d \n",

    nm,r);

//printf("p=%d nm=%d\n",p,nm);fflush(0);

  for(i=0;i<nm;i++){ int n=V.v[m[i]]; Long s=0;

    for(j=0;j<p;j++)s+=R.x[n][j];

//fprintf(outFILE,"%3ld ",-s);	/* sum of degrees == CY-divisor/linebundle */

    for(j=0;j<p;j++)fprintf(outFILE," %2ld",R.x[n][j]);

    //fputs("   I:",outFILE);

    //prnI(e0,IE[m[i]]);

    fputs("\n",outFILE);

    fflush(0);}

  Free_Matrix(&VT); Free_Matrix(&R); Free_Matrix(&G); free(E); free(UT); }









// ===============	Triangulation  <-->  Stanley Reisner	============ //



/*   SR generators:  G=I+2^n  with  I<2^n,  I<=face  and not  G'<G           */

/*   d = dim(Poly) = #(vertices of simplices on *T) = #vertices on simp.facets

 *   OFFSET=0, i.e. IP not on (position 0 of) Inci64

 */



/* Triangles are of the form $T=I+2^n  with  I<2^n,  G not <=I and |I|=dim

 */

void Triang_from_SR(triang *TR,triang *SR){	/* consistency check ... */

  int i=1, p=TR->v=SR->v, j=p/2, d=TR->d=SR->d, s=SR->n, m=0, k, l, r;

  Inci64 *S=SR->I,*T=TR->I, *A,*M,*N; long long binco=TR->v; /* Bino.Coeff */

  assert(p<=64); if(j>d)j=d+1; while(i<j){binco*=(p-i);binco/=++i;}

  if(binco++>2999)binco=2999; A = (Inci64*) malloc(2*binco*sizeof(Inci64));

  assert(A!=NULL); M=A; N=&A[binco];

  for(i=1;i<p;i++)for(j=0;j<i;j++){M[m]= makeN(i)+makeN(j);

    for(k=0;k<s;k++) if(Inci64_LE(S[k],M[m])) break;  /* no edge of triangle */

    if(k==s) assert(++m<binco);}

  for(r=3;r<=d;r++){int x, n=0; for(k=0;k<m;k++)

    for(x=1+MaxBit(M[k],p); x<p; x++) {N[n]=M[k]+makeN(x);

      for(l=0;l<s;l++) if(Inci64_LE(S[l],N[n])) break;

      if(l==s) assert(++n<binco);}

    Inci64 *swapI=M; M=N; N=swapI; m=n;}

  TR->n=m; assert(m<=TR->nmax); for(k=0;k<m;k++) T[k]=M[k]; free(A);}



void StanleyReisner(triang *SR,triang *T){ /* pre-allocate and compute SR(T) */

  Inci64 *S=SR->I,*I=T->I,*A,*M,*N,U=1; long long binco=T->v; /* Binom.Coeff */

  int i=1, p=T->v, j=p/2, d=T->d,nI=T->n, s=0, m=0, k, l, r; 	 assert(p<=64);

  if(j>d)j=d+1; while(i<j){binco*=(p-i);binco/=++i;}if(binco++>2999)binco=2999;



  A = (Inci64*) malloc(2*binco*sizeof(Inci64));assert(A!=NULL);M=A;N=&A[binco];

  for(i=1;i<p;i++)for(j=0;j<i;j++){M[m]=(U<<i)+(U<<j);		/* OFFSET=0 */

    for(k=0;k<nI;k++) if(Inci64_LE(M[m],I[k])) break;/* no primitive collect */

    if(k==nI) S[s++]=M[m]; else assert(++m<binco);}   /* quadratic generator */

  for(r=3;r<=d+1;r++){int x, n=0; for(k=0;k<m;k++)

    for(x=1+MaxBit(M[k],p); x<p; x++) {N[n]=M[k]+(U<<x);

  /*prnI(p,M[k]);printf("=M[%d] x=%d > N[%d]=",k,x,n);prnI(p,N[n]);puts("");*/

      for(l=0;l<nI;l++) if(Inci64_LE(N[n],I[l])) break;

      if(l==nI) {for(l=0;l<s;l++) if(Inci64_LE(S[l],N[n])) break;

        if(l==s) {assert(s<SR->nmax);S[s++]=N[n];}}

      else assert(++n<binco);}

    Inci64 *swapI=M; M=N; N=swapI; m=n;} SR->n=s; SR->v=T->v; SR->d=T->d;

  { int ok=1; triang TeST; TeST.I=A; TeST.nmax=2*binco;   //

    Triang_from_SR(&TeST,SR); if(TeST.n!=T->n) ok=0; else

#ifdef	TRIANG_CHECKSUM

    {Inci64 sumT=0,sumC=0;for(i=0;i<T->n;i++){sumT+=TeST.I[i];sumS+>T->I[i];}

      if(sumT!=sumC) ok=0;}

#else

    {for(i=0;i<T->n;i++) {int sum=0; for(j=0;j<T->n;j++)

      if(TeST.I[i]==T->I[j]) sum++; if(sum!=1) break;} if(i<T->n) ok=0;}

#endif

    if(ok)free(A); else {PRNtriang(T,"Triangulation");PRNtriang(SR,"SR-ideal");

      PRNtriang(&TeST,"Tri(SR) ... test failed !!!"); assert(0);}

  }}









//	INTERSECTION RING (Singular)  /  MORI CONE  /  TRIANGULATIONS



void InterSectionRing(Inci64 *Tri,int *t,PolyPointList *P, int p,

		MORI_Flags *_Flag, FibW *F){

  triang T,SR;

  Inci64 srI[VERT_Nmax];

  T.v=p; T.n=*t; T.I=Tri;

  SR.d=T.d=P->n; SR.I=srI; SR.d=P->n; SR.n=0; SR.v=p; SR.nmax=T.nmax=VERT_Nmax;

  if(Check_Mori(P,p,&T)){

	  if(_Flag->g)

		  PRNtriang(&T,"Triangulation");

	  StanleyReisner(&SR,&T);

	  if(_Flag->g)

		  PRNtriang(&SR,"SR-ideal");

	  if(_Flag->i || _Flag->t || _Flag->c || _Flag->d || _Flag->a || _Flag->b || _Flag->H){

		  if(P->n==4){

			  HyperSurfSingular(P,&T,&SR,_Flag,F,&p);

		  }

		  else{

			  printf("Intersection ring implemented only for dim=4 polytope \n");

		  }

	  }

	  if(_Flag->m)

		  Print_Mori(P,p,*t,Tri);



  }

}









/*   ==============		TRIANGULATIONS		================   */



void Transpose(Matrix M,Matrix MT){int i,j,l=M.v,c=M.d;

  assert((MT.v==c)&&(MT.d==l));

  for(i=0;i<M.v;i++)for(j=0;j<M.d;j++)MT.x[j][i]=M.x[i][j];}

int  Make_G_for_MxG_LT(Matrix M,Matrix *G){  /* MxG lower trian, return rank */

  Matrix MT; Init_Matrix(&MT,M.d,M.v); Transpose(M,MT);

  Init_Matrix(G,M.v,M.v); return Make_G_for_GxMT_UT(MT,*G);}



void GaleTransform(Matrix A,Matrix *B){

  Matrix G; int i,j; Make_G_for_MxG_LT(A,&G);     // Print(A); Print(G);

  Init_Matrix(B,A.v,A.v-A.d);

  for(i=0;i<A.v;i++) for(j=0;j<A.v-A.d;j++) B->x[i][j]=G.x[A.d+j][i];}





/*   GKZ: Gale transform; implement circuit and 2ndary polygon

 *   T[] space for triangulations; t[] #simplices; nt=number of triangulations

 */



int AcuteAngle(Long *L,Long *R){return L[0]*R[0]+L[1]*R[1]>0;} // R.L>0

int ConeAngle(Long *L,Long *R){Long X=L[1]*R[0]-L[0]*R[1]; // 0->0, (0,pi) ->+1

  if(X>0)return 1;if(X<0)return -1;return(AcuteAngle(L,R)>0)?0:-1;}// [pi,)->-1



#define BZangle(a,b)	(ConeAngle(B.x[Z[a]],B.x[Z[b]]))	// Gale-points

/*

 *   2d: secondary POLYGON: return nmt = #(maximal triangulations)

 *   V[]=facet points; rays R [ r<nr ] [ i<nrp[r] ];

 */





#define		ANfan		20	// alloc number of max 2nd-fans

#define		ANtri		20	// alloc number of max triang.



Inci64 FindPolyCircuits(PolyPointList *P,int p,Inci64 F,int f){

  int i,j,k=0,d=P->n,C[VERT_Nmax]; Inci64 X=0; Matrix A,B; Init_Matrix(&A,f,d);

  for(i=0;i<p;i++)if(getN(i,F)){for(j=0;j<d;j++)A.x[k][j]=P->x[i][j];C[k++]=i;}

  assert(f>d);assert(k==f);d=f-d;GaleTransform(A,&B);//Print_CMatrix(B,"Gale");

  for(i=0;i<k;i++) {Long *Y=B.x[i]; for(j=0;j<d;j++)if(Y[j]) break;

    if(j<d) X+=makeN(C[i]);} Free_Matrix(&A); Free_Matrix(&B); return X;}





/*	aux functions for TRI-CIRCUIT 3d-GKZ:	Triang3dSFan()

 */

#define	SameRayBZ(a,b)		(SameRay(B.x[Z[a]],B.x[Z[b]],B.d))



int  SameRay(Long *X,Long *Y,int d){int x=0,y=000; while(d--)

  if(x) {if(y*X[d]!=x*Y[d]) return 0;} else

  if(X[d]) {if(Y[d]){Long g=Fgcd(X[d],Y[d]); x=X[d]/g; y=Y[d]/g;}

    else return 0;} else {if(Y[d]) return 0;}

  if(x) return (x*y>0); else {puts("ZeroVectors in SameRay (forbidden)");

  assert(x); return 000;}}



Long XYZproduct(Long *X,Long *Y,Long *Z){return Z[2]*(X[0]*Y[1]-X[1]*Y[0])

  +Z[0]*(X[1]*Y[2]-X[2]*Y[1])+Z[1]*(X[2]*Y[0]-X[0]*Y[2]);}

void CROSSproduct(Long *X,Long *Y,Long *XxY){XxY[2]=X[0]*Y[1]-X[1]*Y[0];

  XxY[0]=X[1]*Y[2]-X[2]*Y[1]; XxY[1]=X[2]*Y[0]-X[0]*Y[2];}

Long SCALproduct(Long *X,Long *Y){return X[0]*Y[0]+X[1]*Y[1]+X[2]*Y[2];}

//#define BZx(a,b,c)	(XYZproduct(B.x[Z[a]],B.x[Z[b]],B.x[Z[c]]))

//#define BZRx(a,b,c)	(XYZproduct(B.x[Z[*R[a]]],B.x[Z[*R[b]]],B.x[Z[*R[c]]]))



#define		BZR(eqr)	(B.x[Z[*R[eqr]]])

#define		BZRx(a,b,c)	(XYZproduct(BZR(a),BZR(b),BZR(c)))

#define		BZRE(i,j)	(BZR(Eli[i][j]))



void AuxPrintRays(int R[VERT_Nmax][POLY_Dmax],int nrp[VERT_Nmax],int nr){

  int i,j; printf("Rays: "); for(i=0;i<nr;i++) { printf("R[%d]=",i);

    for(j=0;j<nrp[i];j++)printf("%d ",R[i][j]); } puts("");}



int InterSectE(Long *X,Long *Y,Long *U,Long *V){

  Long xyu=XYZproduct(X,Y,U), xyv=XYZproduct(X,Y,V), uvx=XYZproduct(U,V,X),

    uvy=XYZproduct(U,V,Y); if((0<=xyu*xyv)||(0<=uvx*uvy)) return 0;

    if(xyu*uvx<0) return 1; else return -1;}



void MakeVecPrim(Long *X){Long g=1;if(X[0]*X[1])g=Fgcd(X[0],X[1]);

  if(X[2])g=Fgcd(g,X[2]); if(g<0)g=-g; if(g>1){X[0]/=g;X[1]/=g;X[2]/=g;}}



void IntersectEdges(Long *X,Long *Y,Long *U,Long *V,Long *Q){

  Long XY[3],UV[3]; CROSSproduct(X,Y,XY); MakeVecPrim(XY);

  if(SCALproduct(XY,U)<0)CROSSproduct(U,V,UV);else CROSSproduct(V,U,UV);

  MakeVecPrim(UV); CROSSproduct(XY,UV,Q); MakeVecPrim(Q);



#if (TRACE_TRIANGULATION)

printf("X=%ld %ld %ld, Y=%ld %ld %ld : ",X[0],X[1],X[2],Y[0],Y[1],Y[2]);

printf("U=%ld %ld %ld, V=%ld %ld %ld -> Q=%ld %ld %ld ...\n",

U[0],U[1],U[2],V[0],V[1],V[2],Q[0],Q[1],Q[2]);

#endif



assert((SCALproduct(X,Q)>0)||(SCALproduct(Y,Q)>0));

assert((SCALproduct(U,Q)>0)||(SCALproduct(V,Q)>0));

  assert(XYZproduct(X,Y,Q)==0);assert(XYZproduct(U,V,Q)==0);}



void Print_MaxTrian(Inci64 C, Inci64 *CT[ANtri], int nmt, int *nt,int p)

{ int i,j;for(i=0;i<nmt;i++){Inci64 U=0;int N=nt[i];printf("MaxTr(");prnI(p,C);

    for(j=0;j<N;j++){printf("%s",j?",":")={");prnI(p,CT[i][j]);U|=CT[i][j];}

    puts("}");assert(U==C);}

}









int Triang1dSFan(PolyPointList *P,int p,Inci64 I,Inci64 *X,Inci64 *CT[ANtri],

  int *nmt, int *nt)			// P,p,F[c],X,CT[nmf],&nmt[nmf],nt[nmf]

{ int d=P->n,i,j,k=0,tnt=0, F[POLY_Dmax+1]; Inci64 CI=0; Matrix A,B;

  Init_Matrix(&A,d+1,d); for(i=0;i<p;i++) if(getN(i,I)){

    for(j=0;j<d;j++) A.x[k][j]=P->x[i][j]; F[k++]=i;} assert(k==d+1);

  GaleTransform(A,&B); i=j=*nmt=0;

  for(k=0;k<=d;k++) if(*B.x[k]) {CI|=makeN(F[k]); if(*B.x[k]>0) i++; else j++;}

#if (TRACE_TRIANGULATION)

	{

	  Inci64 PC=FindPolyCircuits(P,p,I,d+1);

	  assert(PC==CI);

	}

	prnI(p,CI);

	printf("=C -> ");

	Print_CMatrix(B,"Gale");

#endif

  if(i>1){CT[*nmt]=X; nt[*nmt]=0; for(k=0;k<=d;k++)

    if(*B.x[k]>0) CT[*nmt][nt[*nmt]++]=CI-makeN(F[k]); tnt+=i; (*nmt)++;}

  if(j>1){CT[*nmt]=X+i; nt[*nmt]=0; for(k=0;k<=d;k++)

    if(*B.x[k]<0) CT[*nmt][nt[*nmt]++]=CI-makeN(F[k]); tnt+=j; (*nmt)++;}

#if (TRACE_TRIANGULATION)

    Print_MaxTrian(CI,CT,*nmt,nt,p); //printf("i=%d j=%d nmt=%d\n",i,j,*nmt);

#endif

  assert(*nmt==(i*j+1>i+j)+1); Free_Matrix(&A); Free_Matrix(&B); return tnt;}









int Triang2dSFan(PolyPointList *P,int p,Inci64 FI,Inci64 *X,Inci64 *CT[ANtri],

  int *nmt,int *nt)			// P,p,F[c],X,CT[nmf],&nmt[nmf],nt[nmf]

{ int d=P->n,v=d+2,i,j,k=0,r,F[POLY_Dmax+2],Z[POLY_Dmax+2],z=0, // GKZ polygon:

    R[VERT_Nmax][POLY_Dmax],nrp[VERT_Nmax],nr,neg=0,tnt=0;// num.max.tri

  Matrix A,B; Init_Matrix(&A,v,d); Inci64 U;		 // P.x::F[A]::F[Z[R]]

  for(i=0;i<p;i++) if(getN(i,FI)){for(j=0;j<d;j++) A.x[k][j]=P->x[i][j];

    F[k++]=i;} assert(k==d+2); GaleTransform(A,&B); U=0;	// F_i<d+2

  for(k=0;k<v;k++) if(B.x[k][0]||B.x[k][1]) U += makeN(F[ Z[z++]=k ]);

#if (TRACE_TRIANGULATION)

	 Print_CMatrix(B,"Gale");//Print_LMatrix(A,"bi-circuit"); 	// Z_k<z

#endif

  nr=1; nrp[0]=1; R[0][0]=0;		// Z[z]::T[0] {circuits} in F[]::=facet

  for(i=1;i<z;i++) if(0<=(k=BZangle(i,0))){	// points with positive angle

    if(k==0) r=0; // else if(nr==1) r=1; // insert@ 0 < r <=nr; add@ 0 <= r <nr

    else for(r=1;r<nr;r++) if((k=BZangle(i,*R[r]))<=0) break;

    {

    	int J,L;

    	if(k){for(J=nr++;J>r;J--)

    	{nrp[J]=nrp[J-1];for(L=0;L<nrp[J];L++)R[J][L]=R[J-1][L];}

    	R[r][0]=i;nrp[r]=1;} else R[r][nrp[r]++]=i;

    }// if(k) insert at r; else add to R[r];

  }



  for(i=1;i<z;i++) if(0>(k=BZangle(i,0))){ 	// points with negative angle

    if(neg==0) {r=neg=nr; k=1;}  	// init neg.

    else for(r=neg;r<nr;r++) if((k=BZangle(i,*R[r]))<=0) break;

//printf("i=%d: neg=%d nr=%d  ->  k=%d r=%d\n",i,neg,nr,k,r);

    {

        int J,L;

        if(k){for(J=nr++;J>r;J--)

        {nrp[J]=nrp[J-1];for(L=0;L<nrp[J];L++)R[J][L]=R[J-1][L];}

        R[r][0]=i;nrp[r]=1;} else R[r][nrp[r]++]=i;

        }// if(k) insert at r; else add to R[r];

  }				 		CT[*nmt=0]=X;

#if (TRACE_TRIANGULATION)

	printf("U="); prnI(p,U); printf(" -> rays: ");

	for(r=0;r<nr;r++){printf("R[%d]=",r);

	for(j=0;j<nrp[r];j++)printf("%d ",R[r][j]);}puts("");	// check rays:

	i=0;k=0;for(r=0;r<nr;r++){assert(BZangle(*R[(r+1)%nr],*R[r])>0);

	for(j=0;j<nrp[r];j++){assert(0<=R[r][j]);assert(R[r][j]<z);k+=R[r][j];

	if(j)assert(BZangle(R[r][j-1],R[r][j])==0);}i+=nrp[r];}

	assert(i==z);assert(2*k==z*(z-1));

#endif

  for(r=0;r<nr;r++){int a,s=(r+1)%nr;	// triangulation for cone (R[r],R[s])

    Inci64 PC=0,*CI=CT[*nmt]; nt[*nmt]=0; assert(BZangle(*R[s],*R[r])>0);

//    if((nrp[r]>1)&&(nrp[s]>1))		// maximal triangulations only

    for(j=1;j<nr;j++){int b,y;   if(BZangle(*R[y=(r+j)%nr],*R[r])<=0) break;

      for(i=1;i<nr;i++){int x; if(BZangle(*R[s],*R[x=(nr-i+s)%nr])<=0)break;

	if(BZangle(*R[y],*R[x])>0) for(a=0;a<nrp[x];a++) for(b=0;b<nrp[y];b++)

	PC|=(CI[nt[*nmt]++]=U-makeN(F[Z[R[x][a]]])-makeN(F[Z[R[y][b]]]));}}

    if(PC==U) {tnt+=nt[*nmt]; CT[++(*nmt)]=&X[tnt];} }	assert(*nmt>0);

#if (TRACE_TRIANGULATION)

	Print_MaxTrian(U,CT,*nmt,nt,p);

#endif

  return tnt;}		// nmt = # maximal triangulations



#ifdef	OLD_code	// problem: negative cone -> XYZcone

int ABCline(Long *a,Long *b,Long *c){int i;	// 1:: edge(abc), -1:: acb|cab

  Long A[3],C[3],AA=0,AC=0,CC=0; 		// 0:: no line or a=b or a=c

  for(i=0;i<3;i++){AA+=b[i]*a[i];AC+=b[i]*b[i];CC+=b[i]*c[i];}	// A=ab^2-b(ab)

  for(i=0;i<3;i++){A[i]=a[i]*AC-b[i]*AA; C[i]=c[i]*AC-b[i]*CC;}	// C=cb^2-b(cb)

//	printf("ABCline a=(%ld,%ld,%ld),b=(%ld,%ld,%ld),c=(%ld,%ld,%ld)  ",

//	a[0],a[1],a[2],b[0],b[1],b[2],c[0],c[1],c[2]);printf(

//	"A=(%ld,%ld,%ld),C=(%ld,%ld,%ld)\n",A[0],A[1],A[2],C[0],C[1],C[2]);

  if(0==(AC=SCALproduct(A,C))) return 0; 	// A/b^2, C/b^2 = ortho. proj.

  if(!(AA=SCALproduct(A,A))) return 0; if(!(CC=SCALproduct(C,C))) return 0;

  if(AA*CC!=AC*AC) return 0; return (AC<0) ? 1 : -1;}	// assuming a!=c

#endif



int XYZcone(Long *A,Long *B,Long *C){	//  1: B inside cone <AC>_+

  int i,j,k; Long x,y,z;		// -1: <ABC>_+ is 2d strict conv. cone

  if(XYZproduct(A,B,C)) return 0;	//  0: <XYZ> 3d or non strictly convex

  for(k=0;k<3;k++){i=(k+1)%3; j=(k+2)%3;//

    if((x=B[i]*C[j]-B[j]*C[i])) break;} //  solve    x A + y B + z C = 0

  if(k==3) return 0;                    //  coincident lines => return 0

  if(!(y=C[i]*A[j]-C[j]*A[i])) return 0;

  if(!(z=A[i]*B[j]-A[j]*B[i])) return 0;

  if(x*z>0) return (x*y<0) ? 1 : 0; else return -1; }



int Triang3dSFan(PolyPointList *P,int p,Inci64 FI,Inci64 *X,Inci64 *CT[ANtri],

  int *nmt,int *nt){Matrix A,B; Inci64 C=0; int tmt=0; // #triang in max.tri's

  int R[VERT_Nmax][POLY_Dmax], nrp[VERT_Nmax], r=0;  // Rays, #ray-Pts, #rays

  int d=P->n,i,j,f=0,z=0,F[VERT_Nmax],Z[VERT_Nmax];  // P.x [ F [ Z [ R < z ]]]

  Inci64 Einc[VERT_Nmax]; Long CR[VERT_Nmax][3];	// chamber rays <

  int Tli[VERT_Nmax][3],Eli[VERT_Nmax][2],IEli[VERT_Nmax][9],ien[VERT_Nmax];

  int tn=0,en=0, a,b,c, k,l,iem=0,ies=0,nQuad=0,chi=0,ncr=0;



/*  int ChamberTriangle(Long *Q,int *T){Long a,b,c;   // Q in BZR(T_0,T_1,T_2)

    if(0<=(a=XYZproduct(BZR(T[0]),BZR(T[1]),Q))) if(0<=(b=XYZproduct(BZR(T[1])

      ,BZR(T[2]),Q))) if(0<=(c=XYZproduct(BZR(T[2]),BZR(T[0]),Q)))

    {assert(a*b*c>0); return 1;} return 0;}		 // END of DECLARATIONS */



  Init_Matrix(&A,Inci64_abs(FI),d);for(i=0;i<=p;i++)if(getN(i,FI)){	// GALE

    for(j=0;j<d;j++) A.x[f][j]=P->x[i][j]; F[f++]=i;} 	GaleTransform(A,&B);

  for(i=0;i<f;i++) if(B.x[i][0]||B.x[i][1]||B.x[i][2]) 	C+=makeN(F[Z[z++]=i]);



  for(i=0;i<z;i++) { for(j=0;j<r;j++) if(SameRayBZ(i,R[j][0]))	   // Make RAYS

    {R[j][(nrp[j])++]=i;break;} if(j==r) {nrp[r]=1; R[r++][0]=i;} }



#if (TRACE_TRIANGULATION)

    prnI(p,C);

    printf("=C(Gale) ");

    Print_CMatrix(B,"Gale");

    assert(f==A.v);

    AuxPrintRays(R,nrp,r);

    assert(z>=r);

#endif



//	3d secondary fan A L G O R I T H M

//

//	list triangles / edges / crossing edges; subdivide; check Euler

//	* no multiple crossing: crossing -> 4 chambers, removing 4 min.trian.

//	* multiple crossings: add crossings and split edges; then either make

//	*-* minimal triangles for refined edges; overlaps -> chamber=union. or:

//	*-* compose split edges to minimal polygons

//	general dim: hyperplan wedges split by intersections; combine pieces



  for(c=2;c<r;c++)for(b=1;b<c;b++)for(a=0;a<b;a++)if((k=BZRx(a,b,c))){ //tn TRI

    int ea=1,eb=1,ec=1; Tli[tn][2]=c;

    if(k>0){Tli[tn][0]=a;Tli[tn][1]=b;} else {Tli[tn][0]=b;Tli[tn][1]=a;}

    for(k=0;k<en;k++){if(Einc[k]==(makeN(b)+makeN(c))) ea=0;	    // en EDGES

      if(Einc[k]==(makeN(a)+makeN(c))) eb=0;

      if(Einc[k]==(makeN(a)+makeN(b))) ec=0;}

    if(ea){Eli[en][0]=b;Eli[en][1]=c; Einc[en++]=makeN(b)+makeN(c);}

    if(eb){Eli[en][0]=a;Eli[en][1]=c; Einc[en++]=makeN(a)+makeN(c);}

    if(ec){Eli[en][0]=a;Eli[en][1]=b; Einc[en++]=makeN(a)+makeN(b);} tn++;}



  for(a=0;a<en;a++){for(b=0;b<r;b++) if(!Inci64_LE(makeN(b),Einc[a])) // remove

    if(XYZcone(BZR(Eli[a][0]),BZR(b),BZR(Eli[a][1]))>0)break;if(b<r){ // SPLIT

#if (TRACE_TRIANGULATION)                                         // EDGES

      printf("Ray %d divides edge[%d]=%d%d\n",b,a,Eli[a][0],Eli[a][1]);

#endif

      Einc[a]=Einc[--en]; Eli[a][0]=Eli[en][0]; Eli[a--][1]=Eli[en][1]; }}







  for(l=0;l<en;l++) ien[l]=0; for(l=1;l<en;l++) 	// INTERSECTing EDGES

  for(k=0;k<l;k++) if(0<InterSectE(BZRE(k,0),BZRE(k,1),BZRE(l,0),

    BZRE(l,1))) {IEli[k][ien[k]++]=l;IEli[l][ien[l]++]=k;}



  for(i=0;i<en;i++) if(ien[i]){if(iem<ien[i])iem=ien[i]; ies+=ien[i]; nQuad++;}



  assert(ies%2==0); assert(2*ies<VERT_Nmax);	ies/=2;



  if(iem<2) assert(2*ies==nQuad);		// nQuad redundant => eliminate



#if (TRACE_TRIANGULATION)

// for(l=0;l<r;l++)printf("RZF=%d:%d:%d ",*R[l],Z[*R[l]],F[Z[*R[l]]]);puts("");

   	for(i=0;i<en;i++){

   	printf("E%d=%d%d ",i,Eli[i][0],Eli[i][1]);

   	}

   	printf(" Edges=%d\n",en);

   	for(i=0;i<tn;i++)

   	  printf("T%d=%d%d%d ",i,Tli[i][0],Tli[i][1],Tli[i][2]);

   	printf(" Triangles=%d\n",tn);

	for(i=0;i<en;i++)

	  if(ien[i]){

	    printf("ie[%d]=%d%d->#",i,Eli[i][0],Eli[i][1]);

	    for(l=0;l<ien[i];l++)

	      printf("%d ",IEli[i][l]);

	    puts("");

	  }

    printf("ISecEdge iemax=%d iesum=%d\n",iem,ies);

#endif







	///// ==========	n - G O N  case  /  make chambers :: CR[ncr]



  if(iem>1){int OE[8*VERT_Nmax][2],noe=0,y,q=0; // orient.edges, y = #Rays+#Q's

    Long *Y[VERT_Nmax],YY[VERT_Nmax][3];

    int nse=0;for(i=0;i<en;i++)nse+=1+ien[i]; // printf("SubdivEdge=%d\n",nse);

    for(y=0;y<r;y++) Y[y]=BZR(y); Y[y]=YY[0];	assert(nse<=VERT_Nmax);



#ifdef	CANNOT_HANDLE_CONINCIDENT_INTERSECTION_POINTS	// o.k. for iem==2 ???

{   Inci64 Qinc[VERT_Nmax];

    for(i=0;i<en;i++) switch(ien[i]) {	// edges: DOUBLE; if(ien[i]) SPLIT;

      case 0: OE[noe][0]=OE[noe+1][1]=Eli[i][0];

	   OE[noe][1]=OE[noe+1][0]=Eli[i][1]; noe+=2;			break;

      case 1: if((j=IEli[i][0])>i) { k=q+r;	// if new add Q to YY's

	     IntersectEdges(BZRE(i,0),BZRE(i,1),BZRE(j,0),BZRE(j,1),Y[y]);

	     Qinc[q]=makeN(i)+makeN(j);Y[++y]=YY[++q];}	  // Y[r+q]::E_i & E_j

	   else {for(k=0;k<q;k++) if( Qinc[k]==(makeN(i)+makeN(j)) ) break;

	     assert(k<q); k+=r; }

	   OE[noe][0]=OE[noe+1][1]=Eli[i][0];OE[noe][1]=OE[noe+1][0]=k; noe+=2;

	   OE[noe][0]=OE[noe+1][1]=Eli[i][1];OE[noe][1]=OE[noe+1][0]=k; noe+=2;

//	printf("XXXX	edge %d intersect with %d ... k=%d\n",i,IEli[i][0],k);

	   assert(y==q+r);						break;

      case 2: for(f=0;f<2;f++){if((j=IEli[i][f])>i){k=q+r;// if new add Q to YY

	     IntersectEdges(BZRE(i,0),BZRE(i,1),BZRE(j,0),BZRE(j,1),Y[y]);

	     Qinc[q]=makeN(i)+makeN(j);Y[++y]=YY[++q];}	  // Y[r+q]::E_i & E_j

	   else {for(k=0;k<q;k++) if( Qinc[k]==(makeN(i)+makeN(j)) ) break;

	     assert(k<q); k+=r; }

//	printf("****	edge %d intersect with %d\n",i,j);

	   OE[noe][f]=OE[noe+1][1-f]=k;} 		     // add middle-edge

	   assert( f=XYZcone (Y[OE[noe][0]],Y[OE[noe][1]],Y[Eli[i][1]]) );

           if(f>0){k=OE[noe][0];f=OE[noe][1];}

	   else   {f=OE[noe][0];k=OE[noe][1];} noe+=2;	    // edge=(i0,k,f,i1)

//	printf("edge(%d)=%d %d %d %d\n",i,Eli[i][0],k,f,Eli[i][1]);

	assert(XYZcone(Y[f],Y[k],Y[Eli[i][0]])>0);	    // add end-edges

	OE[noe][0]=OE[noe+1][1]=Eli[i][0]; OE[noe][1]=OE[noe+1][0]=k; noe+=2;

	OE[noe][0]=OE[noe+1][1]=Eli[i][1]; OE[noe][1]=OE[noe+1][0]=f; noe+=2;

									break;

      default: puts("#(intersect.edges)>2 in Triang3dSFan() TO DO"); exit(0);}

    assert(2*nse==noe);	assert(2-y+nse<=VERT_Nmax); // Euler == y-nse+ncr == 2

}

#endif	// => first make intersections, then split edges, then add unsplit



    for(i=0;i<en;i++) if(ien[i]) {	// Y[k]::Y[y]=YY[q]=intersection(Ei,Ej)

      OE[noe][0]=Eli[i][0]; OE[noe++][1]=Eli[i][1];	// intersecting edges

      for(f=0;f<ien[i];f++) if((j=IEli[i][f])>i){	// intersection points

	IntersectEdges(BZRE(i,0),BZRE(i,1),BZRE(j,0),BZRE(j,1),Y[y]);

      	for(k=r;k<y;k++) if(SameRay(Y[k],Y[y],3)>0) break;

      	if(k==y) Y[++y]=YY[++q];}}			// Y[r+q]::Ei & Ej

    for(f=r;f<y;f++){int s=noe; 				// split edges

//for(l=0;l<noe;l++)printf("oe%d=%d%d ",l,OE[l][0],OE[l][1]);puts(" split:");

      for(i=0;i<noe;i++) if(0<XYZcone(Y[OE[i][0]],Y[f],Y[OE[i][1]])) {

      	OE[s][1]=OE[i][1]; OE[s++][0]=OE[i][1]=f;} noe=s;}	// add unsplit:

    for(i=0;i<en;i++)if(!ien[i]){OE[noe][0]=Eli[i][0];OE[noe++][1]=Eli[i][1];}



//for(l=0;l<noe;l++)printf("oe%d=%d%d ",l,OE[l][0],OE[l][1]);puts(" double:");

//for(f=0;f<3;f++){for(k=0;k<y;k++)printf("%3ld",Y[k][f]);puts(" =Y");}exit(0);



#ifdef	FIRST_TRY__TOO_COMPLICATED_BUT_MIGHT_BE_VIABLE

    for(i=0;i<en;i++) if(ien[i]) for(f=1;f<=ien[i];f++){     // add split edge

      if((j=IEli[i][f-1])>i){     		   // if new add Y[k]=Q to YY's

        IntersectEdges(BZRE(i,0),BZRE(i,1),BZRE(j,0),BZRE(j,1),Y[y]);

	for(k=r;k<y;k++) if(SameRay(Y[k],Y[y],3)) break;

	if(k==y){Qinc[q]=makeN(i)+makeN(j);Y[++y]=YY[++q];} //Y[r+q]::E_i & E_j

	else Qinc[k-r]|=makeN(i)+makeN(j);}

      else {for(k=0;k<q;k++) if(Inci64_LE(makeN(i)+makeN(j),Qinc[k])) break;

for(l=k+1;l<q;l++) assert(!Inci64_LE(makeN(i)+makeN(j),Qinc[l]));

	assert(k<q); k+=r;}		// unique ray with Y[k]::E_i & E_j

printf("#ie[%d]::i=%d/j=%d split at k=%d y=%d q=%d\n",f-1,i,j,k,y,q);

for(l=0;l<=noe+f;l++)printf("oe%d=%d%d ",l,OE[l][0],OE[l][1]);puts(" add:");

      if(f==1) {OE[noe][0]=Eli[i][0]; OE[noe][1]=OE[noe+1][1]=k;

assert(0<XYZcone(Y[Eli[i][0]],Y[k],Y[Eli[i][1]]));

	OE[noe+1][0]=Eli[i][1];}		// OE[noe/noe+1]=1st split(E_i)

      else {for(l=noe;l<noe+f;l++)		// split OE[noe<= ... <noe+f]

	if(XYZcone(Y[OE[l][0]],Y[k],Y[OE[l][1]])>0) break;// split OE_l by Y_k

printf("OE[%d]=%d%d k=%d .. l=%d noe=%d f=%d\n",l,OE[l][0],OE[l][1],k,l,noe,f);

	assert(l<noe+f); OE[noe+f][0]=OE[l][0]; OE[noe+f][1]=OE[l][0]=k;

for(l++;l<noe+f;l++) assert(-1==XYZcone(Y[OE[l][0]],Y[k],Y[OE[l][1]]));

	}

      if(f==ien[i]) noe+=f+1;}

    else {OE[noe][0]=Eli[i][0];OE[noe++][1]=Eli[i][1];}	// non-intersecting

#endif



    assert(noe<=4*VERT_Nmax); 	nse=noe;	// double/revers orientation

    for(i=0;i<noe;i++){OE[noe+i][0]=OE[i][1]; OE[noe+i][1]=OE[i][0];} noe*=2;



#if (TRACE_TRIANGULATION)

//	for(f=0;f<3;f++){for(k=0;k<q;k++)printf("%3ld",YY[k][f]);puts(" =YY");}

	for(f=0;f<3;f++){

	  for(k=0;k<y;k++)

	    printf("%3ld",Y[k][f]);

	  puts(" =Y");

	}

	for(k=0;k<noe;k++)

	  printf(" %d%d",OE[k][0],OE[k][1]);

	  puts(" =OE");

//	for(i=0;i<en;i++)if(ien[i]){printf("ie[%d]=%d%d->#",i,Eli[i][0],

//	Eli[i][1]);for(l=0;l<ien[i];l++)printf("%d ",IEli[i][l]); puts("");}

#endif



    for(i=0;i<noe;i++){int face[VERT_Nmax];		//// make CHAMBERS CR[]

      face[0]=OE[i][0]; face[f=1]=OE[i][1];

      do {int fit=0;					// assemble face

        for(j=i+1;j<noe;j++) if(OE[j][0]==face[f]){	// search fitting edges

	  if(XYZproduct(Y[face[f-1]],Y[face[f]],Y[OE[j][1]])<=0)continue;

	  if(fit)if(XYZproduct(Y[face[f]],Y[OE[k][1]],Y[OE[j][1]])<=0)continue;

	  fit=1; k=j;}

	face[++f]=OE[k][1]; OE[k][0]=OE[--noe][0];OE[k][1]=OE[noe][1];} // next

      while(face[f]!=face[0]);				// face complete

      assert(f>1); for(j=0;j<3;j++) CR[ncr][j] = Y[face[0]][j]+Y[face[2]][j];

      if(f==3) for(j=0;j<3;j++) CR[ncr][j] += Y[face[1]][j]; // chamber ray CR

#if (TRACE_TRIANGULATION)

	printf("face[i=%d]=",i);

	for(j=0;j<=f;j++)

	  printf("%d ",face[j]);

	printf(" noe=%d CR=",noe);

	for(j=0;j<3;j++)

	  printf(" %ld",CR[ncr][j]);

	puts("");

#endif

      ncr++;}						//// 	chambers DONE

    assert(ncr-nse+y==2);} 				// check Euler number



	///// ==================		squares and triangles only:



  else {Inci64 IEI[VERT_Nmax]; int nie=0; ncr=0;	// ies==#(squares)

    for(i=0;i<en;i++) if(ien[i]) if((j=IEli[i][0])>i) {Long Q[3]; // make QUAD

      IEI[nie++]=Einc[i]; assert(Einc[i]==makeN(Eli[i][0])+makeN(Eli[i][1]));

      IEI[nie++]=Einc[j]; assert(Einc[j]==makeN(Eli[j][0])+makeN(Eli[j][1]));

      IntersectEdges(BZRE(i,0),BZRE(i,1),BZRE(j,0),BZRE(j,1),Q);

      for(a=0;a<2;a++)for(b=0;b<2;b++){Long *rA=BZRE(i,a),*rB=BZRE(j,b);

        for(c=0;c<3;c++) CR[ncr][c]=Q[c]+rA[c]+rB[c]; ncr++;}

      }	assert(nie==2*ies);	chi=r-en+3*ies; // QUAD DONE, NOW MIN.TRIANG:

    for(i=0;i<tn;i++){a=Tli[i][0]; b=Tli[i][1]; c=Tli[i][2];

      for(k=0;k<r;k++) if(k!=a) if(k!=b) if(k!=c)	// if (k not in (abc))

      	if((0<=BZRx(a,b,k))&&(0<=BZRx(b,c,k))&&(0<=BZRx(c,a,k))) break;

      if(k==r){Inci64 Iabc=makeN(a)+makeN(b)+makeN(c); // & (diag not in (abc))

	for(j=0;j<nie;j++)if(Inci64_LE(IEI[j],Iabc))break;	// then add CR

	if(j==nie) {

	  for(k=0;k<3;k++)

	    CR[ncr][k]=BZR(a)[k]+BZR(b)[k]+BZR(c)[k];

	  ncr++;

	  chi++;

#if (TRACE_TRIANGULATION)

	    printf("Tmin=%d%d%d ",a,b,c);

    }}}

      printf("min.tri=%d #chamb=%d chi=%d\n",ncr-2*nie,ncr,chi);

#else

	  }}}

#endif

    assert(chi==2);} 	// NO MULTIPLE InterSectE DONE, now make MAX.TRIANG.



	///// ======	CR[ncr] = Chamber Rays = Triangulations, MaxTri <= ncr



  CT[*nmt=0]=X; // CT[#][list] -> X, *nmt=#MaxTri nt[num]=#Triang

  for(i=0;i<ncr;i++){Inci64 *S=CT[*nmt], U=0;	// TRIANGulate forall chambers

    Long *Q=CR[i]; nt[*nmt]=0; for(j=0;j<tn;j++)

    {Long a,b,c; int ChamberTriangle = 0;

    	if(0<=(a=XYZproduct(BZR(Tli[j][0]),BZR(Tli[j][1]),Q)))

    		if(0<=(b=XYZproduct(BZR(Tli[j][1]),BZR(Tli[j][2]),Q)))

    			if(0<=(c=XYZproduct(BZR(Tli[j][2]),BZR(Tli[j][0]),Q)))

    			{assert(a*b*c>0); ChamberTriangle = 1;}

    	if(ChamberTriangle){

    		for(a=0;a<nrp[Tli[j][0]];a++)for(b=0;b<nrp[Tli[j][1]];b++)for(c=0;c<

	nrp[Tli[j][2]];c++){S[nt[*nmt]]=C-makeN(F[Z[R[Tli[j][0]][a]]])-makeN(F[

	Z[R[Tli[j][1]][b]]])-makeN(F[Z[R[Tli[j][2]][c]]]);U|=S[nt[*nmt]++];}}}

    if(U==C) {CT[*nmt+1]=&CT[*nmt][nt[*nmt]]; tmt+=nt[(*nmt)++];}

#ifdef	PRINT_NON_MAXIMAL_TRIANGULATIONS

	j=(*nmt)-(U==C);if(U==C)printf("maximal[");else printf("non-max[");

	prnI(p,C);printf(" > ");prnI(p,U);printf("]:");for(a=0;a<nt[j];a++){

	printf(" ");prnI(p,CT[j][a]);}  puts("");

#endif

    } assert(*nmt);   		// check total #triangle in max.triangulations:

  j=0;for(i=0;i<*nmt;i++)j+=nt[i];assert(j==tmt); return tmt;

  }











/* INDUCED TRIANGulation t_i<n with f=\cup(t_i) induced by T_j<N, f < \cup(T_i)

 */



int Induce_Facet_Tri(Inci64 *T,int N,Inci64 f,Inci64 C, Inci64 *iT,int p,int d){

    int i=p,n=1; Inci64 A=f&(~C); (*iT)=(*T)&C; for(i=1;i<N;i++){int j;

    Inci64 X=T[i]&C; for(j=0;j<n;j++) if(Inci64_LE(X,iT[j])) break;

    if(j==n){for(j=0;j<n;j++)if(Inci64_LE(iT[j],X))iT[j--]=iT[--n];iT[n++]=X;}}

  for(i=0;i<n;i++) {iT[i]|=A; assert(Inci64_abs(iT[i])==d);} assert(A+C==f);

#ifdef	TRACE_TRIANGULATION__

	{int a=-1;printf("T=");for(i=0;i<N;i++){prnI(p,T[i]);printf(" ");}

	printf("with f-C=");prnIP(A); printf(" -> iT=");for(i=0;i<n;i++)

	{prnIP(iT[i]);printf(" ");if(Inci64_abs(iT[i])!=d)a=i;} puts("");

	if(a>=0)assert(Inci64_abs(iT[a])==d);}

					static int s; if(++s>13)exit(0);

#endif

	return n;

}



int Compatible_Tri(Inci64 CA,Inci64 CB,int a,Inci64 *A,int b,Inci64 *B,int p){

  Inci64 I=CA&CB, IA[VERT_Nmax], IB[VERT_Nmax]; int i, ia=1, ib=1;

  *IA=(*A)&I; for(i=1;i<a;i++) {int j; Inci64 X=A[i]&I; for(j=0;j<ia;j++)

    if(Inci64_LE(X,IA[j])) break; if(j==ia) {for(j=0;j<ia;j++)

      if(Inci64_LE(IA[j],X)) IA[j--]=IA[--ia]; IA[ia++]=X;}}

  *IB=(*B)&I; for(i=1;i<b;i++) {int j; Inci64 X=B[i]&I; for(j=0;j<ib;j++)

    if(Inci64_LE(X,IB[j])) break; if(j==ib) {for(j=0;j<ib;j++)

      if(Inci64_LE(IB[j],X)) IB[j--]=IB[--ib]; IB[ib++]=X;}}

#if	TRACE_TRIANGULATION

  puts("Compatibility of induced triangulations:");

  prnI(p,CA);printf(" & "); prnI(p,CB); printf(" = "); prnI(p,CA&CB); puts("");

	printf("IA[%d]=",ia);for(i=0;i<ia;i++){prnI(p,IA[i]);printf(" ");}

	printf("[");prnI(p,I);printf("] <- "); prnI(p,CA);printf("=A[%d]=",a);

	for(i=0;i<a;i++){prnI(p,A[i]),printf(" ");}	puts("");

	printf("IB[%d]=",ib);for(i=0;i<ib;i++){prnI(p,IB[i]);printf(" ");}

	printf("[");prnI(p,I);printf("] <- "); prnI(p,CB);printf("=B[%d]=",b);

	for(i=0;i<b;i++){prnI(p,B[i]),printf(" ");}	puts("");

#endif

  if(ia!=ib) return 0; else for(i=0;i<ia;i++) if(IA[i]!=IB[i]) break;

  if(i==ia) return 1; else for(i=0;i<ia;i++) {int sum=0, j;

    for(j=0;j<ib;j++) if(IA[i]==IB[j]) sum++; if(sum!=1) return 0;}

  return 1;}





//	list maximal 2ndary fans of facets with descending dimensions

//	forall compatible max triangulations make (induced) triang of facets



void GKZsubdivide(Inci64 *F,int f,PolyPointList *P,int p,int *Tp,int *ntp,

  int nPS, MORI_Flags *_Flag, FibW *_F ){ int c,d=P->n, d2, i,j,t=00; // d2=dim(2ndaryFan)

  Inci64 T[naT],*X=&T[nPS],C[ANfan],*CT[ANfan][ANtri];	// C=circuit, CT=triang

  int nmt[ANfan], nt[ANfan][ANtri], nmf=0;  // NumMaxTri, NumTriang, NumMax2Fan

  Inci64 MT[VERT_Nmax*POLY_Dmax], *_CT[ANfan];

  int TpMax=0, nMax=00, _nt[ANfan], I[ANfan], comptri=0;



  for(i=0;i<nPS;i++) {

	  if(Tp[c=ntp[i]]>TpMax) TpMax=Tp[ntp[nMax=i]];

	  j=Inci64_abs(F[c]);

	  T[i]=FindPolyCircuits(P,p,F[c],j);

  }

  for(d2=TpMax-d; d2>0; d2--)

	  for(i=0;i<nPS;i++)

		  if(Tp[c=ntp[i]]==d2+d){

			  for(j=0;j<nmf;j++)

				  if(Inci64_LE(T[i],C[j])) break;

			  assert(Inci64_LE(T[i],F[c]));

			  if(j==nmf)

				  switch(d2) {			// TRIANGULATE POLY_CIRCUITS:

				  case 1: t=Triang1dSFan(P,p,F[c],X,CT[nmf],&nmt[nmf],nt[nmf]); break;

				  case 2: t=Triang2dSFan(P,p,F[c],X,CT[nmf],&nmt[nmf],nt[nmf]); break;

				  case 3: t=Triang3dSFan(P,p,F[c],X,CT[nmf],&nmt[nmf],nt[nmf]); break;

				  default: printf("dim(2ndaryFan)=%d: to be done!\n",d2);

				  exit(0);}

#if	(TRACE_TRIANGULATION)

			  if(j<nmf)

				printf("F%d 2ndary dim=%d induced by C%d\n",ntp[i],d2,j);

#endif

			  if(j >= nmf){

#if	(TRACE_TRIANGULATION)

				printf("F%d 2ndary dim=%d new circ.= C%d #maxTri=%d\n",ntp[i],d2,j,nmt[j]);

#endif

				X=&X[t+1]; // CHECK !!

			  	if(nmt[nmf]>1)

				  	for(j=0;j<nmf;j++)

					  	if(nmt[j]>1)

						  	if(Inci64_abs(T[i]&C[j]) > 3){

								printf("Check compatibility of induced triangulation: ");

								prnI(p,T[i]);

							  	printf(" & ");

								prnI(p,C[j]);

								printf(" = ");

								prnI(p,T[i]&C[j]);

								puts("");

							  	comptri=1;

							}

			C[nmf++]=T[i];}

		  }



#ifdef	TRACE_TRIANGULATION__

  if(comptri){puts("Maximal triangulations from secondary fans:");

    for(j=0;j<nmf;j++)Print_MaxTrian(C[j],CT[j],nmt[j],nt[j],p);}

#endif



//	now make maximal triangulations MT_{j<c}



  c=0;

  for(j=0;j<f;j++) // f=#facets

	  if(d==Inci64_abs(F[j])) //d=dim Polytope

		  MT[c++]=F[j]; assert(c+nPS==f);



  t=Init_Multiloop(nmt,I,&c,&nmf);

  do{

    int k, mt=f-nPS, goodtri=1;



//  Inci64 CC; j=0; for(k=1;k<nmf;k++) for(j=0;j<k;j++)

//    if(Inci64_abs(CC=C[j]&C[k])>2) if(!Compatible_Tri(C[j],C[k],

//	nt[j][I[j]],CT[j][I[j]],nt[k][I[k]],CT[k][I[k]],p)) j=k=nmf+1;

//printf("ML(j=%d)=",j);for(k=0;k<nmf;k++)printf("%d",I[k]);puts("");if(j<nmf)



    for(k=0;k<nmf;k++) { _CT[k]=CT[k][I[k]]; _nt[k]=nt[k][I[k]]; } // combine

    if(comptri){ for(k=1;k<nmf;k++) if(nmt[k]>1) 	      // check compatib

      for(j=0;j<k;j++)if(nmt[j]>1)if(Inci64_abs(C[j]&C[k])>3) // triangulations

      if(!Compatible_Tri(C[j],C[k],_nt[j],_CT[j],_nt[k],_CT[k],p))goodtri=0;}

    if(goodtri){

      for(j=0;j<nPS;j++)

        for(k=0;k<nmf;k++)

          if(Inci64_LE(T[j],C[k])){

            mt+=Induce_Facet_Tri(_CT[k],_nt[k],F[ntp[j]],T[j],&MT[mt],p,d);

            break;

          }





      InterSectionRing(MT,&mt,P,p,_Flag,_F);

    }

  t--;

  }

  while(Multiloop(nmt,I,&c,&nmf)); assert(t==0);

}









/*	S u b d i v i d e ():	triangulate and call InterSectionRing

 *

 *   Incidences of ni=*t bounding equations on I[]; subdivide -> T[*t];

 *      assuming 4 dimensions and at most 3 non-vertices, etc.;

 *      triangulate and call "InterSectionRing(T,t,P,p);" for each T[]

 *      if simplicial add <=3 pts by hand else use GKZ for <= 3d secondary fan

 *   Old code: assuming: initial simplex, i.e. facets are 5 tetrahedra;

 *   3 add. points subdivision by points on edges in each step

 *    -> find points on edge = at intersection of >= 3 max-dim. triangles

 */

void Subdivide(PolyPointList *P,int v,Inci64 I[],int p,Inci64 *T,int *t,

		MORI_Flags *_Flag, FibW *F){

  int i,j,ni=*t,d=P->n, ns2=0, nPS=0, nVS=0, 	// ni=E.ne

    C[3][3],c[3]={0,0,0}, ntp[VERT_Nmax], Tv[VERT_Nmax], Tp[VERT_Nmax];



//  if((P->n!=4)||(p>v+3)) {puts("Mori cone only implemented for 3-fold ");

//	puts("hypersurfaces with <= 3 non-vertices"); exit(1);}



  /*

   * Determine the facets that are not simplicial looking for those that

   * have more points than d+1 with d=dim. of the lattice.

   * Gives the number of non-simplicial facets: nPS.

   * Writes a vector, ntp, with nPS entries, whose entries are the numbers

   * corresponding the non simplicial facets: e.g. ntp=(2,4,5) means that

   * there are nPS=3 non simplicial facets, precisely the 2nd, 4th and 5th facet.

   * */

  for(i=0;i<ni;i++) T[i]=I[i];

  /* runs over the facets */

  for(j=0;j<ni;j++) {Tv[j]=0; 			// Tv[] = #vertices on I[]

    /* counts the number of vertexes on the j-facet */

    for(i=0;i<v;i++) if(getN(i,I[j])) Tv[j]++;	// Tp[] = #points on I[]

    /* if the j-facet has more than d vertexes writes an entry +1 in nVS */

    if(Tv[j]>d) nVS++; // nVS = #not-vertex-simplicial

    /* control: facets must have always at least d vertexes */

    else assert(Tv[j]==d);

    Tp[j]=Tv[j]; 				// ntp[] = non-P.Simpl. facets

    /* counts the number of points on the j-facet */

    for(i=v;i<p;i++) if(getN(i,I[j])) Tp[j]++;	// nPS = #not-point-simplicial

    /* if the j-facet is not simplicial +1 to nPS and add entry j at place nPS

     * of vector ntp that keeps track of non-simplicial facets */

    if(Tp[j]>d)ntp[nPS++]=j; 	// C[i][f]= facets containing P.x[v+i], f<c[i]

    /* control that j-facet is in fact a facet i.e. at least d points */

    else assert(Tp[j]==d);}	// ns2= # on triangles  =>  p-v-ns2 on edges



  if(nPS==0){ // no triang. needed

    InterSectionRing(T,t,P,p,_Flag,F);

    return;

  }



  /*

   * C[i][] is a vector whose entries are the numbers corresponding to the

   * facets that contain the non-vertex point x_i.

   */

  for(i=v;i<p;i++){for(j=0;j<ni;j++)		// find non-vertex positions:

    if(getN(i,I[j])) C[i-v][(c[i-v])++]=j;   // C[i][]=facets containing i

    assert(c[i-v]>1); if(c[i-v]==2) ns2++;}	// ns2 = #pts @ codim2=triangle



#if	TRACE_TRIANGULATION

    for(j=0;j<ni;j++)

      printf("vp[%d]=%d%d%s",j,Tv[j],Tp[j],(j==ni-1)?"\n":" ");

    for(i=0;i<p-v;i++)

      if(c[i]>2){

        printf("p%d on edge:",i+v);

        for(j=0;j<c[i];j++)

          printf(" %d",C[i][j]);

        puts("");

      }

      else printf("p%d on triangle(%d,%d)\n",i+v,C[i][0],C[i][1]);

#endif



  GKZsubdivide(I,ni,P,p,Tp,ntp,nPS,_Flag,F); return;

  }







// count solutions in PM of PM->x[i].PN->x[j]+\d(k,j)>=0 forall j

int SectionCount(PolyPointList *PN,int k,PolyPointList *PM){int i,n=0;

  for(i=0;i<PM->np;i++){int j; for(j=0;j<PN->np;j++) {Long e = (k==j);

      int l; for(l=0;l<PM->n;l++) e+=PM->x[i][l]*PN->x[j][l]; if(e<0) break;}

    if(j==PN->np) n++;} return n;}



#if(NewtonMonomCOORD)

#if((NewtonMonomCOORD<1)||(NewtonMonomCOORD>2))

#error		NewtonMonomCOORD has to be 1 or 2

#endif

void NewtonMonomial(Long *X,int d){int num=0,den=0,i,t=NewtonMonomCOORD==1;

  for(i=0;i<d;i++) if(X[i]>0) num++; else if(X[i]) den++;

  if(num) {for(i=0;i<d;i++) if(X[i]>0) {if(t)printf("t_%d",i+1);

    else {char uvw[2]={'u',0}; uvw[0]+=i;printf("%s",uvw);}

    if(X[i]>1)printf("^%d",(int)X[i]);}} else printf("1");

  if(den) {printf("/"); if(den>1)printf("(");

    for(i=0;i<d;i++) if(X[i]<0) {if(t)printf("t_%d",i+1);

    else {char uvw[2]={'u',0}; uvw[0]+=i;printf("%s",uvw);}

    if(X[i]<-1)printf("^%d",(int)-X[i]);} if(den>1)printf(")");}}

#endif





/* Hypersurface divisors Q(charges) permutes the N-lattice points of

 * non-intersecting divisors to the end of the PPL *_P, calls IP_Simplex_Fiber

 * and then prints the charges = linear relations

 */

void HyperSurfDivisorsQ(PolyPointList *_P,VertexNumList*V,EqList *E,FibW *F, MORI_Flags * _Flag){

  int i=V->nv,j, cp=_P->np-1,t=E->ne,d=_P->n, Dh0[VERT_Nmax], Dh2[VERT_Nmax];

  Inci64 I[VERT_Nmax], T[FACE_Nmax]; Long Wsum[VERT_Nmax];

  PolyPointList *DP = (PolyPointList *) malloc(sizeof(PolyPointList));

  EqList *DE = (EqList *) malloc(sizeof(EqList));





 /*************************************************************************

  * RECAST PPL

  *

  * Permutes the points that do not intersect facets

  * to the right of the poly point matrix

  */



  /* initialize the incidence relations I[1]=0,...,I[#facets]=0*/

  for(j=0;j<E->ne;j++)I[j]=0;

  /* i runs over the poly points */

  for(i=0;i<cp;i++){

	int fi=0;

	/* j runs over the facets eqn's */

    for(j=0;j<E->ne;j++)

    	/* add 1 to fi if the i-point lies on the j-facet.

    	 * OPTIM: can't we use INCI arithmetic instead of Eval_Eq_on_V? */

    	if(0==Eval_Eq_on_V(&(E->e[j]),_P->x[i],_P->n)) fi++;

    /* if fi<=1 then permutes the i-point with the (cp-1)-point.

     * And reduces cp to the number of intersecting divisors */

    if(fi<=1) Swap_Vecs(_P->x[i--],_P->x[--cp],_P->n);

  }

  /************************************************************************/



  /* Generates the IP_simplices among point of codim >1 */

  IP_Simplex_Fiber(_P->x,cp,_P->n,F,FIB_Nmax,0);



  if(_Flag->P){

  /************************************************************************

   * PRINT IP

   */



  // fi==2::dual to edges whose length is Dh0

  	/* prints the PPL */

	Print_PPL(_P,"points of P* and IP-simplices");



	for(i=0; i<cp; i++)

		fprintf(outFILE,"-----");



    /* The # IP simplices = the number of weight-relations (dim of matrix W) */

    fprintf(outFILE,"   #IP-simp=%d",F->nw);



    /* If there are more weight-relations than prim.div.cl. print the info

     * that there are more IP-simplexes than the Nr. of independent vectors */

    if(F->nw>cp-_P->n)

    	fprintf(outFILE," > %d=#pts-dim",cp-_P->n);



    fprintf(outFILE,"\n");



    /* Prints the weight matrix */

    for(i=0; i<F->nw; i++){ /* i runs over # weight relations = IP-simplices */

    	int cd=_P->n-cp+1 ;

        /* j runs over the facet intersecting divisors */

    	for(j=0;j<cp;j++)

    		if(!F->W[i][j]) cd++;



    	/* Prints the single rows of the matrix */

    	Aux_IPS_Print_WP(F->W[i],cp,cd); 			// if(F->ZS)



    	/* Prints the quozient group if any */

    	if(F->nz[i])

    		Print_QuotZ(&F->Z[F->n0[i]],&F->M[F->n0[i]],cp,F->nz[i]);

    	fprintf(outFILE,"\n");

    }

    /******************************************************************/

    }



  /********************************************************************

   * FACET INCIDENCES

   */



  /* First creates the Incidences I */

  for(i=0;i<cp;i++) {

	  int fi=0,FI[3];

	  for(j=0;j<E->ne;j++) // INCIdence

		  if(0==Eval_Eq_on_V(&(E->e[j]),_P->x[i],_P->n)) {

			  putN(i,&I[j]);

			  if(fi<3) FI[fi++]=j;

		  }

	  Wsum[i]=0;

	  for(j=0;j<F->nw;j++)Wsum[i]+=F->W[i][j];

	  Dh0[i]=1;

	  if(fi==2){

		  Long g=E->e[FI[1]].a[0]-E->e[FI[0]].a[0];		   // h0(D_i)

		  for(j=1;j<d;j++)g=NNgcd(g,E->e[FI[1]].a[j]-E->e[FI[0]].a[j]);

		  Dh0[i]=g;}

  }

  /*******************************************************************/



  if(_Flag->K){

   /******************************************************************

   * KREUZER polytope

   *

   * prints the PolyPointList matrix in polynomial form if precompiler its on

   */



	  printf("KreuzerPoly=");

	  for(i=0;i<V->nv;i++) {

		  if(i)printf("+");

		  NewtonMonomial(_P->x[i],d);

	  }

	  for(i=V->nv;i<cp;i++){

		  printf("-");

		  NewtonMonomial(_P->x[i],d);

	  }

	  printf(";");

	  if(cp<_P->np-1)

		  printf(" intpts=%d;",_P->np-cp-1);

	  j=cp-_P->n;

	  for(i=0;i<cp;i++)

		  if(Dh0[i]>1){

			  j+=Dh0[i]-1;

		          printf("  multd%d=%d;",i+1,Dh0[i]);

		  }

	  if(_Flag->H == 0)

		  printf("  Pic=%d",j);

	  printf("\n");

  /*****************************************************************/

  }



  if(_Flag->I){

  /*******************************************************************

  * PRINT INCIDENCES

  */





  /* prints the incidences */

  fprintf(outFILE,"Incidence:");

  /* j runs over the facets */

  for(j=0;j<E->ne;j++) {

	  fprintf(outFILE," ");

	  /* prints the incidence of the j-facet*/

	  fprI(cp,I[j]);}

  fprintf(outFILE,"\n");

  /******************************************************************/

  }







#if	(PRINT_MONOMIALS)

  for(j=0;j<E->ne;j++) {for(i=0;i<cp;i++)    		/* print monomials */

    printf(" %4ld",Eval_Eq_on_V(&(E->e[j]),_P->x[i],_P->n));puts(" monomial");}

#endif





// h02(D_i): enumerate $\sum c_kD_k=D_i$ with $c_i=-1$ and 0<=c_k for k!=i

// where "=" means equal charges, i.e. columns in F->W[n]; but: the $c_k$

// are sections $\chi_m$ i.e. $c_k=<m,v_k>=m_i v_ki$ only if $m_i$ are

// admitted to be fractional in case of torsion gradings. Hence the correct

// formula enumerates $m$ with $<m,v_j>+a_j\ge0$ for sections of $\sum a_jD_j$.



//    c=Init_Multiloop(N,L,&j,&J); Dh2[i]=0;	// 0<=I[j]<N[j], 0<=j<J

//    do {K[0]=L[0]=N[0]=J; c--;}

//    while(Multiloop(N,L,&j,&J)); assert(c==0);}



  Make_Dual_Poly(_P,V,E,DP);

  {

	  PairMat PM,DPM; Make_VEPM(_P,V,E,PM);

	  VNL_to_DEL(_P, V, DE); assert(Transpose_PM(PM,DPM,V->nv,E->ne));

	  Complete_Poly(DPM,DE,E->ne,DP);

  }

  for(i=0;i<cp;i++)

	  Dh2[i]=SectionCount(_P,i,DP)-1;

  free(DP);

  free(DE);	//for(i=0;i<cp;i++)printf("%d ",Dh2[i]);puts("=h02");



  Subdivide(_P,V->nv,I,cp,T,&t, _Flag ,F);



}

















