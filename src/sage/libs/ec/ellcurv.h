#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <pari/pari.h>

#define gequalgs(s,y)    (gequalsg((y),(s)))

typedef long long int int64;
#ifdef __alpha__
#define atoll(x) atol((x))
#define gint64(x) itos((x))
#endif
typedef long long unsigned int uint64;
#define ratdoub(inp) rtodbl(gmul(inp,dbltor(1.0)))
#define CL CURVES
struct MYST {double cost; int condcorr; int which;}; struct MYST *pile;

int INITANALRANK; double Manal;

int ELLACC; int TRACE; int VERBOSE; int ROOTNO; int CM;
int PRINT; int DISK; int plimit; int READ; int SPIT; GEN X1VOL;
GEN CURVE; GEN COND; GEN CF; GEN TAMA; int DOADD; FILE *OUTFILE;
GEN TWCOND; GEN TWCF; GEN TWCURVE; GEN TWPROD; GEN TWEFF;
GEN SSCOND; GEN EXNEG; GEN EXPOS; GEN EXEFF; int MDFORCE;
int ISOG; int NI; int AISOG; int PLACE; GEN CURVES; int PCOUNT; int HIRANK;
int ISPRIME; int ISSQFREE; int ISSETZER; int LOOKUP; int OFILE;
int64 PMAX; int PSIZE; int X0NUM; GEN MODDEG; int DOTWIST;

int c6r,c4r,c6s,c4s,numcurves,recsize,residues,c4neg,c6neg;
unsigned char *curvedata; int cdatsize; unsigned char **twdata;
unsigned char *size16bytes; unsigned char *size256bytes;
unsigned char *size4096bytes; unsigned char *leadbyte;
unsigned char *moddegbytes; unsigned char *localsizebytes;
unsigned char headbytes[8],backbytes[4];

double *anarray; double *antwarray; unsigned char *TWFERRY;
unsigned char *auxp; int ansize; int antwsize; int PREVD;

//analrank.c
int analrank(GEN,double*); int qtwistar(int,double*,GEN); GEN mytors(GEN);
//apcompute.c
void computeap(GEN,int); void computetwistarray(int,int);
//arderivs.c
void firstderivs(void);
// arintern.c
void analinitw(void); void recompM(double,int); void recompH(double,int);
double csm(double,int); double cmd(double,int,int); double clg(double,int,int);
double analw(int64,int);
//arith.c
GEN fundamentaltau(GEN*,GEN*); void orderrealroots(GEN*);
void allcubicroots(GEN,GEN*,int); GEN cubicrealroot(GEN,int);
GEN periodvolvec0(GEN,int); GEN periodvolvec(GEN);
GEN volratio(GEN,GEN); GEN myvol(GEN);
//artwists.c
int whichtwists(void); void readtwists(int); int rootnotwist(int);
void do2divismd(int); void outtwist(int,int,double);
void twistdump(int,int,double,GEN);
//arutil.c
int qtunhandle(unsigned char,unsigned char,int*,double*);
void qthandle(int,double,int*); int blocktwists(void);
void createpile(int);
/*int Ordering(struct MYST*, struct MYST*);*/
int Ordering(__const void *a, __const void *b) ;

//degphi.c
void recomp1(double*,double); void recomp2(double*,double);
void InitialiseWeights(void); void getprimes(int64,int);
double csm1(double); double cmd1(double,int); double clg1(double,int);
double csm2(double); double cmd2(double,int); double clg2(double,int);
double WEIGHT(int64); void ADDLL(int64,int,int64,int64,int64,int64,int64);
void degphi(GEN,GEN,int);
//diskio.c
void readfile(char*); void writefile(char*);
//dodisk.c //docurve.c
void dodisk(void);
void mdforce(void);
void docurve(char*,int,int,int,int,int);
GEN sage_docurve(char*,int,int,int,int,int);
//equation.c
GEN avecfromc4c6(GEN,GEN); GEN mineqfromc4c6(GEN,GEN);
GEN qtwist(GEN,GEN); void minimaltwist(void); GEN gettama(GEN);
GEN retminimaltwist(GEN,GEN*);
GEN minimalequation(void); void symsq(void);
//exotic.c
void checkexotic(void);
void exo11A(void); void exo11B(void); void exo11C(void); void exo14(void);
void exo15A(void); void exo15B(void); void exo17A(void); void exo17B(void);
void exo19(void); void exo21A(void); void exo21B(void);
void exo27(void); void exo37(void); void exo43(void); void exo67(void); void exo163(void);
//fixit.c checkit.c
void fixit(void); void checkit(int);
//iisog2.c
int findtwo(GEN,GEN*,GEN*,GEN*); GEN twogetcurve(GEN,GEN);
//iisog3.c
int findthree(GEN,GEN*,GEN*); GEN threegetcurve(GEN,GEN);
//isog.c
void isogenies(void);
void garbagecollect(GEN*,GEN*,GEN*,pari_sp);
//isog2.c //isog23.c //isog24.c
void dotwoisog(GEN,GEN*); void do3twoisog(GEN,GEN);
void do9twoisog(GEN,GEN,GEN);
void do4twoisog(GEN,GEN,GEN,GEN);
//isog3.c //isog5.c //isog52.c //isog713.c
void do3isog(GEN*,GEN*,GEN*); void do5isog(GEN); void do2fiveisog(GEN,GEN);
void do7isog(GEN); void do13isog(GEN);
//isoggen.c
int doisogreal(GEN,GEN*,int); int doisogimag(GEN,GEN*,int);
void shiftperiodi(int,GEN*,GEN,int); void shiftperiodr(int,GEN*,GEN*,int);
int doisoggeneric(GEN,GEN*,int,int);
//isogNprime.c
void isogN11(void); void isogN17(void); void isogN19(void); void isogN37(void); void isogNsz(void);
GEN szcN(GEN); GEN szcP(GEN);
//isogsort.c
void sorteight(GEN,GEN,GEN,GEN,GEN,GEN,int);
void sortsixteen(GEN,GEN,GEN,GEN,GEN,GEN,GEN,GEN,int);
void sorttwelve(GEN,GEN,GEN,GEN,GEN,GEN,GEN,GEN,int);
//isogx0.c
int x0isogeny(void); GEN myznstar(void);
//isogx0branch.c
int whichbent(GEN,GEN); int bentisog(GEN,GEN,GEN);
int realisog(GEN,GEN,GEN); int realisog2(GEN,GEN,GEN,GEN);
int imagisog(GEN,GEN,GEN); int imagisog2(GEN,GEN,GEN,GEN);
//isogx0branch*.c
int x0isog1(GEN,int); int x0isog2(GEN,int); int x0isog3(GEN,int);
int x0isog4(GEN,int); int whichone5(int,int); int x0isog5(GEN,int,int);
int x0isog6(GEN,int,int);
//isogx0getd.c
GEN createcycgenlist(GEN,int); GEN findincycliclist(GEN,GEN,GEN);
GEN getd1(GEN,int,GEN);
//isogx0period.c
GEN linearcomb(GEN,GEN,GEN); int numterms(double,int,double);
double minp(GEN,GEN); GEN computeperiod(GEN,GEN,double);
//lookup.c
int lookup(void); int retlen(GEN); void createnewfile(char*,GEN,GEN,GEN,GEN);
int findcurve(GEN,GEN); void insertcurve(int,GEN,GEN);
int comparecd(unsigned char*,int,unsigned char*,int);
//main.c //readit.c //special.c
int main(int,char**); void doread(void); void readit(int);
GEN specialcurve(char*); void spitout(GEN);
//util.c
int absval(int); GEN lltoi(int64); int64 gint64(GEN); int classsize(void);
int isogdeg(int); void printcurve(GEN); void printcurven(GEN);
int u8(int); int checkconductor(GEN); int gettheint(unsigned char*,int);
void get4bytes(int,unsigned char*); int uwlsb(unsigned char*,int);
void cv(int,GEN); void docount(void); GEN isprimepower(GEN);
GEN icbrt(GEN); int KnownDiffOptimal(GEN); GEN getcurve(char*);
GEN jcurve(char*);
