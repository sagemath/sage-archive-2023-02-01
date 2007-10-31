#ifndef TGB_INTERNAL_H
#define TGB_INTERNAL_H
//!\file tgb_internal.h
/****************************************
*  Computer Algebra System SINGULAR     *
****************************************/
/* $Id: tgb_internal.h,v 1.69 2007/02/27 10:40:57 bricken Exp $ */
/*
 * ABSTRACT: tgb internal .h file
*/
#include <omalloc.h>
#include "p_polys.h"

#include "ideals.h"
#include "ring.h"
#include "febase.h"
#include "structs.h"
#include "polys.h"
#include "stdlib.h"
#include <modulop.h>
//#define USE_NORO 1


#ifdef USE_NORO
#define NORO_CACHE 1
#define NORO_SPARSE_ROWS_PRE 1
#define NORO_NON_POLY 1
#include <algorithm>
#endif
#ifdef NORO_CACHE
//#include <map>
#include <vector>
#endif
//temporarily disabled because SAGE doesn't use boost (malb, 20071131)
//#ifdef HAVE_BOOST_DYNAMIC_BITSET_HPP
//#define  HAVE_BOOST 1
//#endif
//#define HAVE_BOOST 1
//#define USE_STDVECBOOL 1
#ifdef HAVE_BOOST
#include "boost/dynamic_bitset.hpp"
#include <vector>
using boost::dynamic_bitset;
using std::vector;
#endif
#ifdef USE_STDVECBOOL
#include <vector>
using std::vector;
#endif
#include "kutil.h"
#include "kInline.cc"
#include "kstd1.h"
#include "kbuckets.h"




//#define TGB_DEBUG
#define FULLREDUCTIONS
#define HANS_IDEA
//#define HALFREDUCTIONS
//#define HEAD_BIN
//#define HOMOGENEOUS_EXAMPLE
#define REDTAIL_S
#define PAR_N 100
#define PAR_N_F4 5000
#define AC_NEW_MIN 2
#define AC_FLATTEN 1

//#define FIND_DETERMINISTIC
//#define REDTAIL_PROT
//#define QUICK_SPOLY_TEST
class PolySimple{
public:
  PolySimple(poly p){
    impl=p;
  }
  PolySimple(){
    impl=NULL;
  }
  PolySimple(const PolySimple& a){
    //impl=p_Copy(a.impl,currRing);
    impl=a.impl;
  }
  PolySimple& operator=(const PolySimple& p2){
    //p_Delete(&impl,currRing);
    //impl=p_Copy(p2.impl,currRing);
    impl=p2.impl;
    return *this;
  }
  ~PolySimple(){
    //p_Delete(&impl,currRing);
  }
  bool operator< (const PolySimple& other) const{
    return pLmCmp(impl,other.impl)<0;
  }
  bool operator==(const PolySimple& other){
    return pLmEqual(impl,other.impl);
  }
  poly impl;

};
template<class number_type> class DataNoroCacheNode;
/*class MonRedRes{
public:
  poly p;
  number coef;
  BOOLEAN changed;
  int len;
  BOOLEAN onlyBorrowed;
  bool operator<(const MonRedRes& other) const{
    int cmp=p_LmCmp(p,other.p,currRing);
    if ((cmp<0)||((cmp==0)&&((onlyBorrowed)&&(!(other.onlyBorrowed))))){
      return true;
    } else return false;
  }
  DataNoroCacheNode* ref;
  MonRedRes(){
    ref=NULL;
    p=NULL;
  }
};*/
template <class number_type> class MonRedResNP{
public:
  number coef;


  DataNoroCacheNode<number_type>* ref;
  MonRedResNP(){
    ref=NULL;
  }
};
struct sorted_pair_node{
  //criterium, which is stable 0. small lcm 1. small i 2. small j
  wlen_type expected_length;
  poly lcm_of_lm;
  int i;
  int j;
  int deg;


};
#ifdef NORO_CACHE
#ifndef NORO_NON_POLY
class NoroPlaceHolder{
public:
  DataNoroCacheNode* ref;
  number coef;
};
#endif
#endif
//static ideal debug_Ideal;


struct poly_list_node{
  poly p;
  poly_list_node* next;
};

struct int_pair_node{
  int_pair_node* next;
  int a;
  int b;
};
struct monom_poly{
  poly m;
  poly f;
};
struct mp_array_list{
  monom_poly* mp;
  int size;
  mp_array_list* next;
};


struct poly_array_list{
  poly* p;
  int size;
  poly_array_list* next;
};
class slimgb_alg
{
  public:
    slimgb_alg(ideal I, int syz_comp,BOOLEAN F4);
		void introduceDelayedPairs(poly* pa,int s);
    virtual ~slimgb_alg();
    void cleanDegs(int lower, int upper);
#ifndef HAVE_BOOST
#ifdef USE_STDVECBOOL
  vector<vector<bool> > states;
#else
  char** states;
#endif
#else
  vector<dynamic_bitset<> > states;
#endif
  ideal add_later;
  ideal S;
  ring r;
  int* lengths;
  wlen_type* weighted_lengths;
  long* short_Exps;
  kStrategy strat;
  int* T_deg;
  int* T_deg_full;
  poly tmp_lm;
  poly* tmp_pair_lm;
  sorted_pair_node** tmp_spn;
  poly* expandS;
  poly* gcd_of_terms;
  int_pair_node* soon_free;
  sorted_pair_node** apairs;
  #if 0
  BOOLEAN* modifiedS;
  #endif
  #ifdef TGB_RESORT_PAIRS
  bool* replaced;
  #endif
  poly_list_node* to_destroy;
  //for F4
  mp_array_list* F;
  poly_array_list* F_minus;

  //end for F4
#ifdef HEAD_BIN
  struct omBin_s*   HeadBin;
#endif
  unsigned int reduction_steps;
  int n;
  //! array_lengths should be greater equal n;
  int syz_comp;
  int array_lengths;
  int normal_forms;
  int current_degree;
  int Rcounter;
  int last_index;
  int max_pairs;
  int pair_top;
  int easy_product_crit;
  int extended_product_crit;
  int average_length;
  int lastDpBlockStart;
  int lastCleanedDeg;
  BOOLEAN isDifficultField;
  BOOLEAN completed;
  BOOLEAN is_homog;
  BOOLEAN tailReductions;
  BOOLEAN eliminationProblem;
  BOOLEAN F4_mode;
  BOOLEAN nc;
  #ifdef TGB_RESORT_PAIRS
  BOOLEAN used_b;
  #endif
};
class red_object{
 public:
  kBucket_pt bucket;
  poly p;
  unsigned long sev;
  int sugar;
  void flatten();
  void validate();
  wlen_type initial_quality;
  void adjust_coefs(number c_r, number c_ac_r);
  wlen_type guess_quality(slimgb_alg* c);
  int clear_to_poly();
  void canonicalize();
};


enum calc_state
  {
    UNCALCULATED,
    HASTREP//,
    //UNIMPORTANT,
    //SOONTREP
  };
static BOOLEAN pair_cmp(sorted_pair_node* a,sorted_pair_node* b);
template <class len_type, class set_type>  int pos_helper(kStrategy strat, poly p, len_type len, set_type setL, polyset set);
static int add_to_reductors(slimgb_alg* c, poly h, int len, int ecart, BOOLEAN simplified=FALSE);
static int bucket_guess(kBucket* bucket);
static poly redNFTail (poly h,const int sl,kStrategy strat, int len);
static poly redNF2 (poly h,slimgb_alg* c , int &len, number&  m,int n=0);
void free_sorted_pair_node(sorted_pair_node* s, ring r);
static void shorten_tails(slimgb_alg* c, poly monom);
static void replace_pair(int & i, int & j, slimgb_alg* c);
//static sorted_pair_node** add_to_basis(poly h, int i, int j,slimgb_alg* c, int* ip=NULL);
static void do_this_spoly_stuff(int i,int j,slimgb_alg* c);
//ideal t_rep_gb(ring r,ideal arg_I);
static BOOLEAN has_t_rep(const int & arg_i, const int & arg_j, slimgb_alg* state);
static int* make_connections(int from, poly bound, slimgb_alg* c);
static int* make_connections(int from, int to, poly bound, slimgb_alg* c);
void now_t_rep(const int & arg_i, const int & arg_j, slimgb_alg* c);
static void soon_t_rep(const int & arg_i, const int & arg_j, slimgb_alg* c);
static int pLcmDeg(poly a, poly b);
static int simple_posInS (kStrategy strat, poly p,int len, wlen_type wlen);
static BOOLEAN find_next_pair(slimgb_alg* c, BOOLEAN go_higher=TRUE);

static sorted_pair_node* pop_pair(slimgb_alg* c);
static BOOLEAN no_pairs(slimgb_alg* c);
void clean_top_of_pair_list(slimgb_alg* c);
static void super_clean_top_of_pair_list(slimgb_alg* c);
static BOOLEAN state_is(calc_state state, const int & i, const int & j, slimgb_alg* c);
static BOOLEAN pair_better(sorted_pair_node* a,sorted_pair_node* b, slimgb_alg* c=NULL);
static int tgb_pair_better_gen(const void* ap,const void* bp);
static poly redTailShort(poly h, kStrategy strat);
static poly gcd_of_terms(poly p, ring r);
static BOOLEAN extended_product_criterion(poly p1, poly gcd1, poly p2, poly gcd2, slimgb_alg* c);
static poly kBucketGcd(kBucket* b, ring r);
static void multi_reduction(red_object* los, int & losl, slimgb_alg* c);
int slim_nsize(number n, ring r);
sorted_pair_node* quick_pop_pair(slimgb_alg* c);
sorted_pair_node* top_pair(slimgb_alg* c);
sorted_pair_node** add_to_basis_ideal_quotient(poly h, slimgb_alg* c, int* ip);//, BOOLEAN new_pairs=TRUE);
sorted_pair_node**  spn_merge(sorted_pair_node** p, int pn,sorted_pair_node **q, int qn,slimgb_alg* c);
int kFindDivisibleByInS_easy(kStrategy strat,const red_object & obj);
int tgb_pair_better_gen2(const void* ap,const void* bp);
int kFindDivisibleByInS_easy(kStrategy strat,poly p, long sev);
//static int quality(poly p, int len, slimgb_alg* c);
/**
   makes on each red_object in a region a single_step
 **/
class reduction_step{
 public:
  /// we assume hat all occuring red_objects have same lm, and all
  /// occ. lm's in r[l...u] are the same, only reductor does not occur
  virtual void reduce(red_object* r, int l, int u);
  //int reduction_id;
  virtual ~reduction_step();
  slimgb_alg* c;
  int reduction_id;
};
class simple_reducer:public reduction_step{
 public:
  poly p;
  kBucket_pt fill_back;
  int p_len;
  int reducer_deg;
  simple_reducer(poly p, int p_len,int reducer_deg, slimgb_alg* c =NULL){
    this->p=p;
    this->reducer_deg=reducer_deg;
    assume(p_len==pLength(p));
    this->p_len=p_len;
    this->c=c;
  }
  virtual void pre_reduce(red_object* r, int l, int u);
  virtual void reduce(red_object* r, int l, int u);
  ~simple_reducer();


  virtual void do_reduce(red_object & ro);
};

//class sum_canceling_reducer:public reduction_step {
//  void reduce(red_object* r, int l, int u);
//};
struct find_erg{
  poly expand;
  int expand_length;
  int to_reduce_u;
  int to_reduce_l;
  int reduce_by;//index of reductor
  BOOLEAN fromS;//else from los

};

static void multi_reduce_step(find_erg & erg, red_object* r, slimgb_alg* c);
static void finalize_reduction_step(reduction_step* r);

template <class len_type, class set_type>  int pos_helper(kStrategy strat, poly p, len_type len, set_type setL, polyset set){
  //Print("POSHELER:%d",sizeof(wlen_type));
  int length=strat->sl;
  int i;
  int an = 0;
  int en= length;

  if ((len>setL[length])
      || ((len==setL[length]) && (pLmCmp(set[length],p)== -1)))
    return length+1;

  loop
  {
    if (an >= en-1)
    {
      if ((len<setL[an])
          || ((len==setL[an]) && (pLmCmp(set[an],p) == 1))) return an;
      return en;
    }
    i=(an+en) / 2;
    if ((len<setL[i])
        || ((len==setL[i]) && (pLmCmp(set[i],p) == 1))) en=i;
    //else if ((len>setL[i])
    //|| ((len==setL[i]) && (pLmCmp(set[i],p) == -1))) an=i;
    else an=i;
  }

}
#ifdef NORO_CACHE
#define slim_prec_cast(a) (unsigned int) (unsigned long) (a)
#define F4mat_to_number_type(a) (number_type) slim_prec_cast(a)
typedef unsigned short tgb_uint16;
typedef unsigned char tgb_uint8;
typedef unsigned int tgb_uint32;
class NoroCacheNode{
public:
  NoroCacheNode** branches;
  int branches_len;


  NoroCacheNode(){
    branches=NULL;
    branches_len=0;

  }
  NoroCacheNode* setNode(int branch, NoroCacheNode* node){
    if (branch>=branches_len){
      if (branches==NULL){
        branches_len=branch+1;
        branches_len=si_max(branches_len,3);
        branches=(NoroCacheNode**) omalloc(branches_len*sizeof(NoroCacheNode*));
        int i;
        for(i=0;i<branches_len;i++){
          branches[i]=NULL;
        }
      }else{
        int branches_len_old=branches_len;
        branches_len=branch+1;
        branches=(NoroCacheNode**) omrealloc(branches,branches_len*sizeof(NoroCacheNode*));
        int i;
        for(i=branches_len_old;i<branches_len;i++){
          branches[i]=NULL;
        }
      }
    }
    assume(branches[branch]==NULL);
    branches[branch]=node;
    return node;
  }
  NoroCacheNode* getBranch(int branch){
    if (branch<branches_len) return branches[branch];
    return NULL;
  }
  virtual ~NoroCacheNode(){
    int i;
    for(i=0;i<branches_len;i++){
      delete branches[i];
    }
    omfree(branches);
  }
  NoroCacheNode* getOrInsertBranch(int branch){
    if ((branch<branches_len)&&(branches[branch]))
      return branches[branch];
    else{
      return setNode(branch,new NoroCacheNode());
    }
  }
};
class DenseRow{
public:
  number* array;
  int begin;
  int end;
  DenseRow(){
    array=NULL;
  }
  ~DenseRow(){
    omfree(array);
  }
};
template <class number_type> class SparseRow{
public:
  int* idx_array;
  number_type* coef_array;
  int len;
  SparseRow(){
    len=0;
    idx_array=NULL;
    coef_array=NULL;
  }
  SparseRow<number_type>(int n){
    len=n;
    idx_array=(int*) omalloc(n*sizeof(int));
    coef_array=(number_type*) omalloc(n*sizeof(number_type));
  }
  SparseRow<number_type>(int n, const number_type* source){
    len=n;
    idx_array=NULL;
    coef_array=(number_type*) omalloc(n*sizeof(number_type));
    memcpy(coef_array,source,n*sizeof(number_type));
  }
  ~SparseRow<number_type>(){
    omfree(idx_array);
    omfree(coef_array);
  }
};

template <class number_type> class DataNoroCacheNode:public NoroCacheNode{
public:

  int value_len;
  poly value_poly;
  #ifdef NORO_SPARSE_ROWS_PRE
  SparseRow<number_type>* row;
  #else
  DenseRow* row;
  #endif
  int term_index;
  DataNoroCacheNode(poly p, int len){
    value_len=len;
    value_poly=p;
    row=NULL;
    term_index=-1;
  }
  #ifdef NORO_SPARSE_ROWS_PRE
  DataNoroCacheNode(SparseRow<number_type>* row){
    if (row!=NULL)
      value_len=row->len;
    else
      value_len=0;
    value_poly=NULL;
    this->row=row;
    term_index=-1;
  }
  #endif
  ~DataNoroCacheNode(
  ){
    //p_Delete(&value_poly,currRing);
    if (row) delete row;
  }
};
template <class number_type> class TermNoroDataNode{
public:
  DataNoroCacheNode<number_type>* node;
  poly t;
};

template <class number_type> class NoroCache{
public:
  poly temp_term;
#ifndef NORO_NON_POLY
  void evaluatePlaceHolder(number* row,std::vector<NoroPlaceHolder>& place_holders);
  void evaluateRows();
  void evaluateRows(int level, NoroCacheNode* node);
#endif
  void collectIrreducibleMonomials( std::vector<DataNoroCacheNode<number_type>* >& res);
  void collectIrreducibleMonomials(int level,  NoroCacheNode* node, std::vector<DataNoroCacheNode<number_type>* >& res);

#ifdef NORO_RED_ARRAY_RESERVER
  int reserved;
  poly* recursionPolyBuffer;
#endif
  static const int backLinkCode=-222;
  DataNoroCacheNode<number_type>* insert(poly term, poly nf, int len){
    //assume(impl.find(p_Copy(term,currRing))==impl.end());
    //assume(len==pLength(nf));
    assume(npIsOne(p_GetCoeff(term,currRing)));
    if (term==nf){
      term=p_Copy(term,currRing);

      ressources.push_back(term);
      nIrreducibleMonomials++;
      return treeInsertBackLink(term);

    } else{

      if (nf){
        //nf=p_Copy(nf,currRing);
        assume(p_LmCmp(nf,term,currRing)==-1);
        ressources.push_back(nf);
      }
      return treeInsert(term,nf,len);

    }

    //impl[term]=std::pair<PolySimple,int> (nf,len);
  }
  #ifdef NORO_SPARSE_ROWS_PRE
  DataNoroCacheNode<number_type>* insert(poly term, SparseRow<number_type>* srow){
    //assume(impl.find(p_Copy(term,currRing))==impl.end());
    //assume(len==pLength(nf));

      return treeInsert(term,srow);


    //impl[term]=std::pair<PolySimple,int> (nf,len);
  }
  #endif
  DataNoroCacheNode<number_type>* insertAndTransferOwnerShip(poly t, ring r){

    ressources.push_back(t);
    DataNoroCacheNode<number_type>* res=treeInsertBackLink(t);
    res->term_index=nIrreducibleMonomials;
    nIrreducibleMonomials++;
    return res;
  }
  poly lookup(poly term, BOOLEAN& succ, int & len);
  DataNoroCacheNode<number_type>* getCacheReference(poly term);
  NoroCache(){
    buffer=NULL;
#ifdef NORO_RED_ARRAY_RESERVER
    reserved=0;
    recursionPolyBuffer=(poly*)omalloc(1000000*sizeof(poly));
#endif
    nIrreducibleMonomials=0;
    nReducibleMonomials=0;
    temp_term=pOne();
    tempBufferSize=3000;
    tempBuffer=omalloc(tempBufferSize);
  }
  void ensureTempBufferSize(size_t size){
    if (tempBufferSize<size){
      tempBufferSize=2*size;
      omfree(tempBuffer);
      tempBuffer=omalloc(tempBufferSize);
    }
  }
#ifdef NORO_RED_ARRAY_RESERVER
  poly* reserve(int n){
    poly* res=recursionPolyBuffer+reserved;
    reserved+=n;
    return res;
  }
  void free(int n){
    reserved-=n;
  }
#endif
  ~NoroCache(){
    int s=ressources.size();
    int i;
    for(i=0;i<s;i++){
      p_Delete(&ressources[i].impl,currRing);
    }
    p_Delete(&temp_term,currRing);
#ifdef NORO_RED_ARRAY_RESERVER
    omfree(recursionPolyBuffer);
#endif
   omfree(tempBuffer);
  }

  int nIrreducibleMonomials;
  int nReducibleMonomials;
  void* tempBuffer;
  size_t tempBufferSize;
protected:
  DataNoroCacheNode<number_type>* treeInsert(poly term,poly nf,int len){
    int i;
    nReducibleMonomials++;
    int nvars=pVariables;
    NoroCacheNode* parent=&root;
    for(i=1;i<nvars;i++){
      parent=parent->getOrInsertBranch(p_GetExp(term,i,currRing));
    }
    return (DataNoroCacheNode<number_type>*) parent->setNode(p_GetExp(term,nvars,currRing),new DataNoroCacheNode<number_type>(nf,len));
  }
  #ifdef NORO_SPARSE_ROWS_PRE
  DataNoroCacheNode<number_type>* treeInsert(poly term,SparseRow<number_type>* srow){
    int i;
    nReducibleMonomials++;
    int nvars=pVariables;
    NoroCacheNode* parent=&root;
    for(i=1;i<nvars;i++){
      parent=parent->getOrInsertBranch(p_GetExp(term,i,currRing));
    }
    return (DataNoroCacheNode<number_type>*) parent->setNode(p_GetExp(term,nvars,currRing),new DataNoroCacheNode<number_type>(srow));
  }
  #endif
  DataNoroCacheNode<number_type>* treeInsertBackLink(poly term){
    int i;
    int nvars=pVariables;
    NoroCacheNode* parent=&root;
    for(i=1;i<nvars;i++){
      parent=parent->getOrInsertBranch(p_GetExp(term,i,currRing));
    }
    return (DataNoroCacheNode<number_type>*) parent->setNode(p_GetExp(term,nvars,currRing),new DataNoroCacheNode<number_type>(term,backLinkCode));
  }

  //@TODO descruct nodes;
  typedef std::vector<PolySimple> poly_vec;
  poly_vec ressources;
  //typedef std::map<PolySimple,std::pair<PolySimple,int> > cache_map;
  //cache_map impl;
  NoroCacheNode root;
  number* buffer;
};
template<class number_type> SparseRow<number_type> * noro_red_to_non_poly_t(poly p, int &len, NoroCache<number_type>* cache,slimgb_alg* c);
template<class number_type> MonRedResNP<number_type> noro_red_mon_to_non_poly(poly t,  NoroCache<number_type> * cache,slimgb_alg* c)
{
  MonRedResNP<number_type> res_holder;


    DataNoroCacheNode<number_type>* ref=cache->getCacheReference(t);
    if (ref!=NULL){


      res_holder.coef=p_GetCoeff(t,c->r);

      res_holder.ref=ref;
      p_Delete(&t,c->r);
      return res_holder;
    }

  unsigned long sev=p_GetShortExpVector(t,currRing);
  int i=kFindDivisibleByInS_easy(c->strat,t,sev);
  if (i>=0){
    number coef_bak=p_GetCoeff(t,c->r);

    p_SetCoeff(t,npInit(1),c->r);
    assume(npIsOne(p_GetCoeff(c->strat->S[i],c->r)));
    number coefstrat=p_GetCoeff(c->strat->S[i],c->r);


    poly exp_diff=cache->temp_term;
    p_ExpVectorDiff(exp_diff,t,c->strat->S[i],c->r);
    p_SetCoeff(exp_diff,npNeg(npInvers(coefstrat)),c->r);
    p_Setm(exp_diff,c->r);
    assume(c->strat->S[i]!=NULL);

    poly res;
    res=pp_Mult_mm(pNext(c->strat->S[i]),exp_diff,c->r);

    int len=c->strat->lenS[i]-1;
    SparseRow<number_type>* srow;
    srow=noro_red_to_non_poly_t<number_type>(res,len,cache,c);
    ref=cache->insert(t,srow);
    p_Delete(&t,c->r);


    res_holder.coef=coef_bak;
    res_holder.ref=ref;
    return res_holder;

  } else {
    number coef_bak=p_GetCoeff(t,c->r);
    number one=npInit(1);
    p_SetCoeff(t,one,c->r);

    res_holder.ref=cache->insertAndTransferOwnerShip(t,c->r);
    assume(res_holder.ref!=NULL);
    res_holder.coef=coef_bak;

    return res_holder;

  }

}
/*
poly tree_add(poly* a,int begin, int end,ring r){
  int d=end-begin;
  switch(d){
    case 0:
      return NULL;
    case 1:
      return a[begin];
    case 2:
      return p_Add_q(a[begin],a[begin+1],r);
    default:
      int s=d/2;
      return p_Add_q(tree_add(a,begin,begin+s,r),tree_add(a,begin+s,end,r),r);
  }
}
*/
#ifdef __GNUC__
#define LIKELY(expression) (__builtin_expect(!!(expression), 1))
#define UNLIKELY(expression) (__builtin_expect(!!(expression), 0))
#else
#define LIKELY(expression) (expression)
#define UNLIKELY(expression) (expression)
#endif

template<class number_type> SparseRow<number_type>* convert_to_sparse_row(number_type* temp_array,int temp_size,int non_zeros){
SparseRow<number_type>* res=new SparseRow<number_type>(non_zeros);
//int pos=0;
//Print("denseness:%f\n",((double) non_zeros/(double) temp_size));
number_type* it_coef=res->coef_array;
int* it_idx=res->idx_array;
#if 0
for(i=0;i<cache->nIrreducibleMonomials;i++){
  if (!(0==temp_array[i])){

    res->idx_array[pos]=i;
    res->coef_array[pos]=temp_array[i];

    pos++;
    non_zeros--;
    if (non_zeros==0) break;
  }

}
#else
int64* start=(int64*) ((void*)temp_array);
int64* end;
const int multiple=sizeof(int64)/sizeof(number_type);
if (temp_size==0) end=start;

else
{
  int temp_size_rounded=temp_size+(multiple-(temp_size%multiple));
  assume(temp_size_rounded>=temp_size);
  assume(temp_size_rounded%multiple==0);
  assume(temp_size_rounded<temp_size+multiple);
  number_type* nt_end=temp_array+temp_size_rounded;
  end=(int64*)((void*)nt_end);
}
int64* it=start;
while(it!=end){
  if UNLIKELY((*it)!=0){
    int small_i;
    const int temp_index=((number_type*)((void*) it))-temp_array;
    const int bound=temp_index+multiple;
    number_type c;
    for(small_i=temp_index;small_i<bound;small_i++){
      if((c=temp_array[small_i])!=0){
        //res->idx_array[pos]=small_i;
        //res->coef_array[pos]=temp_array[small_i];
        (*(it_idx++))=small_i;
        (*(it_coef++))=c;
        //pos++;
        non_zeros--;

      }
      if UNLIKELY(non_zeros==0) break;
    }

  }
  ++it;
}
#endif
return res;
}
template <class number_type> void add_coef_times_sparse(number_type* const temp_array,
int temp_size,SparseRow<number_type>* row, number coef){
  int j;
  number_type* const coef_array=row->coef_array;
  int* const idx_array=row->idx_array;
  const int len=row->len;
  tgb_uint32 buffer[256];
  const tgb_uint32 prime=npPrimeM;
  const tgb_uint32 c=F4mat_to_number_type(coef);
  assume(!(npIsZero(coef)));
  for(j=0;j<len;j=j+256){
    const int bound=std::min(j+256,len);
    int i;
    int bpos=0;
    for(i=j;i<bound;i++){
      buffer[bpos++]=coef_array[i];
    }
    int bpos_bound=bound-j;
    for(i=0;i<bpos_bound;i++){
       buffer[i]*=c;
     }
    for(i=0;i<bpos_bound;i++){
       buffer[i]=buffer[i]%prime;
    }
    bpos=0;
    for(i=j;i<bound;i++){
      int idx=idx_array[i];
      assume(bpos<256);
      assume(!(npIsZero((number) buffer[bpos])));
      temp_array[idx]=F4mat_to_number_type(npAddM((number) temp_array[idx], (number) buffer[bpos++]));
      assume(idx<temp_size);
    }

  }
}
template <class number_type> void add_coef_times_dense(number_type* const temp_array,
int temp_size,const number_type* row, int len,number coef){
  int j;
  const number_type* const coef_array=row;
  //int* const idx_array=row->idx_array;
  //const int len=temp_size;
  tgb_uint32 buffer[256];
  const tgb_uint32 prime=npPrimeM;
  const tgb_uint32 c=F4mat_to_number_type(coef);
  assume(!(npIsZero(coef)));
  for(j=0;j<len;j=j+256){
    const int bound=std::min(j+256,len);
    int i;
    int bpos=0;
    for(i=j;i<bound;i++){
      buffer[bpos++]=coef_array[i];
    }
    int bpos_bound=bound-j;
    for(i=0;i<bpos_bound;i++){
       buffer[i]*=c;
     }
    for(i=0;i<bpos_bound;i++){
       buffer[i]=buffer[i]%prime;
    }
    bpos=0;
    for(i=j;i<bound;i++){
      //int idx=idx_array[i];
      assume(bpos<256);
      //assume(!(npIsZero((number) buffer[bpos])));
      temp_array[i]=F4mat_to_number_type(npAddM((number) temp_array[i], (number) buffer[bpos++]));
      assume(i<temp_size);
    }

  }
}
template <class number_type> void add_dense(number_type* const temp_array,
int temp_size,const number_type* row, int len){
  //int j;
  const number_type* const coef_array=row;
  //int* const idx_array=row->idx_array;
  //const int len=temp_size;
  tgb_uint32 buffer[256];
  const tgb_uint32 prime=npPrimeM;
  //const tgb_uint32 c=F4mat_to_number_type(coef);

  int i;
  for(i=0;i<len;i++){

      temp_array[i]=F4mat_to_number_type(npAddM((number) temp_array[i], (number) row[i]));
      assume(i<temp_size);
    }

  }
template <class number_type> void sub_dense(number_type* const temp_array,
int temp_size,const number_type* row, int len){
  //int j;
  const number_type* const coef_array=row;
  //int* const idx_array=row->idx_array;
  //const int len=temp_size;
  tgb_uint32 buffer[256];
  const tgb_uint32 prime=npPrimeM;
  //const tgb_uint32 c=F4mat_to_number_type(coef);

  int i;
  for(i=0;i<len;i++){

      temp_array[i]=F4mat_to_number_type(npSubM((number) temp_array[i], (number) row[i]));
      assume(i<temp_size);
    }

  }

template <class number_type> void add_sparse(number_type* const temp_array,int temp_size,SparseRow<number_type>* row){
  int j;

          number_type* const coef_array=row->coef_array;
          int* const idx_array=row->idx_array;
          const int len=row->len;
        for(j=0;j<len;j++){
          int idx=idx_array[j];
          temp_array[idx]=F4mat_to_number_type(   npAddM((number) temp_array[idx],(number) coef_array[j]));
          assume(idx<temp_size);
        }
}
template <class number_type> void sub_sparse(number_type* const temp_array,int temp_size,SparseRow<number_type>* row){
  int j;

          number_type* const coef_array=row->coef_array;
          int* const idx_array=row->idx_array;
          const int len=row->len;
        for(j=0;j<len;j++){
          int idx=idx_array[j];
          temp_array[idx]=F4mat_to_number_type(   npSubM((number) temp_array[idx],(number) coef_array[j]));
          assume(idx<temp_size);
        }
}
template <class number_type> SparseRow<number_type>* noro_red_to_non_poly_dense(MonRedResNP<number_type>* mon, int len,NoroCache<number_type>* cache){
  size_t temp_size_bytes=cache->nIrreducibleMonomials*sizeof(number_type)+8;//use 8bit int for testing
   assume(sizeof(int64)==8);
   cache->ensureTempBufferSize(temp_size_bytes);
   number_type* temp_array=(number_type*) cache->tempBuffer;//omalloc(cache->nIrreducibleMonomials*sizeof(number_type));
   int temp_size=cache->nIrreducibleMonomials;
   memset(temp_array,0,temp_size_bytes);
   number minus_one=npInit(-1);
   int i;
   for(i=0;i<len;i++){
     MonRedResNP<number_type> red=mon[i];
     if ((red.ref)){
       if (red.ref->row){
         SparseRow<number_type>* row=red.ref->row;
         number coef=red.coef;
         int j;
         if (row->idx_array){
         if (!((coef==(number) 1)||(coef==minus_one))){
           add_coef_times_sparse(temp_array,temp_size,row,coef);



         }else{
           if (coef==(number) 1){
              add_sparse(temp_array,temp_size,row);
           } else {

             sub_sparse(temp_array,temp_size,row);
           }
         }
       } else
       //TODO: treat, 1,-1
       if (!((coef==(number) 1)||(coef==minus_one))){
         add_coef_times_dense(temp_array,temp_size,row->coef_array,row->len,coef);
       } else {
         if (coef==(number)1)
           add_dense(temp_array,temp_size,row->coef_array,row->len);
         else{
           assume(coef==minus_one);
           sub_dense(temp_array,temp_size,row->coef_array,row->len);
           //add_coef_times_dense(temp_array,temp_size,row->coef_array,row->len,coef);
         }
       }
       }

       else{
         if (red.ref->value_len==NoroCache<number_type>::backLinkCode){
           temp_array[red.ref->term_index]=F4mat_to_number_type( npAddM((number) temp_array[red.ref->term_index],red.coef));
         } else {
           //PrintS("third case\n");
         }
       }
     }
   }
   int non_zeros=0;
   for(i=0;i<cache->nIrreducibleMonomials;i++){
     //if (!(temp_array[i]==0)){
     //  non_zeros++;
     //}
     assume(((temp_array[i]!=0)==0)|| (((temp_array[i]!=0)==1)));
     non_zeros+=(temp_array[i]!=0);
   }

   if (non_zeros==0){
     //omfree(mon);
     return NULL;
   }
   SparseRow<number_type>* res=new SparseRow<number_type>(temp_size,temp_array);//convert_to_sparse_row(temp_array,temp_size, non_zeros);

   //omfree(temp_array);


   return res;
}
template<class number_type> class CoefIdx{
public:
  number_type coef;
  int idx;
  bool operator<(const CoefIdx<number_type>& other) const{
    return (idx<other.idx);
  }
};
template<class number_type> void write_coef_times_xx_idx_to_buffer(CoefIdx<number_type>* const pairs,int& pos,int* const idx_array, number_type* const coef_array,const int rlen, const number coef){
  int j;
  for(j=0;j<rlen;j++){
    assume(coef_array[j]!=0);
    CoefIdx<number_type> ci;
    ci.coef=F4mat_to_number_type(npMultM((number) coef,(number) coef_array[j]));
    ci.idx=idx_array[j];
    pairs[pos++]=ci;
  }
}
template<class number_type> void write_coef_times_xx_idx_to_buffer_dense(CoefIdx<number_type>* const pairs,int& pos, number_type* const coef_array,const int rlen, const number coef){
  int j;

  for(j=0;j<rlen;j++){
    if (coef_array[j]!=0){
    assume(coef_array[j]!=0);
    CoefIdx<number_type> ci;
    ci.coef=F4mat_to_number_type(npMultM((number) coef,(number) coef_array[j]));
    assume(ci.coef!=0);
    ci.idx=j;
    pairs[pos++]=ci;
  }
  }
}
template<class number_type> void write_coef_idx_to_buffer_dense(CoefIdx<number_type>* const pairs,int& pos, number_type* const coef_array,const int rlen){
  int j;

  for(j=0;j<rlen;j++){
    if (coef_array[j]!=0){
    assume(coef_array[j]!=0);
    CoefIdx<number_type> ci;
    ci.coef=coef_array[j];
    assume(ci.coef!=0);
    ci.idx=j;
    pairs[pos++]=ci;
  }
  }
}

template<class number_type> void write_minus_coef_idx_to_buffer_dense(CoefIdx<number_type>* const pairs,int& pos, number_type* const coef_array,const int rlen){
  int j;

  for(j=0;j<rlen;j++){
    if (coef_array[j]!=0){
    assume(coef_array[j]!=0);
    CoefIdx<number_type> ci;
    ci.coef=F4mat_to_number_type(npNegM((number) coef_array[j]));
    assume(ci.coef!=0);
    ci.idx=j;
    pairs[pos++]=ci;
  }
  }
}
template<class number_type> void write_coef_idx_to_buffer(CoefIdx<number_type>* const pairs,int& pos,int* const idx_array, number_type* const coef_array,const int rlen){
  int j;
  for(j=0;j<rlen;j++){
    assume(coef_array[j]!=0);
    CoefIdx<number_type> ci;
    ci.coef=coef_array[j];
    ci.idx=idx_array[j];
    pairs[pos++]=ci;
  }
}

template<class number_type> void write_minus_coef_idx_to_buffer(CoefIdx<number_type>* const pairs,int& pos,int* const idx_array, number_type* const coef_array,const int rlen){
  int j;
  for(j=0;j<rlen;j++){
    assume(coef_array[j]!=0);
    CoefIdx<number_type> ci;
    ci.coef=F4mat_to_number_type(npNegM(coef_array[j]));
    ci.idx=idx_array[j];
    pairs[pos++]=ci;
  }
}
template <class number_type> SparseRow<number_type>* noro_red_to_non_poly_sparse(MonRedResNP<number_type>* mon, int len,NoroCache<number_type>* cache){
  int i;
int together=0;
for(i=0;i<len;i++){
  MonRedResNP<number_type> red=mon[i];
  if ((red.ref) &&( red.ref->row)){
    together+=red.ref->row->len;
    } else {
      if ((red.ref) &&(red.ref->value_len==NoroCache<number_type>::backLinkCode))
      together++;
  }

}
  //PrintS("here\n");
if (together==0) return 0;
  //PrintS("there\n");
cache->ensureTempBufferSize(together*sizeof(CoefIdx<number_type>));
CoefIdx<number_type>* pairs=(CoefIdx<number_type>*) cache->tempBuffer; //omalloc(together*sizeof(CoefIdx<number_type>));
int pos=0;
int j;
const number one=npInit(1);
const number minus_one=npInit(-1);
for(i=0;i<len;i++){
  MonRedResNP<number_type> red=mon[i];
  if ((red.ref) &&( red.ref->row)){
      //together+=red.ref->row->len;
    int* idx_array=red.ref->row->idx_array;
    number_type* coef_array=red.ref->row->coef_array;
    int rlen=red.ref->row->len;
    number coef=red.coef;
    if (idx_array){
      if ((coef!=one)&&(coef!=minus_one)){
        write_coef_times_xx_idx_to_buffer(pairs,pos,idx_array, coef_array,rlen, coef);
        } else
        {
          if (coef==one){
            write_coef_idx_to_buffer(pairs,pos,idx_array, coef_array,rlen);
          } else {
            assume(coef==minus_one);
            write_minus_coef_idx_to_buffer(pairs,pos,idx_array, coef_array,rlen);
          }
        }
      } else{
        if ((coef!=one)&&(coef!=minus_one)){
          write_coef_times_xx_idx_to_buffer_dense(pairs,pos,coef_array,rlen,coef);
        } else{
          if (coef==one)
            write_coef_idx_to_buffer_dense(pairs,pos,coef_array,rlen);
          else {
            assume(coef==minus_one);
            write_minus_coef_idx_to_buffer_dense(pairs,pos,coef_array,rlen);
          }
        }
      }
    }
    else {
      if ((red.ref) &&(red.ref->value_len==NoroCache<number_type>::backLinkCode)){
        CoefIdx<number_type> ci;
        ci.coef=F4mat_to_number_type(red.coef);
        ci.idx=red.ref->term_index;
        pairs[pos++]=ci;
      }
    }
  }
  assume(pos<=together);
  together=pos;

  std::sort(pairs,pairs+together);

  int act=0;

  assume(pairs[0].coef!=0);
  for(i=1;i<together;i++){
    if (pairs[i].idx!=pairs[act].idx){
      if (pairs[act].coef!=0){
        act=act+1;
      }
      pairs[act]=pairs[i];
    } else{
      pairs[act].coef=F4mat_to_number_type(npAddM((number)pairs[act].coef,(number)pairs[i].coef));
    }
  }

  if (pairs[act].coef==0){

    act--;
  }
  int sparse_row_len=act+1;
  //Print("res len:%d",sparse_row_len);
  if (sparse_row_len==0) {return NULL;}
  SparseRow<number_type>* res=new SparseRow<number_type>(sparse_row_len);
  {
    number_type* coef_array=res->coef_array;
    int* idx_array=res->idx_array;
    for(i=0;i<sparse_row_len;i++){
      idx_array[i]=pairs[i].idx;
      coef_array[i]=pairs[i].coef;
    }
  }
  //omfree(pairs);

  return res;
}
template<class number_type> SparseRow<number_type> * noro_red_to_non_poly_t(poly p, int &len, NoroCache<number_type>* cache,slimgb_alg* c){
  assume(len==pLength(p));
  poly orig_p=p;
  if (p==NULL) {
    len=0;
    return NULL;
  }

  number zero=npInit(0);
  MonRedResNP<number_type>* mon=(MonRedResNP<number_type>*) omalloc(len*sizeof(MonRedResNP<number_type>));
  int i=0;
  double max_density=0.0;
  while(p){

    poly t=p;
    pIter(p);
    pNext(t)=NULL;

#ifndef NDEBUG
    number coef_debug=p_GetCoeff(t,currRing);
#endif
    MonRedResNP<number_type> red=noro_red_mon_to_non_poly(t,cache,c);
    if ((red.ref) && (red.ref->row)){
      double act_density=(double) red.ref->row->len;
      act_density/=(double) cache->nIrreducibleMonomials;
      max_density=std::max(act_density,max_density);
    }
    mon[i]=red;
    i++;
  }

  assume(i==len);
  len=i;
  bool dense=true;
  if (max_density<0.3) dense=false;
  if (dense){
    SparseRow<number_type>* res=noro_red_to_non_poly_dense(mon,len,cache);
    omfree(mon);
    return res;
  } else   {
      SparseRow<number_type>* res=noro_red_to_non_poly_sparse(mon,len,cache);
      omfree(mon);
      return res;
    }
  //in the loop before nIrreducibleMonomials increases, so position here is important

}
#endif
static wlen_type pair_weighted_length(int i, int j, slimgb_alg* c);
wlen_type pELength(poly p, ring r);
int terms_sort_crit(const void* a, const void* b);
//void simplest_gauss_modp(number* a, int nrows,int ncols);
// a: a[0,0],a[0,1]....a[nrows-1,ncols-1]
// assume: field is Zp
#ifdef USE_NORO


template <class number_type > void write_poly_to_row(number_type* row, poly h, poly*terms, int tn, ring r){
  //poly* base=row;
  while(h!=NULL){
    //Print("h:%i\n",h);
    number coef=p_GetCoeff(h,r);
    poly* ptr_to_h=(poly*) bsearch(&h,terms,tn,sizeof(poly),terms_sort_crit);
    assume(ptr_to_h!=NULL);
    int pos=ptr_to_h-terms;
    row[pos]=F4mat_to_number_type(coef);
    //number_type_array[base+pos]=coef;
    pIter(h);
  }
}
template <class number_type > poly row_to_poly(number_type* row, poly* terms, int tn, ring r){
  poly h=NULL;
  int j;
  number_type zero=0;//;npInit(0);
  for(j=tn-1;j>=0;j--){
    if (!(zero==(row[j]))){
      poly t=terms[j];
      t=p_LmInit(t,r);
      p_SetCoeff(t,(number) row[j],r);
      pNext(t)=h;
      h=t;
    }

  }
  return h;
}
template <class number_type > int modP_lastIndexRow(number_type* row,int ncols){
  int lastIndex;
  const number_type zero=0;//npInit(0);
  for(lastIndex=ncols-1;lastIndex>=0;lastIndex--){
    if (!(row[lastIndex]==zero)){
      return lastIndex;
    }
  }
  return -1;
}
template <class number_type> int term_nodes_sort_crit(const void* a, const void* b){
  return -pLmCmp(((TermNoroDataNode<number_type>*) a)->t,((TermNoroDataNode<number_type>*) b)->t);
}

template <class number_type>class ModPMatrixBackSubstProxyOnArray;
template <class number_type > class ModPMatrixProxyOnArray{
public:
  friend class ModPMatrixBackSubstProxyOnArray<number_type>;

  int ncols,nrows;
  ModPMatrixProxyOnArray(number_type* array, int nrows, int ncols){
    this->ncols=ncols;
    this->nrows=nrows;
    rows=(number_type**) omalloc(nrows*sizeof(number_type*));
    startIndices=(int*)omalloc(nrows*sizeof(int));
    int i;
    for(i=0;i<nrows;i++){
      rows[i]=array+(i*ncols);
      updateStartIndex(i,-1);
    }
  }
  ~ModPMatrixProxyOnArray(){
    omfree(rows);
    omfree(startIndices);
  }

  void permRows(int i, int j){
    number_type* h=rows[i];
    rows[i]=rows[j];
    rows[j]=h;
    int hs=startIndices[i];
    startIndices[i]=startIndices[j];
    startIndices[j]=hs;
  }
  void multiplyRow(int row, number_type coef){
    int i;
    number_type* row_array=rows[row];
    for(i=startIndices[row];i<ncols;i++){
      row_array[i]=F4mat_to_number_type(npMult((number) row_array[i],(number) coef));
    }
  }
  void reduceOtherRowsForward(int r){

    //assume rows "under r" have bigger or equal start index
    number_type* row_array=rows[r];
    number_type zero=F4mat_to_number_type(npInit(0));
    int start=startIndices[r];
    number_type coef=row_array[start];
    assume(start<ncols);
    int other_row;
    assume(!(npIsZero((number) row_array[start])));
    if (!(npIsOne((number) coef)))
      multiplyRow(r,F4mat_to_number_type(npInvers((number) coef)));
    assume(npIsOne((number) row_array[start]));
    int lastIndex=modP_lastIndexRow(row_array, ncols);
    number minus_one=npInit(-1);
    for (other_row=r+1;other_row<nrows;other_row++){
      assume(startIndices[other_row]>=start);
      if (startIndices[other_row]==start){
        int i;
        number_type* other_row_array=rows[other_row];
        number coef2=npNeg((number) other_row_array[start]);
        if (coef2==minus_one){
          for(i=start;i<=lastIndex;i++){
            if (row_array[i]!=zero)
              other_row_array[i]=F4mat_to_number_type(npSubM((number) other_row_array[i], (number) row_array[i]));
          }
      }else {
          //assume(FALSE);
          for(i=start;i<=lastIndex;i++){
            if (row_array[i]!=zero)
            other_row_array[i]=F4mat_to_number_type(npAddM(npMult(coef2,(number) row_array[i]),(number) other_row_array[i]));
          }
        }
        updateStartIndex(other_row,start);
        assume(npIsZero((number) other_row_array[start]));
      }
    }
  }
  void updateStartIndex(int row,int lower_bound){
    number_type* row_array=rows[row];
    assume((lower_bound<0)||(npIsZero((number) row_array[lower_bound])));
    int i;
    //number_type zero=npInit(0);
    for(i=lower_bound+1;i<ncols;i++){
      if (!(row_array[i]==0))
        break;
    }
    startIndices[row]=i;
  }
  int getStartIndex(int row){
    return startIndices[row];
  }
  BOOLEAN findPivot(int &r, int &c){
    //row>=r, col>=c

    while(c<ncols){
      int i;
      for(i=r;i<nrows;i++){
        assume(startIndices[i]>=c);
        if (startIndices[i]==c){
          //r=i;
          if (r!=i)
            permRows(r,i);
          return TRUE;
        }
      }
      c++;
    }
    return FALSE;
  }
protected:
  number_type** rows;
  int* startIndices;
};
template <class number_type > class ModPMatrixBackSubstProxyOnArray{
  int *startIndices;
  number_type** rows;
  int *lastReducibleIndices;
  int ncols;
  int nrows;
  int nonZeroUntil;
public:
  void multiplyRow(int row, number_type coef){
    int i;
    number_type* row_array=rows[row];
    for(i=startIndices[row];i<ncols;i++){
      row_array[i]=F4mat_to_number_type(npMult((number) row_array[i],(number) coef));
    }
  }
  ModPMatrixBackSubstProxyOnArray<number_type> (ModPMatrixProxyOnArray<number_type> & p){
//  (number_type* array, int nrows, int ncols, int* startIndices, number_type** rows){
    //we borrow some parameters ;-)
    //we assume, that nobody changes the order of the rows
    this->startIndices=p.startIndices;
    this->rows=p.rows;
    this->ncols=p.ncols;
    this->nrows=p.nrows;
    lastReducibleIndices=(int*) omalloc(nrows*sizeof(int));
    nonZeroUntil=0;
    while(nonZeroUntil<nrows){
      if (startIndices[nonZeroUntil]<ncols){

        nonZeroUntil++;
      } else break;

    }
    if (TEST_OPT_PROT)
      Print("rank:%i\n",nonZeroUntil);
    nonZeroUntil--;
    int i;
    for(i=0;i<=nonZeroUntil;i++){
      assume(startIndices[i]<ncols);
      assume(!(npIsZero((number) rows[i][startIndices[i]])));
      assume(startIndices[i]>=i);
      updateLastReducibleIndex(i,nonZeroUntil+1);
    }
  }
  void updateLastReducibleIndex(int r, int upper_bound){
    number_type* row_array=rows[r];
    if (upper_bound>nonZeroUntil) upper_bound=nonZeroUntil+1;
    int i;
    const number_type zero=0;//npInit(0);
    for(i=upper_bound-1;i>r;i--){
      int start=startIndices[i];
      assume(start<ncols);
      if (!(row_array[start]==zero)){
        lastReducibleIndices[r]=start;
        return;
      }
    }
    lastReducibleIndices[r]=-1;
  }
  void backwardSubstitute(int r){
    int start=startIndices[r];
    assume(start<ncols);
    number_type zero=0;//npInit(0);
    number_type* row_array=rows[r];
    assume((!(npIsZero((number) row_array[start]))));
    assume(start<ncols);
    int other_row;
    if (!(npIsOne((number) row_array[r]))){
      //it should be one, but this safety is not expensive
      multiplyRow(r, F4mat_to_number_type(npInvers((number) row_array[start])));
    }
    int lastIndex=modP_lastIndexRow(row_array, ncols);
    assume(lastIndex<ncols);
    assume(lastIndex>=0);
    for(other_row=r-1;other_row>=0;other_row--){
      assume(lastReducibleIndices[other_row]<=start);
      if (lastReducibleIndices[other_row]==start){
        number_type* other_row_array=rows[other_row];
        number coef=npNeg((number) other_row_array[start]);
        assume(!(npIsZero(coef)));
        int i;
        assume(start>startIndices[other_row]);
        for(i=start;i<=lastIndex;i++){
          if (row_array[i]!=zero)
            other_row_array[i]=F4mat_to_number_type(npAddM(npMult(coef,(number)row_array[i]),(number)other_row_array[i]));
        }
        updateLastReducibleIndex(other_row,r);
      }
    }
  }
  ~ModPMatrixBackSubstProxyOnArray<number_type>(){
    omfree(lastReducibleIndices);
  }
  void backwardSubstitute(){
    int i;
    for(i=nonZeroUntil;i>0;i--){
      backwardSubstitute(i);
    }
  }
};
template <class number_type > void simplest_gauss_modp(number_type* a, int nrows,int ncols){
  //use memmoves for changing rows
  if (TEST_OPT_PROT)
    PrintS("StartGauss\n");
  ModPMatrixProxyOnArray<number_type> mat(a,nrows,ncols);

  int c=0;
  int r=0;
  while(mat.findPivot(r,c)){
    //int pivot=find_pivot()
      mat.reduceOtherRowsForward(r);
    r++;
    c++;
  }
  ModPMatrixBackSubstProxyOnArray<number_type> backmat(mat);
  backmat.backwardSubstitute();
  //backward substitutions
  if (TEST_OPT_PROT)
    PrintS("StopGauss\n");
}
//int term_nodes_sort_crit(const void* a, const void* b);
template <class number_type> void noro_step(poly*p,int &pn,slimgb_alg* c){
  //Print("Input rows %d\n",pn);
  int j;
  if (TEST_OPT_PROT){
    Print("Input rows %d\n",pn);
  }

  NoroCache<number_type> cache;

  SparseRow<number_type> ** srows=(SparseRow<number_type>**) omalloc(pn*sizeof(SparseRow<number_type>*));
  int non_zeros=0;
  for(j=0;j<pn;j++){

    poly h=p[j];
    int h_len=pLength(h);

    //number coef;


    srows[non_zeros]=noro_red_to_non_poly_t<number_type>(h,h_len,&cache,c);
    if (srows[non_zeros]!=NULL) non_zeros++;
  }
  std::vector<DataNoroCacheNode<number_type>*> irr_nodes;
  cache.collectIrreducibleMonomials(irr_nodes);
  //now can build up terms array
  //Print("historic irred Mon%d\n",cache.nIrreducibleMonomials);
  int n=irr_nodes.size();//cache.countIrreducibleMonomials();
  cache.nIrreducibleMonomials=n;
  if (TEST_OPT_PROT){
    Print("Irred Mon:%d\n",n);
    Print("red Mon:%d\n",cache.nReducibleMonomials);
  }
  TermNoroDataNode<number_type>* term_nodes=(TermNoroDataNode<number_type>*) omalloc(n*sizeof(TermNoroDataNode<number_type>));

  for(j=0;j<n;j++){
    assume(irr_nodes[j]!=NULL);
    assume(irr_nodes[j]->value_len==NoroCache<number_type>::backLinkCode);
    term_nodes[j].t=irr_nodes[j]->value_poly;
    assume(term_nodes[j].t!=NULL);
    term_nodes[j].node=irr_nodes[j];
  }


  qsort(term_nodes,n,sizeof(TermNoroDataNode<number_type>),term_nodes_sort_crit<number_type>);
  poly* terms=(poly*) omalloc(n*sizeof(poly));

  int* old_to_new_indices=(int*) omalloc(cache.nIrreducibleMonomials*sizeof(int));
  for(j=0;j<n;j++){
    old_to_new_indices[term_nodes[j].node->term_index]=j;
    term_nodes[j].node->term_index=j;
    terms[j]=term_nodes[j].t;
  }

  //if (TEST_OPT_PROT)
  //  Print("Evaluate Rows \n");
  pn=non_zeros;
  number_type* number_array=(number_type*) omalloc(n*pn*sizeof(number_type));
  memset(number_array,0,sizeof(number_type)*n*pn);
  number zero=npInit(0);

  for(j=0;j<pn;j++){
    int i;
    number_type* row=number_array+n*j;
    /*for(i=0;i<n;i++){
      row[i]=zero;
    }*/

    SparseRow<number_type>* srow=srows[j];

    if (srow){
      int* const idx_array=srow->idx_array;
      number_type* const coef_array=srow->coef_array;
      const int len=srow->len;
      if (srow->idx_array){
        for(i=0;i<len;i++){
         int idx=old_to_new_indices[idx_array[i]];
         row[idx]=F4mat_to_number_type(coef_array[i]);
        }
      }
      else {
        for(i=0;i<len;i++){
          row[old_to_new_indices[i]]=F4mat_to_number_type(coef_array[i]);
        }
      }
      delete srow;
    }
  }

  static int export_n=0;
  //export_mat(number_array,pn,n,"mat%i.py",++export_n);
  simplest_gauss_modp(number_array,pn,n);

  int p_pos=0;
  for(j=0;j<pn;j++){
    poly h=row_to_poly(number_array+j*n,terms,n,c->r);
    if(h!=NULL){
      p[p_pos++]=h;
    }
  }
  pn=p_pos;
  omfree(terms);
  omfree(term_nodes);
  omfree(number_array);
  #ifdef NORO_NON_POLY
  omfree(srows);
  omfree(old_to_new_indices);
  #endif
  //don't forget the rank

}

template <class number_type> void NoroCache<number_type>::collectIrreducibleMonomials( std::vector<DataNoroCacheNode<number_type> *>& res){
  int i;
  for(i=0;i<root.branches_len;i++){
    collectIrreducibleMonomials(1,root.branches[i],res);
  }
}
template <class number_type> void NoroCache<number_type>::collectIrreducibleMonomials(int level, NoroCacheNode* node, std::vector<DataNoroCacheNode<number_type>*>& res){
  assume(level>=0);
  if (node==NULL) return;
  if (level<pVariables){
    int i,sum;
    for(i=0;i<node->branches_len;i++){
      collectIrreducibleMonomials(level+1,node->branches[i],res);
    }
  } else {
    DataNoroCacheNode<number_type>* dn=(DataNoroCacheNode<number_type>*) node;
    if (dn->value_len==backLinkCode){
      res.push_back(dn);
    }
  }
}

template<class number_type> DataNoroCacheNode<number_type>* NoroCache<number_type>::getCacheReference(poly term){
  int i;
  NoroCacheNode* parent=&root;
  for(i=1;i<pVariables;i++){
    parent=parent->getBranch(p_GetExp(term,i,currRing));
    if (!(parent)){
      return NULL;
    }
  }
  DataNoroCacheNode<number_type>* res_holder=(DataNoroCacheNode<number_type>*) parent->getBranch(p_GetExp(term,i,currRing));
  return res_holder;
}
template<class number_type> poly NoroCache<number_type>::lookup(poly term, BOOLEAN& succ, int & len){
  int i;
  NoroCacheNode* parent=&root;
  for(i=1;i<pVariables;i++){
    parent=parent->getBranch(p_GetExp(term,i,currRing));
    if (!(parent)){
      succ=FALSE;
      return NULL;
    }
  }
  DataNoroCacheNode<number_type>* res_holder=(DataNoroCacheNode<number_type>*) parent->getBranch(p_GetExp(term,i,currRing));
  if (res_holder){
    succ=TRUE;
    if ((res_holder->value_len==backLinkCode)){
      len=1;
      return term;
    }
    len=res_holder->value_len;
    return res_holder->value_poly;
  } else {
    succ=FALSE;
    return NULL;
  }
}
#endif

#endif
