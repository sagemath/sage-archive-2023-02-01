#include "polynomial.h"

#include <sstream>

#include "printer.h"
#include "timer.h"
#include "linalg.h"

//static Timer polynomialTimer("Polynomial subtraction",10);
//static Timer polynomialTimer1("Polynomial monomial multiplication",1);

//-----------------------------------------
// Polynomial
//-----------------------------------------

Polynomial::Polynomial(const Term &t):
  marked(t),
  sugar(0),
  theRing(t.getRing())
{
  terms[t.m]=t.c;
  isMarkedBool=true;
}


Polynomial::Polynomial(PolynomialRing const &r):
  theRing(r),
  sugar(0),
  marked(r)
{
	isMarkedBool=false;
}

void Polynomial::mark(TermOrder const &termOrder)
{
  TermMap::iterator i=terms.begin();

  if(i!=terms.end())
    {
      Term best=Term(i->second,i->first);

      for(;i!=terms.end();i++)
	if(termOrder(best.m.exponent,i->first.exponent))best=Term(i->second,i->first);
      marked=best;
    }
  isMarkedBool=true;
}


void Polynomial::mark(Monomial const &monomial)
{
  assert(!terms.empty());
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(i->first.exponent==monomial.exponent)
      {
	marked=Term(i->second,monomial);
	isMarkedBool=true;
	return;
      }
  fprintf(Stderr,"Monomial ");
  AsciiPrinter(Stderr).printMonomial(monomial);
  fprintf(Stderr," not found in ");
  AsciiPrinter(Stderr).printPolynomial(*this);
  fprintf(Stderr,"\n");
  assert(0);
}


bool Polynomial::checkMarking(TermOrder const &termOrder)const
{
	if(!isMarked())return false;
	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      //      if(!termOrder(i->first.exponent,marked.m.exponent))return false;
      if(termOrder(marked.m.exponent,i->first.exponent))return false;
    }
  return true;
}


bool Polynomial::isHomogeneous(IntegerVector const &v)const
{
  return degree(v)==-degree(-v);
}


void Polynomial::scaleMarkedCoefficientToOne()
{
	assert(isMarked());
  FieldElement a=marked.c.inverse();
  Monomial k(theRing);
  k.exponent=IntegerVector(getNumberOfVariables());
  Term s(a,k);
  *this*=s;
  marked.c=marked.c.one();
}


int Polynomial::getNumberOfVariables()const
{
  TermMap::const_iterator i=terms.begin();
  if(i==terms.end())return 0;

  return i->first.exponent.size();
}


void Polynomial::changeNumberOfVariables(int n)
{
  PolynomialRing newRing=theRing;//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!111111111111!!!!!!
  Polynomial q(newRing);

  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      v.resize(n);
      FieldElement e=i->second;
      q+=Term(e,Monomial(theRing,v));
    }
  if(isMarked())marked.m.exponent.resize(n);
  terms=q.terms;
}


void Polynomial::madd(const Term &m, const Polynomial &p)
{
  if(p.terms.empty())return; //added May 7 2005
  int sugar2=p.getSugar()+m.m.exponent.sum();
  if(sugar2>sugar)sugar=sugar2;
  TermMap::iterator i=terms.lower_bound(Monomial(theRing,p.terms.begin()->first.exponent+m.m.exponent));

  for(TermMap::const_iterator j=p.terms.begin();j!=p.terms.end();j++)
    {
      while(i!=terms.end() && TermMapCompare()(i->first,Monomial(theRing,j->first.exponent+m.m.exponent)))i++;
      if(i==terms.end())
	{
	  terms.insert(i,TermMap::value_type(Monomial(theRing,j->first.exponent+m.m.exponent),j->second*m.c));
	}
      else
	{
	  if(!TermMapCompare()(Monomial(theRing,j->first.exponent+m.m.exponent),i->first))
	    { // they must be equal
	      FieldElement c=i->second+j->second*m.c;
	      if(c.isZero())
		{
		  TermMap::iterator oldI=i;
		  i++;
		  terms.erase(oldI);
		}
	      else
		{
		  i->second=c;
		}
	    }
	  else
	    {
	      terms.insert(i,TermMap::value_type(Monomial(theRing,j->first.exponent+m.m.exponent),j->second*m.c));
	    }
	}
    }
}

void Polynomial::operator+=(const Polynomial &p)
{
  if(p.terms.empty())return; //added May 7 2005
  if(p.getSugar()>sugar)sugar=p.getSugar();
  // fast addition
  //  TimerScope ts(&polynomialTimer);

  TermMap::iterator i=terms.lower_bound(p.terms.begin()->first);

  for(TermMap::const_iterator j=p.terms.begin();j!=p.terms.end();j++)
    {
      while(i!=terms.end() && (terms.value_comp()(*i,*j)))i++;
      if(i==terms.end())
	{
	  terms.insert(i,TermMap::value_type(j->first,j->second));
	}
      else
	{
	  if(!terms.value_comp()(*j,*i))
	    { // they must be equal
	      FieldElement c=i->second+j->second;
	      if(c.isZero())
		{
		  TermMap::iterator oldI=i;
		  i++;
		  terms.erase(oldI);
		}
	      else
		{
		  i->second=c;
		}
	    }
	  else
	    {
	      terms.insert(i,TermMap::value_type(j->first,j->second));
	    }
	}
    }
  // slow addition
  /*  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      if(terms.count(i->first)==1)
        {
          terms[i->first]=terms[i->first]+i->second;
          if(terms[i->first].isZero())terms.erase(i->first);
        }
      else
        terms[i->first]=i->second;
    }
  */
}

void Polynomial::operator-=(const Polynomial &p)
{
  if(p.terms.empty())return; //added May 7 2005
  if(p.getSugar()>sugar)sugar=p.getSugar();
  // fast subtraction
  //TimerScope ts(&polynomialTimer);

  TermMap::iterator i=terms.lower_bound(p.terms.begin()->first);

  for(TermMap::const_iterator j=p.terms.begin();j!=p.terms.end();j++)
    {
      while(i!=terms.end() && (terms.value_comp()(*i,*j)))i++;
      if(i==terms.end())
	{
	  terms.insert(i,TermMap::value_type(j->first,-j->second));
	}
      else
	{
	  if(!terms.value_comp()(*j,*i))
	    { // they must be equal
	      FieldElement c=i->second-j->second;
	      if(c.isZero())
		{
		  //		  if(i->fist.exponent==marked.exponent)marked=Term);
		  TermMap::iterator oldI=i;
		  i++;
		  terms.erase(oldI);
		}
	      else
		{
		  i->second=c;
		}
	    }
	  else
	    {
	      terms.insert(i,TermMap::value_type(j->first,-j->second));
	    }
	}
    }

  // slow subtraction
  /*  for(TermMap::const_iterator i=p.terms.begin();i!=p.terms.end();i++)
    {
      if(terms.count(i->first)==1)
        {
          terms[i->first]=terms[i->first]-i->second;
          if(terms[i->first].isZero())terms.erase(i->first);
        }
      else
        terms[i->first]=-i->second;
    }
  */
}


Polynomial operator+(const Polynomial &p, const Polynomial &q)
{
  Polynomial r(p);
  r+=q;
  return r;
}


Polynomial operator-(const Polynomial &p, const Polynomial &q)
{
  Polynomial r(p);
  r-=q;
  return r;
}


void Polynomial::operator*=(const Term &t)
{
  sugar+=t.m.exponent.sum();
  // faster multiplication
  Polynomial p(theRing);
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      FieldElement prod=i->second;
      prod*=t.c;
      p.terms.insert(p.terms.end(),TermMap::value_type(Monomial(theRing,i->first.exponent+t.m.exponent),prod));
    }
  terms=p.terms;

  // slow multiplication
  /*  Polynomial p;
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      Term T(i->second,i->first);
      T*=t;
      p+=T;
    }
  terms=p.terms;
  */
}


void Polynomial::operator*=(const Monomial &m)
{
  sugar+=m.exponent.sum();
  Polynomial p(theRing);
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      FieldElement prod=i->second;
      p.terms.insert(p.terms.end(),TermMap::value_type(Monomial(theRing,i->first.exponent+m.exponent),prod));
    }
  terms=p.terms;
}


void Polynomial::operator*=(FieldElement const &c)
{
  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      i->second*=c;
    }
}

void Polynomial::operator*=(Polynomial const &p)
{
  Polynomial r(theRing);

  for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
    {
      r+=p*Term(i->second,i->first);
    }
  *this=r;
}

Polynomial operator*(const Polynomial &p, Term const &t)
{
  Polynomial r(p);
  r*=t;

  return r;
}


Polynomial operator*(const Polynomial &p, Monomial const &m)
{
  Polynomial r(p);
  r*=m;

  return r;
}


Polynomial operator*(const Polynomial &p, FieldElement const &c)
{
  Polynomial r(p);
  r*=c;

  return r;
}


Polynomial operator*(const Polynomial &p, const Polynomial &q)
{
  Polynomial r(p);
  r*=q;

  return r;
}


bool Polynomial::isZero()const
{
  return terms.begin()==terms.end();
}

int Polynomial::numberOfTerms()const
{
  return terms.size();
}

IntegerVector Polynomial::exponentsSum()const
{
  IntegerVector sum(numberOfVariablesInRing());

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    sum+=i->first.exponent;

  return sum;
}

bool Polynomial::isMonomial()const
{
  return terms.size()==1;//could it be faster to compare begin and end iterators?
}

int Polynomial::totalDegree()const
{
  int d=-1;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(i->first.exponent.sum()>d)d=i->first.exponent.sum();

  return d;
}


int64 Polynomial::degree(IntegerVector const &w)const
{
  bool first=true;
  int64 d=0;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      if(first||dotLong(i->first.exponent,w)>d)d=dotLong(i->first.exponent,w);
      first=false;
    }

  assert(!first);

  return d;
}


Polynomial Polynomial::homogenization(PolynomialRing const &newRing, IntegerVector const *w)const   //does not compute sugar
{
  int degree;
  Polynomial ret(newRing);

  if(w)
    degree=this->degree(*w);
  else
    degree=totalDegree();

  IntegerVector m;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      int d;
      if(w)
	d=dot(v,*w);
      else
	d=v.sum();
      IntegerVector a(v.size());

      a.grow(a.size()+1);
      v.grow(v.size()+1);
      a[a.size()-1]=degree-d;

      v+=a;
      ret+=Term(i->second,Monomial(newRing,v));
      if(isMarked())if(marked.m.exponent==i->first.exponent)m=v;
    }

  if(isMarked())
    {
	  IntegerVector v=m;
      if(m.size()==newRing.getNumberOfVariables())
      {
        ret.mark(Monomial(newRing,v));
      }
	  //      assert(m.size()==newRing.getNumberOfVariables());
//      ret.mark(Monomial(newRing,m));

    }

  return ret;
}


Polynomial Polynomial::torusAct(FieldVector const &w)const
{
  Polynomial ret=*this;
  int n=theRing.getNumberOfVariables();
  assert(w.size()==n);
  for(TermMap::iterator i=ret.terms.begin();i!=ret.terms.end();i++)
    {
      FieldElement c=theRing.getField().zHomomorphism(1);
      for(int j=0;j<n;j++)
	for(int k=0;k<i->first.exponent[j];k++)
	  c*=w[j];
      i->second=i->second*c;
    }

  if(isMarked())
    {
      ret.mark(marked.m);
    }
  return ret;
}


Polynomial Polynomial::derivative()const
{
  Polynomial ret(theRing);

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      assert(v.size()==1);
      if(v[0]!=0)
	{
	  v[0]--;
	  ret+=Term(i->second*theRing.getField().zHomomorphism(v[0]+1),Monomial(theRing,v));
	}
    }

  return ret;
}


Polynomial Polynomial::deHomogenization()const
{
  Polynomial ret=*this;
  int n=numberOfVariablesInRing();
  assert(n>0);
  ret.changeNumberOfVariables(n-1);
  return ret;
}


Polynomial Polynomial::deHomogenizationInSameRing()const
{
  Polynomial ret(getRing());
  int n=getRing().getNumberOfVariables();
  assert(n>0);
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      v[n-1]=0;
      FieldElement e=i->second;
      ret+=Term(e,Monomial(theRing,v));
    }
  if(isMarked())
  {
	  ret.marked.m.exponent=marked.m.exponent.subvector(0,n-1);
  ret.isMarkedBool=true;
  }
	  return ret;
}


int Polynomial::numberOfVariablesInRing()const
{
  assert(terms.size()!=0);
  return terms.begin()->first.exponent.size();
}


void Polynomial::saturate(int variableNum)//does not compute sugar
{
  if(!terms.empty())
    {
      IntegerVector smallest=terms.begin()->first.exponent;

      for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
	smallest=min(smallest,i->first.exponent);

      if(variableNum!=-1)
	{
	  for(int j=0;j<smallest.size();j++)if(j!=variableNum)smallest[j]=0;
	}

      Polynomial p(theRing);
      for(TermMap::iterator i=terms.begin();i!=terms.end();i++)
	p.terms.insert(p.terms.end(),TermMap::value_type(Monomial(theRing,i->first.exponent-smallest),i->second));
      terms=p.terms;

      if(isMarked())
	{
	  marked.m.exponent-=smallest;
	}
    }
}


void Polynomial::computeInitialSugar()
{
  sugar=totalDegree();
}

int Polynomial::getSugar()const
{
  return sugar;
}


bool Polynomial::isMarked()const
{
	return isMarkedBool;
//  return marked.m.exponent.size()!=0;
}


bool Polynomial::isValid(int numberOfVariables)const
{
  if(!terms.empty())
    {
      if(numberOfVariables==-1)numberOfVariables=numberOfVariablesInRing();
      for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
	{
	  if(i->first.exponent.size()!=numberOfVariables)
	    {
	      fprintf(Stderr,"Polynomial::isValid failed!!!!\n");
	      return false;
	    }
	}
      if(isMarked())
	{
	  assert(marked.m.exponent.size()==numberOfVariables);
	}
    }
  return true;
}


FieldElement Polynomial::evaluate(const FieldElement &x)const
{
  FieldElement r=x-x;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      assert(v.size()==1);
      FieldElement s=theRing.getField().zHomomorphism(1);
      for(int j=0;j<v[0];j++)
	{
	  s=s*x;
	}
      r=r+i->second*s;
    }
  return r;
}

int Polynomial::maximalIndexOfVariableInSupport()const
{
  int ret=-1;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      for(int j=ret+1;j<v.size();j++)
	if(v[j])ret=j;
    }
  return ret;
}


Polynomial Polynomial::embeddedInto(PolynomialRing const &r2, list<int> const *chosenVariables)const
{
  Polynomial q(r2);

  if(chosenVariables)
    {
      assert(chosenVariables->size()==r2.getNumberOfVariables());
    }

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      IntegerVector v=i->first.exponent;
      if(chosenVariables)
	v=v.subvector(*chosenVariables);
      else
	v.resize(r2.getNumberOfVariables());
      FieldElement e=i->second;
      q+=Term(e,Monomial(r2,v));
    }
  if(isMarked())
    {
      IntegerVector m=marked.m.exponent;
      //AsciiPrinter(Stderr)<<*this<<"\n";
      if(chosenVariables)
	m=m.subvector(*chosenVariables);
      else
	m.resize(r2.getNumberOfVariables());
      q.mark(Monomial(r2,m));
    }
  //  q.marked.m.exponent.resize(r2.getNumberOfVariables());
  return q;
}


int Polynomial::numberOfVariablesInUseConsecutive()const
{
  int ret=0;

  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    {
      int a=i->first.exponent.indexOfLargestNonzeroEntry()+1;
      if(a>ret)ret=a;
    }

  return ret;
}


string Polynomial::toString(bool latex/*, bool mathMode*/)const
{
  stringstream s;
  /*
  if(latex && !mathMode)
    s << "$";
  */
  bool first=true;

  if(terms.empty())
    {
      s << "0";
      return s.str();
    }
  // If the polynomial has a marked term it is written first
  //   printString("_");
  IntegerVector e=getMarked().m.exponent;
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(e==i->first.exponent)
      {
	s << i->second.toString(i->first.exponent.isZero(),!first,latex);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))s<<"*";
	s << i->first.toString(false,false,latex);
	first=false;
      }
  //    printString("_");
  for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
    if(e!=i->first.exponent)
      {
	s << i->second.toString(i->first.exponent.isZero(),!first,latex);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))s<<"*";
	s << i->first.toString(false,false,latex);
	/*	printFieldElement(i->second,i->first.exponent.isZero(),!first);
	if((!i->first.exponent.isZero())&&(!i->second.isOne())&&(!(-(i->second)).isOne()))printString("*");
	printMonomial(i->first,false,false);*/
	first=false;
      }
  /*  if(latex && !mathMode)
    s << "$";
  */
   return s.str();
}


bool Polynomial::checkExponentVectors()const
{
	int n=getRing().getNumberOfVariables();
	for(TermMap::const_iterator i=terms.begin();i!=terms.end();i++)
		if(i->first.exponent.v.size()!=n)
		{
		//	log1
			{
			//	AsciiPrinter(Stderr)<<"Exponent vector length does not match ring\n"<<*this;
				assert(0);
			}
			return false;
		}
	return true;
}

//-----------------------------------------
// PolynomialSet
//-----------------------------------------

void PolynomialSet::saturate(int variableNum)
{
  for(iterator i=begin();i!=end();i++)
    i->saturate(variableNum);
}

PolynomialSet PolynomialSet::polynomialRingIntersection(PolynomialRing const &newRing)const
{
  PolynomialSet ret(newRing);

  for(PolynomialSet::const_iterator i=begin();i!=end();i++)
    if(i->maximalIndexOfVariableInSupport()<newRing.getNumberOfVariables())
      {
	ret.push_back(i->embeddedInto(newRing));
      }

  return ret;
}

void PolynomialSet::changeNumberOfVariables(int n)
{
  int numberOfVariables=0;

  if(n==-1)
    {
      for(PolynomialSet::const_iterator i=begin();i!=end();i++)
	{
	  if(i->getNumberOfVariables()>numberOfVariables)numberOfVariables=i->getNumberOfVariables();
	}
    }
  else
    {
      numberOfVariables=n;
    }
  //  fprintf(Stderr,"numberOFva%i\n",numberOfVariables);

  for(iterator i=begin();i!=end();i++)
    {
      i->changeNumberOfVariables(numberOfVariables);
    }
}


void PolynomialSet::markAndScale(TermOrder const &termOrder)
{
  for(iterator i=begin();i!=end();i++)
    i->mark(termOrder);
  scaleMarkedCoefficientsToOne();
}


void PolynomialSet::scaleMarkedCoefficientsToOne()
{
  for(iterator i=begin();i!=end();i++)
    i->scaleMarkedCoefficientToOne();
}


bool PolynomialSet::checkMarkings(TermOrder const &termOrder)const
{
  for(const_iterator i=begin();i!=end();i++)
    if(!i->checkMarking(termOrder))return false;

  return true;
}


bool PolynomialSet::containsInClosedGroebnerCone(IntegerVector const &v)const
{
  for(const_iterator i=begin();i!=end();i++)
    {
      if(i->degree(v)>dotLong(v,i->getMarked().m.exponent))return false;
    }
  return true;
}


bool PolynomialSet::isHomogeneous(IntegerVector const &v)const
{
  for(const_iterator i=begin();i!=end();i++)
    {
      if(!i->isHomogeneous(v))return false;
    }
  return true;
}


void PolynomialSet::unionPolynomial(const Polynomial &p)
{
  const_iterator j;
  for(j=begin();j!=end();j++)
    if((p-*j).isZero())break;
  if(j==end())push_back(p);
}


void PolynomialSet::unionSet(const PolynomialSet &s)
{
  for(const_iterator i=s.begin();i!=s.end();i++)
      unionPolynomial(*i);
}



bool PolynomialCompare::operator()(const Polynomial &a, const Polynomial &b)const
{
  if(a.terms.size()<b.terms.size())return true;
  if(b.terms.size()<a.terms.size())return false;

  TermMap::const_iterator i=a.terms.begin();
  for(TermMap::const_iterator j=b.terms.begin();j!=b.terms.end();j++)
    {
      if(LexicographicTermOrder()(i->first.exponent,j->first.exponent))return true;
      if(LexicographicTermOrder()(j->first.exponent,i->first.exponent))return false;
      i++;
    }
  if((a-b).isZero())return false;
  if((a+b).isZero())return false; // we need some way of comparing field elements

  AsciiPrinter(Stderr).printPolynomial(a);
  fprintf(Stderr,"\n");
  AsciiPrinter(Stderr).printPolynomial(b);

  assert("Polynomial compare must be improved to handle this case"==0);
  return false;
}


PolynomialCompareMarkedTerms::PolynomialCompareMarkedTerms(TermOrder const &termOrder_):
  termOrder(termOrder_)
{
}


bool PolynomialCompareMarkedTerms::operator()(const Polynomial &a, const Polynomial &b)const
{
  return termOrder(a.getMarked().m.exponent,b.getMarked().m.exponent);
}


PolynomialCompareMarkedTermsReverse::PolynomialCompareMarkedTermsReverse(TermOrder const &termOrder_):
  termOrder(termOrder_)
{
}


bool PolynomialCompareMarkedTermsReverse::operator()(const Polynomial &a, const Polynomial &b)const
{
  return termOrder(b.getMarked().m.exponent,a.getMarked().m.exponent);
}


int PolynomialSet::totalDegree()const
{
  int d=-1;
  for(const_iterator i=begin();i!=end();i++)
    if(d<i->totalDegree())d=i->totalDegree();

  return d;
}


void PolynomialSet::sort_()
{
  sort(PolynomialCompare());
}

int PolynomialSet::numberOfVariablesInRing()const
{
  return theRing.getNumberOfVariables();
  //  assert(size()!=0);
  //  return begin()->numberOfVariablesInRing();
}


bool operator==(PolynomialSet const &a, PolynomialSet const &b)
{
  return b.isEqualTo(a);
}


bool PolynomialSet::isEqualTo(PolynomialSet const &a)const
{
  if(a.size()!=size())return false;

  PolynomialSet::const_iterator j=begin();

  for(PolynomialSet::const_iterator i=a.begin();i!=a.end();i++)
    {
      if(!(*i-*j).isZero())return false;
      j++;
    }
  return true;
}


IntegerVector PolynomialSet::exponentsSum()const
{
  IntegerVector sum(numberOfVariablesInRing());

  for(const_iterator i=begin();i!=end();i++)
    sum+=i->exponentsSum();

  return sum;
}


PolynomialSet PolynomialSet::markedTermIdeal()const
{
  PolynomialSet LT(theRing);

  for(const_iterator i=begin();i!=end();i++)
    {
      LT.push_back(Polynomial(i->getMarked()));
    }

  return LT;
}


void PolynomialSet::computeInitialSugar()
{
  for(iterator i=begin();i!=end();i++)
    i->computeInitialSugar();
}


bool PolynomialSet::isMarked()const
{
  for(const_iterator i=begin();i!=end();i++)
    if(!i->isMarked())return false;

  return true;
}


PolynomialSet PolynomialSet::torusAct(FieldVector const &w)const
{
  PolynomialSet ret(theRing);
  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->torusAct(w));

  return ret;
}


PolynomialSet PolynomialSet::homogenization(PolynomialRing const &newRing, IntegerVector const *w)const
{
  PolynomialSet ret(newRing);
  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->homogenization(newRing,w));

  return ret;
}


PolynomialSet PolynomialSet::deHomogenization()const
{
  PolynomialSet ret(theRing);

  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->deHomogenization());

  return ret;
}


PolynomialSet PolynomialSet::deHomogenizationInSameRing()const
{
  PolynomialSet ret(theRing);

  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->deHomogenizationInSameRing());

  return ret;
}


bool PolynomialSet::isValid()const
{
  if(size()!=0)
    {
      int n=numberOfVariablesInRing();
      for(const_iterator i=begin();i!=end();i++)
	if(!i->isValid(n))return false;
    }

  return true;
}


PolynomialSet PolynomialSet::embeddedInto(PolynomialRing const &r2, list<int> const *chosenVariables)const
{
  PolynomialSet ret(r2);
  for(const_iterator i=begin();i!=end();i++)
    ret.push_back(i->embeddedInto(r2,chosenVariables));

  return ret;
}


int PolynomialSet::numberOfVariablesInUseConsecutive()const
{
  int ret=0;

  for(const_iterator i=begin();i!=end();i++)
    {
      int a=i->numberOfVariablesInUseConsecutive();
      if(a>ret)ret=a;
    }

  return ret;
}


int PolynomialSet::totalNumberOfTerms()const
{
  int ret=0;
  for(const_iterator i=begin();i!=end();i++)
    ret+=i->numberOfTerms();
  return ret;
}
