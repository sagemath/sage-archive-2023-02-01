#include <iostream>
#include "parser.h"
#include "printer.h"
#include "polynomial.h"
#include "division.h"
#include "buchberger.h"
#include "wallideal.h"
#include "lp.h"
#include "reversesearch.h"
#include "breadthfirstsearch.h"
#include "termorder.h"
#include "ep_standard.h"
#include "ep_xfig.h"
#include "gfanapplication.h"
#include "polyhedralfan.h"
#include "minkowskidual.h"
#include "log.h"

class MinkowskiEnumerationTarget:public EnumerationTarget
{
  SymmetryGroup const &s;
 public:
  IntegerVectorList weightVectors;
  MinkowskiEnumerationTarget(SymmetryGroup const &s_):
    s(s_)
  {
  }
  virtual void beginEnumeration(const PolynomialSet &groebnerBasis)
  {
  }
  virtual void endEnumeration()
  {
  }
  virtual bool basis(const PolynomialSet &groebnerBasis)
  {
    //    log0 cerr<<"ETST"<<endl;
    PolyhedralCone c=coneFromMarkedBasis(groebnerBasis);
    weightVectors.push_back(PolyhedralFan::stableRay(c,&s));
    return true;
  } /* return false to break enumaration */
};


class MinkowskiApplication : public GFanApplication
{
  SimpleOption optionSymmetry;
  SimpleOption optionDisableSymmetryTest;
  SimpleOption optionIgnoreCones;
  SimpleOption optionPartOne;
  SimpleOption optionPartTwo;
public:
  const char *helpText()
  {
    return "This is a program for computing the normal fan of the Minkowski sum of the Newton polytopes of a list of polynomials.\n";
  }
  MinkowskiApplication():
    optionSymmetry("--symmetry",
		   "Tells the program to read in generators for a group of symmetries (subgroup of $S_n$) after having read in the ideal. The program checks that the ideal stays fixed when permuting the variables with respect to elements in the group. The program uses breadth first search to compute the set of reduced Groebner bases up to symmetry with respect to the specified subgroup.\n"),
    optionDisableSymmetryTest("--disableSymmetryTest","When using --symmetry this option will disable the check that the group read off from the input actually is a symmetry group with respect to the input ideal.\n"),
    optionIgnoreCones("--nocones","Tell the program to not list cones in the output."),
    optionPartOne("--partOne","Do only the fist part of the computation i.e. compute a weight vector for each vertex.\n"),
    optionPartTwo("--partTwo","Do only the second part of the computation .\n")
  {
    registerOptions();
    optionPartOne.hide();
    optionPartTwo.hide();
  }

  char *name()
  {
    return "_minkowskisum";
  }

  int main()
  {
    LexicographicTermOrder myOrder;

    PolynomialSet g=FileParser(Stdin).parsePolynomialSetWithRing();
    g.sort_();
    g.markAndScale(myOrder);


    SymmetryGroup s(g.numberOfVariablesInRing());

    IntegerVectorList generators;
    if(optionSymmetry.getValue())
      {
	generators=FileParser(Stdin).parseIntegerVectorList();
	if(!optionDisableSymmetryTest.getValue())
	  {
	    /*	    for(IntegerVectorList::iterator i=generators.begin();i!=generators.end();i++)
	      {
		PolynomialSet G=SymmetryGroup::permutePolynomialSet(g,*i);
		G.sort_();
		assert(g==G);
	      }
	    */
	    assertSymmetriesMatch(generators, g,0,true);//eldMatrix const *torusActions);
	  }
	s.computeClosure(generators);

	s.createByteTable();
      }
    //    EnumerationTargetCollector ep;

    IntegerVectorList weightVectors;
    if(optionPartTwo.getValue())
      {
	weightVectors=FileParser(Stdin).parseIntegerVectorList();
      }
    else
      {
	MinkowskiEnumerationTarget ep(s);
	log1 fprintf(Stderr,"Traversing normal fan.\n");
	BreadthFirstSearch rs(s,true);
	rs.setEnumerationTarget(&ep);
	rs.enumerate(g);
	weightVectors=ep.weightVectors;
      }
    if(optionPartOne.getValue())
      {
	AsciiPrinter P(Stdout);
	P.printPolynomialSet(g);
	if(optionSymmetry.getValue())
	  P.printVectorList(generators);
	P.printVectorList(weightVectors);
      }
    else
      {
	log1 fprintf(Stderr,"Computing polyhedral cones.\n");
	int n=g.getRing().getNumberOfVariables();
	PolyhedralFan f(n);
	int counter=0;
	for(IntegerVectorList::const_iterator i=weightVectors.begin();i!=weightVectors.end();i++)
	  {
      {
	static int t;
	if(!((t++)%10))log1 fprintf(Stderr,"%i\n",t);
      }

      //            if(counter++==35)break;
	    //log0 fprintf(Stderr,"Weight vektor:\n");
	    //log0 AsciiPrinter(Stderr).printVector(*i);
	    // log0 fprintf(Stderr,"1");
	    WeightTermOrder T(*i);
	    g.markAndScale(T);
	    // log0 fprintf(Stderr,"2");
	    PolyhedralCone c=coneFromMarkedBasis(g);
	    // log0 fprintf(Stderr,"3");
	    //inequalities=
	      /*wallInequalities(g);
	    fprintf(stderr,"%i ",inequalities.size());
	    inequalities=normalizedWithSumsAndDuplicatesRemoved(inequalities);
	    fprintf(stderr,"%i\n",inequalities.size());
	    IntegerVectorList empty;
	    PolyhedralCone c(inequalities,empty,n);
	      */

	    c.canonicalize();
	    //log0 fprintf(Stderr,"4");
	    f.insert(c);
	    //log0 fprintf(Stderr,"5\n");
	    //static int i; // XXX redefinition of i -- an error with GCC 4.7.0
	    //log0 fprintf(Stderr,"inserted:%i\n",++i);
	  }
	log1 fprintf(Stderr,"Resolving symmetries.\n");

	/*	{
	  SymmetricComplex c=dualMinkowskiMixed(g,s,f);
	  return 0;
	  }*/


	AsciiPrinter P(Stdout);
	f.printWithIndices(&P,
			   FPF_conesExpanded|
			   (optionIgnoreCones.getValue()?0:FPF_cones|FPF_maximalCones)|
			   (optionSymmetry.getValue()?FPF_group|FPF_conesCompressed:0),
			   &s);
	//	f.printWithIndices(&P,false, &s,optionSymmetry.getValue(),optionIgnoreCones.getValue());
      }

    return 0;
  }
};

static MinkowskiApplication theApplication;

