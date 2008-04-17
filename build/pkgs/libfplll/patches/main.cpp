/*
  Copyright 2005, 2006, 2007 David Cadé, Damien Stehlé.

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2 of the License, or (at your
  option) any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; see the file COPYING.  If not, write to the Free
  Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  02111-1307, USA.

  This program implements ideas from the paper "Floating-point LLL Revisited",
  by Phong Nguyen and Damien Stehlé, in the Proceedings of Eurocrypt'2005,
  Springer-Verlag; and was partly inspired by Shoup's NTL library:
  http://www.shoup.net/ntl/

*/

#include <iostream>
#include <cstring>
#include "fplll.h"

using namespace std;


#define ZMPZ 0
#define ZINT 1

#define FDPE 0
#define FDOUBLE 1
#define FMPFR 2

#define MPROVED 0
#define MHEURISTIC 1
#define MFASTEARLY 3
#define MFAST 2
#define MWRAPPER 4
#define MHEUREARLY 5

int main (int argc, char ** argv)
{
  /*arg parsing (rewrite)*/
  int r=0;
  int c=0;

  int pr=0;

  double eta=0.51;
  double delta=0.99;

  char z=ZMPZ;
  char f=FDPE;
  char m=MWRAPPER;

  int ac=1;
  while (ac<argc)
    {
      if (argv[ac][0]!='-')
	{
	  cerr<<"Parse error : option expected\n";
	  abort();
	}
      else
	{
	  /* TODO  :better bound checking, name complete check */
	  switch (argv[ac][1])
	    {
	    case 'p':
	      ++ac;
	      pr=atoi(argv[ac]);
	      break;
	    case 'r':
	      ++ac;
	      r=atoi(argv[ac]);
	      break;
	    case 'c':
	      ++ac;
	      c=atoi(argv[ac]);
	      break;
	    case 'e':
	      ++ac;
	      eta=atof(argv[ac]);
	      break;
	    case 'd':
	      ++ac;
	      delta=atof(argv[ac]);
	      break;
	    case 'z':
	      ++ac;
	      if (strcmp("mpz",argv[ac])==0)
		z=ZMPZ;
	      else if (strcmp("int",argv[ac])==0)
		z=ZINT;
	      else
		{
		  cerr<<"Parse error in -z switch : int or mpz expected\n";
		  abort();
		}
	      break;
	    case 'f':
	      ++ac;
	      if (strcmp("mpfr",argv[ac])==0)
		f=FMPFR;
	      else if (strcmp("dpe",argv[ac])==0)
		f=FDPE;
	      else if (strcmp("double",argv[ac])==0)
		f=FDOUBLE;
	      else
		{
		  cerr<<"Parse error in -f switch : mpfr, dpe or double expected\n";
		  abort();
		}
	      break;
	    case 'm':
	      ++ac;
	      if (strcmp("proved",argv[ac])==0)
		m=MPROVED;
	      else if (strcmp("heuristic",argv[ac])==0)
		m=MHEURISTIC;
	      else if (strcmp("fast",argv[ac])==0)
		m=MFAST;
	      else if (strcmp("fastearly",argv[ac])==0)
		m=MFASTEARLY;
	      else if (strcmp("heuristicearly",argv[ac])==0)
		m=MHEUREARLY;
	      else if (strcmp("wrapper",argv[ac])==0)
		m=MWRAPPER;
	      else
		{
		  cerr<<"Parse error in -m switch : proved, heuristic, fast, heuristicearly, fastearly, or wrapper expected\n";
		  abort();
		}
	      break;
	    }
	  ++ac;
	}
    }

#ifdef DEBUG
  printf("Arguments read c=%d r=%d\n", c, r);
#endif

  /* here we instantiate according to the desired method*/
  switch (m)
    {
    case MWRAPPER:
      {
	ZZ_mat<mpz_t> * zzmat= new ZZ_mat<mpz_t>(r,c);
	if (zzmat->read()!=0) {cerr.flush();return 1;}

	wrapper* w=new wrapper ( zzmat,0,eta,delta);
	w->LLL();
	w->GetBase()->print();
	delete w;
	delete zzmat;
      }
      break;

    case MPROVED:
      switch (z)
	{
	case ZINT:
	  {
	    ZZ_mat<long int> * zzmat= new ZZ_mat<long int>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FMPFR:
		{
		  proved<long int,mpfr_t>* p=new proved<long int,mpfr_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDPE:
		{
		  proved<long int,dpe_t>* p=new proved<long int,dpe_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDOUBLE:
		{
		  proved<long int,double>* p=new proved<long int,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      }
	    delete zzmat;
	  }
	  break;
	case ZMPZ:
	  {
	    ZZ_mat<mpz_t> * zzmat= new ZZ_mat<mpz_t>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FMPFR:
		{
		  proved<mpz_t,mpfr_t>* p=new proved<mpz_t,mpfr_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDPE:
		{
		  proved<mpz_t,dpe_t>* p=new proved<mpz_t,dpe_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDOUBLE:
		{
		  proved<mpz_t,double>* p=new proved<mpz_t,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      }
	    delete zzmat;
	  }
	  break;
	}
      break;
    case MHEURISTIC:
       switch (z)
	{
	case ZINT:
	  {
	    ZZ_mat<long int> * zzmat= new ZZ_mat<long int>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FMPFR:
		{
		  heuristic<long int,mpfr_t>* p=new heuristic<long int,mpfr_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDPE:
		{
		  heuristic<long int,dpe_t>* p=new heuristic<long int,dpe_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDOUBLE:
		{
		  heuristic<long int,double>* p=new heuristic<long int,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      }
	    delete zzmat;
	  }
	  break;
	case ZMPZ:
	  {
	    ZZ_mat<mpz_t> * zzmat= new ZZ_mat<mpz_t>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}
	    switch (f)
	      {
	      case FMPFR:
		{
		  heuristic<mpz_t,mpfr_t>* p=new heuristic<mpz_t,mpfr_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDPE:
		{
		  heuristic<mpz_t,dpe_t>* p=new heuristic<mpz_t,dpe_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDOUBLE:
		{
		  heuristic<mpz_t,double>* p=new heuristic<mpz_t,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      }
	    delete zzmat;
	  }
	  break;
	}
       break;
    case MHEUREARLY:
       switch (z)
	{
	case ZINT:
	  {
	    ZZ_mat<long int> * zzmat= new ZZ_mat<long int>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FMPFR:
		{
		  heuristic_early_red<long int,mpfr_t>* p=new heuristic_early_red<long int,mpfr_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDPE:
		{
		  heuristic_early_red<long int,dpe_t>* p=new heuristic_early_red<long int,dpe_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDOUBLE:
		{
		  heuristic_early_red<long int,double>* p=new heuristic_early_red<long int,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      }
	    delete zzmat;
	  }
	  break;
	case ZMPZ:
	  {
	    ZZ_mat<mpz_t> * zzmat= new ZZ_mat<mpz_t>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}
	    switch (f)
	      {
	      case FMPFR:
		{
		  heuristic_early_red<mpz_t,mpfr_t>* p=new heuristic_early_red<mpz_t,mpfr_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDPE:
		{
		  heuristic_early_red<mpz_t,dpe_t>* p=new heuristic_early_red<mpz_t,dpe_t> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      case FDOUBLE:
		{
		  heuristic_early_red<mpz_t,double>* p=new heuristic_early_red<mpz_t,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      }
	    delete zzmat;
	  }
	  break;
	}
       break;
    case MFAST:
      switch (z)
	{
	case ZINT:
	  {
	    ZZ_mat<long int> * zzmat= new ZZ_mat<long int>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FDOUBLE:
		{
		  fast<long int,double>* p=new fast<long int,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      default:
		cerr << "Double required \n";
		break;
	      }
	    delete zzmat;
	  }
	  break;
	case ZMPZ:
	  {
	    ZZ_mat<mpz_t> * zzmat= new ZZ_mat<mpz_t>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FDOUBLE:
		{
		  fast<mpz_t,double>* p=new fast<mpz_t,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      default:
		cerr<<"Double required\n";
		break;
	      }
	    delete zzmat;
	  }
	  break;
	}
      break;
    case MFASTEARLY:
      switch (z)
	{
	case ZINT:
	  {
	    ZZ_mat<long int> * zzmat= new ZZ_mat<long int>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FDOUBLE:
		{
		  fast_early_red<long int,double>* p=new fast_early_red<long int,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      default:
		cerr << "Double required \n";
		break;
	      }
	    delete zzmat;
	  }
	  break;
	case ZMPZ:
	  {
	    ZZ_mat<mpz_t> * zzmat= new ZZ_mat<mpz_t>(r,c);
	    if (zzmat->read()!=0) {cerr.flush();return 1;}

	    switch (f)
	      {
	      case FDOUBLE:
		{
		  fast_early_red<mpz_t,double>* p=new fast_early_red<mpz_t,double> (zzmat,pr,eta,delta);
		  p->LLL();
		  p->GetBase()->print();
		  delete p;
		}
		break;
	      default:
		cerr<<"Double required\n";
		break;
	      }
	    delete zzmat;
	  }
	  break;
	}
      break;
    }
}

