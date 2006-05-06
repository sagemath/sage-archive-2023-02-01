/* emacs edit mode for this file is -*- C++ -*- */
/* $Id: ftmpl_inst.cc,v 1.9 1998/06/30 16:25:47 schmidt Exp $ */

//{{{ docu
//
// ftmpl_inst.cc - Factory's template instantiations.
//
// For a detailed description how to instantiate Factory's
// template classes and functions and how to add new
// instantiations see the `README' file.
//
//}}}

#include <factoryconf.h>

#ifdef macintosh
#include <:templates:ftmpl_array.cc>
#include <:templates:ftmpl_factor.cc>
#include <:templates:ftmpl_list.cc>
#include <:templates:ftmpl_functions.h>
#include <:templates:ftmpl_matrix.cc>
#else
#include <templates/ftmpl_array.cc>
#include <templates/ftmpl_factor.cc>
#include <templates/ftmpl_list.cc>
#include <templates/ftmpl_functions.h>
#include <templates/ftmpl_matrix.cc>
#endif

#include <factory.h>


//{{{ explicit template class instantiations
template class Factor<CanonicalForm>;
template class List<CFFactor>;
template class ListItem<CFFactor>;
template class ListIterator<CFFactor>;
template class List<CanonicalForm>;
template class ListItem<CanonicalForm>;
template class ListIterator<CanonicalForm>;
template class Array<CanonicalForm>;
template class List<MapPair>;
template class ListItem<MapPair>;
template class ListIterator<MapPair>;
template class Matrix<CanonicalForm>;
template class SubMatrix<CanonicalForm>;
template class Array<REvaluation>;
//}}}

//{{{ explicit template function instantiations
#ifndef NOSTREAMIO
template ostream & operator << ( ostream &, const List<CanonicalForm> & );
template ostream & operator << ( ostream &, const List<CFFactor> & );
template ostream & operator << ( ostream &, const List<MapPair> & );
template ostream & operator << ( ostream &, const Array<CanonicalForm> & );
template ostream & operator << ( ostream &, const Factor<CanonicalForm> & );
template ostream & operator << ( ostream &, const Matrix<CanonicalForm> & );
template ostream & operator << ( ostream &, const Array<REvaluation> & );
#endif /* NOSTREAMIO */

template int operator == ( const Factor<CanonicalForm> &, const Factor<CanonicalForm> & );

template List<CFFactor> Union ( const List<CFFactor> &, const List<CFFactor> & );

#if ! defined(WINNT) || defined(__GNUC__)
template CanonicalForm tmax ( const CanonicalForm &, const CanonicalForm & );
template CanonicalForm tmin ( const CanonicalForm &, const CanonicalForm & );

template Variable tmax ( const Variable &, const Variable & );
template Variable tmin ( const Variable &, const Variable & );

template int tmax ( const int &, const int & );
template int tmin ( const int &, const int & );
template int tabs ( const int & );
#endif
//}}}

//
// place here your own template stuff, not yet instantiated by factory
//
