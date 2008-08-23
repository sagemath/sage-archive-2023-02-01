#include "ccobject.h"
#include <ginac/ginac.h>

#include <string>

using namespace GiNaC;


const symbol & get_symbol(const std::string & s)
{
    static std::map<std::string, symbol> directory;
    std::map<std::string, symbol>::iterator i = directory.find(s);
    if (i != directory.end())
        return i->second;
    else
        return directory.insert(std::make_pair(s, symbol(s))).first->second;
}

double GEx_to_double(ex& e, int* success) {
  /*
     Convert an expression to a double if possible.

     INPUT:
         e -- an expression
         success -- pointer to int
     OUTPUT:
         double -- the answer,
         *AND* sets success to true on success, and false on failure
         if e is not coercible to a double
  */
  ex f = e.evalf();
  if (is_a<numeric>(f)) {
    *success = true;
    return (ex_to<numeric>(f)).to_double();
  } else {
    *success = false;
    return 0;
  }
}

#define ASSIGN_WRAP(x,y) x = y

#define ADD_WRAP(x,y) (x)+(y)
#define SUB_WRAP(x,y) (x)-(y)
#define MUL_WRAP(x,y) (x)*(y)
#define DIV_WRAP(x,y) (x)/(y)

#define LT_WRAP(x,y)  (x)<(y)
#define EQ_WRAP(x,y)  (x)==(y)
#define GT_WRAP(x,y)  (x)>(y)
#define LE_WRAP(x,y)  (x)<=(y)
#define NE_WRAP(x,y)  (x)!=(y)
#define GE_WRAP(x,y)  (x)>=(y)

