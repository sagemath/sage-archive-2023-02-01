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

#define ADD_WRAP(x,y) (x)+(y)
#define MUL_WRAP(x,y) (x)*(y)
