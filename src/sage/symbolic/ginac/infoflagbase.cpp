/*
 * File:   infoflagbase.cpp
 * Author: Ralf Stephan <ralf@ark.in-berlin.de>
 *
 * Created on August 3, 2015, 9:31 AM
 */

#include <stdexcept>
#include "infoflagbase.h"

namespace GiNaC {

infoflagbase::infoflagbase() {
        if (index == nullptr)
                init_index();
}

infoflagbase::infoflagbase(const infoflagbase& orig) {
}

infoflagbase::~infoflagbase() {
}

//------------------------------------------
constexpr unsigned const infoflagbase::flags[];
unsigned* infoflagbase::index = nullptr;
unsigned infoflagbase::max;

void infoflagbase::init_index()
{
        max = 0;
        for (unsigned i : flags)
                if (max < i)
                        max = i;

        index = new unsigned(max+1);
        unsigned ctr = 0;
        for (unsigned i : flags)
                index[i] = ctr++; 
}

bool infoflagbase::get(unsigned flag) const
{
        if (flag > max)
                throw(std::runtime_error("requested wrong info flag"));
        return bits[index[flag]];
}


} // namespace GiNaC
