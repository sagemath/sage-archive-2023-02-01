// Logic for info flags of function expressions
//
// (c) 2016,2017  Ralf Stephan <ralf@ark.in-berlin.de>
// Distributed under GPL2, see http://www.gnu.org

#include <unordered_map>

#include "inifcns.h"
#include "function.h"
#include "operators.h"
#include "utils.h"

namespace GiNaC {

static bool exp_info(const function& f, unsigned inf)
{
        const ex& arg = f.op(0);
        switch (inf) {
        case info_flags::nonzero:
                return true;
        case info_flags::real:
        case info_flags::positive:
                return arg.info(info_flags::real);
        }
        return false;
}

static bool log_info(const function& f, unsigned inf)
{
        const ex& arg = f.op(0);
        switch (inf) {
        case info_flags::real:
                return arg.info(info_flags::positive);
        case info_flags::positive:
                return arg.info(info_flags::real) and
                        (arg-_ex1).info(info_flags::positive);
        case info_flags::negative:
                return arg.info(info_flags::real) and
                        arg.info(info_flags::positive) and
                        (arg-_ex1).info(info_flags::negative);
        }
        return false;
}

static bool trig_info(const function& f, unsigned inf)
{
        const ex& arg = f.op(0);
        switch (inf) {
        case info_flags::real:
                return arg.info(info_flags::real);
        }
        return false;
}

static bool sin_info(const function& f, unsigned inf)
{
        return trig_info(f, inf);
}

static bool cos_info(const function& f, unsigned inf)
{
        return trig_info(f, inf);
}

static bool tan_info(const function& f, unsigned inf)
{
        return trig_info(f, inf);
}

static bool cot_info(const function& f, unsigned inf)
{
        return trig_info(f, inf);
}

static bool sec_info(const function& f, unsigned inf)
{
        return trig_info(f, inf);
}

static bool csc_info(const function& f, unsigned inf)
{
        return trig_info(f, inf);
}

static bool atan_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::negative:
        case info_flags::positive:
        case info_flags::nonzero:
                return f.op(0).info(inf);
        }
        return trig_info(f, inf);
}

static bool acot_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::positive:
        case info_flags::nonzero:
                return true;
        }
        return trig_info(f, inf);
}

static bool sinh_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::negative:
        case info_flags::positive:
        case info_flags::nonzero:
                return f.op(0).info(inf);
        }
        return trig_info(f, inf);
}

static bool cosh_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::positive:
        case info_flags::nonzero:
                return f.op(0).info(info_flags::real);
        }
        return trig_info(f, inf);
}

static bool tanh_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::negative:
        case info_flags::positive:
        case info_flags::nonzero:
                return f.op(0).info(inf);
        }
        return trig_info(f, inf);
}

static bool coth_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::positive:
        case info_flags::negative:
                return f.op(0).info(inf);
        case info_flags::nonzero:
                return true;
        }
        return trig_info(f, inf);
}

static bool sech_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::positive:
        case info_flags::nonnegative:
        case info_flags::nonzero:
                return f.op(0).info(info_flags::real);
        }
        return trig_info(f, inf);
}

static bool csch_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::negative:
        case info_flags::positive:
                return f.op(0).info(inf);
        case info_flags::nonzero:
                return true;
        }
        return trig_info(f, inf);
}

static bool asinh_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::positive:
        case info_flags::negative:
        case info_flags::nonzero:
                return f.op(0).info(inf);
        }
        return trig_info(f, inf);
}

static bool acsch_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::nonnegative:
        case info_flags::negative:
        case info_flags::positive:
                return f.op(0).info(inf);
        case info_flags::nonzero:
                return true;
        }
        return trig_info(f, inf);
}

static bool gamma_info(const function& f, unsigned inf)
{
        const ex& arg = f.op(0);
        switch (inf) {
        case info_flags::real:
        case info_flags::positive:
        case info_flags::nonzero:
                return arg.info(info_flags::positive);
        case info_flags::integer:
                return arg.info(info_flags::integer)
                   and arg.info(info_flags::positive);
        case info_flags::even:
                return arg.info(info_flags::integer)
                   and (arg+_ex_2).info(info_flags::positive);
        }
        return false;
}

static bool zeta_info(const function& f, unsigned inf)
{
        const ex& arg = f.op(0);
        switch (inf) {
        case info_flags::real:
                return arg.info(inf);
        }
        return false;
}

static bool abs_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::real:
        case info_flags::nonnegative:
                return true;
        case info_flags::integer:
        case info_flags::nonzero:
        case info_flags::even:
                return f.op(0).info(inf);
        case info_flags::positive:
                return f.op(0).info(info_flags::nonzero);
        }
        return false;
}

static bool real_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::real:
                return true;
        case info_flags::nonnegative:
        case info_flags::negative:
        case info_flags::positive:
        case info_flags::integer:
        case info_flags::even:
                return f.op(0).info(inf);
        }
        return false;
}

static bool imag_info(const function& f, unsigned inf)
{
        switch (inf) {
        case info_flags::real:
                return true;
        }
        return false;
}

static bool factorial_info(const function& f, unsigned inf)
{
        const ex& arg = f.op(0);
        switch (inf) {
        case info_flags::real:
        case info_flags::nonnegative:
        case info_flags::integer:
                return arg.info(inf);
        case info_flags::even:
                return (arg.info(info_flags::integer)
                        and (arg+_ex_1).info(info_flags::positive))
                       or (arg.info(info_flags::even)
                        and arg.info(info_flags::positive));
        }
        return false;
}

static bool binomial_info(const function& f, unsigned inf)
{
        const ex& arg = f.op(0);
        switch (inf) {
        case info_flags::real:
        case info_flags::nonnegative:
                return arg.info(inf);
        case info_flags::integer:
                return arg.info(inf) and f.op(1).info(inf);
        }
        return false;
}


bool function::info(unsigned inf) const
{
        using ifun_t = decltype(exp_info);
        static std::unordered_map<unsigned int,ifun_t*> funcmap {{
                {exp_SERIAL::serial, &exp_info},
                {log_SERIAL::serial, &log_info},
                {sin_SERIAL::serial, &sin_info},
                {cos_SERIAL::serial, &cos_info},
                {tan_SERIAL::serial, &tan_info},
                {cot_SERIAL::serial, &cot_info},
                {sec_SERIAL::serial, &sec_info},
                {csc_SERIAL::serial, &csc_info},
                {atan_SERIAL::serial, &atan_info},
                {acot_SERIAL::serial, &acot_info},
                {sinh_SERIAL::serial, &sinh_info},
                {cosh_SERIAL::serial, &cosh_info},
                {tanh_SERIAL::serial, &tanh_info},
                {coth_SERIAL::serial, &coth_info},
                {sech_SERIAL::serial, &sech_info},
                {csch_SERIAL::serial, &csch_info},
                {asinh_SERIAL::serial, &asinh_info},
                {acsch_SERIAL::serial, &acsch_info},
                {gamma_SERIAL::serial, &gamma_info},
                {zeta1_SERIAL::serial, &zeta_info},
                {abs_SERIAL::serial, &abs_info},
                {real_part_function_SERIAL::serial, &real_info},
                {imag_part_function_SERIAL::serial, &imag_info},
                {binomial_SERIAL::serial, &binomial_info},
                {factorial_SERIAL::serial, &factorial_info},
        }};

        auto search = funcmap.find(serial);
        if (search == funcmap.end())
                return iflags.get(inf);
        return (*(search->second))(*this, inf);
}

}
