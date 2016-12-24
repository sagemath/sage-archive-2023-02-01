// Algorithms for indefinite and definite summation
//
// (c) 2016  Ralf Stephan <ralf@ark.in-berlin.de>
// Distributed under GPL2, see http://www.gnu.org
//
// Ref.:       1. W. Koepf, Algorithms for m-fold Hypergeometric Summation,
//                  Journal of Symbolic Computation (1995) 20, 399-417

#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

#include "inifcns.h"
#include "ex.h"
#include "mul.h"
#include "add.h"
#include "numeric.h"
#include "function.h"
#include "symbol.h"
#include "utils.h"
#include "operators.h"
#include "relational.h"
#include "normal.h"
#include "upoly.h"
#include "mpoly.h"
#include "matrix.h"

namespace GiNaC {

class gosper_domain_error : public std::domain_error {
        public:
                gosper_domain_error() : std::domain_error("") {}
};

static bool is_rational_linear(const ex& the_ex)
{
        if (is_exactly_a<symbol>(the_ex))
                return true;
        if (is_exactly_a<numeric>(the_ex))
                return (ex_to<numeric>(the_ex).is_mpz()
                                or ex_to<numeric>(the_ex).is_mpq());

        if (is_exactly_a<mul>(the_ex)) {
                const mul& m = ex_to<mul>(the_ex);
                for (unsigned int i=0; i<m.nops(); i++) {
                        if (not is_exactly_a<symbol>(m.op(i))
                            and not is_exactly_a<numeric>(m.op(i)))
                                return false;
                }
                const ex& oc = m.op(m.nops());
                return (is_exactly_a<numeric>(oc)
                       and (ex_to<numeric>(oc).is_mpz()
                            or ex_to<numeric>(oc).is_mpq()));
        }
        if (is_exactly_a<add>(the_ex)) {
                const add& a = ex_to<add>(the_ex);
                for (unsigned int i=0; i<a.nops(); i++) {
                        if (not is_rational_linear(a.op(i)))
                                return false;
                }
                const ex& oc = a.op(a.nops());
                return (is_exactly_a<numeric>(oc)
                       and (ex_to<numeric>(oc).is_mpz()
                            or ex_to<numeric>(oc).is_mpq()));
        }
        return false;
}

static ex factorial_to_gamma(const function& f)
{
        return tgamma(f.op(0) + _ex1);
}

static ex gamma_to_gamma(const function& f) { return ex(f); }

static ex binomial_to_gamma(const function& f)
{
        const ex& a = f.op(0);
        const ex& k = f.op(1);
        if (is_exactly_a<numeric>(a)) {
                numeric anum = ex_to<numeric>(a);
                if (anum.info(info_flags::integer)
                    and anum.info(info_flags::negative))
                        return pow(_ex_1, k) * 
                                (tgamma(k - a) / (tgamma(k+1) * (anum-*_num1_p).factorial()));
        }
        ex t = (k - a).expand();
        if (is_exactly_a<numeric>(t)
            and ex_to<numeric>(t).info(info_flags::integer)
            and ex_to<numeric>(t).info(info_flags::negative))
                return _ex0;

        return tgamma(a+1) / (tgamma(k+1) * tgamma(a-k+1));
}

static ex rising_factorial_to_gamma(const function& f)
{
        return tgamma(f.op(0) + f.op(1)) / tgamma(f.op(0));
}

static ex falling_factorial_to_gamma(const function& f)
{
        return tgamma(f.op(0) + _ex1) / tgamma(f.op(0) - f.op(1) + _ex_1);
}

using tgfun_t = decltype(gamma_to_gamma);
static std::unordered_map<unsigned int,tgfun_t*> funcmap {{
        {factorial_SERIAL::serial, &factorial_to_gamma},
        {tgamma_SERIAL::serial, &gamma_to_gamma},
        {binomial_SERIAL::serial, &binomial_to_gamma},
        {rising_factorial_SERIAL::serial, &rising_factorial_to_gamma},
        {falling_factorial_SERIAL::serial, &falling_factorial_to_gamma},
}};

static bool has_suitable_form(ex the_ex)
{
        if (is_rational_linear(the_ex))
                return true;
        if (is_exactly_a<power>(the_ex)) {
                power pow = ex_to<power>(the_ex);
                const ex& expo = pow.op(1);
                if (is_exactly_a<numeric>(expo)
                    and expo.info(info_flags::integer))
                        return has_suitable_form(pow.op(0));
                return (is_rational_linear(pow.op(0))
                     and is_rational_linear(pow.op(1)));
        }
        if (is_exactly_a<function>(the_ex)) {
                function f = ex_to<function>(the_ex);
                if (funcmap.find(f.get_serial()) == funcmap.end())
                        return false;
                for (unsigned int i=0; i<f.nops(); i++)
                        if (not is_rational_linear(f.op(i)))
                                return false;
                return true;
        }
        if (is_exactly_a<mul>(the_ex)) {
                const mul& m = ex_to<mul>(the_ex);
                for (unsigned int i=0; i<m.nops(); i++) {
                        if (not has_suitable_form(m.op(i)))
                                return false;
                }
                const ex& oc = m.op(m.nops());
                return (is_exactly_a<numeric>(oc)
                       and (ex_to<numeric>(oc).is_mpz()
                            or ex_to<numeric>(oc).is_mpq()));
        }
        return false;
}

static ex to_gamma(const ex& the_ex)
{
        if (is_rational_linear(the_ex))
                return the_ex;
        if (is_exactly_a<power>(the_ex)) {
                power pow = ex_to<power>(the_ex);
                const ex& expo = pow.op(1);
                if (is_exactly_a<numeric>(expo)
                    and expo.info(info_flags::integer))
                        return power(to_gamma(pow.op(0)), expo);
                return the_ex;
        }
        if (is_exactly_a<function>(the_ex)) {
                function f = ex_to<function>(the_ex);
                auto search = funcmap.find(f.get_serial());
                if (search == funcmap.end())
                        return the_ex;
                return (*search->second)(f);
        }
        if (is_exactly_a<mul>(the_ex)) {
                const mul& m = ex_to<mul>(the_ex);
                exvector vec;
                for (unsigned int i=0; i<m.nops(); i++)
                        vec.push_back(to_gamma(m.op(i)));
                vec.push_back(m.op(m.nops()));
                return mul(vec);
        }
        throw std::runtime_error("can't happen in to_gamma");
}

using ex_int_pair = std::pair<ex, int>;
static void simp_gamma_pair(ex& res, ex_int_pair& np, ex_int_pair& dp,
                int m, int d)
{
        np.second -= m;
        dp.second -= m;
        if (d == 0)
                return;
        if (d>0) {
                ex base = dp.first;
                for (int k=0; k<d; ++k) {
                        if (m == 1)
                                res *= base;
                        else
                                res *= power(base, m);
                        base = base + _ex1;
                }
        }
        else {
                ex base = np.first;
                for (int k=0; k<-d; ++k) {
                        res *= power(base, -m);
                        base = base + _ex1;
                }
        }
}

static ex simplify_gamma(ex the_ex)
// not a general purpose algorithm because e.g. we know exponents are integer
{
        if (not is_exactly_a<mul>(the_ex))
                return the_ex;

        // The main data struct holds gamma arguments and their exponents
        std::vector<ex_int_pair> nums, dens;
        ex res = _ex1;
        const mul& m = ex_to<mul>(the_ex);
        for (unsigned int i=0; i<m.nops(); i++) {
                const ex& term = m.op(i);
                if (is_exactly_a<function>(term)) {
                        function f = ex_to<function>(term);
                        if (f.get_serial() == tgamma_SERIAL::serial) {
                                nums.push_back({f.op(0), 1});
                                continue;
                        }
                }
                else if (is_exactly_a<power>(term)) {
                        const ex& basis = term.op(0);
                        numeric e = ex_to<numeric>(term.op(1));
                        if (is_exactly_a<function>(basis)) {
                                function f = ex_to<function>(basis);
                                if (f.get_serial() == tgamma_SERIAL::serial) {
                                        if (e > 0)
                                                nums.push_back({f.op(0), e.to_int()});
                                        else
                                                dens.push_back({f.op(0), -(e.to_int()) });
                                        continue;
                                }
                        }
                }
                res = res * term;
        }

        ex simplified = _ex1;
        bool anything_changed = true;
        while (anything_changed) {
                anything_changed = false;
                for (auto& np : nums) {
                        if (np.second == 0)
                                continue;
                        for (auto& dp : dens) {
                                if (dp.second == 0)
                                        continue;
                                ex diff = np.first - dp.first;
                                int imin = std::min(np.second, dp.second);
                                if (not is_exactly_a<numeric>(diff))
                                        continue;
                                numeric dnum = ex_to<numeric>(diff);
                                if (not dnum.info(info_flags::integer)
                                    or dnum>numeric(4096)) // arbitrary cutoff
                                        continue;
                                int d = dnum.to_int();
                                simp_gamma_pair(simplified, np, dp, imin, d);
                                anything_changed = true;
                        }
                }
        }
        ex remaining_gammas = _ex1;
        for (const auto& p : nums)
                if (p.second == 1)
                        remaining_gammas = remaining_gammas * tgamma(p.first);
                else if (p.second != 0)
                        remaining_gammas = remaining_gammas * power(tgamma(p.first), numeric(p.second));
        for (const auto& p : dens)
                if (p.second != 0)
                        remaining_gammas = remaining_gammas * power(tgamma(p.first), numeric(-p.second));
        return res * simplified * remaining_gammas;
}

static ex combine_powers(const ex& the_ex)
{
        if (not is_exactly_a<mul>(the_ex))
                return the_ex;

        // The main data struct holds gamma arguments and their exponents
        exmap factors;
        const mul& m = ex_to<mul>(the_ex);
        ex res = _ex1;
        for (unsigned int i=0; i<m.nops(); i++) {
                const ex& term = m.op(i);
                if (is_exactly_a<function>(term)) {
                        res = res * term;
                        continue;
                }
                if (is_exactly_a<power>(term)) {
                        const power& p = ex_to<power>(term);
                        ex basis = p.op(0);
                        ex expo = p.op(1);
                        if (is_exactly_a<numeric>(expo)
                            and is_exactly_a<power>(basis)
                            and ex_to<numeric>(expo) == *_num_1_p) {
                                const power& bbasis = ex_to<power>(basis);
                                expo = bbasis.op(1) * _ex_1;
                                basis = bbasis.op(0);
                        }
                        auto search = factors.find(basis);
                        if (search == factors.end())
                                factors[basis] = expo;
                        else
                                search->second += expo;
                }
                else {
                        auto search = factors.find(term);
                        if (search == factors.end())
                                factors[term] = _ex1;
                        else
                                search->second += _ex1;
                }
        }

        for (auto& f : factors)
                res *= power(f.first, f.second);
        return res;
}

ex hypersimp(const ex& e, const ex& k)
// See Algorithm 2.1 in the Koepf reference
{
        ex f = e.expand();
        ex g = (f.subs(k == k-_ex_1)) / f;
        ex gf;
        bool red = factor(g, gf);
        ex& gr = gf;
        if (not red)
                gr = g;
        if (not has_suitable_form(gr))
                throw gosper_domain_error();

        return combine_powers(simplify_gamma(to_gamma(gr)));
}

static numeric get_root(const ex& f, const ex& sym)
{
        if (is_exactly_a<symbol>(f)) {
                if (f.is_equal(sym))
                        return *_num0_p;
                else
                        throw std::runtime_error("can't happen in get_root()");
        }
        if (is_exactly_a<power>(f))
                return get_root((ex_to<power>(f)).op(0), sym);
        numeric c, t;
        if (is_exactly_a<add>(f)) {
                const add& a = ex_to<add>(f);
                if (a.nops() > 2)
                        throw std::runtime_error("can't happen in get_root()");
                const ex& a0 = a.op(0);
                const ex& a1 = a.op(1);
                if (is_exactly_a<numeric>(a1) and is_exactly_a<mul>(a0)) {
                        c = ex_to<numeric>(a1);
                        mul m = ex_to<mul>(a0);
                        if (m.nops()>2 or not m.op(0).is_equal(sym)
                                        or not is_exactly_a<numeric>(m.op(1)))
                                throw std::runtime_error("can't happen in get_root()");
                        t = ex_to<numeric>(m.op(1));
                }
                else if (is_exactly_a<numeric>(a1) and is_exactly_a<symbol>(a0)) {
                        c = ex_to<numeric>(a1);
                        t = *_num1_p;
                }
                else
                        throw std::runtime_error("can't happen in get_root()");
        }
        else
                throw std::runtime_error("can't happen in get_root()");
        return -c/t;
}

static ex diagonal_poly(const exvector& syms, const ex& var)
// Return sum(0<=i<n, sym_i * var^i)
{
        unsigned int n = syms.size();
        ex res = _ex0;
        for (unsigned int i=0; i<n; ++i)
                res += power(var, i) * syms[n-i-1];
        return res;
}

static ex binomial_poly(const exvector& syms, const ex& var)
// Return sum(0<=i<n, sym_i * (var+1)^i) already expanded
{
        unsigned int n = syms.size();
        ex res = _ex0;
        for (unsigned int row=0; row<n; ++row) {
                const ex& v = power(var, row);
                for (unsigned int col=0; col<n; ++col) {
                        if (row+col < n)
                                res += v * binomial(row+col, col) * syms[n-col-row-1];
                }
        }
        return res;
}

static matrix solve_system(ex mpoly,
                const exvector& syms, const ex& msym)
// Solve mpoly==0 by converting to a linear system
{
        mpoly = mpoly.expand();
        if (not is_exactly_a<add>(mpoly))
                throw gosper_domain_error();
        ex_int_map sym_idx;
        const size_t nc = syms.size(), nr = mpoly.degree(msym) + 1;
        for (size_t i=0; i<nc; ++i)
                sym_idx[syms[i]] = i;
        matrix mat(nr, nc), vars(nc, 1), rhs(nr, 1);
        const symbol& main = ex_to<symbol>(msym);
        const add& a = ex_to<add>(mpoly);
        for (size_t i=0; i<a.nops(); ++i) {
                const ex& term = a.op(i);
                int s = -1;
                int e = 0;
                ex coeff(_ex0);
                if (is_exactly_a<mul>(term)) {
                        const mul& m = ex_to<mul>(term);
                        coeff = *_num1_p;
                        for (size_t j=0; j<m.nops(); ++j) {
                                const ex& mterm = m.op(j);
                                if (is_exactly_a<symbol>(mterm)) {
                                        const symbol& sy = ex_to<symbol>(mterm);
                                        if (sy.is_equal(main)) {
                                                e = 1;
                                                continue;
                                        }
                                        auto search = sym_idx.find(sy);
                                        if (search == sym_idx.end())
                                                coeff *= sy;
                                        else
                                                s = search->second;
                                        continue;
                                }
                                if (is_exactly_a<numeric>(mterm)) {
                                        coeff *= ex_to<numeric>(mterm);
                                        continue;
                                }
                                if (is_exactly_a<power>(mterm)) {
                                        const power& p = ex_to<power>(mterm);
                                        if (not p.op(0).is_equal(msym))
                                                throw std::runtime_error("can't happen in solve_system()");
                                        if (not is_exactly_a<numeric>(p.op(1)))
                                                throw std::runtime_error("can't happen in solve_system()");
                                        const numeric& expo = ex_to<numeric>(p.op(1));
                                        if (not expo.is_mpz() or expo < 0)
                                                throw std::runtime_error("can't happen in solve_system()");
                                        e = expo.to_int();
                                        continue;
                                }
                                throw std::runtime_error("can't happen in solve_system()");
                        }
                }
                else if (is_exactly_a<power>(term)) {
                        coeff = *_num1_p;
                        const power& p = ex_to<power>(term);
                        if (not p.op(0).is_equal(msym))
                                throw std::runtime_error("can't happen in solve_system()");
                        if (not is_exactly_a<numeric>(p.op(1)))
                                throw std::runtime_error("can't happen in solve_system()");
                        const numeric& expo = ex_to<numeric>(p.op(1));
                        if (not expo.is_mpz() or expo < 0)
                                throw std::runtime_error("can't happen in solve_system()");
                        e = expo.to_int();
                        }
                else if (is_exactly_a<symbol>(term)) {
                        const symbol& sy = ex_to<symbol>(term);
                        if (sy.is_equal(main)) {
                                e = 1;
                        }
                        else {
                                auto search = sym_idx.find(sy);
                                if (search == sym_idx.end())
                                        coeff = sy;
                                else {
                                        coeff = *_num1_p;
                                        s = search->second;
                                }
                        }
                }
                else if (is_exactly_a<numeric>(term)) {
                        coeff = ex_to<numeric>(term);
                }
                else
                        throw std::runtime_error("can't happen in solve_system()");
                if (s >= 0)
                        mat(e, s) = coeff;
                else
                        rhs(e, 0) = -coeff;
        }
        for (size_t i=0; i<nc; ++i)
                vars(i, 0) = syms[i];
        matrix res = mat.solve(vars, rhs, solve_algo::automatic);
        return res;
}

ex gosper_term(ex e, ex n)
{
        ex the_ex = hypersimp(e, n);
        ex num = the_ex.numer().expand();
        ex den = the_ex.denom().expand();
        ex cn = num.lcoeff(n);
        ex cd = den.lcoeff(n);
        numeric ldq = ex_to<numeric>(cn) / ex_to<numeric>(cd);
        ex A = num / cn;
        ex B = den / cd;
        symbol h;
        ex res = resultant(A, B.subs(n == n+h), n);
        ex prod;
        bool factored = factorpoly(res.expand(), prod);
        ex C = _ex1;
        if (factored and is_exactly_a<mul>(prod)) {
                const mul& m = ex_to<mul>(prod);
                std::set<int> roots;
                for (unsigned int i=0; i<m.nops(); i++) {
                        ex f = m.op(i);
                        if (is_exactly_a<numeric>(f))
                                continue;
                        if (is_exactly_a<numeric>(f.subs(h == 0))) {
                                numeric root = get_root(f, h);
                                if (not root.info(info_flags::integer)
                                    or root < 0)
                                        continue;
                                roots.insert(root.to_int());
                        }
                }
                for (int root : roots) {
                        ex d = gcd(A, B.subs(n == n+ex(root)).expand());
                        A = quo(A, d, n, false);
                        B = quo(B, d.subs(n == n-ex(root)), false);
                        for (long j=1; j<root+1; ++j)
                                C *= d.subs(n == n-ex(j));
                }
        }
        else {
                if (not is_exactly_a<numeric>(res)
                    and is_exactly_a<numeric>(res.subs(h == 0))) {
                        numeric root;
                        if (factored and is_exactly_a<power>(prod))
                                root = get_root(ex_to<power>(prod).op(0), h);
                        else
                                root = get_root(res, h);
                        if (root.info(info_flags::integer) and root >= 0) {
                                ex d = gcdpoly(A, B.subs(n == n+root).expand());
                                A = quo(A, d, n, false);
                                B = quo(B, d.subs(n == n-root), false);
                                for (long j=1; j<root.to_long()+1; ++j)
                                        C *= d.subs(n == n-ex(j));
                        }
                }
        }
        A *= ldq;
        B = B.subs(n == n-1);
        int N = A.degree(n);
        int M = B.degree(n);
        int K = C.degree(n);
        std::unordered_set<int> D;
        if (N != M or not A.lcoeff(n).is_equal(B.lcoeff(n)))
                D.insert(K - std::max(M,N));
        else if (N == 0) {
                D.insert(K - N + 1);
                D.insert(0);
        }
        else {
                D.insert(K - N + 1);
                ex t = (B.coeff(n,N-1) - A.coeff(n,N-1)) / A.lcoeff(n);
                if (is_exactly_a<numeric>(t)
                    and ex_to<numeric>(t).info(info_flags::integer)
                    and ex_to<numeric>(t) >= *_num0_p)
                                D.insert(ex_to<numeric>(t).to_int());
        }
        if (D.empty())
                throw gosper_domain_error();
        int d = *std::max_element(D.begin(), D.end());
        exvector syms;
        for (int i=0; i<d+1; ++i)
                syms.push_back((new symbol)->setflag(status_flags::dynallocated));
        ex xshifted = binomial_poly(syms, n);
        ex x = diagonal_poly(syms, n);
        ex H = A*xshifted - B*x -C;
        auto solution = solve_system(H, syms, n);
        for (size_t i=0; i<solution.rows(); ++i) {
                ex sym = syms[i];
                ex val = ex_to<numeric>(solution(i,0));
                x = x.subs(sym == val);
        }
        for (size_t i=0; i<syms.size(); ++i) {
                ex sym = syms[i];
                x = x.subs(sym == _ex0);
        }
        return B*x / C;
}

ex gosper_sum_definite(ex f, ex s, ex a, ex b, int* success)
{
        try {
                ex g = gosper_term(f, s);
                ex t = (f*(g + _ex1)).subs(s==b) - (f*g).expand().subs(s==a);
                *success = 1;
                ex res;
                bool changed = factor(t, res);
                if (changed)
                        return res;
                else
                        return t;
        }
        catch (gosper_domain_error) {
                *success = 0;
                return _ex0;
        }
}

ex gosper_sum_indefinite(ex f, ex s, int* success)
{
        try {
                ex t = f*gosper_term(f, s);
                *success = 1;
                ex res;
                bool changed = factor(t, res);
                if (changed)
                        return res;
                else
                        return t;
        }
        catch (gosper_domain_error) {
                *success = 0;
                return _ex0;
        }
}

}
