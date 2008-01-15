\\Load the modular equations. There are of two types; Canonical or Atkin.
\\For each prime ell the smallest in size was chosen

global(modular_eqn,modular_eqn_type);
/* modular_eqn is a vector of modular equations (canonical or Atkin type).
 *  The file sea_init.gp contains such eqns for all prime ell 2 <= ell <= MAXL
 * modular_eqn_type is a vector of types (canonical or Atkin) of the
 *  corresponding modular equation */
modular_eqn=read("sea_mod.gp");
modular_eqn_type =["C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "C", "A", "A", "A", "C", "C", "A", "C", "A", "A", "A", "C", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A", "A"];
\\We chose Atkin modular equation for the prime number ell when its size is
\\significantly smaller than the corresponding canonical modular equation

\\install(FpXQ_powers, "GLGG","bkinit");
\\bkinit(f, n, h, p) computes the list of powers L = [1, f, f^2,..., f^n] mod h in Fp[x]

\\install(FpX_FpXQV_compo, "GGGG","brent_kung");
\\brent_kung(g, L, h, p) computes f(g) mod h from g and L. To compute x^(p^r)
\\quickly from x^(p^(r-1))

install("FpXQ_pow", "GGGG");
\\FpXQ_pow(f, n, g, p ) computes f^n in Fp[x]/(g)
install("FpX_gcd", "GGG");
\\FpXQ_gcd(f, g) computes gcd(f,g) in Fp[x]
install("FpXQ_matrix_pow",GLLGG)
\\FpXQ_matrix_pow(a,n,m,T,p) compute the matrix of the map | Fp[x] -> Fp[x]/(T)
\\restricted to nxm entries                                | P(x)  -> P(a)
install("CM_CardEFp", GG)
\\CM_CardEFp(E, p) = #E(Fp) provided E has CM by a principal order, 0 otherwise
install("FpM_gauss", GGG)
install("FpM_ker", GG)
install("FpM_rank", lGG)
install("ZV_search",lGG)
install("ZV_sort_uniq",G)
install("Fp_sqrt", GG)
