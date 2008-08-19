/* emacs edit mode for this file is -*- C++ -*- */
/* $Id: fac_distrib.cc,v 1.11 2008/08/05 15:16:56 Singular Exp $ */

#include <config.h>

#include "assert.h"
#include "debug.h"

#include "cf_defs.h"
#include "canonicalform.h"
#include "cf_algorithm.h"
#include "cf_iter.h"
#include "fac_util.h"

bool
nonDivisors ( CanonicalForm omega, CanonicalForm delta, const CFArray & F, CFArray & d )
{
    DEBINCLEVEL( cerr, "nonDivisors" );
    CanonicalForm q, r;
    int k = F.size();
    d = CFArray( 0, k );
    d[0] = delta * omega;
    for ( int i = 1; i <= k; i++ )
    {
        q = abs(F[i]);
        for ( int j = i-1; j >= 0; j-- )
        {
            r = d[j];
            do
            {
                r = gcd( r, q );
                q = q / r;
            } while ( !r.isOne() );
            if ( q == 1 )
            {
                DEBDECLEVEL( cerr, "nonDivisors" );
                return false;
            }
        }
        d[i] = q;
    }
    DEBDECLEVEL( cerr, "nonDivisors" );
    return true;
}

bool
checkEvaluation ( const CanonicalForm & U, const CanonicalForm & lcU, const CanonicalForm & omega, const CFFList & F, const Evaluation & A, CanonicalForm & delta )
{
    CanonicalForm Vn, U0 = A( U );
    CFFListIterator I;
    int j;
    CFArray FF = CFArray( 1, F.length() );
    CFArray D;
    Vn = A( lcU );
    if ( Vn.isZero() )
        return false;
    delta = content( U0 );
    for ( I = F, j = 1; I.hasItem(); I++, j++ )
        FF[j] = A( I.getItem().factor() );
    return nonDivisors( omega, delta, FF, D );
}

bool
distributeLeadingCoeffs ( CanonicalForm & U, CFArray & G, CFArray & lcG, const CFFList & F, const CFArray & D, CanonicalForm & delta, CanonicalForm & omega, const Evaluation & A, int r )
{
    DEBINCLEVEL( cerr, "distributeLeadingCoeffs" );
    CanonicalForm ut, gt, d, ft;
    CanonicalForm dd;
    CFFListIterator I;
    int m, j, i;
    lcG = CFArray( 1, r );
    for ( j = 1; j <= r; j ++ )
        lcG[j] = 1;

    for ( I = F, i = 1; I.hasItem(); I++, i++ )
    {
        ft = I.getItem().factor();
        m = I.getItem().exp();
        DEBOUTLN( cerr, "trying to distribute " << ft );
        DEBOUTLN( cerr, "which is tested with " << D[i] );
        DEBOUTLN( cerr, "and contained to the power of " << m );
        j = 1;
        while ( m > 0 && j <= r )
        {
            ut = lc( G[j] );
            DEBOUTLN( cerr, "checking with " << ut );
            while ( m > 0 && fdivides( D[i], ut ) )
            {
                DEBOUTLN( cerr, "match found" );
                m--; ut /= D[i];
                lcG[j] *= ft;
            }
            j++;
        }
        if (m != 0)
        {
            DEBDECLEVEL( cerr, "distributeLeadingCoeffs" );
            return false;
        }
    }
    DEBOUTLN( cerr, "the leading coeffs before omega and delta correction: " << lcG );
    if ( !omega.isOne() )
    {
        for ( j = 1; j <= r; j++ )
        {
//            G[j] *= omega;
            lcG[j] *= omega;
            if(lc( G[j] ).isZero()) return false;
            G[j] = G[j] * ( A( lcG[j] ) / lc( G[j] ) );
        }
        U *= power( omega, r-1 );
    }
    if ( !delta.isOne() )
    {
        for ( j = 1; j <= r; j++ )
        {
            lcG[j] *= delta;
            if(lc( G[j] ).isZero()) return false;
            G[j] = G[j] * ( A( lcG[j] ) / lc( G[j] ) );
        }
        U *= power( delta, r );
    }
    DEBDECLEVEL( cerr, "distributeLeadingCoeffs" );
    return true;
}


static void
gfbAdjoin ( const CanonicalForm & F, CFList & L )
{
    if ( F.isOne() )
        return;
    if ( L.isEmpty() )
    {
        L.append( F );
        return;
    }
    CanonicalForm h, f = F;
    CFListIterator i, j;
    for ( i = L; i.hasItem() && ! f.isOne(); )
    {
        h = gcd( f, i.getItem() );
        if ( h.isOne() )
        {
            i++;
            continue;
        }
        while ( fdivides( h, f ) )
            f /= h;
        CFList D( h );
        gfbAdjoin( i.getItem() / h, D );
        for ( j = D; j.hasItem(); j++ )
            i.append( j.getItem() );
        i.remove( true );
    }
    if ( ! f.isOne() )
        L.append( f );
}


CFList
gcdFreeBasis ( const CFList L )
{
    CFListIterator i;
    CFList R;
    for ( i = L; i.hasItem(); i++ )
        gfbAdjoin( i.getItem(), R );
    return R;
}

bool
Univar2Bivar(const CanonicalForm & U, CFArray & G, const Evaluation & A, const modpk & bound, const Variable & x )
{
    CanonicalForm l = LC( U, Variable(1) );
    int n = G.size();
    CFArray lcG(1,n);
    for ( int i = 1; i <= n; i++ )
    {
        G[i] *= A(l)/lc(G[i]);
        lcG[i] = l;
    }
    return Hensel( U * power( l, n-1 ), G, lcG, A, bound, x );
}

bool
Hensel2 ( const CanonicalForm & U, CFArray & G, const Evaluation & A, const modpk & bound, const Variable & x )
{
    int i,n = G.size(); // n is number of factors of U
    CFArray TrueLcs(1, n);
    for (i=1; i <= n; i++)
        TrueLcs[i] = 1;
    Variable y;
    CanonicalForm lcU = LC ( U, Variable(1) );
    while (! lcU.inCoeffDomain())
    {
        y = lcU.mvar(); // should make a more intelligent choice
        CanonicalForm BivariateU = A ( U, 2, y.level() - 1 );
        CFArray BivariateFactors = G;
        CFArray lcFactors(1,n);
        Univar2Bivar(BivariateU, BivariateFactors, A, bound, y);
        for (i = 1; i <= n; i++)
        {
            BivariateFactors[i] /= content(BivariateFactors[i]);
            lcFactors[i] = LC(BivariateFactors[i],Variable(1));
        }
//        GFB = GcdFreeBasis(lcFactors); // this is not right... should probably make everything squarefree
//        Hensel2(lcU, GFB, A, bound, y);
    }
    for (i = 1; i <= n; i++)
        G[i] *= A(TrueLcs[i])/lc(G[i]);
    return Hensel(U, G, TrueLcs, A, bound, x);
}
