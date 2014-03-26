"""
Module of Supersingular Points

The module of divisors on the modular curve `X_0(N)` over `F_p` supported at supersingular points.

AUTHORS:

- William Stein

- David Kohel

- Iftikhar Burhanuddin

EXAMPLES::

    sage: x = SupersingularModule(389)
    sage: m = x.T(2).matrix()
    sage: a = m.change_ring(GF(97))
    sage: D = a.decomposition()
    sage: D[:3]
    [
    (Vector space of degree 33 and dimension 1 over Finite Field of size 97
    Basis matrix:
    [ 0  0  0  1 96 96  1 96 96  0  2 96 96  0  1  0  1  2 95  0  1  1  0  1  0 95  0 96 95  1 96  0  2], True),
    (Vector space of degree 33 and dimension 1 over Finite Field of size 97
    Basis matrix:
    [ 0  1 96 75 16 81 22 17 17  0  0 80 80  1 16 40 74  0  0 96 81 23 57 74  0  0  0 24  0 23 73  0  0], True),
    (Vector space of degree 33 and dimension 1 over Finite Field of size 97
    Basis matrix:
    [ 0  1 96 90 90  7  7  6 91  0  0 91  6 13  7  0 91  0  0 84 90  6  0  6  0  0  0 90  0 91  7  0  0], True)
    ]
    sage: len(D)
    9

We compute a Hecke operator on a space of huge dimension!::

    sage: X = SupersingularModule(next_prime(10000))
    sage: t = X.T(2).matrix()            # long time (21s on sage.math, 2011)
    sage: t.nrows()                      # long time
    835

TESTS::

    sage: X = SupersingularModule(389)
    sage: T = X.T(2).matrix().change_ring(QQ)
    sage: d = T.decomposition()
    sage: len(d)
    6
    sage: [a[0].dimension() for a in d]
    [1, 1, 2, 3, 6, 20]
    sage: loads(dumps(X)) == X
    True
    sage: loads(dumps(d)) == d
    True
"""

#*****************************************************************************
#       Copyright (C) 2004,2006 William Stein <wstein@gmail.com>
#       Copyright (C) 2006 David Kohel <kohel@maths.usyd.edu.au>
#       Copyright (C) 2006 Iftikhar Burhanuddin <burhanud@usc.edu>
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************

import math

import sage.modular.hecke.all as hecke
import sage.rings.all as rings
from sage.matrix.matrix_space import MatrixSpace
from sage.modular.arithgroup.all import Gamma0
from sage.databases.db_class_polynomials import HilbertClassPolynomialDatabase
from sage.databases.db_modular_polynomials \
     import ClassicalModularPolynomialDatabase

from sage.misc.misc import verbose


def Phi2_quad(J3, ssJ1, ssJ2):
    r"""
    This function returns a certain quadratic polynomial over a finite
    field in indeterminate J3.

    The roots of the polynomial along with ssJ1 are the
    neighboring/2-isogenous supersingular j-invariants of ssJ2.

    INPUT:

    - ``J3`` -- indeterminate of a univariate polynomial ring defined over a finite
      field with p^2 elements where p is a prime number

    - ``ssJ2``, ``ssJ2`` -- supersingular j-invariants over the finite field

    OUTPUT:

    - polynomial -- defined over the finite field

    EXAMPLES:
    The following code snippet produces a factor of the modular polynomial
    `\Phi_{2}(x,j_{in})`, where `j_{in}` is a supersingular j-invariant
    defined over the finite field with `37^2` elements::

        sage: F = GF(37^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: poly = sage.modular.ssmod.ssmod.Phi_polys(2,X,j_in)
        sage: poly.roots()
        [(8, 1), (27*a + 23, 1), (10*a + 20, 1)]
        sage: sage.modular.ssmod.ssmod.Phi2_quad(X, F(8), j_in)
        x^2 + 31*x + 31

    .. note::

        Given a root (j1,j2) to the polynomial `Phi_2(J1,J2)`, the pairs
        (j2,j3) not equal to (j2,j1) which solve `Phi_2(j2,j3)` are roots of
        the quadratic equation:

        J3^2 + (-j2^2 + 1488*j2 + (j1 - 162000))*J3 + (-j1 + 1488)*j2^2 +
        (1488*j1 + 40773375)*j2 + j1^2 - 162000*j1 + 8748000000

        This will be of use to extend the 2-isogeny graph, once the initial
        three roots are determined for `Phi_2(J1,J2)`.

    AUTHORS:

    - David Kohel -- kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin -- burhanud@usc.edu
    """
    ssJ1_pow2 = ssJ1**2
    ssJ2_pow2 = ssJ2**2

    return J3.parent()([(-ssJ1 + 1488)*ssJ2_pow2+ (1488*ssJ1 +
    40773375)*ssJ2 + ssJ1_pow2 - 162000*ssJ1 + 8748000000,
    -ssJ2_pow2 + 1488*ssJ2 + (ssJ1 - 162000),
    1])



def Phi_polys(L, x, j):
    r"""
    This function returns a certain polynomial of degree `L+1` in the
    indeterminate x over a finite field.

    The roots of the **modular** polynomial `\Phi(L, x, j)` are the
    `L`-isogenous supersingular j-invariants of j.

    INPUT:

    - ``L`` -- integer, either 2,3,5,7 or 11

    - ``x`` -- indeterminate of a univariate polynomial ring defined over a
      finite field with p^2 elements, where p is a prime number

    - ``j`` -- supersingular j-invariant over the finite field

    OUTPUT:

    - polynomial -- defined over the finite field

    EXAMPLES:
    The following code snippet produces the modular polynomial
    `\Phi_{L}(x,j_{in})`, where `j_{in}` is a supersingular j-invariant
    defined over the finite field with `7^2` elements::

        sage: F = GF(7^2, 'a')
        sage: X = PolynomialRing(F, 'x').gen()
        sage: j_in = supersingular_j(F)
        sage: sage.modular.ssmod.ssmod.Phi_polys(2,X,j_in)
        x^3 + 3*x^2 + 3*x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(3,X,j_in)
        x^4 + 4*x^3 + 6*x^2 + 4*x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(5,X,j_in)
        x^6 + 6*x^5 + x^4 + 6*x^3 + x^2 + 6*x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(7,X,j_in)
        x^8 + x^7 + x + 1
        sage: sage.modular.ssmod.ssmod.Phi_polys(11,X,j_in)
        x^12 + 5*x^11 + 3*x^10 + 3*x^9 + 5*x^8 + x^7 + x^5 + 5*x^4 + 3*x^3 + 3*x^2 + 5*x + 1

    AUTHORS:

    - William Stein -- wstein@gmail.com

    - David Kohel -- kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin -- burhanud@usc.edu
    """
    if not(L in [2,3,5,7,11]):
        raise ValueError, "L should be either 2,3,5,7 or 11. For other values use ClassicalModularPolynomialDatabase()."

    j_tmp = 1
    j_pow = [j_tmp]
    #degree of polynomial = L+1
    for i in range(int(L+2)):
        j_tmp = j_tmp*j
        j_pow += [j_tmp]

    if L == 2:
        return x.parent()([j_pow[3] - 162000*j_pow[2] + 8748000000*j_pow[1] -
                           157464000000000,
                           1488*j_pow[2] + 40773375*j_pow[1] + 8748000000,
                           - (j_pow[2] - 1488*j_pow[1] + 162000),
                           1])
    elif L == 3:
        return x.parent()([1855425871872000000000*j_pow[1] +
                           452984832000000*j_pow[2] + 36864000*j_pow[3] + j_pow[4],
                           1855425871872000000000 - 770845966336000000*j_pow[1] +
                           8900222976000*j_pow[2] - 1069956*j_pow[3],
                           452984832000000 + 8900222976000*j_pow[1] + 2587918086*j_pow[2] + 2232*j_pow[3],
                           36864000 - 1069956*j_pow[1] + 2232*j_pow[2] - j_pow[3],
                           1])
    elif L == 5:
        return x.parent()([
            141359947154721358697753474691071362751004672000 +
            53274330803424425450420160273356509151232000*j_pow[1] +
            6692500042627997708487149415015068467200*j_pow[2] +
            280244777828439527804321565297868800*j_pow[3] +
            1284733132841424456253440*j_pow[4] + 1963211489280*j_pow[5] +
            j_pow[6],
            53274330803424425450420160273356509151232000 -
            264073457076620596259715790247978782949376*j_pow[1] +
            36554736583949629295706472332656640000*j_pow[2] -
            192457934618928299655108231168000*j_pow[3] +
            128541798906828816384000*j_pow[4] - 246683410950*j_pow[5],
            6692500042627997708487149415015068467200 +
            36554736583949629295706472332656640000*j_pow[1] +
            5110941777552418083110765199360000*j_pow[2] +
            26898488858380731577417728000*j_pow[3] +
            383083609779811215375*j_pow[4] + 2028551200*j_pow[5],
            280244777828439527804321565297868800 -
            192457934618928299655108231168000*j_pow[1] +
            26898488858380731577417728000*j_pow[2] -
            441206965512914835246100*j_pow[3] +
            107878928185336800*j_pow[4] - 4550940*j_pow[5],
            1284733132841424456253440 + 128541798906828816384000*j_pow[1]
            + 383083609779811215375*j_pow[2] + 107878928185336800*j_pow[3]
            + 1665999364600*j_pow[4] + 3720*j_pow[5],
            1963211489280 - 246683410950*j_pow[1] + 2028551200*j_pow[2] - 4550940*j_pow[3]
            + 3720*j_pow[4] - j_pow[5],
            1])
    elif L == 7:
        return x.parent()([1464765079488386840337633731737402825128271675392000000000000000000*j_pow[2]
            +    13483958224762213714698012883865296529472356352000000000000000*j_pow[3]
            +    41375720005635744770247248526572116368162816000000000000*j_pow[4]
            +    42320664241971721884753245384947305283584000000000*j_pow[5] +
            3643255017844740441130401792000000*j_pow[6] +
            104545516658688000*j_pow[7] + j_pow[8],
            1221349308261453750252370983314569119494710493184000000000000000000*j_pow[1]
            -    838538082798149465723818021032241603179964268544000000000000000*j_pow[2]
            -    129686683986501811181602978946723823397619367936000000000000*j_pow[3]
            +    553293497305121712634517214392820316998991872000000000*j_pow[4]
            -    40689839325168186578698294668599003971584000000*j_pow[5] +
            1038063543615451121419229773824000*j_pow[6] -
            34993297342013192*j_pow[7],
            1464765079488386840337633731737402825128271675392000000000000000000
            -    838538082798149465723818021032241603179964268544000000000000000*j_pow[1]
            -    46666007311089950798495647194817495401448341504000000000000*j_pow[2]
            +    72269669689202948469186346100000679630099972096000000000*j_pow[3]
            +    308718989330868920558541707287296140145328128000000*j_pow[4] +
            11269804827778129625111322263056523132928000*j_pow[5] +
            10685207605419433304631062899228*j_pow[6] +
            720168419610864*j_pow[7],
            13483958224762213714698012883865296529472356352000000000000000 -
            129686683986501811181602978946723823397619367936000000000000*j_pow[1]
            +    72269669689202948469186346100000679630099972096000000000*j_pow[2]
            -    5397554444336630396660447092290576395211374592000000*j_pow[3] +
            17972351380696034759035751584170427941396480000*j_pow[4] -
            901645312135695263877115693740562092344*j_pow[5] +
            16125487429368412743622133040*j_pow[6] - 4079701128594*j_pow[7],
            41375720005635744770247248526572116368162816000000000000 +
            553293497305121712634517214392820316998991872000000000*j_pow[1] +
            308718989330868920558541707287296140145328128000000*j_pow[2] +
            17972351380696034759035751584170427941396480000*j_pow[3] +
            88037255060655710247136461896264828390470*j_pow[4] +
            14066810691825882583305340438456800*j_pow[5] +
            4460942463213898353207432*j_pow[6] + 9437674400*j_pow[7],
            42320664241971721884753245384947305283584000000000 -
            40689839325168186578698294668599003971584000000*j_pow[1] +
            11269804827778129625111322263056523132928000*j_pow[2] -
            901645312135695263877115693740562092344*j_pow[3] +
            14066810691825882583305340438456800*j_pow[4] -
            18300817137706889881369818348*j_pow[5] +
            177089350028475373552*j_pow[6] - 10246068*j_pow[7],
            3643255017844740441130401792000000 +
            1038063543615451121419229773824000*j_pow[1] +
            10685207605419433304631062899228*j_pow[2] +
            16125487429368412743622133040*j_pow[3] +
            4460942463213898353207432*j_pow[4] +
            177089350028475373552*j_pow[5] + 312598931380281*j_pow[6] +
            5208*j_pow[7],
            104545516658688000 - 34993297342013192*j_pow[1] +
            720168419610864*j_pow[2] - 4079701128594*j_pow[3] +
            9437674400*j_pow[4] - 10246068*j_pow[5] + 5208*j_pow[6] -
            j_pow[7],
            1])
    elif L == 11:
        return x.parent()([3924233450945276549086964624087200490995247233706746270899364206426701740619416867392454656000000000000000000000000000000000000
            -    3708476896661234261166595138586620846782660237574536888784393380944856551532392652692520960000000000000000000000000000000000*j_pow[1]
            +    1509199706449264373105244249368970977209959173066491449939153900434037998316228131684352000000000000000000000000000000000*j_pow[2]
            -    337500037290942764495395868386562971754016116785390841072048221617443316658082155384012800000000000000000000000000000*j_pow[3]
            +    43714682637171236021367604966833305309923746974850894665325331604362303109715777067941888000000000000000000000000*j_pow[4]
            -    3111357148902865912417988391836350251682805385917571877568422664218078901010004935966720000000000000000000000*j_pow[5]
            +    95356266594731795079493309965756674711058734831164489212811553129058773080352804044800000000000000000000*j_pow[6]
            +    618840723107761889896363016885251574078635388443306832549992828319945330157158400000000000000000*j_pow[7]
            +    1338586400912357073420399795635643400599836918986297982928179335149920452608000000000000*j_pow[8]
            +    965122546660349298406724063940884252743873633176129290337528305418240000000000*j_pow[9]
            +    29298331981110197366602526090413106879319244800000000*j_pow[10]
            +    296470902355240575283200000*j_pow[11]
            +    j_pow[12],
            -    3708476896661234261166595138586620846782660237574536888784393380944856551532392652692520960000000000000000000000000000000000
            +    6950986496704390042399105433049126860396103535300642728895074819467726754375236055025582080000000000000000000000000000000*j_pow[1]
            -    4175190947377089941611452135383204997172948465221368432119554418845446929655566146994176000000000000000000000000000000*j_pow[2]
            +    493751729222149651035457063068642305508233453469401395944974296438196687728770695603159040000000000000000000000000*j_pow[3]
            +    59659609577030961637541110289112021078091104767187787822549078869394205439302452893450240000000000000000000000*j_pow[4]
            -    7840379248214196729643062796493269425081859930100141304047932909346022483171510017064960000000000000000000*j_pow[5]
            -    95333447356443287210404497374050404132491763274506548619337189691919811046970438451200000000000000000*j_pow[6]
            -    24155957253764418975307742823129586187061243620756339515602571075061236992294518784000000000000*j_pow[7]
            +    66806304467998310581793391194791115184805127528413091235284315294143736709120000000000*j_pow[8]
            -    1458178254597295207839980786768623018650234306932331393013634952069120000000*j_pow[9]
            +    33446467926379842030532687838341039552110187929600000*j_pow[10]
            -    374642006356701393515817612*j_pow[11],
                 1509199706449264373105244249368970977209959173066491449939153900434037998316228131684352000000000000000000000000000000000
            -    4175190947377089941611452135383204997172948465221368432119554418845446929655566146994176000000000000000000000000000000*j_pow[1]
            -    301851634381591833346238394387907563828793379391119445614595161272769455527698270716428288000000000000000000000000*j_pow[2]
            +    1038677201789914991362090465961377302769147065985487222285672689158918175716097236444119040000000000000000000000*j_pow[3]
            +    378494977797549959360178068152933818044335078157093771639955480261351930169113765048483840000000000000000000*j_pow[4]
            +    9718148718139346647384449201643833517488848029697396574289278515913329360524510494720000000000000000000*j_pow[5]
            +    30494044246550310117871895628421273379173050630568397072391110688366558535804457582592000000000000*j_pow[6]
            +    44681231489418997440503069818655052635806384532381152777755381649015689662976491520000000000*j_pow[7]
            +    171790435018380416903247878610824648919543398246401012395341432490921925017600000000*j_pow[8]
            +    804436418307995738740132598166893365099468842089705900525050627891200000*j_pow[9]
            +    1587728122949690904187089204116332301200302760915266*j_pow[10]
            +    27209811658056645815522600*j_pow[11],
            -    337500037290942764495395868386562971754016116785390841072048221617443316658082155384012800000000000000000000000000000
            +    493751729222149651035457063068642305508233453469401395944974296438196687728770695603159040000000000000000000000000*j_pow[1]
            +    1038677201789914991362090465961377302769147065985487222285672689158918175716097236444119040000000000000000000000*j_pow[2]
            -    925461466455522523607980072366478440235575959511945288268604770825451300845059605937520640000000000000000000*j_pow[3]
            -    51038778870467375317174627414281203016789153392265449880353463871004348816411677478092800000000000000000*j_pow[4]
            -    1328993907465108152135763886999825071444084099881098607565574716140191426369978927939584000000000000*j_pow[5]
            -    7211912299746007510535159486199919697482960389278446632552985263875183091897870581760000000000*j_pow[6]
            -    22093249696627933419655226823604057638897222562682635800269909178325710985117040640000000*j_pow[7]
            +    79513247125057906492841989395207442300133781750924860449090230806481243648000000*j_pow[8]
            -    199188452917764242987050083089364860927274115441197382331866126825820*j_pow[9]
            +    14131378888778142661582693947549844785863493325800*j_pow[10]
            -    529134841844639613861795*j_pow[11],
                 43714682637171236021367604966833305309923746974850894665325331604362303109715777067941888000000000000000000000000
            +    59659609577030961637541110289112021078091104767187787822549078869394205439302452893450240000000000000000000000*j_pow[1]
            +    378494977797549959360178068152933818044335078157093771639955480261351930169113765048483840000000000000000000*j_pow[2]
            -    51038778870467375317174627414281203016789153392265449880353463871004348816411677478092800000000000000000*j_pow[3]
            +    15043423165563966645618284609730360176005265392518745580151910727157028699006028388237312000000000000*j_pow[4]
            -    177994641867075262695184980920462608604060357466681128822395417442867019643767352197120000000000*j_pow[5]
            +    1938738373821740121470446368665797412833082873875468530371642913339302678999680942080000000*j_pow[6]
            +    2973119672716212219456471881112888569835575578534065127175856819648732682854604800000*j_pow[7]
            +    8498500708725193890718329655230574962816784139443636591086906768989729050095*j_pow[8]
            +    22148485195925584385790489089697473918894904664093860668378292000*j_pow[9]
            +    35372414460361796790312007060191890803134127320*j_pow[10]
            +    4297837238774928467520*j_pow[11],
            -    3111357148902865912417988391836350251682805385917571877568422664218078901010004935966720000000000000000000000
            -    7840379248214196729643062796493269425081859930100141304047932909346022483171510017064960000000000000000000*j_pow[1]
            +    9718148718139346647384449201643833517488848029697396574289278515913329360524510494720000000000000000000*j_pow[2]
            -    1328993907465108152135763886999825071444084099881098607565574716140191426369978927939584000000000000*j_pow[3]
            -    177994641867075262695184980920462608604060357466681128822395417442867019643767352197120000000000*j_pow[4]
            -    15057297311708922526580514410563848478334693758624999774108600968667487260827388477440000000*j_pow[5]
            +    224080399886627495149771654692369177094059649940825305182078225594292057242702643200000*j_pow[6]
            -    75948585201267973403627533631138995089882647284307484579413691458563029509971992*j_pow[7]
            +    208334210762751500564946204497082337222910461284651050215872586641463200*j_pow[8]
            -    994774826102691960922410649494629085486856242714439003812180*j_pow[9]
            +    28890545335855949285086003898461917345026160*j_pow[10]
            -    17899526272883039048*j_pow[11],
                 95356266594731795079493309965756674711058734831164489212811553129058773080352804044800000000000000000000
            -    95333447356443287210404497374050404132491763274506548619337189691919811046970438451200000000000000000*j_pow[1]
            +    30494044246550310117871895628421273379173050630568397072391110688366558535804457582592000000000000*j_pow[2]
            -    7211912299746007510535159486199919697482960389278446632552985263875183091897870581760000000000*j_pow[3]
            +    1938738373821740121470446368665797412833082873875468530371642913339302678999680942080000000*j_pow[4]
            +    224080399886627495149771654692369177094059649940825305182078225594292057242702643200000*j_pow[5]
            +    1168150167526575837857761510359647773943258089269992605255478096499695783789300124*j_pow[6]
            +    247900233561939294388612799857476424364856251769094880288086537904279396400*j_pow[7]
            +    987807801334019988631500819088661487281712947788833523552559299560*j_pow[8]
            +    14690460927260804690751501000083244161647396386205851440*j_pow[9]
            +    7848482999227584325448694633580010490867*j_pow[10]
            +    42570393135641712*j_pow[11],
                 618840723107761889896363016885251574078635388443306832549992828319945330157158400000000000000000
            -    24155957253764418975307742823129586187061243620756339515602571075061236992294518784000000000000*j_pow[1]
            +    44681231489418997440503069818655052635806384532381152777755381649015689662976491520000000000*j_pow[2]
            -    22093249696627933419655226823604057638897222562682635800269909178325710985117040640000000*j_pow[3]
            +    2973119672716212219456471881112888569835575578534065127175856819648732682854604800000*j_pow[4]
            -    75948585201267973403627533631138995089882647284307484579413691458563029509971992*j_pow[5]
            +    247900233561939294388612799857476424364856251769094880288086537904279396400*j_pow[6]
            -    64999046469909490143435875140651300541119093852394968074094803537810*j_pow[7]
            +    636861023141767565580039581191818069063579259290464688398880*j_pow[8]
            -    51135193038502008150804190472844550800569441050500*j_pow[9]
            +    645470833566425875717489618904152240*j_pow[10]
            -    61058988656490*j_pow[11],
                 1338586400912357073420399795635643400599836918986297982928179335149920452608000000000000
            +    66806304467998310581793391194791115184805127528413091235284315294143736709120000000000*j_pow[1]
            +    171790435018380416903247878610824648919543398246401012395341432490921925017600000000*j_pow[2]
            +    79513247125057906492841989395207442300133781750924860449090230806481243648000000*j_pow[3]
            +    8498500708725193890718329655230574962816784139443636591086906768989729050095*j_pow[4]
            +    208334210762751500564946204497082337222910461284651050215872586641463200*j_pow[5]
            +    987807801334019988631500819088661487281712947788833523552559299560*j_pow[6]
            +    636861023141767565580039581191818069063579259290464688398880*j_pow[7]
            +    29211180544704743418963619709378403797452606969172658*j_pow[8]
            +    24228593349948582884094197811518266845689352*j_pow[9]
            +    12407796387712093514736413264496*j_pow[10]
            +    53686822816*j_pow[11],
                 965122546660349298406724063940884252743873633176129290337528305418240000000000
            -    1458178254597295207839980786768623018650234306932331393013634952069120000000*j_pow[1]
            +    804436418307995738740132598166893365099468842089705900525050627891200000*j_pow[2]
            -    199188452917764242987050083089364860927274115441197382331866126825820*j_pow[3]
            +    22148485195925584385790489089697473918894904664093860668378292000*j_pow[4]
            -    994774826102691960922410649494629085486856242714439003812180*j_pow[5]
            +    14690460927260804690751501000083244161647396386205851440*j_pow[6]
            -    51135193038502008150804190472844550800569441050500*j_pow[7]
            +    24228593349948582884094197811518266845689352*j_pow[8]
            -    573388748843683532691009051194955437*j_pow[9]
            +    30134971854812981978547264*j_pow[10] - 28278756*j_pow[11],
                 29298331981110197366602526090413106879319244800000000
            +    33446467926379842030532687838341039552110187929600000*j_pow[1]
            +    1587728122949690904187089204116332301200302760915266*j_pow[2]
            +    14131378888778142661582693947549844785863493325800*j_pow[3]
            +    35372414460361796790312007060191890803134127320*j_pow[4]
            +    28890545335855949285086003898461917345026160*j_pow[5]
            +    7848482999227584325448694633580010490867*j_pow[6]
            +    645470833566425875717489618904152240*j_pow[7]
            +    12407796387712093514736413264496*j_pow[8]
            +    30134971854812981978547264*j_pow[9]
            +    1608331026427734378*j_pow[10]
            +    8184*j_pow[11],
                 296470902355240575283200000
            -    374642006356701393515817612*j_pow[1]
            +    27209811658056645815522600*j_pow[2]
            -    529134841844639613861795*j_pow[3]
            +    4297837238774928467520*j_pow[4]
            -    17899526272883039048*j_pow[5]
            +    42570393135641712*j_pow[6]
            -    61058988656490*j_pow[7]
            +    53686822816*j_pow[8]
            -    28278756*j_pow[9]
            +    8184*j_pow[10]
            -    j_pow[11],
                 1])


def dimension_supersingular_module(prime, level=1):
    r"""
    This function returns the dimension of the Supersingular module, which is
    equal to the dimension of the space of modular forms of weight `2`
    and conductor equal to prime times level.

    INPUT:

    - ``prime`` -- integer, prime

    - ``level`` -- integer, positive

    OUTPUT:
       dimension -- integer, nonnegative

    EXAMPLES:
    The code below computes the dimensions of Supersingular modules
    with level=1 and prime = 7, 15073 and 83401::

        sage: dimension_supersingular_module(7)
        1

        sage: dimension_supersingular_module(15073)
        1256

        sage: dimension_supersingular_module(83401)
        6950

    NOTES:
    The case of level > 1 has not been implemented yet.

    AUTHORS:

    - David Kohel -- kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin - burhanud@usc.edu
    """
    if not(rings.Integer(prime).is_prime()):
        raise ValueError, "%s is not a prime"%prime

    if level == 1:
        return Gamma0(prime).dimension_modular_forms(2)

    #list of genus(X_0(level)) equal to zero
    #elif (level in [ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 13, 16, 18, 25]):
    #compute basis

    else:
        raise NotImplementedError


def supersingular_D(prime):
    r"""
    This function returns a fundamental discriminant `D` of an
    imaginary quadratic field, where the given prime does not split
    (see Silverman's Advanced Topics in the Arithmetic of Elliptic
    Curves, page 184, exercise 2.30(d).)

    INPUT:

    - prime -- integer, prime

    OUTPUT:
        D -- integer, negative

    EXAMPLES:

    These examples return *supersingular discriminants* for 7,
    15073 and 83401::

        sage: supersingular_D(7)
        -4

        sage: supersingular_D(15073)
        -15

        sage: supersingular_D(83401)
        -7

    AUTHORS:

    - David Kohel - kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin - burhanud@usc.edu
    """
    if not(rings.Integer(prime).is_prime()):
        raise ValueError, "%s is not a prime"%prime

    #Making picking D more intelligent
    D = -1
    while True:
        Dmod4 = rings.Mod(D,4)
        if Dmod4 in (0,1) and (rings.kronecker(D,prime) != 1):
            return D
        D = D - 1

def supersingular_j(FF):
    r"""
    This function returns a supersingular j-invariant over the finite
    field FF.

    INPUT:

    - ``FF``  -- finite field with p^2 elements, where p is a prime number

    OUTPUT:
       finite field element -- a supersingular j-invariant
       defined over the finite field FF

    EXAMPLES:

    The following examples calculate supersingular j-invariants for a
    few finite fields::

        sage: supersingular_j(GF(7^2, 'a'))
        6

    Observe that in this example the j-invariant is not defined over
    the prime field::

        sage: supersingular_j(GF(15073^2,'a'))  # optional - database
        10630*a + 6033

        sage: supersingular_j(GF(83401^2, 'a'))
        67977

    AUTHORS:

    - David Kohel -- kohel@maths.usyd.edu.au

    - Iftikhar Burhanuddin -- burhanud@usc.edu
    """
    if not(FF.is_field()) or not(FF.is_finite()):
        raise ValueError, "%s is not a finite field"%FF
    prime = FF.characteristic()
    if not(rings.Integer(prime).is_prime()):
        raise ValueError, "%s is not a prime"%prime
    if not(rings.Integer(FF.cardinality())) == rings.Integer(prime**2):
        raise ValueError, "%s is not a quadratic extension"%FF
    if rings.kronecker(-1, prime) != 1:
        j_invss = 1728                 #(2^2 * 3)^3
    elif rings.kronecker(-2, prime) != 1:
        j_invss = 8000                 #(2^2 * 5)^3
    elif rings.kronecker(-3, prime) != 1:
        j_invss = 0                    #0^3
    elif rings.kronecker(-7, prime) != 1:
        j_invss = 16581375             #(3 * 5 * 17)^3
    elif rings.kronecker(-11, prime) != 1:
        j_invss = -32768               #-(2^5)^3
    elif rings.kronecker(-19, prime) != 1:
        j_invss = -884736              #-(2^5 * 3)^3
    elif rings.kronecker(-43, prime) != 1:
        j_invss = -884736000           #-(2^6 * 3 * 5)^3
    elif rings.kronecker(-67, prime) != 1:
        j_invss = -147197952000        #-(2^5 * 3 * 5 * 11)^3
    elif rings.kronecker(-163, prime) != 1:
        j_invss = -262537412640768000  #-(2^6 * 3 * 5 * 23 * 29)^3
    else:
        D = supersingular_D(prime)
        DBCP = HilbertClassPolynomialDatabase()
        hc_poly = FF['x'](DBCP[D])
        root_hc_poly_list = list(hc_poly.roots())

        j_invss = root_hc_poly_list[0][0]
    return FF(j_invss)

class SupersingularModule(hecke.HeckeModule_free_module):
    r"""
    The module of supersingular points in a given characteristic, with
    given level structure.

    The characteristic must not divide the level.

    NOTE: Currently, only level 1 is implemented.

    EXAMPLES::

        sage: S = SupersingularModule(17)
        sage: S
        Module of supersingular points on X_0(1)/F_17 over Integer Ring
        sage: S = SupersingularModule(16)
        Traceback (most recent call last):
        ...
        ValueError: the argument prime must be a prime number
        sage: S = SupersingularModule(prime=17, level=34)
        Traceback (most recent call last):
        ...
        ValueError: the argument level must be coprime to the argument prime
        sage: S = SupersingularModule(prime=17, level=5)
        Traceback (most recent call last):
        ...
        NotImplementedError: supersingular modules of level > 1 not yet implemented
    """
    def __init__(self, prime=2, level=1, base_ring=rings.IntegerRing()):
        r"""
        Create a supersingular module.

        EXAMPLE::

            sage: SupersingularModule(3)
            Module of supersingular points on X_0(1)/F_3 over Integer Ring
        """

        if not prime.is_prime():
            raise ValueError, "the argument prime must be a prime number"
        if prime.divides(level):
            raise ValueError, "the argument level must be coprime to the argument prime"
        if level != 1:
            raise NotImplementedError, "supersingular modules of level > 1 not yet implemented"
        self.__prime = prime
        from sage.rings.all import FiniteField
        self.__finite_field = FiniteField(prime**2,'a')
        self.__level = level
        self.__hecke_matrices = {}
        hecke.HeckeModule_free_module.__init__(
            self, base_ring, prime*level, weight=2)

    def _repr_(self):
        """
        String representation of self

        EXAMPLES::

            sage: SupersingularModule(11)._repr_()
            'Module of supersingular points on X_0(1)/F_11 over Integer Ring'
        """

        return "Module of supersingular points on X_0(%s)/F_%s over %s"%(
            self.__level, self.__prime, self.base_ring())

    def __cmp__(self, other):
        r"""
        Compare self to other.

        EXAMPLES::

            sage: SupersingularModule(37) == ModularForms(37, 2)
            False
            sage: SupersingularModule(37) == SupersingularModule(37, base_ring=Qp(7))
            False
            sage: SupersingularModule(37) == SupersingularModule(37)
            True
        """
        if not isinstance(other, SupersingularModule):
            return cmp(type(self), type(other))
        else:
            return cmp( (self.__level, self.__prime, self.base_ring()), (other.__level, other.__prime, other.base_ring()))

    def free_module(self):
        """
        EXAMPLES::

            sage: X = SupersingularModule(37)
            sage: X.free_module()
            Ambient free module of rank 3 over the principal ideal domain Integer Ring

        This illustrates the fix at :trac:`4306`::

            sage: X = SupersingularModule(389)
            sage: X.basis()
            Traceback (most recent call last):
            ...
            NotImplementedError
        """
        return rings.ZZ**self.dimension()

    def dimension(self):
        r"""
        Return the dimension of the space of modular forms of weight 2
        and level equal to the level associated to self.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT:
            integer -- dimension, nonnegative

        EXAMPLES::

            sage: S = SupersingularModule(7)
            sage: S.dimension()
            1

            sage: S = SupersingularModule(15073)
            sage: S.dimension()
            1256

            sage: S = SupersingularModule(83401)
            sage: S.dimension()
            6950

        NOTES:
           The case of level > 1 has not yet been implemented.

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        try:
            return self.__dimension
        except AttributeError:
            pass
        if self.__level == 1:
            G = Gamma0(self.__prime)
            self.__dimension = G.dimension_modular_forms(2)
        else:
            raise NotImplementedError
        return self.__dimension

    rank = dimension

    def level(self):
        r"""
        This function returns the level associated to self.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT:
            integer -- the level, positive

        EXAMPLES::

            sage: S = SupersingularModule(15073)
            sage: S.level()
            1

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        return self.__level

    def prime(self):
        r"""
        This function returns the characteristic of the finite field
        associated to self.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT:

        - integer -- characteristic, positive

        EXAMPLES::

            sage: S = SupersingularModule(19)
            sage: S.prime()
            19

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        return self.__prime

    def weight(self):
        r"""
        This function returns the weight associated to self.

        INPUT:

        - ``self`` -- SupersingularModule object

        OUTPUT:
            integer -- weight, positive

        EXAMPLES::

            sage: S = SupersingularModule(19)
            sage: S.weight()
            2

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        return 2

    def supersingular_points(self):
        r"""
        This function computes the supersingular j-invariants over the
        finite field associated to self.

        INPUT:

        -  ``self`` -- SupersingularModule object

        OUTPUT: list_j, dict_j -- list_j is the list of supersingular
            j-invariants, dict_j is a dictionary with these
            j-invariants as keys and their indexes as values. The
            latter is used to speed up j-invariant look-up. The
            indexes are based on the order of their *discovery*.

        EXAMPLES:

        The following examples calculate supersingular j-invariants
        over finite fields with characteristic 7, 11 and 37::

            sage: S = SupersingularModule(7)
            sage: S.supersingular_points()
            ([6], {6: 0})

            sage: S = SupersingularModule(11)
            sage: S.supersingular_points()
            ([1, 0], {0: 1, 1: 0})

            sage: S = SupersingularModule(37)
            sage: S.supersingular_points()
            ([8, 27*a + 23, 10*a + 20], {8: 0, 10*a + 20: 2, 27*a + 23: 1})

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        try:
            return (self._ss_points_dic, self._ss_points)
        except AttributeError:
            pass
        Fp2 = self.__finite_field
        level = self.__level
        prime = Fp2.characteristic()
        X = Fp2['x'].gen()
        jinv = supersingular_j(Fp2)

        dim = dimension_supersingular_module(prime, level)

        pos = int(0)
        #using list to keep track of explored nodes using pos
        ss_points = [jinv]

        #using  to keep track of index of the previous node
        ss_points_pre = [-1]

        #using dictionary for fast j-invariant look-up
        ss_points_dic = {jinv:pos}

        T2_matrix = MatrixSpace(rings.Integers(), dim, sparse=True)(0)

        while pos < len(ss_points):
            if pos == 0:
                neighbors = Phi_polys(2,X,ss_points[pos]).roots()
            else:
                j_prev = ss_points_pre[pos]
                # TODO: These are quadratic polynomials -- maybe we should use the
                # quadratic formula and fast square root finding (??)
                neighbors = Phi2_quad(X, ss_points[j_prev], ss_points[pos]).roots()

            for (xj,ej) in neighbors:
                if xj not in ss_points_dic:
                    j = len(ss_points)
                    ss_points += [xj]
                    ss_points_pre += [pos]
                    ss_points_dic[xj] = j
                else:
                    j = ss_points_dic[xj]
                T2_matrix[pos, j] += ej
            # end for
            if pos != 0:
                # also record the root from j_prev
                T2_matrix[pos, j_prev] += 1
            pos += int(1)

        self.__hecke_matrices[2] = T2_matrix
        return (ss_points, ss_points_dic)


    def upper_bound_on_elliptic_factors(self, p=None, ellmax=2):
        r"""
        Return an upper bound (provably correct) on the number of
        elliptic curves of conductor equal to the level of this
        supersingular module.

        INPUT:

        - ``p`` - (default: 997) prime to work modulo

        ALGORITHM: Currently we only use `T_2`.  Function will be
        extended to use more Hecke operators later.

        The prime p is replaced by the smallest prime that doesn't
        divide the level.

        EXAMPLE::

            sage: SupersingularModule(37).upper_bound_on_elliptic_factors()
            2

        (There are 4 elliptic curves of conductor 37, but only 2 isogeny
        classes.)
        """
        # NOTE: The heuristic runtime is *very* roughly `p^2/(2\cdot 10^6)`.
        #ellmax -- (default: 2) use Hecke operators T_ell with ell <= ellmax
        if p is None:
            p = 997

        while self.level() % p == 0:
             p = rings.next_prime(p)

        ell = 2
        t = self.hecke_matrix(ell).change_ring(rings.GF(p))

        # TODO: temporarily try using sparse=False
        # turn this off when sparse rank is optimized.
        t = t.dense_matrix()

        B = 2*math.sqrt(ell)
        bnd = 0
        lower = -int(math.floor(B))
        upper = int(math.floor(B))+1
        for a in range(lower, upper):
            tm = verbose("computing T_%s - %s"%(ell, a))
            if a == lower:
                c = a
            else:
                c = 1
            for i in range(t.nrows()):
                t[i,i] += c
            tm = verbose("computing kernel",tm)
            #dim = t.kernel().dimension()
            dim = t.nrows() - t.rank()
            bnd += dim
            verbose('got dimension = %s; new bound = %s'%(dim, bnd), tm)
        return bnd

    def hecke_matrix(self,L):
        r"""
        This function returns the `L^{\text{th}}` Hecke matrix.

        INPUT:

        - ``self`` -- SupersingularModule object

        - ``L`` -- integer, positive

        OUTPUT:
            matrix -- sparse integer matrix

        EXAMPLES:
        This example computes the action of the Hecke operator `T_2`
        on the module of supersingular points on `X_0(1)/F_{37}`::

            sage: S = SupersingularModule(37)
            sage: M = S.hecke_matrix(2)
            sage: M
            [1 1 1]
            [1 0 2]
            [1 2 0]

        This example computes the action of the Hecke operator `T_3`
        on the module of supersingular points on `X_0(1)/F_{67}`::

            sage: S = SupersingularModule(67)
            sage: M = S.hecke_matrix(3)
            sage: M
            [0 0 0 0 2 2]
            [0 0 1 1 1 1]
            [0 1 0 2 0 1]
            [0 1 2 0 1 0]
            [1 1 0 1 0 1]
            [1 1 1 0 1 0]

        .. note::

            The first list --- list_j --- returned by the supersingular_points
            function are the rows *and* column indexes of the above hecke
            matrices and its ordering should be kept in mind when interpreting
            these matrices.

        AUTHORS:

        - David Kohel -- kohel@maths.usyd.edu.au

        - Iftikhar Burhanuddin -- burhanud@usc.edu
        """
        if L in self.__hecke_matrices:
            return self.__hecke_matrices[L]
        SS, II = self.supersingular_points()
        if L == 2:
            # since T_2 gets computed as a side effect of computing the supersingular points
            return self.__hecke_matrices[2]
        Fp2 = self.__finite_field
        h = len(SS)
        R = self.base_ring()
        T_L = MatrixSpace(R,h)(0)
        S, X = Fp2['x'].objgen()

        if L in [3,5,7,11]:
            for i in range(len(SS)):
                ss_i = SS[i]
                phi_L_in_x = Phi_polys(L, X, ss_i)
                rts = phi_L_in_x.roots()
                for r in rts:
                    T_L[i,int(II[r[0]])] = r[1]
        else:
            DBMP = ClassicalModularPolynomialDatabase()
            phi_L = DBMP[L]
            M, (x,y) =Fp2['x,y'].objgens()
            phi_L = phi_L(x,y)

            # As an optimization, we compute the coefficients of y and evaluate
            # them, since univariate polynomial evaluation is much faster than
            # multivariate evaluation (in Sage :-( ).
            uni_coeff = [phi_L(x,0).univariate_polynomial()] + \
                              [phi_L.coefficient(y**i).univariate_polynomial() for
                                          i in range(1,phi_L.degree(y)+1)]
            for i in range(len(SS)):
                ss_i = SS[i]
                ## We would do the eval below, but it is too slow (right now).
                #phi_L_in_x = phi_L(X, ss_i)

                phi_L_in_x = S([f(ss_i) for f in uni_coeff])
                rts = phi_L_in_x.roots()
                for r in rts:
                    T_L[i,int(II[r[0]])] = r[1]

        self.__hecke_matrices[L] = T_L
        return T_L
