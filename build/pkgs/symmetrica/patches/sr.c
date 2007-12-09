/* SYMMETRICA sr.c */

#include "def.h"
#include "macro.h"


static OP mahh_h = NULL;
static INT add_schur_schur_co();
static INT scan_schur_co();
INT splitpart();


INT schur_ende() /* AK 100692 */
{
    INT erg = OK;
    INT thp_ende(),thm_ende(),tmh_ende(),tep_ende(),teh_ende(),tem_ende();
    INT mps_ende(),mem_ende(),mpp_ende(),mmm_ende(),mss_ende();
    INT mes_ende(),tme_ende();

    erg += thp_ende();
    erg += thm_ende();
    erg += tmh_ende();
    erg += tme_ende();
    erg += tep_ende();
    erg += teh_ende();
    erg += tem_ende();
    erg += mem_ende();
    erg += mpp_ende();
    erg += mps_ende();
    erg += mmm_ende();
    erg += mss_ende();
    erg += mes_ende();
    if (mahh_h != NULL) {
        FREEALL(mahh_h);
        mahh_h=NULL;
        }
    ENDR("schur_ende");
}

INT schnitt_schur(a,b,c) OP a,b,c;
/* gemeinsame bestandteile der schurfunktionen a und b nach c*/
/* AK 310789 V1.0 */ /* AK 181289 V1.1 */ /* AK 200891 V1.3 */
{
    OP zeigera, zeigerb;
    INT erg=OK,e;
    CTO(SCHUR,"schnitt_schur(1)",a);
    CTO(SCHUR,"schnitt_schur(2)",b);

    CE3(a,b,c,schnitt_schur);

    zeigera=a; zeigerb=b;
    while(zeigera != NULL && zeigerb !=NULL)
    {
        e =  comp(S_S_S(zeigera),S_S_S(zeigerb));
        if (e == (INT)0)
        {
            OP neu = callocobject();
            erg += m_pa_s(S_S_S(zeigera),neu); /* AK 201192 */
            if (ge(S_S_K(zeigerb),S_S_K(zeigera)))
                erg +=    copy(S_S_K(zeigera),S_S_K(neu));
            else
                erg += copy(S_S_K(zeigerb),S_S_K(neu));

            insert(neu,c,add_koeff,comp);
            zeigera=S_S_N(zeigera);
            zeigerb=S_S_N(zeigerb);
        }
        else if (e < (INT)0) zeigera=S_S_N(zeigera);
        else if (e > (INT)0) zeigerb=S_S_N(zeigerb);
    };
    ENDR("schnitt_schur");
}


INT einsp_symfunc(p) OP p;
/*
return TRUE if constant and coeff is eins
*/
/* AK 181103 */
{
OP z;
FORALL(z,p,{
    if (S_PA_LI(S_MO_S(z)) != 0)
        {
        if (not NULLP(S_MO_K(z))) return FALSE;
        }
    else
        {
        if (not EINSP(S_MO_K(z))) return FALSE;
        }
    });
return TRUE;
}

INT tex_schur(poly) OP poly;
/* AK 101187 */ /* zur ausgabe eines Schurpolynoms */
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */
/* AK 070291 V1.2 prints to texout */ /* AK 200891 V1.3 */
/* AK 021199      works for SCHUR,MONOMIAL,HOMSYM,ELMSYM,POWSYM */
{
    OP zeiger = poly;

    fprintf(texout,"\\ ");
    if (EMPTYP(poly)) return(OK);
    while (zeiger != NULL)
    {
        if (not einsp (S_S_K(zeiger)))
            /* der koeffizient wird nur geschrieben wenn er
            ungleich 1 ist */
            if (listp(S_S_K(zeiger))) { /* AK 130397 */
                fprintf(texout,"(");
                tex(S_S_K(zeiger));
                fprintf(texout,")");
                        }
            else    {
                tex(S_S_K(zeiger));
                }

        if (S_O_K(zeiger) == SCHUR)
            fprintf(texout,"\\ $S_{ ");
        else if (S_O_K(zeiger) == MONOMIAL)
                    fprintf(texout,"\\ $m_{ ");
        else if (S_O_K(zeiger) == HOM_SYM)
                    fprintf(texout,"\\ $h_{ ");
        else if (S_O_K(zeiger) == POW_SYM)
                    fprintf(texout,"\\ $p_{ ");
        else if (S_O_K(zeiger) == ELM_SYM)
                    fprintf(texout,"\\ $e_{ ");

        fprint(texout,S_S_S(zeiger));
        fprintf(texout," } $\\ ");
        zeiger = S_S_N(zeiger);
        if (zeiger != NULL) fprintf(texout," $+$ ");
        texposition += 15;
        if (texposition >tex_row_length) {
            fprintf(texout,"\n");
            texposition = 0;
        }
    };
    fprintf(texout,"\\ ");
    texposition += 3;
    return(OK);
}

INT compute_skewschur_with_alphabet_det(a,b,c) OP a,b,c;
/* skewschurpolyomial with det */
/* AK 090790 V1.1 */ /* AK 250291 V1.2 */ /* AK 200891 V1.3 */
{
    INT erg = OK,i,j,gli,kli;
    OP d,h;
    CTO(SKEWPARTITION,"compute_skewschur_with_alphabet_det",a);
    CTO(INTEGER,"compute_skewschur_with_alphabet_det",b);
    d = callocobject();
    h = callocobject();
    gli = S_SPA_GLI(a);
    kli = S_SPA_KLI(a); /* alt gli */
    erg += m_ilih_m(gli,gli,d);
    for (i=(INT)0; i<gli; i++)
        for (j=(INT)0; j<gli; j++)
            {
            if (i >= (gli - kli) )
                m_i_i(S_SPA_GII(a,j)+j-i-
                      S_SPA_KII(a,i-gli+kli)
                      ,h);
            else
                m_i_i(S_SPA_GII(a,j)+j-i,h);
            erg += compute_complete_with_alphabet(h,b,S_M_IJ(d,i,j));
            }
    erg += det_mat_imm(d,c);
    erg += freeall(d);
    erg += freeall(h);  /* AK 160893 */
    ENDR("compute_skewschur_with_alphabet_det");
}





INT compute_schur_with_alphabet_det(a,b,c) OP a,b,c;
/* schurpolyomial with det */
/* AK 090790 V1.1 */ /* AK 200891 V1.3 */
{
    INT erg = OK;
    CTO(PARTITION,"compute_schur_with_alphabet_det(1)",a);
    CTO(INTEGER,"compute_schur_with_alphabet_det(2)",b);
    CE3(a,b,c,compute_schur_with_alphabet_det);
	{
	    OP d,h;
	    INT i,j;
	    d = callocobject();
	    h = callocobject();
	    erg += m_ilih_m(S_PA_LI(a),S_PA_LI(a),d);
	    for (i=(INT)0; i<S_PA_LI(a); i++)
		for (j=(INT)0; j<S_PA_LI(a); j++)
		    {
		M_I_I(S_PA_II(a,j)+j-i,h);
		erg += compute_complete_with_alphabet(h,b,S_M_IJ(d,i,j));
		    }
	    erg += det_mat_imm(d,c);
	    erg += freeall(d);
	    erg += freeall(h);   /* fehlte AK 301092 */
	  }
    ENDR("compute_schur_with_alphabet_det");
}



INT compute_schur(part,res) OP part,res;
/* AK 161187 */ /* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 260291 V1.2 */
/* AK 200891 V1.3 */
{
    OP l;
    INT erg = OK;
    CTO(PARTITION,"compute_schur(1)",part);
    l=callocobject();
    erg += weight( part,l);
    erg += compute_schur_with_alphabet(part,l,res);
    erg += freeall(l);
    ENDR("compute_schur");
}

INT t_POLYNOM_MONOMIAL(a,b) OP a,b;
/* assumes a is symmetric */
/* AK 080591 V1.2 */ /* AK 200891 V1.3 */
    {
    OP c,d,e;
    INT erg = OK;
    CTO(POLYNOM,"t_POLYNOM_MONOMIAL(1)",a);
    CE2(a,b,t_POLYNOM_MONOMIAL);
    erg += init(MONOMIAL,b);
    c = callocobject();
    erg += copy(a,c);
    e = callocobject();
    while (not NULLP(c))
        {
        d = callocobject();
        erg += m_v_mon(S_PO_S(c),d);
        erg += copy(S_PO_K(c),S_S_K(d));
        erg += compute_monomial_with_alphabet(
                S_S_S(d),S_V_L(S_PO_S(c)),e);
        MULT_APPLY(S_PO_K(c),e);
        erg += sub(c,e,c);
        INSERT_LIST(d,b,add_koeff,comp_monommonomial);
        }
    erg += freeall(e);
    erg += freeall(c);
    ENDR("t_POLYNOM_MONOMIAL");
    }

INT t_POLYNOM_POWER(a,b) OP a,b;
/* assumes a is symmetric */
/* AK 080591 V1.2 */        /* AK 200891 V1.3 */
    {
    OP c,d,e,f=NULL;
    OP pa,z;
    INT pain,i;
    INT erg = OK;
    CTO(POLYNOM,"t_POLYNOM_POWER",a);
    CE2(a,b,t_POLYNOM_POWER);
    if (consp_polynom(a))
        {
        erg += m_scalar_powsym(S_PO_K(a),b);
        goto endr_ende;
        }
    FREESELF(b);
    c = callocobject();
    erg += copy(a,c);
    e = callocobject();
    while (not NULLP(c))
        {
        d = callocobject();
        z = c; pa = c;
        pain = (INT)0;
        while (z != NULL)
            {
            erg += m_v_pa(S_PO_S(z),d);
            if ((i=indexofpart(d)) > pain)
                {
                pain = i;
                pa = z;
                }
            z = S_PO_N(z);
            }
        /* pain ist index der lex kleinsten partition
                   pa der zugehoerige POLYNOM zeiger */
        erg += m_v_ps(S_PO_S(pa),d);
        erg += copy(S_PO_K(pa),S_S_K(d));
        erg += compute_power_with_alphabet(
                S_S_S(d),S_V_L(S_PO_S(pa)),e);
        z = e;
        while (z != NULL)
            {
            if (EQ(S_PO_S(z),S_PO_S(pa)))
            /* find coeff of the leading monom */
                {
                f = callocobject();
                erg += copy(S_PO_K(z),f);
                erg += invers_apply(f);
                MULT_APPLY(f,e);
                break;
                }
            z = S_PO_N(z);
            }
        erg += mult_apply(S_PO_K(pa),e);
        erg += sub(c,e,c);
        erg += mult_apply(f,d); /* AK 020394 */
        erg += freeall(f);
        insert(d,b,NULL,NULL);
        }
    FREEALL2(e,c);
    ENDR("t_POLYNOM_POWER");
    }

INT t_POLYNOM_SCHUR(a,b) OP a,b;
/* assumes a is symmetric */
/* AK 020394 */
{
    OP c;
    INT erg = OK;
    CTO(POLYNOM,"t_POLYNOM_SCHUR(1)",a);
    CE2(a,b,t_POLYNOM_SCHUR);
    if (consp_polynom(a))
        {
        erg += m_scalar_schur(S_PO_K(a),b);
        goto endr_ende;
        }
    init(SCHUR,b); /* AK 080502 */
    c = callocobject();
    erg += t_POLYNOM_POWER(a,c);
    erg += t_POWSYM_SCHUR(c,b);
    erg += freeall(c);
    ENDR("t_POLYNOM_SCHUR");
}

INT t_POLYNOM_ELMSYM(a,b) OP a,b;
/* assumes a is symmetric */
/* AK 120995 faster version */
/* AK 240603 a,b may be equal */
{
    OP c,d,e,f,g;
    INT erg = OK;
    CTO(POLYNOM,"t_POLYNOM_ELMSYM(1)",a);
    CE2(a,b,t_POLYNOM_ELMSYM);

    erg += init(ELMSYM,b);
    if (NULLP(a)) goto ee;

    d = callocobject();
    e = callocobject();
    f = callocobject();
    erg += numberofvariables(a,d);
    erg += copy(a,f);
    while (not NULLP(f))
        {
        c = callocobject();
        g = callocobject();
        erg += m_v_pa(S_PO_S(f),c);
        erg += conjugate(c,c);
        erg += compute_elmsym_with_alphabet(c,d,e);
        erg += b_skn_e(c,callocobject(),NULL,g);
        erg += copy(S_PO_K(f),S_S_K(g));
        insert(g,b,NULL,NULL);
        erg += mult_apply(S_PO_K(f),e);
        erg += sub(f,e,f);
        }

    FREEALL(d);
    FREEALL(e);
    FREEALL(f);
ee:
    ENDR("t_POLYNOM_ELMSYM");
}

static INT c_m_w_a_vp(a,b,c) OP a,b,c;
/* AK 200891 V1.3 */
    {
    OP e,f,g;
    INT erg = OK,i;
    e = CALLOCOBJECT();
    erg += first_permutation(b,e);
    f = CALLOCOBJECT();
    m_l_v(b,f);
    for (i=(INT)0;i<S_I_I(b);i++)
        if (i < S_PA_LI(a)) M_I_I(S_PA_II(a,i),S_V_I(f,i));
        else M_I_I((INT)0,S_V_I(f,i));
    /* f is vector */
    do
        {
        g = CALLOCOBJECT();
        b_skn_po(CALLOCOBJECT(),CALLOCOBJECT(),NULL,g);
        M_I_I((INT)1,S_PO_K(g));
        operate_perm_vector(e,f,S_PO_S(g));
        insert(g,c,NULL,NULL);
        } while(next_apply(e));

    /* nur koeff mit 1 */
    g = c;
    while (g != NULL)
        {
        if (not einsp(S_PO_K(g)))
            m_i_i((INT)1,S_PO_K(g));
        g = S_PO_N(g);
        }

    FREEALL(f);
    FREEALL(e);
    ENDR("c_m_w_a_vp");
    }


INT compute_monomial_with_alphabet(number,l,res) OP number,res,l;
/* l = length of alphabet */
/* AK 090790 V1.1 */ /* AK 090591 V1.2 */ /* AK 200891 V1.3 */
/* AK 060498 V2.0 */
    {
    INT erg = OK;
    OP c,z;


    CTO(INTEGER,"compute_monomial_with_alphabet(2)",l);
    CE3(number,l,res,compute_monomial_with_alphabet);

    erg += init(POLYNOM,res);
    if (S_O_K(number) == PARTITION)
        {
        if (S_PA_K(number) != VECTOR)
            {
            OP c = callocobject();
            erg += t_VECTOR_EXPONENT(number,c);
            erg += compute_monomial_with_alphabet(c,l,res);
            erg += freeall(c);
            goto endr_ende;
            }
        /* number is VECTOR partition */
        if (GR(S_PA_L(number),l))
            goto endr_ende;
        erg += c_m_w_a_vp(number,l,res);
        goto endr_ende;
        }
    else if (S_O_K(number) == INTEGER)
        {
        c = callocobject();
        erg += m_i_pa(number,c);
        erg += compute_monomial_with_alphabet(c,l,res);
        erg += freeall(c);
        goto endr_ende;
        }
    else if (S_O_K(number) == MONOMIAL)
    {
    erg += init(POLYNOM,res);
    if (S_L_S(number) == NULL)
        goto endr_ende;
    c = callocobject();
    z = number;
    while (z != NULL)
        {
        erg += compute_monomial_with_alphabet(S_S_S(z),l,c);
        erg += mult_apply(S_S_K(z),c);
        erg += add_apply(c,res);
        z = S_S_N(z);
        }
    erg += freeall(c);
    goto endr_ende;
    }
    else
      return WTT("compute_monomial_with_alphabet",number,l);

    ENDR("compute_monomial_with_alphabet");
    }


INT compute_complete_with_alphabet(number,l,res) OP number,res,l;
/* AK 090790 V1.1 */ /* AK 260291 V1.2 */ /* AK 200891 V1.3 */
/* number may be INTEGER,PARTITION or HOM_SYM */
    {
    OP b,z;
    INT erg=OK,i;
    if (not EMPTYP(res))
        erg += freeself(res);

    if (S_O_K(number) == INTEGER) /* AK 080992 */
    {
    if (S_I_I(number) == (INT)0)
        {
        M_I_I((INT)1,res);
        goto endr_ende;
        }
    else if (S_I_I(number) < (INT)0 )
        {
        M_I_I((INT)0,res);
        goto endr_ende;
        }

    b = callocobject();
    erg += m_i_pa(number,b);
    erg += compute_schur_with_alphabet(b,l,res);
    erg += freeall(b);
    }
    else if (S_O_K(number) == PARTITION) /* AK 080992 */
    {
    if (S_PA_K(number) != VECTOR)
        return ERROR;
    m_i_i((INT)1,res);
    b = callocobject();
    for (i=(INT)0;i<S_PA_LI(number);i++)
        {
        erg += compute_complete_with_alphabet(S_PA_I(number,i),l,b);
        erg += mult_apply(b,res);
        erg += freeself(b);
        }
    erg += freeall(b);
    }
    else if (S_O_K(number) == HOM_SYM)
    {
    m_i_i((INT)0,res);
    b = callocobject();
    z = number;
    while (z != NULL)
        {
        erg += compute_complete_with_alphabet(S_S_S(z),l,b);
        erg += mult_apply(S_S_K(z),b);
        erg += add_apply(b,res);
        z = S_S_N(z);
        erg += freeself(b);
        }
    erg += freeall(b);
    }
    else
        erg += ERROR;
    ENDR("compute_complete_with_alphabet");
    }


INT compute_elmsym_with_alphabet(label,l,result) OP l,label,result;
/* AK 120391 V1.2 */ /* AK 200891 V1.3 */
{
    INT erg = OK;
    INT i;
    OP zw,z;
    CTTTO(INTEGER,ELMSYM,PARTITION,"compute_elmsym_with_alphabet(1)",label);
    CTO(INTEGER,"compute_elmsym_with_alphabet(2)",l);
    CE3(label,l,result,compute_elmsym_with_alphabet);

    if (S_O_K(label) == INTEGER)
        {
/*
        erg += init(POLYNOM,result);
        zw = callocobject();
        erg += last_partition(label,zw);
        erg += compute_monomial_with_alphabet(zw,l,result);
        erg += freeall(zw);
*/
        if (S_I_I(label) > S_I_I(l)) {
            erg += init(POLYNOM,result);
            goto ende;
            }
        zw = CALLOCOBJECT();
        erg += m_il_nv(S_I_I(l),zw);
        for (i=0;i<S_I_I(label);i++) M_I_I(1,S_V_I(zw,S_V_LI(zw)-1-i));
        erg += lehmercode(zw,zw);
        erg += m_perm_schubert_monom_summe(zw,result);
        FREEALL(zw);
        }
    else if (S_O_K(label) == PARTITION)
        {
        zw = callocobject();
        erg += m_scalar_polynom(cons_eins,result);
        for (i=(INT)0; i<S_PA_LI(label); i++)
            {
            erg += compute_elmsym_with_alphabet(S_PA_I(label,i),
                                      l,zw);
            erg += mult_apply(zw,result);
            }
        erg += freeall(zw);
        }
    else if (S_O_K(label) == ELM_SYM)
        {
        zw = callocobject();
        m_i_i((INT)0,result);
        z = label;
        while (z != NULL)
            {
            erg += compute_elmsym_with_alphabet(S_S_S(z),l,zw);
            erg += mult_apply(S_S_K(z),zw);
            erg += add_apply(zw,result);
            z = S_S_N(z);
            erg += freeself(zw);
            }
        erg += freeall(zw);
        }


ende:
    ENDR("compute_elmsym_with_alphabet");
}


INT compute_power_with_alphabet(label,l,result) OP l,label,result;
/* AK 120391 V1.2 */ /* AK 200891 V1.3 */
{
    INT erg = OK;
    INT i;
    OP zw,z;
    CTTTO(INTEGER,PARTITION,POWSYM,"compute_power_with_alphabet(1)",label);
    CTO(INTEGER,"compute_power_with_alphabet(2)",l);
    FREESELF(result);

    if (S_O_K(label) == INTEGER)
        {
        erg += init(POLYNOM,result);
        for (i=(INT)0;i<S_I_I(l); i++)
            {
            zw = callocobject();
            erg += m_iindex_iexponent_monom(i,S_I_I(label),zw);
            insert(zw,result,NULL,NULL);
            }
        }
    else if (S_O_K(label) == PARTITION)
        {
        zw = callocobject();
        erg += m_scalar_polynom(cons_eins,result);
        for (i=(INT)0; i<S_PA_LI(label); i++)
            {
            erg += compute_power_with_alphabet(S_PA_I(label,i),
                                      l,zw);
            erg += mult_apply(zw,result);
            }
        erg += freeall(zw);
        }
    else if (S_O_K(label) == POW_SYM)
        {
        zw = callocobject();
        erg += init(POLYNOM,result);
        z = label;
        while (z != NULL)
            {
            erg += compute_power_with_alphabet(S_S_S(z),l,zw);
            erg += mult_apply(S_S_K(z),zw);
            erg += add_apply(zw,result);
            z = S_S_N(z);
            erg += freeself(zw);
            }
        erg += freeall(zw);
        }
    else    {
        printobjectkind(label);
        erg = error("compute_power_with_alphabet:wrong kind of label");
        }
    ENDR("compute_power_with_alphabet");
}


INT compute_schur_with_alphabet(part,l,res) OP part,res,l;
/* AK 101187 */
/* AK 161187 l ist die laenge des alphabets */
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 150591 V1.2 */
/* AK 200891 V1.3 */
/* AK 161195 much faster with m_umriss_tableaux etc. */
{
    OP e;
    INT erg = OK;
    CE3(part,l,res,compute_schur_with_alphabet);
    CTTTO(PARTITION,SCHUR,HASHTABLE,"compute_schur_with_alphabet(1)",part);

    FREESELF(res);

    if (S_O_K(part) == PARTITION)
        {
        if (GR(S_PA_L(part),l))
            { m_scalar_polynom(cons_null,res); goto ende; }
        if (S_PA_LI(part) == 0)
            { m_scalar_polynom(cons_eins,res); goto ende; }
        e=callocobject();
        erg += m_umriss_tableaux(part,l,e);
        erg += m_tableaux_polynom(e,res);
        erg += freeall(e);
        goto ende;
        }
    else /* SCHUR, HASHTABLE */
        {
        OP z,e;
        init(POLYNOM,res);
        FORALL(z,part, {
            e = CALLOCOBJECT();
            erg += compute_schur_with_alphabet(S_MO_S(z),l,e);
            MULT_APPLY(S_MO_K(z),e);
            erg += insert(e,res,add_koeff,comp_monomvector_monomvector);
            });
        goto ende;
        }
ende:
    CTO(POLYNOM,"compute_schur_with_alphabet(e3)",res);
    ENDR("compute_schur_with_alphabet");
}


INT m_pa_s(part,schur) OP part, schur;
/* AK 090790 V1.1 */ /* AK 260291 V1.2 */ /* AK 050891 V1.3 */
/* AK 151298 V2.0 */
/* part  and schur may be equal */
{
    INT erg = OK;
    OP d;
    COP("m_pa_s",schur);
    CTO(PARTITION,"m_pa_s(1)",part);

    d = CALLOCOBJECT();
    erg += copy_partition(part,d);
    erg += b_pa_s(d,schur);
    ENDR("m_pa_s");
}

INT m_pa_mon(part,schur) OP part, schur;
/* AK 090790 V1.1 */ /* AK 260291 V1.2 */ /* AK 050891 V1.3 */
/* AK 151298 V2.0 */
/* part  and schur may be equal */
{
    INT erg = OK;
    OP d;
    COP("m_pa_mon",schur);
    CTO(PARTITION,"m_pa_mon(1)",part);

    d = CALLOCOBJECT();
    erg += copy_partition(part,d);
    erg += b_pa_s(d,schur);
    C_O_K(schur,MONOMIAL);
    ENDR("m_pa_mon");
}

INT m_pa_e(part,schur) OP part, schur;
/* AK 090790 V1.1 */ /* AK 260291 V1.2 */ /* AK 050891 V1.3 */
/* AK 151298 V2.0 */
/* part  and schur may be equal */
{
    INT erg = OK;
    OP d;
    COP("m_pa_e",schur);
    CTO(PARTITION,"m_pa_e(1)",part);

    d = CALLOCOBJECT();
    erg += copy_partition(part,d);
    erg += b_pa_s(d,schur);
    C_O_K(schur,ELMSYM);
    ENDR("m_pa_e");
}

INT m_pa_h(part,schur) OP part, schur;
/* AK 090790 V1.1 */ /* AK 260291 V1.2 */ /* AK 050891 V1.3 */
/* AK 151298 V2.0 */
/* part  and schur may be equal */
{
    INT erg = OK;
    OP d;
    COP("m_pa_h",schur);
    CTO(PARTITION,"m_pa_h(1)",part);

    d = CALLOCOBJECT();
    erg += copy_partition(part,d);
    erg += b_pa_s(d,schur);
    C_O_K(schur,HOMSYM);
    ENDR("m_pa_h");
}

INT m_pa_ps(part,schur) OP part, schur;
/* AK 090790 V1.1 */ /* AK 260291 V1.2 */ /* AK 050891 V1.3 */
/* AK 151298 V2.0 */
/* part  and schur may be equal */
{
    INT erg = OK;
    OP d;
    COP("m_pa_ps",schur);
    CTO(PARTITION,"m_pa_ps(1)",part);

    d = CALLOCOBJECT();
    erg += copy_partition(part,d);
    erg += b_pa_s(d,schur);
    C_O_K(schur,POWSYM);
    ENDR("m_pa_ps");
}

INT b_pa_mon(part,schur) OP part, schur;
/* AK 140687 */
/* AK 110789 V1.0 */ /* AK 241189 V1.1 */ /* AK 200891 V1.3 */
/* AK 151298 V2.0 */
/* input: partition p
   output: monomial symmetric  function m_p with coefficient 1 */
/* input must be a partition object of vector type */
{
    INT erg=OK;
    OP d;
    COP("m_pa_mon",schur);
    CTO(PARTITION,"b_pa_mon(1)",part);
    d = CALLOCOBJECT();
    erg += copy_partition(part,d);
    erg += b_pa_s(d,schur);
    C_O_K(schur,MONOMIAL);
    ENDR("b_pa_mon");
}

INT b_pa_s(part,schur) OP part, schur;
/* AK 140687 */
/* AK 110789 V1.0 */ /* AK 241189 V1.1 */ /* AK 200891 V1.3 */
/* AK 151298 V2.0 */
/* input: partition p
   output: schur function s_p with coefficient 1 */
/* input must be a partition object of vector type */
{
    INT erg=OK;
    CTO(PARTITION,"b_pa_s(1)",part);
    SYMCHECK(part == schur,"b_pa_s:identic objects");
    SYMCHECK(S_PA_K(part) != VECTOR,"b_pa_s:partition must be of kind vector");
    erg += b_skn_s(part,CALLOCOBJECT(),NULL,schur);
    M_I_I(1,S_S_K(schur));
    ENDR("b_pa_s");
}



INT m_v_s(vec,schur) OP vec, schur;
/* AK 110187 */
/* AK 110789 V1.0 */ /* AK 170590 V1.1 */ /* AK 200891 V1.3 */
{
    INT erg=OK;
    CE2(vec,schur,m_v_s);
    erg += b_skn_s(    CALLOCOBJECT(), CALLOCOBJECT(),
        NULL, schur);
    erg += m_v_pa(vec,S_S_S(schur));
    M_I_I((INT)1,S_S_K(schur));
    ENDR("m_v_s");
}


INT add_monomial_monomial(a,b,c) OP a, b, c;
{
    return add_schur_schur_co(a,b,c,MONOMIAL);
}
INT add_elmsym_elmsym(a,b,c) OP a, b, c;
{
    return add_schur_schur_co(a,b,c,ELM_SYM);
}
INT add_powsym_powsym(a,b,c) OP a, b, c;
{
    return add_schur_schur_co(a,b,c,POW_SYM);
}



INT add_homsym_homsym(a,b,c) OP a, b, c;
/* AK 200891 V1.3 */
{
    return add_schur_schur_co(a,b,c,HOM_SYM);
}



INT add_schur_schur(a,b,c) OP a, b, c;
/* AK 200891 V1.3 */
{
    INT erg = OK;
    CTTO(HASHTABLE,SCHUR,"add_schur_schur(1)",a);
    CTTO(HASHTABLE,SCHUR,"add_schur_schur(2)",b);
    CTTTO(EMPTY,HASHTABLE,SCHUR,"add_schur_schur(3)",c);
    erg += add_schur_schur_co(a,b,c,SCHUR);
    ENDR("add_schur_schur");
}



static INT add_schur_schur_co(a,b,c,typ) OP a, b, c;OBJECTKIND typ;
/* AK 110789 V1.0 */ /* AK 201289 V1.1 */ /* AK 050891 V1.3 */
{
    INT erg = OK;
    CTTTTTO(SCHUR,POWSYM,ELMSYM,HOMSYM,MONOMIAL,"add_schur_schur_co(1)",a);
    erg += add_polynom_polynom(a,b,c);
    if (listp(c))
        C_O_K(c,typ);
    ENDR("internal routine:add_schur_schur_co");
}



INT m_skn_s(self,koeff,n,ergebnis) OP self,koeff,n,ergebnis;
/* AK 200891 V1.3 */
    {
    INT erg = OK;
    CTO(PARTITION,"m_skn_s(1)",self);
    COP("m_skn_s(4)",ergebnis);
    erg += m_skn_po(self,koeff,n,ergebnis);
    C_O_K(ergebnis,SCHUR);
    ENDR("m_skn_s");
    }


INT b_skn_s(a,b,c,d) OP a,b,c,d;
/* AK 110789 V1.0 */ /* AK 130391 V1.2 */
/* AK 130891 V1.3 */
{
    INT erg = OK;
    CTTO(EMPTY,PARTITION,"b_skn_s(1)",a);
    COP("b_skn_s(4)",d);

    erg += b_sn_l(CALLOCOBJECT(),c,d);
    C_O_K(d,SCHUR);
    erg += b_sk_mo(a,b,S_L_S(d));
    ENDR("b_skn_s");
}


INT objectread_schur(filename,poly) OP poly; FILE *filename;
/* AK 291086 */ /* AK 110789 V1.0 */ /* AK 221289 V1.1 */
/* AK 200891 V1.3 */
{
    char antwort[2];
    INT erg = OK; /* AK 040893 */
    COP("objectread_schur(2)",poly);

    erg += b_skn_s(callocobject(), callocobject(), callocobject(), poly);
    erg += objectread(filename,PARTITION,S_S_S(poly));
    erg += objectread(filename,INTEGER,S_S_K(poly));
    fscanf(filename,"%s",antwort);
    if (antwort[0] == 'j')
        erg += objectread(filename,SCHUR,S_S_N(poly));
    else if (antwort[0] == 'n')     {
        SYM_free(S_S_N(poly));
        C_S_N(poly,NULL);
    }
    else
        error("objectread_schur:wrong data");
    ENDR("objectread_schur");
}

INT objectwrite_schur(filename,poly) FILE *filename; OP poly;
/* AK 291086 */ /* AK 110789 V1.0 */ /* AK 221289 V1.1 */
/* AK 200891 V1.3 */
{
    INT erg = OK;
    COP("objectwrite_schur(2)",poly);
    erg += objectwrite(filename,S_S_S(poly));
    erg += objectwrite(filename,S_S_K(poly));
    if (not lastp(poly))
    {
        fprintf(filename,"j\n");
        erg += objectwrite(filename,S_S_N(poly));
    }
    else fprintf(filename,"n\n");
    ENDR("objectwrite_schur");
}

INT scan_monomial(a) OP a;
{
    printeingabe("Input of a monomial symmetric function");
    return scan_schur_co(a,MONOMIAL);
}
INT scan_elmsym(a) OP a;
{
    printeingabe("Input of a elementary symmetric function");
    return scan_schur_co(a,ELM_SYM);
}
INT scan_powsym(a) OP a;
{
    printeingabe("Input of a powersum symmetric function");
    return scan_schur_co(a,POW_SYM);
}
INT scan_homsym(a) OP a;
{
    printeingabe("Input of a complete symmetric function");
    return scan_schur_co(a,HOM_SYM);
}
INT scan_schur(a) OP a;
{
    printeingabe("Input of a Schur function");
    return scan_schur_co(a,SCHUR);
}
INT sscan_schur(t,a) char *t; OP a;
/* AK 171296 */
{
    INT i,n=1,erg = OK;
    OP c,d,e;
    char *v;
    int SYM_isdigit();

    COP("sscan_schur(2)",a);

    c = callocobject();
    d = callocobject();
    e = callocobject();m_i_i(1L,e);
    erg += init(SCHUR,a);
    v = t;
again:
    if (*v == '\0') goto sse;
        while (*v == ' ') v++;
    if (*v == '[') {
        i = sscan(v,PARTITION,c);
        if (i != OK) goto sse;
        while (*v != ']') v++;
        v++;
        erg += m_skn_s(c,e,NULL,d);
        erg += add_apply(d,a);
        m_i_i(1L,e);
        goto again;
        }
    else if (*v == '+') {
        v++; n=1;
        goto again;
        }
    else if (*v == '-') {
        v++; n= -1;
        goto again;
        }
    else if (SYM_isdigit(*v))
        {
        i = sscan(v,INTEGER,e);
        if (i != OK) goto sse;
        while (SYM_isdigit(*v)) v++;
        v++;
        if (n == -1) addinvers_apply(e);
        n=1;
        goto again;
        }

sse:
    erg += freeall(c);
    erg += freeall(d);
    erg += freeall(e);
    ENDR("sscan_schur");

}

INT sscan_homsym(t,a) char *t; OP a;
/* AK 050901 */
{
    INT i,n=1,erg = OK;
    OP c,d,e;
    char *v;
    int SYM_isdigit();
    COP("sscan_homsym(1)",t);
    COP("sscan_homsym(2)",a);

    c = callocobject();
    d = callocobject();
    e = callocobject();m_i_i(1L,e);
    erg += init(HOMSYM,a);
        v = t;
again:
    if (*v == '\0') goto sse;
        while (*v == ' ') v++;
    if (*v == '[') {
        i = sscan(v,PARTITION,c);
        if (i != OK) goto sse;
        while (*v != ']') v++;
        v++;
        erg += m_skn_h(c,e,NULL,d);
        erg += add_apply(d,a);
        m_i_i(1L,e);
        goto again;
        }
    else if (*v == '+') {
        v++; n=1;
        goto again;
        }
    else if (*v == '-') {
        v++; n= -1;
        goto again;
        }
    else if (SYM_isdigit(*v))
        {
        i = sscan(v,INTEGER,e);
        if (i != OK) goto sse;
        while (SYM_isdigit(*v)) v++;
        v++;
        if (n == -1) addinvers_apply(e);
        n=1;
        goto again;
        }

sse:
    erg += freeall(c);
    erg += freeall(d);
    erg += freeall(e);
    ENDR("sscan_homsym");
}

INT sscan_elmsym(t,a) char *t; OP a;
/* AK 050901 */
{
    INT i,n=1,erg = OK;
    OP c,d,e;
    char *v;
    int SYM_isdigit();
    COP("sscan_elmsym(1)",t);
    CTO(EMPTY,"sscan_elmsym(2)",a);

    c = callocobject();
    d = callocobject();
    e = callocobject();M_I_I(1L,e);
    erg += init_elmsym(a);
    v = t;
again:
    if (*v == '\0') goto sse;
        while (*v == ' ') v++;
    if (*v == '[') {
        i = sscan(v,PARTITION,c);
        if (i != OK) goto sse;
        while (*v != ']') v++;
        v++;
        erg += m_skn_e(c,e,NULL,d);
        erg += add_apply(d,a);
        m_i_i(1L,e);
        goto again;
        }
    else if (*v == '+') {
        v++; n=1;
        goto again;
        }
    else if (*v == '-') {
        v++; n= -1;
        goto again;
        }
    else if (SYM_isdigit(*v))
        {
        i = sscan(v,INTEGER,e);
        if (i != OK) goto sse;
        while (SYM_isdigit(*v)) v++;
        v++;
        if (n == -1) addinvers_apply(e);
        n=1;
        goto again;
        }

sse:
    FREEALL(c);
    FREEALL(d);
    FREEALL(e);
    ENDR("sscan_elmsym");
}


static INT scan_schur_co(a,typ) OP a; OBJECTKIND typ;
/* AK  for the input of a symmetric function (Schur etc.) */
/* AK 110789 V1.0 */ /* AK 221289 V1.1 */ /* AK 050891 V1.3 */

{
    char antwort[2];
    OBJECTKIND kind;
    INT erg=OK;
    OP d;

    COP("scan_schur_co(1)",a);

    erg += b_skn_s( callocobject(), callocobject(), NULL, a);
    C_O_K(a,typ);
    erg += printeingabe("Input of a partition type monom");
    erg += scan(PARTITION,S_S_S(a));
    erg += printeingabe("Input of coefficent");
    kind = scanobjectkind();
    erg += scan(kind,S_S_K(a));
    erg += printeingabe("one more monom y/n");
    scanf("%s",antwort);
    if (antwort[0]  == 'y')
        {
        d = callocobject();
        erg += scan_schur_co(d,typ);
        erg += insert(d,a,NULL,NULL); /* AK 091195 reihenfolge beliebig */
        }
    ENDR("scan_schur internal routine");
}


OP s_s_s(a) OP a;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_mo_s(s_l_s(a))); }

OP s_s_k(a) OP a;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_mo_k(s_l_s(a))); }

OP s_s_n(a) OP a;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_l_n(a)); }

OP s_s_si(a,i) OP a; INT i;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_pa_i(s_mo_s(s_l_s(a)),i)); }

OP s_s_sl(a) OP a;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_pa_l(s_mo_s(s_l_s(a)))); }

INT s_s_ki(a) OP a;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_mo_ki(s_l_s(a))); }

INT s_s_sii(a,i) OP a; INT i;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_pa_ii(s_mo_s(s_l_s(a)),i)); }


INT s_s_sli(a) OP a;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{ return(s_pa_li(s_mo_s(s_l_s(a)))); }

INT c_s_n(a,b) OP a,b;
/* AK 110789 V1.0 */ /* AK 181289 V1.1 */ /* AK 020791 V1.2 */
/* AK 200891 V1.3 */
{
    OBJECTSELF c;
    c = s_o_s(a);
    c.ob_list->l_next = b;
    return(OK);
}


INT test_schur()
/* AK 181289 V1.1 */ /* AK 020791 V1.2 */ /* AK 200891 V1.3 */
{
    OP a = callocobject();
    OP b = callocobject();
    OP c = callocobject();

    printeingabe("test_schur:scan(a)");
    scan(SCHUR,a);
    println(a);
    printeingabe("test_schur:copy(a,b)");
    copy(a,b);
    println(b);
    printeingabe("test_schur:add(a,b,b)");
    add(a,b,b);
    println(b);
    printeingabe("test_schur:mult(a,b,b)");
    mult(a,b,b);
    println(b);
    printeingabe("test_schur:addinvers(b,a)");
    addinvers(b,a);
    println(a);
    printeingabe("test_schur:mult_apply(b,a)");
    mult_apply(b,a);
    println(a);

    freeall(a);
    freeall(b);
    freeall(c);
    return(OK);
}






INT comp_colex_schurmonom(a,b) OP a,b;
/* AK 091189 */ /* AK V1.1 201189 */
/* AK 200891 V1.3 */
{
    INT erg = OK;
    CTO(MONOM,"comp_colex_schurmonom",a);
    CTO(MONOM,"comp_colex_schurmonom",b);
    return(comp_colex_part(S_MO_S(a),S_MO_S(b)));
    ENDR("comp_colex_schurmonom");
}


INT comp_colex_part(a,b) OP a,b;
/* a,b partitions colex order */
/* AK V1.1 151189 */ /* AK 200891 V1.3 */
{
    INT i = S_PA_LI(a)-(INT)1;
    INT j = S_PA_LI(b)-(INT)1;
    INT erg;

    if (S_O_K(a) != PARTITION)
        error("comp_colex_part:kind != PARTITION");
    if (S_O_K(b) != PARTITION)
        error("comp_colex_part:kind != PARTITION");


    for (;(i >= (INT)0) || (j>=(INT)0); i--,j--)
    {
        if (i<(INT)0) return((INT)1);
        if (j<(INT)0) return((INT)-1);
        erg = S_PA_II(a,i) - S_PA_II(b,j);
        if (erg <(INT)0) return((INT)1);
        if (erg >(INT)0) return((INT)-1);
    }
    return((INT)0);
}




INT hall_littlewood_tafel(a,b) OP a,b;
/* AK 191289 a ist grad der sn b wird tafel */ /* AK 201289 V1.1 */
/* AK 200891 V1.3 */
{
    INT i,j;
    OP c = callocobject();
    OP d = callocobject();
    OP z,zz;
    INT erg = OK;
    CTO(INTEGER,"hall_littlewood_tafel",a);
    erg += makevectorofpart(a,c);
    erg += m_ilih_nm(S_V_LI(c),S_V_LI(c),b);



    for (i=(INT)0;i<S_V_LI(c);i++)
    {
        erg += hall_littlewood(S_V_I(c,i),d);
        z = d;
        while (z != NULL) {
            zz = S_MO_S(S_L_S(z)); /* partition */
            for (j=0;j<S_V_LI(c);j++)
                if (EQ(zz,S_V_I(c,j))) break;
            erg += copy(S_MO_K(S_L_S(z)),S_M_IJ(b,i,j));
            /* koef sind polynom */
            z = S_L_N(z);
        }
    }

    erg += freeall(c);
    erg += freeall(d);
    ENDR("hall_littlewood_tafel");
}


INT hall_littlewood_alt(a,b) OP a,b;
/* AK 191289 a ist partition
b wird das zugehoerige hall littlewood polynom */
/* mittels d_ij = langsam */ /* schneller morris mit skew schur */
/* AK 201289 V1.1 */
/* AK 200891 V1.3 */
{
    INT i,j;
    OP c = callocobject();

    if (not EMPTYP(b))
        freeself(b);
    init_hall_littlewood(a,c);

    for (i = 0;i<S_PA_LI(a);i++)
        for (j = i+1;j<S_PA_LI(a);j++)
            hall_littlewood_dij(c,c,i,j);

    reorder_hall_littlewood(c,b);
    return freeall(c);
}


INT init_hall_littlewood(a,b) OP a,b;
/* AK 200891 V1.3 */
{
    b_skn_s(callocobject(),callocobject(),NULL,b);
    copy_partition(a,S_S_S(b));
    m_skn_po(callocobject(),callocobject(),NULL,S_S_K(b));
    m_il_v((INT)1,S_PO_S(S_S_K(b)));
    M_I_I((INT)0,S_PO_SI(S_S_K(b),(INT)0));
    M_I_I((INT)1,S_PO_K(S_S_K(b)));
    return(OK);
}


INT reorder_vector(a,b) OP a,b;
/* AK 280901 */
/* sorts the vector to a partition */
/* return +1 , -1, 0 */
/* return 0 means zero schur function */
/* a is integervector */
{
    INT erg = OK;
    CTO(INTEGERVECTOR,"reorder_vector(1)",a);
    CTO(EMPTY,"reorder_vector(2)",b);
    erg += copy_integervector(a,b);
    return reorder_vector_apply(b);
    ENDR("reorder_vector");
}

INT reorder_vector_apply(b) OP b;
/* AK 041201 */
/* reorders a schur function partition
   return 0,+1,-1 */
{
    INT erg = OK;
    INT i,res=1,t;
    CTO(INTEGERVECTOR,"reorder_vector_apply",b);
    i=1;
    while(i<S_V_LI(b))
        {
        if (i==0) i++;
        if (i==1)
            if (S_V_II(b,0) < 0) { return 0; }
        if (S_V_II(b,i) == S_V_II(b,i-1)-1) {
            return 0; }
        else if (S_V_II(b,i) < S_V_II(b,i-1)) {
            res = res * (-1);
            INC_INTEGER(S_V_I(b,i));
            DEC_INTEGER(S_V_I(b,i-1));
            t=S_V_II(b,i);
            M_I_I(S_V_II(b,i-1),S_V_I(b,i));
            M_I_I(t,S_V_I(b,i-1));
            i--;
            }
        else i++;
        }

    /* nach links packen */
    for (i=0;i<S_V_LI(b);i++)
    if (S_V_II(b,i) != 0) break;
    for (t=0;i<S_V_LI(b);i++,t++)
    M_I_I(S_V_II(b,i),S_V_I(b,t));
    M_I_I(t,S_V_L(b));
    return res;
    ENDR("reorder_vector_apply");
}



INT reorder_schur(a,b) OP a,b;
{
    INT erg = OK;
    INT i;

    OP d,z,m;
    CTTO(SCHUR,HASHTABLE,"reorder_schur(1)",a);
    if ((S_O_K(b) != SCHUR) && (S_O_K(b) != HASHTABLE))
        {
        CE2(a,b,reorder_schur);
        }
    if (S_O_K(b) == EMPTY)
        {
        erg += init(SCHUR,b);
        }
    FORALL(z,a, {
        d = CALLOCOBJECT();
        i = reorder_vector(S_PA_S(S_MO_S(z)),d);
        if (i == 0) { FREEALL(d); }
        else {
            m = CALLOCOBJECT();
            b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),m);
            b_ks_pa(VECTOR,d,S_MO_S(m));
            if (i==1) COPY(S_MO_K(z),S_MO_K(m));
            else ADDINVERS(S_MO_K(z),S_MO_K(m));

            if (S_O_K(b) == SCHUR)
                INSERT_LIST(m,b,add_koeff,comp_monomschur);
            else
                insert_scalar_hashtable(m,b,add_koeff,eq_monomsymfunc,hash_monompartition);
            }
        });


    ENDR("reorder_schur");
}

INT reorder_hall_littlewood(a,b) OP a,b;
/* AK 191289
es werden die partitionen wieder neu sortiert
in die ansteigende form */
/* a is SCHUR */
/* AK 200891 V1.3 */
{
    OP z = a,d,e;
    OP zzz;
    INT i,j;
    INT erg = OK;
    CTO(SCHUR,"reorder_hall_littlewood",a);
    CE2(a,b,reorder_hall_littlewood);

    init(SCHUR,b);
    while (z != NULL)
    {
        d = CALLOCOBJECT();
        erg += copy_monom(S_L_S(z),d);
        zzz = S_MO_S(d); /* zzz ist partition */
re_again:
        for (i=1;i<S_PA_LI(zzz);i++)
        {

            if (S_PA_II(zzz,(INT)0) < (INT)0)
            {
                FREEALL(d);
                goto re_while_ende;
            }
            else if (S_PA_II(zzz,i) == S_PA_II(zzz,i-1) -1)
            {
                FREEALL(d);
                goto re_while_ende;
            }
            else if (S_PA_II(zzz,i) < S_PA_II(zzz,i-1) )
            {
                ADDINVERS_APPLY(S_MO_K(d));
                INC_INTEGER(S_PA_I(zzz,i));
                DEC_INTEGER(S_PA_I(zzz,i-1));
                SWAP(S_PA_I(zzz,i),S_PA_I(zzz,i-1));
                goto re_again;
            }
        }
        for (i=(INT)0;i<S_PA_LI(zzz);i++)
            if (S_PA_II(zzz,i)>(INT)0) break;
        /* noch nach links schieben */

        for (j=i;j<S_PA_LI(zzz);j++)
            M_I_I(S_PA_II(zzz,j),S_PA_I(zzz,j-i));

        if ((S_PA_LI(zzz)-i) == (INT)0) /* AK 300197 */
            {
            j = S_PA_II(zzz,(INT)0);
            erg += m_il_v((INT)1,S_PA_S(zzz));
            C_O_K(S_PA_S(zzz),INTEGERVECTOR);
            M_I_I(j,S_PA_I(zzz,(INT)0));
            }
        else if ((S_PA_LI(zzz)-i) == (INT)1) /* AK 121093 */
            {
            j = S_PA_II(zzz,(INT)0);
            erg += m_il_v((INT)1,S_PA_S(zzz));
            C_O_K(S_PA_S(zzz),INTEGERVECTOR);
            M_I_I(j,S_PA_I(zzz,(INT)0));
            }
        else     {
            M_I_I(S_PA_LI(zzz)-i,S_PA_L(zzz));
            }

        e = CALLOCOBJECT();
        erg += b_sn_s(d,NULL,e);
        INSERT_LIST(e,b,add_koeff,comp_monomvector_monomvector);
re_while_ende:
        z = S_L_N(z);
    }
    if (S_L_S(b) == NULL)
        FREESELF(b);
    ENDR("reorder_hall_littlewood");
}



INT hall_littlewood_dij(a,b,i,j) OP a,b; INT i,j;
/* AK 181289
bei der berechnung von hall littlewood polynomen benoetigt
man die anwendung vomuliplikation mit
(1 + t * d_ij + t^2 * d_ij^2 ... )
eingabe: a = hall_littlewood polynom
     i<j indices
ausgabe: b neues hall_littlewood polynom
*/
/* hall littlewood polynome sind schurpolynome mit polynomen in t als koeff*/
/* AK 201289 V1.1 */
/* AK 200891 V1.3 */
{
    INT k,tt;
    OP sp = callocobject();
    OP z,zz,zzz;

    copy_list(a,sp); /* funktioniert auch beim
                aufruf mit gleichen variablen */
    copy_list(sp,b); /* die multiplikation mit 1 */
    for (k=1;;k++)
    {
        tt = (INT)0;
        z = sp;
        while (z != NULL)
        {
            zz = S_L_S(z);
            zzz = S_MO_S(zz);
            if (j <= S_PA_LI(zzz)) /* index j zulaessig */
            if (S_PA_II(zzz,i) >= (k-i) ) {
                OP d = callocobject();
                OP e = callocobject();
                tt = (INT)1;
                copy(zz,d);
                M_I_I(S_PA_II(zzz,i)-k,S_PA_I(S_MO_S(d),i));
                M_I_I(S_PA_II(zzz,j)+k,S_PA_I(S_MO_S(d),j));
                b_skn_po(callocobject(),callocobject(),NULL,e);
                m_il_v((INT)1,S_PO_S(e));
                M_I_I(k,S_PO_SI(e,(INT)0));
                M_I_I((INT)1,S_PO_K(e)); /* e = t^k */
                mult(e,S_MO_K(d),S_MO_K(d));
            insert(d,b,add_koeff,comp_monomvector_monomvector);
                /* add(d,b,b);*/
                freeall(e);
            }
            z = S_L_N(z);
        }
        if (tt == (INT)0) break; /* ende */
    }
    freeall(sp);
    return(OK);
}


INT tex_hall_littlewood(a) OP a;
/* AK 191289 tex ausgabe */
/* AK 200891 V1.3 */
{
    return tex(a);
}





INT hall_littlewood(a,b) OP a,b;
/* AK 221289 V1.1 die zweite methode, siehe morris 1963 */
/* Math. Zeit 81 112-123 (1963) */
/* AK 200891 V1.3 */
{
    OP c,d,e,ff,g,z;
    INT erg = OK; /* AK 180893 */
    INT i;
    CE2(a,b,hall_littlewood);
    CTO(PARTITION,"hall_littlewood(1)",a);
    if (S_PA_LI(a) == (INT)1)
        {
        erg += b_skn_s(callocobject(),callocobject(),NULL,b);
        erg += copy(a,S_S_S(b));
        erg += b_skn_po(callocobject(),callocobject(),NULL,S_S_K(b));
        M_I_I((INT)1,S_PO_K(S_S_K(b)));
        erg += m_il_v((INT)1,S_PO_S(S_S_K(b)));
        M_I_I((INT)0,S_PO_SI(S_S_K(b),(INT)0));
        goto hl_ende;
        }
    /* wenn die laenge groesser 1 ist */

    erg += init(SCHUR,b);

    c = callocobject(); d = callocobject(); e = callocobject(); g = callocobject();

    erg += copy_partition(a,c);
    erg += dec_partition(c);
    erg += hall_littlewood(c,d);
    erg += weight_partition(c,e);
    erg += copy(d,c);
    z = c;
    while (z != NULL)
        {
        erg += inc_partition(S_S_S(z));
        M_I_I(    S_PA_II(a,S_PA_LI(a)-(INT)1), S_S_SI(z,S_S_SLI(z)-(INT)1));
        z = S_S_N(z);
        }

    ff = callocobject();
    erg += reorder_hall_littlewood(c,ff);
    insert(ff,b,NULL,NULL);

    erg += copy(d,c);
    for (i=(INT)1;i<=S_I_I(e); i++)
        {
        erg += m_i_pa(e,g);
        M_I_I(i,S_PA_I(g,(INT)0));
        z = c; /* c ist das ergebnis der rekursion */
        erg += init(SCHUR,d);
        while (z != NULL)
            {
            ff = callocobject();
            erg += part_part_skewschur(S_S_S(z),g,ff);
            if (not NULLP(ff))
                {
                MULT_APPLY(S_S_K(z),ff);
                INSERT_LIST(ff,d,add_koeff,comp_monomschur);
                }
            else
                erg += freeall(ff);
            z = S_S_N(z);
            }
    /* d ist nun die liste mit den expansion der skewpartition */
        z = d;
    /* nun noch die multiplikation mit t^i */
        erg += b_skn_po(callocobject(),callocobject(),NULL,g);
        M_I_I((INT)1,S_PO_K(g));
        erg += m_il_v((INT)1,S_PO_S(g));
        M_I_I(i,S_PO_SI(g,(INT)0));
        while (z != NULL)
            {
            erg += inc_partition(S_S_S(z));
            M_I_I(S_PA_II(a,S_PA_LI(a)-(INT)1)+i, S_S_SI(z,S_S_SLI(z)-(INT)1));
            erg += mult_apply(g,S_S_K(z));
            z = S_S_N(z);
            }

        ff = callocobject();
        erg += reorder_hall_littlewood(d,ff);
        erg += insert(ff,b,NULL,NULL);
        }
    erg += freeall(e);
    erg += freeall(d);
    erg += freeall(c);
    erg += freeall(g);
hl_ende:
    ENDR("hall_littlewood");
}



INT copy_monomial(a,b) OP a,b;
/* AK 270901 */
{
    INT erg = OK;
    CTO(MONOMIAL,"copy_monomial(1)",a);
    erg += transformlist(a,b,copy_monom);
    ENDR("copy_monomial");
}

INT copy_schur(a,b) OP a,b;
/* AK 270901 */
{
    INT erg = OK;
    CTO(SCHUR,"copy_schur(1)",a);
    erg += transformlist(a,b,copy_monom);
    ENDR("copy_schur");
}

INT copy_homsym(a,b) OP a,b;
/* AK 270901 */
{
    INT erg = OK;
    CTO(HOMSYM,"copy_homsym(1)",a);
    erg += transformlist(a,b,copy_monom);
    ENDR("copy_homsym");
}

INT copy_elmsym(a,b) OP a,b;
/* AK 270901 */
{
    INT erg = OK;
    CTO(ELMSYM,"copy_elmsym(1)",a);
    erg += transformlist(a,b,copy_monom);
    ENDR("copy_elmsym");
}

INT copy_powsym(a,b) OP a,b;
/* AK 270901 */
{
    INT erg = OK;
    CTO(POWSYM,"copy_powsym(1)",a);
    erg += transformlist(a,b,copy_monom);
    ENDR("copy_powsym");
}

INT add_apply_symfunc_symfunc(a,b) OP a,b;
/* AK 200891 V1.3 */
    {
    OP c = callocobject();
    copy_polynom(a,c);
    return(insert(c,b,add_koeff,comp_monomvector_monomvector));
    }



INT add_apply_symfunc(a,b) OP a,b;
/* result ist von typ b, falls beides sym func */
{
    OP c;
    INT erg = OK;

    if (S_O_K(a) == S_O_K(b))
        erg += add_apply_symfunc_symfunc(a,b);
    else {
        c = CALLOCOBJECT();
        SWAP(b,c);
        add(c,a,b);
        FREEALL(c);
        }
    ENDR("add_apply_symfunc");
}


INT dimension_schur(a,b) OP a,b;
/* AK 020890 V1.1 */ /* AK 200891 V1.3 */
/* AK 260198 V2.0 */
/*
   input: schur ( may be SCHUR or HASHTABLE )
   output: dimension of corresponding representation of sn
*/
{
    OP z,res;
    INT erg = OK;
    CTTO(HASHTABLE,SCHUR,"dimension_schur(1)",a);
    CE2(a,b,dimension_schur);

    /* b is freed */

    res = CALLOCOBJECT();
    M_I_I(0,b);
    FORALL(z,a,
        {
        erg += dimension(S_MO_S(z),res);
        MULT_APPLY(S_MO_K(z),res);
        ADD_APPLY(res,b);
        } );

    FREEALL(res);
    ENDR("dimension_schur");
}




INT add_staircase_part(a,n,b) OP a,n,b;
/* adds the vector 0,1,...,n-1 to the partition a */
/* AK 050990 V1.1 */
/* AK 200891 V1.3 */
    {
    OP c = callocobject();
    INT i,j;
    m_l_v(n,c);
    for (i=S_V_LI(c)-(INT)1,j=S_PA_LI(a)-(INT)1;i>=(INT)0;i--,j--)
        if (j>=(INT)0) M_I_I(S_PA_II(a,j)+i,S_V_I(c,i));
        else M_I_I(i,S_V_I(c,i));

    b_ks_pa(VECTOR,c,b);
    return OK;
    }



INT mod_part(a,b,c) OP a,b,c;
/* the single parts of partition a mod b gives c */
/* AK 050990 V1.1 */
/* AK 200891 V1.3 */
    {
    INT i;
    if (a != c) copy(a,c);
    for (i=0;i<S_PA_LI(c);i++)
        M_I_I(S_PA_II(c,i) % S_I_I(b), S_PA_I(c,i));
    return OK;
    }



INT p_root_schur(a,n,p,b) OP a,n,p,b;
/* a ist schur n ist integer p ist integer b wird schur */
/* AK 050990 V1.1 */
/* wie folgt wird gerechnet: addiere 0,..,n-2,n-1 zu den
partitionen, die die laenger als n sind werden gestrichen,
dann die partition modulo p, dann wieder zu aufsteigenden
partitionen um sortieren */
/*dies entspricht dem einsetzen von p-ten einheitswurzeln */
/* AK 200891 V1.3 */
    {
    OP z,c,d;
    if (a == b ) { c = callocobject(); copy(a,c);
        p_root_schur(c,n,p,b); freeall(c); return OK; }

    z = a;
    if (not EMPTYP(b) )
        freeself(b);
    b_sn_s(NULL,NULL,b);

    while (z != NULL)
        {
        if (S_S_SLI(z) <= S_I_I(n))
            {
            c = callocobject();
            d = callocobject();
            p_root_part(S_S_S(z),n,p,c);
            b_skn_s(c,callocobject(),NULL,d);
            copy(S_S_K(z),S_S_K(d));
            insert(d,b,NULL,NULL);
            }
        z = S_S_N(z);
        }
    reorder_hall_littlewood(b,b);
    return OK;
    }



INT p_root_part(a,n,p,b) OP a,n,p,b;
/* a ist part, n ist integer, p ist integer */
/* AK 050990 V1.1 */
/* AK 200891 V1.3 */
    {
    INT i;
    OP c = callocobject();
    m_l_v(n,c);
    for (i=(INT)0; i<S_V_LI(c); i++) M_I_I(i,S_V_I(c,i));
    add_staircase_part(a,n,b);
    mod_part(b,p,b);
    sub(S_PA_S(b),c,S_PA_S(b));
    freeall(c); return OK;
    }



static INT m_int_qelm(a,b) INT a; OP b;
{
    INT erg = OK;
    COP("m_int_qelm(1)",b);
    erg += b_skn_e(callocobject(),callocobject(),NULL,b);
    erg += m_i_i(1L , S_S_K(b));
    erg += m_int_pa(a,S_S_S(b));
    ENDR("m_int_qelm");
}

static INT m_int_int_qelm(a,b,d) INT a,b; OP d;
{
    OP c = callocobject();
    INT erg = OK;
    COP("m_int_int_qelm(3)",d);
    SYMCHECK( (a>b) , "m_int_int_qelm: para1 > para2");

    erg += b_ks_pa(VECTOR,callocobject(),c);
    erg += m_il_v(2L,S_PA_S(c));
    C_O_K(S_PA_S(c),INTEGERVECTOR);
    erg += m_i_i(a,S_PA_I(c,0));
    erg += m_i_i(b,S_PA_I(c,1));
    erg += m_part_qelm(c,d);
    erg += freeall(c);
    ENDR("m_int_int_qelm");

}

INT m_part_qelm(a,b) OP a,b;
/* AK 060995 */
/* computes q polynomial as elmsym */
{
    INT i,j;
    OP c,d,e;
    INT erg = OK;

    CTO(PARTITION,"m_part_qelm",a);
    if (S_PA_LI(a) == 1)
        {
        erg += m_int_qelm(S_PA_II(a,0),b);
        }
    else if (S_PA_LI(a) == 2)
        {
        c = callocobject();
        erg += m_int_qelm(S_PA_II(a,0),c);
        d = callocobject();
        erg += m_int_qelm(S_PA_II(a,1),d);
        erg += mult(c,d,b);
        e = callocobject();
        for (i=1;i<=S_PA_II(a,0);i++)
            {
            erg += m_int_qelm(S_PA_II(a,0)-i,c);
            erg += m_int_qelm(S_PA_II(a,1)+i,d);
            erg += mult(c,d,e);
            erg += mult_apply(cons_zwei,e);
            if (i%2 == 1)
                erg += mult_apply(cons_negeins,e);
            erg += add_apply(e,b);
            }
        erg += freeall(c);
        erg += freeall(d);
        erg += freeall(e);
        }
    else if  (S_PA_LI(a) %2 == 0)
        {
        c = callocobject();
        erg += m_ilih_m(S_PA_LI(a), S_PA_LI(a),c);
        for (i=0;i<S_M_HI(c);i++)
            for (j=i+1;j<S_M_LI(c);j++)
                {
                m_int_int_qelm(S_PA_II(a,S_PA_LI(a)-1-j),
                    S_PA_II(a,S_PA_LI(a)-1-i),S_M_IJ(c,i,j));
                }
        pfaffian_matrix(c,b);
        erg += freeall(c);
        }
    else    {
        d = callocobject();
        b_ks_pa(VECTOR,callocobject(),d);
        m_il_nv(S_PA_LI(a)+1,S_PA_S(d));
        C_O_K(S_PA_S(d),INTEGERVECTOR);
        for (i=0;i<S_PA_LI(a);i++)
            M_I_I(S_PA_II(a,i),S_PA_I(d,i+1));

        c = callocobject();
        erg += m_ilih_m(S_PA_LI(d), S_PA_LI(d),c);
        for (i=0;i<S_M_HI(c);i++)
            for (j=i+1;j<S_M_LI(c);j++)
                {
                m_int_int_qelm(S_PA_II(d,S_PA_LI(d)-1-j),
                    S_PA_II(d,S_PA_LI(d)-1-i),S_M_IJ(c,i,j));
                }
        pfaffian_matrix(c,b);
        erg += freeall(c);
        freeall(d);
        }
    ENDR("m_part_qelm");
}

INT m_i_powsym(a,b) INT a; OP b;
/* changes a INT into a POWSYMpolynomial with this INT as
koeffizent and labeled by the part with one zero part */
{
    OP c;
    INT erg = OK;
    COP("m_i_powsym(2)",b);
    c = callocobject();
    erg += m_i_i(a,c);
    erg += m_scalar_powsym(c,b);
    erg += freeall(c);
    ENDR("m_i_elmsym");
}

INT m_i_elmsym(a,b) INT a; OP b;
/* changes a INT into a ELMSYMpolynomial with this INT as
koeffizent and labeled by the part with one zero part */
/* AK 060995 */
{
    OP c;
    INT erg = OK;
    COP("m_i_elmsym(2)",b);
    NEW_INTEGER(c,a);
    erg += m_scalar_elmsym(c,b);
    FREEALL(c);
    ENDR("m_i_elmsym");
}

INT m_i_schur(a,b) INT a; OP b;
/* changes a INT into a SCHURpolynomial with this INT as
koeffizent and labeled by the part with one zero part */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 060498 V2.0 */
{
    OP c = callocobject();
    INT erg = OK;
    erg += m_i_i(a,c);
    erg += m_scalar_schur(c,b);
    FREEALL(c);
    ENDR("m_i_schur");
}



INT m_scalar_powsym(a,b)  OP a,b;
/* changes a scalar into a POWSYMpolynomial with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 060995 */
/* AK 060498 V2.0 */
{
    INT erg = OK;
    CE2(a,b,m_scalar_powsym);
    erg += b_skn_ps(callocobject(),callocobject(),NULL,b);
    COPY(a,S_S_K(b));
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("m_scalar_powsym");
}

INT m_scalar_homsym(a,b)  OP a,b;
/* changes a scalar into a HOMSYMpolynomial with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 060995 */
/* AK 060498 V2.0 */
{
    INT erg = OK;
    CE2(a,b,m_scalar_homsym);
    erg += b_skn_h(callocobject(),callocobject(),NULL,b);
    COPY(a,S_S_K(b));
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("m_scalar_homsym");
}

INT m_scalar_elmsym(a,b)  OP a,b;
/* changes a scalar into a ELMSYMpolynomial with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 060995 */
/* AK 060498 V2.0 */
{
    INT erg = OK;
    CE2(a,b,m_scalar_elmsym);
    erg += b_skn_e(callocobject(),callocobject(),NULL,b);
    COPY(a,S_S_K(b));
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("m_scalar_elmsym");
}

INT m_scalar_schur(a,b)  OP a,b;
/* changes a scalar into a SCHURpolynomial with this scalar as
coefficent and labeled by the part of zero length */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 030498 V2.0 */ /* AK 210704 V3.0 */
{
    INT erg = OK;
    CE2(a,b,m_scalar_schur);
        {
        erg += b_skn_s(CALLOCOBJECT(),CALLOCOBJECT(),NULL,b);
        COPY(a,S_S_K(b));
        erg += first_partition(cons_null,S_S_S(b));
        }
    ENDR("m_scalar_schur");
}

INT b_scalar_schur(a,b)  OP a,b;
/* changes a scalar into a SCHURpolynomial with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 030498 V2.0 */
{
    INT erg = OK;
    if (a == b) {
        erg += error("b_scalar_schur:identical objects");
        goto endr_ende;
        }
    erg += b_skn_s(CALLOCOBJECT(),a,NULL,b);
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("b_scalar_schur");
}

INT b_scalar_powsym(a,b)  OP a,b;
/* changes a scalar into a POWSYM with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 030498 V2.0 */
{
    INT erg = OK;
    if (a == b) {
        erg += error("b_scalar_powsym:identical objects");
        goto endr_ende;
        }
    erg += b_skn_ps(CALLOCOBJECT(),a,NULL,b);
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("b_scalar_powsym");
}

INT b_scalar_homsym(a,b)  OP a,b;
/* changes a scalar into a HOMSYM with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 030498 V2.0 */
{
    INT erg = OK;
    if (a == b) {
        erg += error("b_scalar_homsym:identical objects");
        goto endr_ende;
        }
    erg += b_skn_h(CALLOCOBJECT(),a,NULL,b);
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("b_scalar_homsym");
}

INT b_scalar_elmsym(a,b)  OP a,b;
/* changes a scalar into a ELMSYM with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 030498 V2.0 */
{
    INT erg = OK;
    if (a == b) {
        erg += error("b_scalar_elmsym:identical objects");
        goto endr_ende;
        }
    erg += b_skn_e(CALLOCOBJECT(),a,NULL,b);
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("b_scalar_elmsym");
}

INT b_scalar_monomial(a,b)  OP a,b;
/* changes a scalar into a MONOMIAL with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 030498 V2.0 */
{
    INT erg = OK;
    if (a == b) {
        erg += error("b_scalar_monomial:identical objects");
        goto endr_ende;
        }
    erg += b_skn_mon(CALLOCOBJECT(),a,NULL,b);
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("b_scalar_monomial");
}


INT m_scalar_monomial(a,b)  OP a,b;
/* changes a scalar into a MONOMIAL with this scalar as
koeffizent and labeled by the part with zero length self part */
/* AK 181290 V1.1 */ /* AK 200891 V1.3 */
/* AK 030498 V2.0 */
{
    INT erg = OK;
    CE2(a,b,m_scalar_monomial);
    erg += b_skn_mon(callocobject(),callocobject(),NULL,b);
    COPY(a,S_S_K(b));
    erg += first_partition(cons_null,S_S_S(b));
    ENDR("m_scalar_monomial");
}


INT mult_schur_powsym(a,b,c) OP a,b,c;
/* AK 031001 */
{
    OP d;
    INT erg = OK;
    CTTO(SCHUR,PARTITION,"mult_schur_powsym(1)",a);
    CTTO(POWSYM,PARTITION,"mult_schur_powsym(2)",b);
    CTO(EMPTY,"mult_schur_powsym(3)",c);

    d = CALLOCOBJECT();
    erg += t_SCHUR_POWSYM(a,d);
    erg += mult_powsym_powsym(d,b,c);
    FREEALL(d);
    ENDR("mult_schur_powsym");
}

INT mult_powsym_monomial(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTO(POWSYM,PARTITION,"mult_powsym_monomial(1)",a);
    CTTO(MONOMIAL,PARTITION,"mult_powsym_monomial(2)",b);
    CTO(EMPTY,"mult_powsym_monomial(3)",c);

    d = callocobject();
    erg += t_POWSYM_MONOMIAL(a,d);
    erg += mult_monomial_monomial(d,b,c);
    erg += freeall(d);
    ENDR("mult_powsym_monomial");
}

INT mult_powsym_homsym(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTO(POWSYM,PARTITION,"mult_powsym_homsym",a);
    CTTO(HOMSYM,PARTITION,"mult_powsym_homsym",b);
    CTO(EMPTY,"mult_powsym_homsym",c);

    d = callocobject();
    erg += t_POWSYM_HOMSYM(a,d);
    erg += mult_homsym_homsym(d,b,c);
    erg += freeall(d);
    ENDR("mult_powsym_homsym");
}

INT mult_elmsym_homsym(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTTO(HASHTABLE,ELMSYM,PARTITION,"mult_elmsym_homsym(1)",a);
    CTTTO(HASHTABLE,HOMSYM,PARTITION,"mult_elmsym_homsym(2)",b);
    CTTTO(HOMSYM,HASHTABLE,EMPTY,"mult_elmsym_homsym(3)",c);

    NEW_HASHTABLE(d);
    erg += t_ELMSYM_HOMSYM(a,d);
    erg += mult_homsym_homsym(d,b,c);
    FREEALL(d);
    ENDR("mult_elmsym_homsym");
}

INT mult_monomial_homsym(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTO(MONOMIAL,PARTITION,"mult_monomial_homsym(1)",a);
    CTTO(HOMSYM,PARTITION,"mult_monomial_homsym(2)",b);
    CTO(EMPTY,"mult_monomial_homsym",c);

    d = callocobject();
    erg += t_MONOMIAL_HOMSYM(a,d);
    erg += mult_homsym_homsym(d,b,c);
    erg += freeall(d);
    ENDR("mult_monomial_homsym");
}


INT mult_powsym_elmsym(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTO(POWSYM,PARTITION,"mult_powsym_elmsym",a);
    CTTO(ELMSYM,PARTITION,"mult_powsym_elmsym",b);
    CTO(EMPTY,"mult_powsym_elmsym",c);

    d = callocobject();
    erg += t_POWSYM_ELMSYM(a,d);
    erg += mult_elmsym_elmsym(d,b,c);
    erg += freeall(d);
    ENDR("mult_powsym_elmsym");
}

INT mult_monomial_elmsym(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTO(MONOMIAL,PARTITION,"mult_monomial_elmsym",a);
    CTTO(ELMSYM,PARTITION,"mult_monomial_elmsym",b);
    CTO(EMPTY,"mult_monomial_elmsym",c);

    d = CALLOCOBJECT();init_hashtable(d);
    erg += t_MONOMIAL_ELMSYM(a,d);
    erg += mult_elmsym_elmsym(d,b,c);
    erg += freeall(d);
    ENDR("mult_monomial_elmsym");
}






INT mult_elmsym_powsym(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTO(ELMSYM,PARTITION,"mult_elmsym_powsym",a);
    CTTO(POWSYM,PARTITION,"mult_elmsym_powsym",b);
    CTO(EMPTY,"mult_elmsym_powsym",c);

    d = callocobject();
    erg += t_ELMSYM_POWSYM(a,d);
    erg += mult_powsym_powsym(d,b,c);
    erg += freeall(d);
    ENDR("mult_elmsym_powsym");
}

INT mult_monomial_powsym(a,b,c) OP a,b,c;
/* AK 081001 */
{
    OP d;
    INT erg = OK;
    CTTO(MONOMIAL,PARTITION,"mult_monomial_powsym",a);
    CTTO(POWSYM,PARTITION,"mult_monomial_powsym",b);
    CTO(EMPTY,"mult_monomial_powsym",c);

    d = callocobject();
    erg += t_MONOMIAL_POWSYM(a,d);
    erg += mult_powsym_powsym(d,b,c);
    erg += freeall(d);
    ENDR("mult_monomial_powsym");
}

#define ADDINVERS_APPLY_SF(a) \
    { \
    OP z; \
    FORALL(z,a,ADDINVERS_APPLY(S_MO_K(z))); \
    }

INT addinvers_apply_schur(a) OP a;
/* AK 290402 */
{
    INT erg = OK;
    CTTO(SCHUR,HASHTABLE,"addinvers_apply_schur(1)",a);
    ADDINVERS_APPLY_SF(a);
    ENDR("addinvers_apply_schur");
}

INT addinvers_apply_homsym(a) OP a;
/* AK 290402 */
{
    INT erg = OK;
    CTTO(HOMSYM,HASHTABLE,"addinvers_apply_homsym(1)",a);
    ADDINVERS_APPLY_SF(a);
    ENDR("addinvers_apply_homsym");
}

INT addinvers_apply_powsym(a) OP a;
/* AK 290402 */
{
    INT erg = OK;
    CTTO(POWSYM,HASHTABLE,"addinvers_apply_powsym(1)",a);
    ADDINVERS_APPLY_SF(a);
    ENDR("addinvers_apply_powsym");
}

INT addinvers_apply_elmsym(a) OP a;
/* AK 290402 */
{
    INT erg = OK;
    CTTO(ELMSYM,HASHTABLE,"addinvers_apply_elmsym(1)",a);
    ADDINVERS_APPLY_SF(a);
    ENDR("addinvers_apply_elmsym");
}

INT addinvers_apply_monomial(a) OP a;
/* AK 290402 */
{
    INT erg = OK;
    CTTO(MONOMIAL,HASHTABLE,"addinvers_apply_monomial(1)",a);
    ADDINVERS_APPLY_SF(a);
    ENDR("addinvers_apply_monomial");
}



#define ADD_SF(a,b,c,typ,cf,addf) \
    switch((int)S_O_K(b)) {\
        case typ: erg += addf(a,b,c); goto ende;\
        default: {\
            OP add_sf_c;\
            add_sf_c = CALLOCOBJECT();\
            cf(b,add_sf_c); \
            erg += addf(a,add_sf_c,c); \
            FREEALL(add_sf_c);\
            goto ende; }\
        }\
ende:

INT add_schur(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(SCHUR,"add_schur(1)",a);
    CTO(EMPTY,"add_schur(3)",c);
    ADD_SF(a,b,c,SCHUR,cast_schur,add_schur_schur);
    CTO(SCHUR,"add_schur(e3)",c);
    ENDR("add_schur");
}

INT add_monomial(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(MONOMIAL,"add_monomial(1)",a);
    CTO(EMPTY,"add_monomial(3)",c);
    ADD_SF(a,b,c,MONOMIAL,cast_monomial,add_monomial_monomial);
    CTO(MONOMIAL,"add_monomial(e3)",c);
    ENDR("add_monomial");
}
INT add_homsym(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(HOMSYM,"add_homsym(1)",a);
    CTO(EMPTY,"add_homsym(3)",c);
    ADD_SF(a,b,c,HOMSYM,cast_homsym,add_homsym_homsym);
    CTO(HOMSYM,"add_homsym(e3)",c);
    ENDR("add_homsym");
}
INT add_powsym(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(POWSYM,"add_powsym(1)",a);
    CTO(EMPTY,"add_powsym(3)",c);
    ADD_SF(a,b,c,POWSYM,cast_powsym,add_powsym_powsym);
    CTO(POWSYM,"add_powsym(e3)",c);
    ENDR("add_powsym");
}
INT add_elmsym(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(ELMSYM,"add_elmsym(1)",a);
    CTO(EMPTY,"add_elmsym(3)",c);
    ADD_SF(a,b,c,ELMSYM,cast_elmsym,add_elmsym_elmsym);
    CTO(ELMSYM,"add_elmsym(e3)",c);
    ENDR("add_elmsym");
}





INT mult_apply_monomial_monomial(a,b) OP a,b;
/* platzhalter */
/* b = b*a */
{
    INT erg = OK;
    OP c;
    CTTO(HASHTABLE,MONOMIAL,"mult_apply_monomial_monomial",a);
    CTTO(HASHTABLE,MONOMIAL,"mult_apply_monomial_monomial",b);

    c = CALLOCOBJECT();
    if (S_O_K(b) == HASHTABLE) erg += init_hashtable(c);
    erg += mult_monomial_monomial(a,b,c);
    SWAP(c,b);
    FREEALL(c);

    ENDR("mult_apply_monomial_monomial");
}
INT mapp_hashtable_hashtable_();

INT mult_apply_powsym_powsym(a,b) OP a,b;
/* AK 060901 */
/* a and b of equal type */
/* b := a * b */
/* AK 160603 case of empty lists */
{
    INT erg = OK;
    OP z,y;
    CTTO(HASHTABLE,POWSYM,"mult_apply_powsym_powsym(1)",a);
    CTTO(HASHTABLE,POWSYM,"mult_apply_powsym_powsym(2)",b);

    if (S_O_K(a) == HASHTABLE) {
        erg += mapp_hashtable_hashtable_(a,b,cons_eins);
        goto ende;
        }

    if (S_O_K(b) == HASHTABLE) {
        erg += mapp_hashtable_hashtable_(a,b,cons_eins);
        goto ende;
        }

    if (S_L_S(b) == NULL)  goto endr_ende;  /* AK 160603 */
    if (S_L_S(a) == NULL)  { init(POWSYM,b); goto endr_ende; }  /* AK 160603 */

    z = a;
    if (S_S_N(z) == NULL) /* only one summand in a */
        {
        y = b;
        while (y != NULL)
            {
            MULT_APPLY(S_S_K(z),S_S_K(y));
            erg += append_apply_part(S_S_S(y),S_S_S(z)); /* y = y+z */
            y = S_S_N(y);
            }
        goto endr_ende;
        }


    /* the function b must be copied */
    z = CALLOCOBJECT();
    *z = *b;
    C_O_K(b,EMPTY);
    erg += mult_powsym_powsym(a,z,b);
    FREEALL(z);
ende:
    ENDR("mult_apply_powsym_powsym");
}


INT mult_apply_elmsym_elmsym(a,b) OP a,b;
/* b := a * b */
/* AK 161001 */
{
    INT erg = OK;
    OP z,y;
    CTTTO(HASHTABLE,PARTITION,ELMSYM,"mult_apply_elmsym_elmsym(1)",a);
    CTTO(ELMSYM,HASHTABLE,"mult_apply_elmsym_elmsym(2)",b);
    if (NULLP(b)) goto ende;
    if (NULLP(a)) {
        if (S_O_K(b) == HASHTABLE) CLEAR_HASHTABLE(b);
        else erg += init(ELMSYM,b);
        goto ende;
        }

    if (
        (S_O_K(b) == HASHTABLE) ||
        (S_O_K(a) == HASHTABLE)
       )
        {
        z = CALLOCOBJECT();
        *z = *b;
        C_O_K(b,EMPTY);
        init(S_O_K(z),b);
        erg += mult_elmsym_elmsym(a,z,b);
        FREEALL(z);
        goto endr_ende;
        }
    if (S_O_K(a) == PARTITION)
        {
        y = b;
        while (y != NULL)
            {
            erg += append_apply_part(S_S_S(y),a); /* y = y+a */
            y = S_S_N(y);
            }
        goto endr_ende;
        }
    /* a is ELMSYM */
    z = a;
    if (S_S_N(z) == NULL) /* only one summand in a */
        {
        y = b;
        while (y != NULL)
            {
            erg += mult_apply(S_S_K(z),S_S_K(y));
            erg += append_apply_part(S_S_S(y),S_S_S(z)); /* y = y+z */
            y = S_S_N(y);
            }
        goto endr_ende;
        }
    /* the function b must be copied */
    z = CALLOCOBJECT();
    *z = *b;
    C_O_K(b,EMPTY);
    erg += mult_elmsym_elmsym(a,z,b);
    erg += freeall(z);
ende:
    ENDR("mult_apply_elmsym_elmsym");
}

INT mult_apply_homsym_schur(a,b) OP a,b;
/* b := a * b */
/* AK 221001 */
{
    INT erg = OK;
    OP c;

    CTTTO(HASHTABLE,PARTITION,HOMSYM,"mult_apply_homsym_schur(1)",a);
    CTTTO(HASHTABLE,PARTITION,SCHUR,"mult_apply_homsym_schur(2)",b);
    c = CALLOCOBJECT();
    *c = *b;
    C_O_K(b,EMPTY);

    if (S_O_K(c) == HASHTABLE) init_hashtable(b);
    else init_schur(b);

    erg += mult_homsym_schur(a,c,b);
    FREEALL(c);
    ENDR("mult_apply_homsym_schur");
}

INT mahh_partition_hashtable_(a,b,f) OP a,b;OP f;
/* AK 240102 */
{
    INT erg = OK;
    OP z,m;
    if (mahh_h==NULL) {
        NEW_HASHTABLE(mahh_h);
        }

    TCTO(PARTITION,"mahh_partition_hashtable_(1)",a);
    TCTO(HASHTABLE,"mahh_partition_hashtable_(2)",b);
    TCTO(HASHTABLE,"mahh_partition_hashtable_(i3)",mahh_h);
    FORALL(z,b, {
        m = CALLOCOBJECT();
        *m = *z;
        C_O_K(z,EMPTY);
        append_apply_part(S_MO_S(m),a); /* m = m + a */
        MULT_APPLY(f,S_MO_K(m));
        insert_scalar_hashtable(m,mahh_h,add_koeff,eq_monomsymfunc,hash_monompartition);
        DEC_INTEGER(S_V_I(b,S_V_LI(b)));
        });
    TCTO(HASHTABLE,"mahh_partition_hashtable_(i2)",b);
    SYMCHECK(WEIGHT_HASHTABLE(b) != 0,"mahh_partition_hashtable_:i1");
    SWAP(mahh_h,b);

    ENDR("mahh_partition_hashtable_");
}


INT mahh_hashtable_hashtable_(a,b,f) OP a,b,f;
/* AK 240102 */
{
    INT erg = OK;
    OP c;
    INT mhh_partition_partition_();
    CTTO(HOMSYM,HASHTABLE,"mahh_hashtable_hashtable_(1)",a);
    CTO(HASHTABLE,"mahh_hashtable_hashtable_(2)",b);
    NEW_HASHTABLE(c);
    SWAP(b,c);
    M_FORALL_MONOMIALS_IN_AB(a,c,b,f,mhh_partition_partition_);
    FREEALL(c);
    ENDR("mahh_hashtable_hashtable_");
}

INT mapp_hashtable_hashtable_(a,b,f) OP a,b,f;
/* AK 240102 */
{
    INT erg = OK;
    OP c;
    INT mpp_partition_partition_();
    CTTO(POWSYM,HASHTABLE,"mapp_hashtable_hashtable_(1)",a);
    CTTO(HASHTABLE,POWSYM,"mapp_hashtable_hashtable_(2)",b);
    NEW_HASHTABLE(c);
    SWAP(b,c);
    M_FORALL_MONOMIALS_IN_AB(a,c,b,f,mpp_partition_partition_);
    FREEALL(c);
    ENDR("mapp_hashtable_hashtable_");
}


INT mahh_integer_homsym_(a,b,f) OP a,b,f;
/* AK 290102 */
{
    INT erg = OK;
    OP z;
    CTO(INTEGER,"mahh_integer_homsym_(1)",a);
    CTO(HOMSYM,"mahh_integer_homsym_(2)",b);
    FORALL(z,b,{
        erg += append_apply_part(S_MO_S(z),a);
        if (not EINSP(f)) {
             MULT_APPLY(f,S_MO_K(z));
             }
        });
    ENDR("mahh_integer_homsym_");
}

INT mult_apply_homsym_homsym(a,b) OP a,b;
/* b := a * b */
/* AK 161001 */
{
    INT erg = OK;
    OP z,y;
    CTTTTO(INTEGER,HASHTABLE,PARTITION,HOMSYM,"mult_apply_homsym_homsym(1)",a);
    CTTO(HASHTABLE,HOMSYM,"mult_apply_homsym_homsym(2)",b);

    if (S_O_K(a) == INTEGER)
        {
        OP c;
        if (S_O_K(b) == HOMSYM)
             {
             erg += mahh_integer_homsym_(a,b,cons_eins);
             goto ende;
             }
        else {
             c = CALLOCOBJECT();
             m_i_pa(a,c);
             erg += mult_apply_homsym_homsym(c,b);
             FREEALL(c);
             goto ende;
             }
        }
    if (S_O_K(b) == HASHTABLE)
        {
        if (S_O_K(a) == PARTITION)  erg += mahh_partition_hashtable_(a,b,cons_eins);
        else erg+=mahh_hashtable_hashtable_(a,b,cons_eins);
        goto ende;
        }
    if (S_O_K(a) == HASHTABLE)
        {
        z = CALLOCOBJECT();
        init_homsym(z);
        SWAP(z,b);
        erg += mult_homsym_homsym(a,z,b);
        FREEALL(z);
        goto endr_ende;
        }
    /* no more hashtables */
    if (S_L_S(b) == NULL) { /* null */
        goto endr_ende;
        }
    if (S_O_K(a) == PARTITION)
        {
        y = b;
        while (y != NULL)
            {
            erg += append_apply_part(S_S_S(y),a); /* y = y+a */
            y = S_S_N(y);
            }
    goto endr_ende;
        }
    /* a is HOMSYM */

    if (NULLP(a)) {
        erg += init(HOMSYM,b);
        goto endr_ende;
        }

    z = a;
    if (S_S_N(z) == NULL) /* only one summand in a */
        {
        y = b;
        while (y != NULL)
            {
            MULT_APPLY(S_S_K(z),S_S_K(y));
            erg += append_apply_part(S_S_S(y),S_S_S(z)); /* y = y+z */
            y = S_S_N(y);
            }
    goto endr_ende;
        }

    /* the function b must be copied */
    z = callocobject();
    *z = *b;
    C_O_K(b,EMPTY);
    erg += mult_homsym_homsym(a,z,b);
    erg += freeall(z);
ende:
    ENDR("mult_apply_homsym_homsym");
}

INT mult_apply_homsym(a,b) OP a,b;
/* AK 311001 */
{
    INT erg = OK;
    CTO(HOMSYM,"mult_apply_homsym(1)",a);
    if (S_O_K(b) == HOMSYM)
        erg += mult_apply_homsym_homsym(a,b);
    else if (S_O_K(b) == SCHUR)
        erg += mult_apply_homsym_schur(a,b);
    else
        {
        OP c = CALLOCOBJECT();
        *c = *b;
        C_O_K(b,EMPTY);
        erg += mult_homsym(a,c,b);
        FREEALL(c);
        }
    ENDR("mult_apply_homsym");
}

INT mult_apply_powsym(a,b) OP a,b;
/* AK 311001 */
{
    INT erg = OK;
    CTO(POWSYM,"mult_apply_powsym(1)",a);
    if (S_O_K(b) == POWSYM)
        erg += mult_apply_powsym_powsym(a,b);
    else
        {
        OP c = CALLOCOBJECT();
        *c = *b;
        C_O_K(b,EMPTY);
        erg += mult_powsym(a,c,b);
        FREEALL(c);
        }
    ENDR("mult_apply_powsym");
}

INT mult_apply_elmsym(a,b) OP a,b;
/* AK 311001 */
{
    INT erg = OK;
    CTO(ELMSYM,"mult_apply_elmsym(1)",a);
    if (S_O_K(b) == ELMSYM)
        erg += mult_apply_elmsym_elmsym(a,b);
    else
        {
        OP c = CALLOCOBJECT();
        *c = *b;
        C_O_K(b,EMPTY);
        erg += mult_elmsym(a,c,b);
        FREEALL(c);
        }
    ENDR("mult_apply_elmsym");
}

INT mult_apply_schur(a,b) OP a,b;
/* AK 311001 */
/* b = a*b */
{
    INT erg = OK;
    CTO(SCHUR,"mult_apply_schur(1)",a);
    if (S_O_K(b) == SCHUR)
        erg += mult_apply_schur_schur(a,b);
    else
        {
        OP c = CALLOCOBJECT();
        *c = *b;
        C_O_K(b,EMPTY);
        erg += mult_schur(a,c,b);
        FREEALL(c);
        }
    ENDR("mult_apply_schur");
}

INT mult_apply_monomial(a,b) OP a,b;
/* AK 270901 */
/* a is monomial
   b is unknown
   they are different */
{
    INT erg = OK;
    OP c;
    CTO(MONOMIAL,"mult_apply_monomial(1)",a);
    c = callocobject();
    erg += swap(b,c);
    erg += mult_monomial(a,c,b);
    erg += freeall(c);
    ENDR("mult_apply_monomial");
}


#define MULT_SF_SCALAR(a,b,c) \
    {\
    OP z;\
    if (NULLP(b)) init(S_O_K(a),c);\
    else {\
         COPY(a,c);\
         FORALL(z,c,MULT_APPLY(b,S_MO_K(z)));\
         }\
    }

INT mult_schur_scalar(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(SCHUR,"mult_schur_scalar(1)",a);
    CTO(EMPTY,"mult_schur_scalar(3)",c);
    MULT_SF_SCALAR(a,b,c);
    ENDR("mult_schur_scalar");
}
INT mult_monomial_scalar(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(MONOMIAL,"mult_monomial_scalar(1)",a);
    CTO(EMPTY,"mult_monomial_scalar(3)",c);
    MULT_SF_SCALAR(a,b,c);
    ENDR("mult_monomial_scalar");
}
INT mult_homsym_scalar(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(HOMSYM,"mult_homsym_scalar(1)",a);
    CTO(EMPTY,"mult_homsym_scalar(3)",c);
    MULT_SF_SCALAR(a,b,c);
    ENDR("mult_homsym_scalar");
}
INT mult_powsym_scalar(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(POWSYM,"mult_powsym_scalar(1)",a);
    CTO(EMPTY,"mult_powsym_scalar(3)",c);
    MULT_SF_SCALAR(a,b,c);
    ENDR("mult_powsym_scalar");
}
INT mult_elmsym_scalar(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(ELMSYM,"mult_elmsym_scalar(1)",a);
    CTO(EMPTY,"mult_elmsym_scalar(3)",c);
    MULT_SF_SCALAR(a,b,c);
    ENDR("mult_elmsym_scalar");
}

#define MULT_SF(a,b,c,scalarf,schurf,homsymf,monomialf,elmsymf,powsymf,t)\
    switch(S_O_K(b)) \
        {\
        case INTEGER:\
        case LONGINT:\
        case FF:\
        case CYCLOTOMIC:\
        case SQ_RADICAL:\
        case POLYNOM:\
        case BRUCH: erg += scalarf(a,b,c); goto endr_ende;\
        case SCHUR: erg += schurf(a,b,c); goto endr_ende;\
        case HOMSYM: erg += homsymf(a,b,c); goto endr_ende;\
        case MONOMIAL: erg += monomialf(a,b,c); goto endr_ende;\
        case ELMSYM: erg += elmsymf(a,b,c); goto endr_ende;\
        case POWSYM: erg += powsymf(a,b,c); goto endr_ende;\
        default: erg += WTO(t,b); goto endr_ende;\
        }

INT mult_schur(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(SCHUR,"mult_schur(1)",a);
    CTO(EMPTY,"mult_schur(3)",c);
    MULT_SF(a,b,c,mult_schur_scalar,mult_schur_schur,mult_schur_homsym,
                  mult_schur_monomial,mult_schur_elmsym,mult_schur_powsym,
                  "mult_schur(2)");
    ENDR("mult_schur");
}

INT mult_homsym(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(HOMSYM,"mult_homsym(1)",a);
    CTO(EMPTY,"mult_homsym(3)",c);
    MULT_SF(a,b,c,mult_homsym_scalar,mult_homsym_schur,mult_homsym_homsym,
                  mult_homsym_monomial,mult_homsym_elmsym,mult_homsym_powsym,
                  "mult_homsym(2)");
    ENDR("mult_homsym");
}
INT mult_powsym(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(POWSYM,"mult_powsym(1)",a);
    CTO(EMPTY,"mult_powsym(3)",c);
    MULT_SF(a,b,c,mult_powsym_scalar,mult_powsym_schur,mult_powsym_homsym,
                  mult_powsym_monomial,mult_powsym_elmsym,mult_powsym_powsym,
                  "mult_powsym(2)");
    ENDR("mult_powsym");
}
INT mult_elmsym(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(ELMSYM,"mult_elmsym(1)",a);
    CTO(EMPTY,"mult_elmsym(3)",c);
    MULT_SF(a,b,c,mult_elmsym_scalar,mult_elmsym_schur,mult_elmsym_homsym,
                  mult_elmsym_monomial,mult_elmsym_elmsym,mult_elmsym_powsym,
                  "mult_elmsym(2)");
    ENDR("mult_elmsym");
}
INT mult_monomial(a,b,c) OP a,b,c;
/* AK 080102 */
{
    INT erg = OK;
    CTO(MONOMIAL,"mult_monomial(1)",a);
    CTO(EMPTY,"mult_monomial(3)",c);
    MULT_SF(a,b,c,mult_monomial_scalar,mult_monomial_schur,mult_monomial_homsym,
                  mult_monomial_monomial,mult_monomial_elmsym,mult_monomial_powsym,
                  "mult_monomial(2)");
    ENDR("mult_monomial");
}




INT number_nat_matrices(a,b,c) OP a,b,c;
/* AK 180901 */
/* number of matrices of rowsum a
                         coulmnsum b
   of the non negative numbers */
/* computed using the kostka number of skew shape */
{
    INT erg = OK;
    INT i,w;
    OP inner,outer,spa;
    CTTTO(PARTITION,VECTOR,INTEGERVECTOR,"number_nat_matrices",a);
    CTTTO(PARTITION,VECTOR,INTEGERVECTOR,"number_nat_matrices",b);
    if (S_O_K(a) == VECTOR) {
        for (i=0;i<S_V_LI(a);i++)
            if (S_O_K(S_V_I(a,i)) != INTEGER) {
                error("number_nat_matrices:no integer vector");
                goto endr_ende;
                }
        }
    if (S_O_K(b) == VECTOR) {
        for (i=0;i<S_V_LI(b);i++)
            if (S_O_K(S_V_I(b,i)) != INTEGER) {
                error("number_nat_matrices:no integer vector");
                goto endr_ende;
                }
        }
    if (S_O_K(a) == PARTITION) {
        if (S_PA_K(a) != VECTOR)
            {
            erg += error("number_nat_matrices:only for vector type partitions");
            goto endr_ende;
            }
        erg += number_nat_matrices(S_PA_S(a),b,c);
        goto endr_ende;
        }
    if (S_O_K(b) == PARTITION) {
        if (S_PA_K(b) != VECTOR)
            {
            erg += error("number_nat_matrices:only for vector type partitions");
            goto endr_ende;
            }
        erg += number_nat_matrices(a,S_PA_S(b),c);
        goto endr_ende;
        }

    for (i=0,w=0;i<S_V_LI(a);i++)
        w += S_V_II(a,i);
    inner = callocobject();
    outer = callocobject();

    erg += m_il_v(S_V_LI(a),outer);
    erg += m_il_v(S_V_LI(a)-1,inner);

    for (i=S_V_LI(a)-1;i>0;i--)
        {
        M_I_I(w,S_V_I(outer,i));
        w = w - S_V_II(a,i);
        M_I_I(w,S_V_I(inner,i-1));
        }
    M_I_I(w,S_V_I(outer,i));
    erg += m_v_pa(inner,inner);
    erg += m_v_pa(outer,outer);

    spa = callocobject();
    erg += m_gk_spa(outer,inner,spa);
    erg += freeall(inner);
    erg += freeall(outer);

    erg += kostka_number_skewpartition(b,spa,c);
    erg += freeall(spa);
    ENDR("number_nat_matrices");
}

INT t_ELMSYM_MONOMIAL(a,b) OP a,b;
/* AK 270901 */
/* using multiplication e_I * m_0 -> \sum m_J */
/* fastest up to now */
{
    INT erg = OK;
    OP m;
    CTTTTO(HASHTABLE,INTEGER,PARTITION,ELMSYM,"t_ELMSYM_MONOMIAL",a);
    TCE2(a,b,t_ELMSYM_MONOMIAL,MONOMIAL);

    m=CALLOCOBJECT();

    erg += first_partition(cons_null,m);
    erg += m_pa_mon(m,m);
    erg += mult_elmsym_monomial(a,m,b);

    FREEALL(m);
    ENDR("t_ELMSYM_MONOMIAL");
}

static INT all_01_matrices_rek_160802(a,c,d,i,b) OP a,b,c,d; INT i;
/* AK 160802 */
{
    INT erg=OK;
    if (i>=S_V_LI(a)) {
        INC(b);
        COPY(d,S_V_I(b,S_V_LI(b)-1));
        }
    else {
        OP e;
        INT j;
        erg = ERROR;
        e = callocobject();
        first_subset(S_V_L(c),S_V_I(a,i), e);
        do {
            for (j=0;j<S_V_LI(c);j++)
                {
                if ((S_V_II(e,j) == 1) &&
                        (S_V_II(c,j) == 0) )
                                goto stop;
                }
                /* the composition is ok */
            for (j=0;j<S_V_LI(c);j++)
                {
                if (S_V_II(e,j) == 1)
                    {
                    M_I_I(1,S_M_IJ(d,i,j));
                    DEC_INTEGER(S_V_I(c,j));
                    }
                }
                /* now recursion */
            erg = all_01_matrices_rek_160802(a,c,d,i+1,b);
            for (j=0;j<S_V_LI(c);j++)
                {
                if (S_V_II(e,j) == 1)
                    {
                    M_I_I(0,S_M_IJ(d,i,j));
                    INC_INTEGER(S_V_I(c,j));
                    }
                }
stop:   ;
            } while(next_apply(e));
        freeall(e);
    }
    ENDR("internal routine:all_01_matrices_rek_160802");
}

INT all_01_matrices(a,b,c) OP a,b,c;
/* AK 160802 */
/* generates a vector of 01 matrices */
{
    INT erg = OK;
    CTTTO(PARTITION,VECTOR,INTEGERVECTOR,"all_01_matrices(1)",a);
    CTTTO(PARTITION,VECTOR,INTEGERVECTOR,"all_01_matrices(2)",b);

    CE3(a,b,c,all_01_matrices);
    if (S_O_K(a) == PARTITION) a = S_PA_S(a);
    if (S_O_K(b) == PARTITION) b = S_PA_S(b);
    if (S_O_K(a) == VECTOR) {
        INT i;
        for (i=0;i<S_V_LI(a);i++)
            if (S_O_K(S_V_I(a,i)) != INTEGER) {
                error("all_01_matrices:no integer vector");
                goto endr_ende;
                }
        }
    if (S_O_K(b) == VECTOR) {
        INT i;
        for (i=0;i<S_V_LI(b);i++)
            if (S_O_K(S_V_I(b,i)) != INTEGER) {
                error("all_01_matrices:no integer vector");
                goto endr_ende;
                }
        }

    {
    OP d;
    d = CALLOCOBJECT();
    erg += m_il_v(0,c);
    erg += m_lh_nm(S_V_L(b),S_V_L(a),d);
    C_O_K(d,INTEGERMATRIX);
    erg += all_01_matrices_rek_160802(a,b,d,0,c);
    FREEALL(d);
    }
    ENDR("all_01_matrices");
}

INT number_01_matrices(a,b,c) OP a,b,c;
/* number of 01 matrices of rowsum a
                         coulmnsum b  */
/* computed using the kostka number of skew shape */
/* AK 180901 */
{
    INT erg = OK;
    CTTTO(PARTITION,VECTOR,INTEGERVECTOR,"number_01_matrices(1)",a);
    CTTTO(PARTITION,VECTOR,INTEGERVECTOR,"number_01_matrices(1)",b);

    if (S_O_K(a) == PARTITION) a = S_PA_S(a);
    if (S_O_K(b) == PARTITION) b = S_PA_S(b);
    if (S_O_K(a) == VECTOR) {
        INT i;
        for (i=0;i<S_V_LI(a);i++)
            if (S_O_K(S_V_I(a,i)) != INTEGER) {
                error("number_01_matrices:no integer vector");
                goto endr_ende;
                }
        }
    if (S_O_K(b) == VECTOR) {
        INT i;
        for (i=0;i<S_V_LI(b);i++)
            if (S_O_K(S_V_I(b,i)) != INTEGER) {
                error("number_nat_matrices:no integer vector");
                goto endr_ende;
                }
        }

    {
    INT i,w;
    OP inner,outer,spa;
    for (i=0,w=0;i<S_V_LI(a);i++)
        w += S_V_II(a,i);
    inner = callocobject();
    outer = callocobject();

    erg += m_il_v(S_V_LI(a),outer);
    erg += m_il_v(S_V_LI(a)-1,inner);

    for (i=S_V_LI(a)-1;i>0;i--)
        {
        M_I_I(w,S_V_I(outer,i));
        w = w - S_V_II(a,i);
        M_I_I(w,S_V_I(inner,i-1));
        }
    M_I_I(w,S_V_I(outer,i));
    erg += m_v_pa(inner,inner);
    erg += m_v_pa(outer,outer);

    spa = callocobject();
    erg += m_gk_spa(outer,inner,spa);
    FREEALL(inner);
    FREEALL(outer);
    erg += conjugate(spa,spa);

    erg += kostka_number(b,spa,c);
    FREEALL(spa);
    }
    ENDR("number_01_matrices");
}

INT t_SCHUR_SCHUR(a,b) OP a,b; { return copy(a,b); }
INT t_ELMSYM_ELMSYM(a,b) OP a,b; { return copy(a,b); }
INT t_POWSYM_POWSYM(a,b) OP a,b; { return copy(a,b); }
INT t_HOMSYM_HOMSYM(a,b) OP a,b; { return copy(a,b); }
INT t_MONOMIAL_MONOMIAL(a,b) OP a,b; { return copy(a,b); }

#define CAST_SF(a,b,scalarf,partf,schurf,homsymf,powsymf,monomialf,elmsymf,t,t2)\
    COP(t,a);\
    COP(t2,b);\
    switch(S_O_K(a)) \
        {\
        case INTEGER:\
        case LONGINT:\
        case FF:\
        case CYCLOTOMIC:\
        case SQ_RADICAL:\
        case POLYNOM:\
        case BRUCH: erg += scalarf(a,b); goto ende;\
        case PARTITION: erg += partf(a,b); goto ende;\
        case SCHUR: erg += schurf(a,b); goto ende;\
        case HOMSYM: erg += homsymf(a,b); goto ende;\
        case MONOMIAL: erg += monomialf(a,b); goto ende;\
        case ELMSYM: erg += elmsymf(a,b); goto ende;\
        case POWSYM: erg += powsymf(a,b); goto ende;\
        default: erg += WTO(t,a); goto ende;\
        }\
ende:

#define CAST_APPLY_SF(a,scalarf,partf,schurf,homsymf,powsymf,monomialf,elmsymf,t)\
CAST_SF(a,a,scalarf,partf,schurf,homsymf,powsymf,monomialf,elmsymf,t,t)

INT cast_apply_schur(a) OP a;
/* AK 080102 */
{
    INT erg = OK;
    CAST_APPLY_SF(a,m_scalar_schur,m_pa_s,t_SCHUR_SCHUR,t_HOMSYM_SCHUR,
                                 t_POWSYM_SCHUR,t_MONOMIAL_SCHUR,
                                 t_ELMSYM_SCHUR, "cast_apply_schur(1)");
    CTO(SCHUR,"cast_apply_schur(e1)",a);
    ENDR("cast_apply_schur");
}
INT cast_apply_elmsym(a) OP a;
/* AK 080102 */
{
    INT erg = OK;
    CAST_APPLY_SF(a,m_scalar_elmsym,m_pa_e,t_SCHUR_ELMSYM,t_HOMSYM_ELMSYM,
                                 t_POWSYM_ELMSYM,t_MONOMIAL_ELMSYM,
                                 t_ELMSYM_ELMSYM, "cast_apply_elmsym(1)");
    CTO(ELMSYM,"cast_apply_elmsym(e1)",a);
    ENDR("cast_apply_elmsym");
}
INT cast_apply_homsym(a) OP a;
/* AK 080102 */
{
    INT erg = OK;
    CAST_APPLY_SF(a,m_scalar_homsym,m_pa_h,t_SCHUR_HOMSYM,t_HOMSYM_HOMSYM,
                                 t_POWSYM_HOMSYM,t_MONOMIAL_HOMSYM,
                                 t_ELMSYM_HOMSYM, "cast_apply_homsym(1)");
    CTO(HOMSYM,"cast_apply_homsym(e1)",a);
    ENDR("cast_apply_homsym");
}
INT cast_apply_powsym(a) OP a;
/* AK 080102 */
{
    INT erg = OK;
    CAST_APPLY_SF(a,m_scalar_powsym,m_pa_ps,t_SCHUR_POWSYM,t_HOMSYM_POWSYM,
                                 t_POWSYM_POWSYM,t_MONOMIAL_POWSYM,
                                 t_ELMSYM_POWSYM, "cast_apply_powsym(1)");
    CTO(POWSYM,"cast_apply_powsym(e1)",a);
    ENDR("cast_apply_powsym");
}
INT cast_apply_monomial(a) OP a;
/* AK 080102 */
{
    INT erg = OK;
    CAST_APPLY_SF(a,m_scalar_monomial,m_pa_mon,t_SCHUR_MONOMIAL,t_HOMSYM_MONOMIAL,
                                 t_POWSYM_MONOMIAL,t_MONOMIAL_MONOMIAL,
                                 t_ELMSYM_MONOMIAL, "cast_apply_monomial(1)");
    CTO(MONOMIAL,"cast_apply_monomial(e1)",a);
    ENDR("cast_apply_monomial");
}

INT cast_schur(a,b) OP a,b;
/* AK 080102 */
{
    INT erg = OK;
    CAST_SF(a,b,m_scalar_schur,m_pa_s,t_SCHUR_SCHUR,t_HOMSYM_SCHUR,
                                 t_POWSYM_SCHUR,t_MONOMIAL_SCHUR,
                                 t_ELMSYM_SCHUR, "cast_schur(1)","cast_schur(2)");
    CTO(SCHUR,"cast_schur(e2)",b);
    ENDR("cast_schur");
}
INT cast_elmsym(a,b) OP a,b;
/* AK 080102 */
{
    INT erg = OK;
    CAST_SF(a,b,m_scalar_elmsym,m_pa_e,t_SCHUR_ELMSYM,t_HOMSYM_ELMSYM,
                                 t_POWSYM_ELMSYM,t_MONOMIAL_ELMSYM,
                                 t_ELMSYM_ELMSYM, "cast_elmsym(1)","cast_elmsym(2)");
    CTO(ELMSYM,"cast_elmsym(e2)",b);
    ENDR("cast_elmsym");
}
INT cast_homsym(a,b) OP a,b;
/* AK 080102 */
{
    INT erg = OK;
    CAST_SF(a,b,m_scalar_homsym,m_pa_h,t_SCHUR_HOMSYM,t_HOMSYM_HOMSYM,
                                 t_POWSYM_HOMSYM,t_MONOMIAL_HOMSYM,
                                 t_ELMSYM_HOMSYM, "cast_homsym(1)","cast_homsym(2)");
    CTO(HOMSYM,"cast_homsym(e2)",b);
    ENDR("cast_homsym");
}
INT cast_powsym(a,b) OP a,b;
/* AK 080102 */
{
    INT erg = OK;
    CAST_SF(a,b,m_scalar_powsym,m_pa_ps,t_SCHUR_POWSYM,t_HOMSYM_POWSYM,
                                 t_POWSYM_POWSYM,t_MONOMIAL_POWSYM,
                                 t_ELMSYM_POWSYM, "cast_powsym(1)", "cast_powsym(2)");
    CTO(POWSYM,"cast_powsym(e2)",b);
    ENDR("cast_powsym");
}
INT cast_monomial(a,b) OP a,b;
/* AK 080102 */
{
    INT erg = OK;
    CAST_SF(a,b,m_scalar_monomial,m_pa_mon,t_SCHUR_MONOMIAL,t_HOMSYM_MONOMIAL,
                                 t_POWSYM_MONOMIAL,t_MONOMIAL_MONOMIAL,
                                 t_ELMSYM_MONOMIAL, "cast_monomial(1)","cast_monomial(2)");
    CTO(MONOMIAL,"cast_monomial(e2)",b);
    ENDR("cast_monomial");
}

INT frobenius_elmsym(a,b) OP a,b;
/* result is n basis of elmsym */
{
    INT erg = OK;
    OP z;
    CTTTO(PARTITION,ELMSYM,HASHTABLE,"frobenius_elmsym(1)",a);
    CTTTO(EMPTY,ELMSYM,HASHTABLE,"frobenius_elmsym(2)",b);
    if (S_O_K(b) == EMPTY) erg += init(ELMSYM,b);
    if (S_O_K(a) == PARTITION) {
        erg += t_HOMSYM_ELMSYM(a,b);
        goto ende;
        }
    else if (S_O_K(a) == HASHTABLE) {
        erg += t_HOMSYM_ELMSYM(a,b);
        goto ende;
        }
    else { /* ELMSYM */
        z = a; while (z != NULL) { C_O_K(z,HOMSYM); z = S_L_N(z); }
        erg += t_HOMSYM_ELMSYM(a,b);
        z = a; while (z != NULL) { C_O_K(z,HOMSYM); z = S_L_N(z); }
        goto ende;
        }
ende:
    ENDR("frobenius_elmsym");
}

INT frobenius_schur(a,b) OP a,b;
/* result is n basis of schur */
{
    INT erg = OK;
    CTTTO(PARTITION,SCHUR,HASHTABLE,"frobenius_schur(1)",a);
    CTTTO(EMPTY,SCHUR,HASHTABLE,"frobenius_schur(2)",b);
    if (S_O_K(b) == EMPTY) erg += init(SCHUR,b);
    if (S_O_K(a) == PARTITION) {
        OP d = CALLOCOBJECT();
        erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
        M_I_I(1,S_MO_K(d));
        erg += conjugate_partition(a,S_MO_S(d));
        INSERT_SCHURMONOM_(d,b);
        goto ende;
        }
    else { /* SCHUR */ /* HASHTABLE */
        erg += conjugate_schur(a,b);
        goto ende;
        }
ende:
    ENDR("frobenius_schur");
}


INT frobenius_powsym(a,b) OP a,b;
/* result is n basis of powsym */
{
    INT erg = OK;
    INT sig,i;
    OP z,d;
    CTTTO(PARTITION,POWSYM,HASHTABLE,"frobenius_powsym(1)",a);
    CTTTO(EMPTY,POWSYM,HASHTABLE,"frobenius_powsym(2)",b);
    if (S_O_K(b) == EMPTY) erg += init(POWSYM,b);
    if (S_O_K(a) == PARTITION) {
        d = CALLOCOBJECT();
        erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
        for (sig=1,i=0;i<S_PA_LI(a);i++)
            if ( (S_PA_II(a,i) % 2) == 0 ) sig *= (-1);
        M_I_I(sig,S_MO_K(d));
        erg += copy_partition(a,S_MO_S(d));
        INSERT_POWSYMMONOM_(d,b);
        goto ende;
        }
    else { /* POWSYM */ /* HASHTABLE */
        FORALL(z,a, {
            d = CALLOCOBJECT();
            erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
            for (sig=1,i=0;i<S_PA_LI(S_MO_S(z));i++)
                if ( (S_PA_II(S_MO_S(z),i) % 2) == 0 ) sig *= (-1);
            if (sig == 1)
                COPY(S_MO_K(z),S_MO_K(d));
            else
                ADDINVERS(S_MO_K(z),S_MO_K(d));
            erg += copy_partition(S_MO_S(z),S_MO_S(d));
            INSERT_POWSYMMONOM_(d,b);
            });
        goto ende;
        }
ende:
    ENDR("frobenius_powsym");
}


INT frobenius_homsym(a,b) OP a,b;
/* result is n basis of homsym */
{
    INT erg = OK;
    OP z;
    CTTTO(PARTITION,HOMSYM,HASHTABLE,"frobenius_homsym(1)",a);
    CTTTO(EMPTY,HOMSYM,HASHTABLE,"frobenius_homsym(2)",b);
    if (S_O_K(b) == EMPTY) erg += init(HOMSYM,b);
    if (S_O_K(a) == PARTITION) {
        erg += t_ELMSYM_HOMSYM(a,b);
        goto ende;
        }
    else if (S_O_K(a) == HASHTABLE) {
        erg += t_ELMSYM_HOMSYM(a,b);
        goto ende;
        }
    else { /* ELMSYM */
        z = a; while (z != NULL) { C_O_K(z,ELMSYM); z = S_L_N(z); }
        erg += t_ELMSYM_HOMSYM(a,b);
        z = a; while (z != NULL) { C_O_K(z,ELMSYM); z = S_L_N(z); }
        goto ende;
        }
ende:
    ENDR("frobenius_homsym");
}

INT frobenius_monomial(a,b) OP a,b;
/* result is n basis of monomial */
{
    INT erg = OK;
    OP z,y;
    CTTTO(PARTITION,MONOMIAL,HASHTABLE,"frobenius_monomial(1)",a);
    CTTTO(EMPTY,MONOMIAL,HASHTABLE,"frobenius_monomial(2)",b);

    /* via powsym */
    z = CALLOCOBJECT();
    t_MONOMIAL_POWSYM(a,z);
    y = CALLOCOBJECT();
    erg += frobenius_powsym(z,y);
    FREEALL(z);
    erg += t_POWSYM_MONOMIAL(y,b);
    FREEALL(y);
    ENDR("frobenius_monomial");
}






INT scalarproduct_schur_schur(a,b,c) OP a,b,c;
/* <S_I,S_J> = delta_I,J */
{
    OP za,zb;
    OP d;
    INT res;
    INT erg = OK;
    CTTO(PARTITION,SCHUR,"scalarproduct_schur_schur(1)",a);
    CTTO(PARTITION,SCHUR,"scalarproduct_schur_schur(2)",b);
    CTO(EMPTY,"scalarproduct_schur_schur(3)",c);
    if (S_O_K(a) == PARTITION) {
        d = CALLOCOBJECT();
        erg += m_pa_s(a,d);
        erg += scalarproduct_schur_schur(d,b,c);
        FREEALL(d);
        goto ende;
        }
    if (S_O_K(b) == PARTITION) {
        d = CALLOCOBJECT();
        erg += m_pa_s(b,d);
        erg += scalarproduct_schur_schur(a,d,c);
        FREEALL(d);
        goto ende;
        }


    d = CALLOCOBJECT();
    za = a;
    zb = b;
    M_I_I((INT)0,c);
    do {
        if (za == NULL) goto preende;
        if (zb == NULL) goto preende;
        res = comp(S_S_S(za),S_S_S(zb));
        if (res == (INT)0)
            {
            FREESELF(d);
            MULT(S_S_K(za),S_S_K(zb),d);
            ADD_APPLY(d,c);
            za = S_S_N(za);
            zb = S_S_N(zb);
            }
        else if (res < (INT)0)
            {
            za = S_S_N(za);
            }
        else
            {
            zb = S_S_N(zb);
            }
    } while(1);
preende:
    FREEALL(d);
ende:
    ENDR("scalarproduct_schur_schur");
}

INT scalarproduct_powsym_powsym(a,b,c) OP a,b,c;
/* <P_I,P_J> = ordcen() delta_I,J */
{
    OP za,zb;
    OP d,e;
    INT res;
    INT erg = OK;
    CTTO(PARTITION,POWSYM,"scalarproduct_powsym_powsym(1)",a);
    CTTO(PARTITION,POWSYM,"scalarproduct_powsym_powsym(2)",b);
    CTO(EMPTY,"scalarproduct_powsym_powsym(3)",c);
    if (S_O_K(a) == PARTITION) {
        d = CALLOCOBJECT();
        erg += m_pa_ps(a,d);
        erg += scalarproduct_powsym_powsym(d,b,c);
        FREEALL(d);
        goto ende;
        }
    if (S_O_K(b) == PARTITION) {
        d = CALLOCOBJECT();
        erg += m_pa_ps(b,d);
        erg += scalarproduct_powsym_powsym(a,d,c);
        FREEALL(d);
        goto ende;
        }

    d = CALLOCOBJECT();
    e = CALLOCOBJECT();
    za = a;
    zb = b;
    M_I_I((INT)0,c);
    do {
        if (za == NULL) goto preende;
        if (zb == NULL) goto preende;
        res = comp(S_S_S(za),S_S_S(zb));
        if (res == (INT)0)
            {
            FREESELF(d);
            MULT(S_S_K(za),S_S_K(zb),d);
            ordcen(S_S_S(za),e);
            MULT_APPLY(e,d);
            ADD_APPLY(d,c);
            za = S_S_N(za);
            zb = S_S_N(zb);
            }
        else if (res < (INT)0)
            {
            za = S_S_N(za);
            }
        else
            {
            zb = S_S_N(zb);
            }
    } while(1);
preende:
    FREEALL(d);
    FREEALL(e);
ende:
    ENDR("scalarproduct_powsym_powsym");
}



INT scalarproduct_homsym_monomial(a,b,c) OP a,b,c;
/* <H_I,M_J> = delta_I,J */
{
    OP za,zb;
    OP d;
    INT res;
    INT erg = OK;
    CTTO(PARTITION,HOMSYM,"scalarproduct_homsym_monomial(1)",a);
    CTTO(PARTITION,MONOMIAL,"scalarproduct_homsym_monomial(2)",b);
    CTO(EMPTY,"scalarproduct_homsym_monomial(3)",c);
    if (S_O_K(a) == PARTITION) {
        d = CALLOCOBJECT();
        erg += m_pa_h(a,d);
        erg += scalarproduct_homsym_monomial(d,b,c);
        FREEALL(d);
        goto ende;
        }
    if (S_O_K(b) == PARTITION) {
        d = CALLOCOBJECT();
        erg += m_pa_mon(b,d);
        erg += scalarproduct_homsym_monomial(a,d,c);
        FREEALL(d);
        goto ende;
        }


    d = CALLOCOBJECT();
    za = a;
    zb = b;
    M_I_I((INT)0,c);
    do {
        if (za == NULL) goto preende;
        if (zb == NULL) goto preende;
        res = comp(S_S_S(za),S_S_S(zb));
        if (res == (INT)0)
            {
            FREESELF(d);
            MULT(S_S_K(za),S_S_K(zb),d);
            ADD_APPLY(d,c);
            za = S_S_N(za);
            zb = S_S_N(zb);
            }
        else if (res < (INT)0)
            {
            za = S_S_N(za);
            }
        else
            {
            zb = S_S_N(zb);
            }
    } while(1);
preende:
    FREEALL(d);
ende:
    ENDR("scalarproduct_homsym_monomial");
}

INT scalarproduct_schur(a,b,c) OP a,b,c;
{
    INT erg = OK;
    CTO(SCHUR,"scalarproduct_schur(1)",a);
    CTO(EMPTY,"scalarproduct_schur(3)",c);
    if (S_O_K(b) == SCHUR)
        {
        erg += scalarproduct_schur_schur(a,b,c);
        }
    else if (S_O_K(b) == HOMSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_HOMSYM_SCHUR(b,d);
        erg += scalarproduct_schur_schur(a,d,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == ELMSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_ELMSYM_SCHUR(b,d);
        erg += scalarproduct_schur_schur(a,d,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == POWSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_POWSYM_SCHUR(b,d);
        erg += scalarproduct_schur_schur(a,d,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == MONOMIAL)
        {
        OP d;
        d = CALLOCOBJECT();
        t_MONOMIAL_SCHUR(b,d);
        erg += scalarproduct_schur_schur(a,d,c);
        FREEALL(d);
        }
    else
        WTO("scalarproduct_schur(2)",b);
    ENDR("scalarproduct_schur");
}

INT scalarproduct_powsym(a,b,c) OP a,b,c;
{
    INT erg = OK;
    CTO(POWSYM,"scalarproduct_powsym(1)",a);
    CTO(EMPTY,"scalarproduct_powsym(3)",c);
    if (S_O_K(b) == POWSYM)
        {
        erg += scalarproduct_powsym_powsym(a,b,c);
        goto sppende;
        }
    else if (S_O_K(b) == HOMSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_HOMSYM_POWSYM(b,d);
        erg += scalarproduct_powsym_powsym(a,d,c);
        FREEALL(d);
        goto sppende;
        }
    else if (S_O_K(b) == ELMSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_ELMSYM_POWSYM(b,d);
        erg += scalarproduct_powsym_powsym(a,d,c);
        FREEALL(d);
        goto sppende;
        }
    else if (S_O_K(b) == SCHUR)
        {
        OP d;
        d = CALLOCOBJECT();
        t_SCHUR_POWSYM(b,d);
        erg += scalarproduct_powsym_powsym(a,d,c);
        FREEALL(d);
        goto sppende;
        }
    else if (S_O_K(b) == MONOMIAL)
        {
        OP d;
        d = CALLOCOBJECT();
        t_MONOMIAL_POWSYM(b,d);
        erg += scalarproduct_powsym_powsym(a,d,c);
        FREEALL(d);
        goto sppende;
        }
    else
        WTO("scalarproduct_powsym(2)",b);
sppende:
    ENDR("scalarproduct_powsym");
}


INT scalarproduct_homsym(a,b,c) OP a,b,c;
{
    INT erg = OK;
    CTO(HOMSYM,"scalarproduct_homsym(1)",a);
    CTO(EMPTY,"scalarproduct_homsym(3)",c);
    if (S_O_K(b) == MONOMIAL)
        {
        erg += scalarproduct_homsym_monomial(a,b,c);
        }
    else if (S_O_K(b) == HOMSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_HOMSYM_MONOMIAL(b,d);
        erg += scalarproduct_homsym_monomial(a,d,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == ELMSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_ELMSYM_MONOMIAL(b,d);
        erg += scalarproduct_homsym_monomial(a,d,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == POWSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_POWSYM_MONOMIAL(b,d);
        erg += scalarproduct_homsym_monomial(a,d,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == SCHUR)
        {
        OP d;
        d = CALLOCOBJECT();
        t_SCHUR_MONOMIAL(b,d);
        erg += scalarproduct_homsym_monomial(a,d,c);
        FREEALL(d);
        }
    else
        WTO("scalarproduct_homsym(2)",b);
    ENDR("scalarproduct_homsym");
}

INT scalarproduct_monomial(a,b,c) OP a,b,c;
{
    INT erg = OK;
    CTO(MONOMIAL,"scalarproduct_monomial(1)",a);
    CTO(EMPTY,"scalarproduct_monomial(3)",c);
    if (S_O_K(b) == HOMSYM)
        {
        erg += scalarproduct_homsym_monomial(b,a,c);
        }
    else if (S_O_K(b) == MONOMIAL)
        {
        OP d;
        d = CALLOCOBJECT();
        t_MONOMIAL_HOMSYM(b,d);
        erg += scalarproduct_homsym_monomial(d,a,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == ELMSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_ELMSYM_HOMSYM(b,d);
        erg += scalarproduct_homsym_monomial(d,a,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == POWSYM)
        {
        OP d;
        d = CALLOCOBJECT();
        t_POWSYM_HOMSYM(b,d);
        erg += scalarproduct_homsym_monomial(d,a,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == SCHUR)
        {
        OP d;
        d = CALLOCOBJECT();
        t_SCHUR_HOMSYM(b,d);
        erg += scalarproduct_homsym_monomial(d,a,c);
        FREEALL(d);
        }
    else
        WTO("scalarproduct_monomial(2)",b);
    ENDR("scalarproduct_monomial");
}

INT scalarproduct_elmsym(a,b,c) OP a,b,c;
{
    INT erg = OK;
    OP e;
    CTO(ELMSYM,"scalarproduct_elmsym(1)",a);
    CTO(EMPTY,"scalarproduct_elmsym(3)",c);
    e = CALLOCOBJECT();
    if (S_O_K(b) == HOMSYM)
        {
        OP d;
        t_ELMSYM_SCHUR(a,e);
        d = CALLOCOBJECT();
        t_HOMSYM_SCHUR(b,d);
        erg += scalarproduct_schur_schur(d,e,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == MONOMIAL)
        {
        OP d;
        t_ELMSYM_SCHUR(a,e);
        d = CALLOCOBJECT();
        t_MONOMIAL_SCHUR(b,d);
        erg += scalarproduct_schur_schur(d,e,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == ELMSYM)
        {
        OP d;
        t_ELMSYM_SCHUR(a,e);
        d = CALLOCOBJECT();
        t_ELMSYM_SCHUR(b,d);
        erg += scalarproduct_schur_schur(d,e,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == POWSYM)
        {
        OP d;
        t_ELMSYM_SCHUR(a,e);
        d = CALLOCOBJECT();
        t_POWSYM_SCHUR(b,d);
        erg += scalarproduct_schur_schur(d,e,c);
        FREEALL(d);
        }
    else if (S_O_K(b) == SCHUR)
        {
        t_ELMSYM_SCHUR(a,e);
        erg += scalarproduct_schur_schur(b,e,c);
        }
    else
        WTO("scalarproduct_elmsym(2)",b);
    FREEALL(e);
    ENDR("scalarproduct_elmsym");
}




static INT co_1611(a,b) OP a,b;
/* transform polynom to vector of monoms */
{
    OP z;
    INT i=0,j;
    z = a;
    while (z != NULL)
        {
        i += S_PO_KI(z);
        z = S_PO_N(z);
        }

    m_il_v(i,b);
    z = a;
    for (i=0;i<S_V_LI(b);)
        {
        if (z == NULL) error("");
        for (j=0;j<S_PO_KI(z);j++)
            {
            m_skn_po(S_PO_S(z),cons_eins,NULL,S_V_I(b,i));
            i++;
            }
        z = S_PO_N(z);
        }
    return OK;
}

INT plet_schur_schur_pol(parta,partb, len, res)
        OP parta,partb, len, res;
/* plethysm of schur polynomials S_b(S_a(A_len)) */
/* res is polynom */
/* parameters may be equal */
{
    OP c ,e;
    INT erg = OK;
    CTO(PARTITION,"plet_schur_schur_pol(1)",parta);
    CTO(PARTITION,"plet_schur_schur_pol(2)",partb);
    CTO(INTEGER,"plet_schur_schur_pol(3)",len);
    c = CALLOCOBJECT();
    e = CALLOCOBJECT();
    erg += compute_schur_with_alphabet(parta,len,c);
    erg += co_1611(c,e);
    erg += compute_schur_with_alphabet(partb,S_V_L(e),c);
    erg += eval_polynom(c,e,res);
    erg += freeall(c);
    erg += freeall(e);
    CTO(POLYNOM,"plet_schur_schur_pol(4-end)",res);
    ENDR("plet_schur_schur_pol");
}



static OP ll=NULL;
static INT co111(a) OP a;
{
        return S_PA_LI(S_MO_S(a)) <= S_I_I(ll);
}
static INT co222(a) OP a;
{
        return S_PA_LI(S_MO_S(a)) == S_I_I(ll);
}

INT put_length_limit_schur(i,a) OP a; INT i;
/* i ist the max length for part */
{
    INT erg = OK;
    CTTO(SCHUR,HASHTABLE,"put_length_limit_schur(2)",a);
    ll = CALLOCOBJECT(); M_I_I(i,ll);
    if (S_O_K(a) == SCHUR)
        erg += filter_apply_list(a,co111);
    else /*HASHTABLE */
        erg += filter_apply_hashtable(a,co111);
    FREEALL(ll); ll = NULL;
    ENDR("put_length_limit_schur");
}
INT put_exactlength_limit_schur(i,a) OP a; INT i;
/* i ist the exact length for part */
{
    INT erg = OK;
    CTTO(SCHUR,HASHTABLE,"put_length_limit_schur(2)",a);
    ll = CALLOCOBJECT(); M_I_I(i,ll);
    if (S_O_K(a) == SCHUR)
        erg += filter_apply_list(a,co222);
    else /*HASHTABLE */
        erg += filter_apply_hashtable(a,co222);
    FREEALL(ll); ll = NULL;
    ENDR("put_exactlength_limit_schur");
}


static INT m_int_Qelm(a,b) INT a; OP b;
{
    INT erg = OK;
    OP c,d,e;

    COP("m_int_Qelm(2)",b);
    c = callocobject();
    d = callocobject();
    if (a == (INT)0)
        {
        m_i_i(1,b);
        goto ee;
        }
    m_i_i(a,c);
    init(MONOMIAL,b);
    first_partition(c,d);
    do {
    e = callocobject();
    erg += m_skn_mon(d,callocobject(),NULL,e);
    hoch(cons_zwei,S_PA_L(d),S_S_K(e));
    insert(e,b,NULL,NULL);
    } while (next(d,d));
ee:
    freeall(c);
    freeall(d);
    ENDR("m_int_Qelm");
}
static INT m_int_int_Qelm(a,b,d) INT a,b; OP d;
{
    OP c;
    INT erg = OK;
    COP("m_int_int_Qelm(3)",d);
    c = callocobject();
    if (a>b) error("::");
    erg += b_ks_pa(VECTOR,callocobject(),c);
    erg += m_il_v(2L,S_PA_S(c));
    erg += m_i_i(a,S_PA_I(c,0));
    erg += m_i_i(b,S_PA_I(c,1));
    erg += m_part_Qschur(c,d);
    erg += freeall(c);
    ENDR("m_int_int_qelm");

}

INT m_part_Qschur(a,b) OP a,b;
/* AK 291295 */
/* computes Q schur polynomial as monomial sym */
{
    INT i,j;
    OP c,d,e;
    INT erg = OK;

    CTO(PARTITION,"m_part_Qschur",a);
    if (S_PA_LI(a) == 1)
        {
        erg += m_int_Qelm(S_PA_II(a,0),b);
        }
    else if (S_PA_LI(a) == 2)
        {
        c = callocobject();
        erg += m_int_Qelm(S_PA_II(a,0),c);
        d = callocobject();
        erg += m_int_Qelm(S_PA_II(a,1),d);
        erg += mult(c,d,b);
        e = callocobject();
        for (i=1;i<=S_PA_II(a,0);i++)
            {
            erg += m_int_Qelm(S_PA_II(a,0)-i,c);
            erg += m_int_Qelm(S_PA_II(a,1)+i,d);
            erg += mult(c,d,e);
            erg += mult_apply(cons_zwei,e);
            if (i%2 == 1)
                erg += mult_apply(cons_negeins,e);
            erg += add_apply(e,b);
            }
        erg += freeall(c);
        erg += freeall(d);
        erg += freeall(e);
        }
    else if  (S_PA_LI(a) %2 == 0)
        {
        c = callocobject();
        erg += m_ilih_m(S_PA_LI(a), S_PA_LI(a),c);
        for (i=0;i<S_M_HI(c);i++)
            for (j=i+1;j<S_M_LI(c);j++)
                {
                m_int_int_Qelm(S_PA_II(a,S_PA_LI(a)-1-j),
                    S_PA_II(a,S_PA_LI(a)-1-i),S_M_IJ(c,i,j));
                }
        pfaffian_matrix(c,b);
        erg += freeall(c);
        }
    else    {
        d = callocobject();
        b_ks_pa(VECTOR,callocobject(),d);
        m_il_nv(S_PA_LI(a)+1,S_PA_S(d));
        for (i=0;i<S_PA_LI(a);i++)
            M_I_I(S_PA_II(a,i),S_PA_I(d,i+1));

        c = callocobject();
        erg += m_ilih_m(S_PA_LI(d), S_PA_LI(d),c);
        for (i=0;i<S_M_HI(c);i++)
            for (j=i+1;j<S_M_LI(c);j++)
                {
                m_int_int_Qelm(S_PA_II(d,S_PA_LI(d)-1-j),
                    S_PA_II(d,S_PA_LI(d)-1-i),S_M_IJ(c,i,j));
                }
        pfaffian_matrix(c,b);
        erg += freeall(c);
        freeall(d);
        }
    ENDR("m_part_Qschur");
}


#ifdef UNDEF
static dec_step();
static dec_step_2();
static dec_step_co();

/* fehlerhafte werte 180901 */
static OP dd,ee,ff;
INT old_number_01_matrices(a,b,c) OP a,b,c;
/* zeilen und spaltensumme in vector schreibweise als VECTOR */
/* a = zeilen, b = spalten */
{
    INT erg = OK;
    OP d,e;
    if (S_O_K(a) == PARTITION) a = S_PA_S(a);
    if (S_O_K(b) == PARTITION) b = S_PA_S(b);
    CTTO(VECTOR,INTEGERVECTOR,"old_number_01_matrices",a);
    CTTO(VECTOR,INTEGERVECTOR,"old_number_01_matrices",b);
    CE3(a,b,c,old_number_01_matrices);

    dd = callocobject();
    ee = callocobject();
    ff = callocobject();
    d = callocobject();
    e = callocobject();
    erg += m_ilih_nm(S_V_LI(a)+1,S_V_LI(b)+1,ff);
    erg += m_ks_pa(VECTOR,b,ee);
    erg += t_VECTOR_EXPONENT(ee,e);
    erg += freeself(ee);
    erg += dec_step(e,b,d);
    if (S_L_S(d) != NULL)
        erg += copy(S_S_K(d),c);
    else
        erg += m_i_i(0,c);
    erg += freeall(d);
    erg += freeall(dd);
    erg += freeall(ee);
    erg += freeall(ff);
    erg += freeall(e);
    ENDR("number_01_matrices");
}


static my_binom(a,b,c) OP a,b,c;
{
    OP z;
    INT ai=S_I_I(a), bi=S_I_I(b);
    INT j=1;
    z = s_m_ij(ff,ai,bi);
    if (not nullp(z))
        copy(z,c);
    else {
    m_i_i(1,c);
    m_i_i(1,ee);
    while (ai > bi)
        {
        M_I_I(ai,dd); MULT_APPLY_INTEGER(dd,c);
        ai--;
        M_I_I(j,dd);  MULT_APPLY_INTEGER(dd,ee);
        j++;
        }
    ganzdiv_apply(ee,c);
        copy(c,z);
        }
}

static dec_step(a,b,c) OP a,b,c;
{
    INT erg = OK;
    CE3(a,b,c,dec_step);
    if (S_O_K(b) == VECTOR || S_O_K(b) == INTEGERVECTOR)
        {
        INT i;
        dec_step_2(a,S_V_I(b,0),c);
        for (i=1;i<S_V_LI(b);i++)
            {
            dec_step_2(c,S_V_I(b,i),c);
            }
        }
    else
        dec_step_2(a,b,c);
    ENDR("dec_step");
}

static INT my_comp(a,b) OP a,b; {
    return (INT)memcmp(
        (char *) S_V_S(S_PA_S(S_MO_S(a))),
        (char *) S_V_S(S_PA_S(S_MO_S(b))),
        sizeof(struct object) * S_PA_LI(S_MO_S(a)));
    }

static dec_step_2(a,b,c) OP a,b,c;
{
    INT erg = OK;
    CE3(a,b,c,dec_step_2);
    if ( (S_O_K(a) == PARTITION) && (S_O_K(b) == INTEGER) )
        {
        dec_step_co(a,b,c,cons_eins);
        }
    else if ( (S_O_K(a) == SCHUR) && (S_O_K(b) == INTEGER) )
        {
        OP z;
        z = a;
        init(SCHUR,c);
        if (S_L_S(z) != NULL)
            do {
                OP d = callocobject();
                dec_step_co(S_S_S(z),b,d, S_S_K(z) );
                insert(d,c,add_koeff,my_comp);
                z = S_S_N(z);
            } while (z != NULL);
        }
    ENDR("dec_step");
}


static dec_step_co(a,b,c,faktor) OP a,b,c; OP faktor;
{
    OP e,f,g,d = callocobject();
    INT i,j;
    OP anzahl,indexvec;
    init(BINTREE,c);
    g = callocobject();
    for (i=0,j=0;i<S_PA_LI(a);i++)
        if (S_PA_II(a,i)>0) j++;
    /* j ist die anzahl der teile != 0 */
    anzahl = callocobject();
    indexvec = callocobject();
    M_I_I(j,anzahl);
    b_l_v(anzahl,indexvec);
    for (i=0,j=0;i<S_PA_LI(a);i++)
        if (S_PA_II(a,i)>0)
            { M_I_I(i,S_V_I(indexvec,j)); j++; }
    /* der j-te eintrag in indexvec
       ist der index des j-ten eintrag != 0 in a */

    first_composition(b,S_V_L(indexvec),d);
    do {
        for (i=0;i<S_V_LI(indexvec);i++)
            if ( S_PA_II(a,S_V_II(indexvec,i) ) < S_V_II(d,i) )
                goto next_step;
        e = callocobject();
        init(MONOM,e);
        copy_partition(a,S_MO_S(e));
        for (i=0;i<S_V_LI(indexvec);i++)
            M_I_I(S_PA_II(S_MO_S(e),S_V_II(indexvec,i))
                  -S_V_II(d,i),
                S_PA_I(S_MO_S(e),S_V_II(indexvec,i))
                );


        copy(faktor, S_MO_K(e));
        for (i=0;i<S_V_LI(d);i++)
            if (S_V_II(d,i)  != 0)
            {
            if (S_V_II(indexvec,i) > 0)
                M_I_I(S_PA_II(S_MO_S(e),S_V_II(indexvec,i)-1)+S_V_II(d,i),
                    S_PA_I(S_MO_S(e),S_V_II(indexvec,i)-1));
            if (S_PA_II(a,S_V_II(indexvec,i)) != S_V_II(d,i) )
                                {
                my_binom(S_PA_I(a,S_V_II(indexvec,i)), S_V_I(d,i), g);
                                mult_apply(g,S_MO_K(e));
                }
            }
        insert(e,c,add_koeff,my_comp);
next_step:              ;
    } while(next_apply(d));
    t_BINTREE_SCHUR(c,c); println(c);
    freeall(d);
    freeall(indexvec);
    freeall(g);
}

#endif


#ifdef UNDEF
INT class_mult_schurmonom(OP a, OP b, OP c)
{
    INT erg = OK;
    OP e,d;
    CE3(a,b,c,class_mult_schurmonom);
    CTO(MONOM,"class_mult_schurmonom(1)",a);
    CTO(MONOM,"class_mult_schurmonom",b);
    CTO(PARTITION,"class_mult_schurmonom-selfpart",S_MO_S(a));
    CTO(PARTITION,"class_mult_schurmonom-selfpart",S_MO_S(b));

    e = callocobject();
    d = callocobject();
    erg += weight(S_MO_S(a),e);
    erg += weight(S_MO_S(b),d);
    if (neq(e,d))
        {
        erg += error("class_mult_schurmonom:different weights of partitions");
        goto ee;
        }
    C2R(S_MO_S(a),S_MO_S(b),"class_mult_part",c);
    erg += init(SCHUR,c);
    erg += first_partition(e,d);
    do {
        erg += c_ijk_sn(S_MO_S(a),S_MO_S(b),d,e);
        if (not NULLP(e))
            {
            OP f,g;
            f = callocobject();
            g = callocobject();
            erg += copy_partition(d,g);
            erg += b_skn_s(g,callocobject(),NULL,f);
            erg += mult(S_MO_K(a),S_MO_K(b),S_S_K(f));
            erg += mult_apply(e,S_S_K(f));
            insert(f,c,NULL,NULL);
            }

    } while(next(d,d));
    S2R(S_MO_S(a),S_MO_S(b),"class_mult_part",c);
ee:
    erg += freeall(e);
    erg += freeall(d);
    ENDR("class_mult_schurmonom");
}
#endif

INT class_mult_schur(OP a, OP b, OP c)
{
#ifdef UNDEF
    INT erg = OK;
    OP z1,z2;
    CTO(SCHUR,"class_mult_schur",a);
    CTO(SCHUR,"class_mult_schur",b);
    CE3(a,b,c,class_mult_schur);
    erg += init(SCHUR,c);
    z1 = a;
    while (z1 != NULL)
        {
        z2 = b;
        while (z2 != NULL)
            {
            OP e;
            e = callocobject();
            if (le(S_S_S(z1), S_S_S(z2)))
                erg += class_mult_schurmonom(S_L_S(z1),S_L_S(z2),e);
            else
                erg += class_mult_schurmonom(S_L_S(z2),S_L_S(z1),e);
            insert(e,c,NULL,NULL);
            z2 = S_S_N(z2);
            };
        z1 = S_S_N(z1);
        }
    ENDR("class_mult_schur");
#endif
    return class_mult(a,b,c);
}




INT init_elmsym(a) OP a;
{
    INT erg = OK;
    CTO(EMPTY,"init_elmsym",a);
    erg += b_sn_e(NULL,NULL,a);
    ENDR("init_elmsym");
}
INT init_homsym(a) OP a;
{
    INT erg = OK;
    CTO(EMPTY,"init_homsym",a);
    erg += b_sn_h(NULL,NULL,a);
    ENDR("init_homsym");
}
INT init_powsym(a) OP a;
{
    INT erg = OK;
    CTO(EMPTY,"init_powsym",a);
    erg += b_sn_ps(NULL,NULL,a);
    ENDR("init_powsym");
}
INT init_schur(a) OP a;
{
    INT erg = OK;
    CTO(EMPTY,"init_schur",a);
    erg += b_sn_s(NULL,NULL,a);
    ENDR("init_schur");
}

INT init_monomial(a) OP a;
{
    INT erg = OK;
    CTO(EMPTY,"init_monomial",a);
    erg += b_sn_mon(NULL,NULL,a);
    ENDR("init_monomial");
}

INT conjugate_schur(a,b) OP a,b;
/* AK 111001 */
{
    INT erg=OK,t=0;
    OP z;
    CTTO(HASHTABLE,SCHUR,"conjugate_schur(1)",a);
    CTTTO(EMPTY,SCHUR,HASHTABLE,"conjugate_schur(2)",b);

    if (S_O_K(b) == EMPTY) {
        init_hashtable(b); t=1; }

    FORALL(z,a,{
        OP d = CALLOCOBJECT();
        erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
        COPY(S_MO_K(z),S_MO_K(d));
        erg += conjugate_partition(S_MO_S(z),S_MO_S(d));

        if (S_O_K(b) == SCHUR)
            insert_list(d,b,NULL,comp_monomschur);
        else
            insert_hashtable(d,b,NULL,NULL,hash_monompartition);
        } );

    if (t==1) t_HASHTABLE_SCHUR(b,b);

    ENDR("conjugate_schur");
}

INT conjugate_elmsym(a,b) OP a,b;
/* AK 111001 */
{
    INT erg=OK,t=0;
    OP z;
    CTTO(HASHTABLE,ELMSYM,"conjugate_elmsym(1)",a);
    CTTTO(EMPTY,ELMSYM,HASHTABLE,"conjugate_elmsym(2)",b);

    if (S_O_K(b) == EMPTY)
        {
        init_hashtable(b);
        t=1;
        }

    FORALL(z,a,{
        OP d = CALLOCOBJECT();
        erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
        COPY(S_MO_K(z),S_MO_K(d));
        erg += conjugate_partition(S_MO_S(z),S_MO_S(d));

        if (S_O_K(b) == ELMSYM)
            insert_list(d,b,NULL,comp_monomelmsym);
        else
            insert_hashtable(d,b,NULL,NULL,hash_monompartition);
        } );

    if (t==1) t_HASHTABLE_ELMSYM(b,b);

    ENDR("conjugate_elmsym");
}

INT conjugate_homsym(a,b) OP a,b;
/* AK 111001 */
{
    INT erg=OK,t=0;
    OP z;
    CTTO(HASHTABLE,HOMSYM,"conjugate_homsym(1)",a);
    CTTTO(EMPTY,HOMSYM,HASHTABLE,"conjugate_homsym(2)",b);

    if (S_O_K(b) == EMPTY)
        {
        init_hashtable(b);
        t=1;
        }

    FORALL(z,a,{
        OP d = CALLOCOBJECT();
        erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
        COPY(S_MO_K(z),S_MO_K(d));
        erg += conjugate_partition(S_MO_S(z),S_MO_S(d));

        if (S_O_K(b) == HOMSYM)
            insert_list(d,b,NULL,comp_monomhomsym);
        else
            insert_hashtable(d,b,NULL,NULL,hash_monompartition);
        } );

    if (t==1) t_HASHTABLE_HOMSYM(b,b);

    ENDR("conjugate_homsym");
}

INT conjugate_powsym(a,b) OP a,b;
/* AK 111001 */
{
    INT erg=OK,t=0;
    OP z;
    CTTO(HASHTABLE,POWSYM,"conjugate_powsym(1)",a);
    CTTTO(EMPTY,POWSYM,HASHTABLE,"conjugate_powsym(2)",b);

    if (S_O_K(b) == EMPTY)
        {
        init_hashtable(b);
        t=1;
        }

    FORALL(z,a,{
        OP d = CALLOCOBJECT();
        erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
        COPY(S_MO_K(z),S_MO_K(d));
        erg += conjugate_partition(S_MO_S(z),S_MO_S(d));

        if (S_O_K(b) == POWSYM)
            insert_list(d,b,NULL,comp_monompowsym);
        else
            insert_hashtable(d,b,NULL,NULL,hash_monompartition);
        } );

    if (t==1) t_HASHTABLE_POWSYM(b,b);

    ENDR("conjugate_powsym");
}

INT conjugate_monomial(a,b) OP a,b;
/* AK 111001 */
{
    INT erg=OK,t=0;
    OP z;
    CTTO(HASHTABLE,MONOMIAL,"conjugate_monomial(1)",a);
    CTTTO(EMPTY,MONOMIAL,HASHTABLE,"conjugate_monomial(2)",b);

    if (S_O_K(b) == EMPTY)
        {
        init_hashtable(b);
        t=1;
        }

    FORALL(z,a,{
        OP d = CALLOCOBJECT();
        erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),d);
        COPY(S_MO_K(z),S_MO_K(d));
        erg += conjugate_partition(S_MO_S(z),S_MO_S(d));

        if (S_O_K(b) == MONOMIAL)
            insert_list(d,b,NULL,comp_monommonomial);
        else
            insert_hashtable(d,b,NULL,NULL,hash_monompartition);
        } );

    if (t==1) t_HASHTABLE_MONOMIAL(b,b);

    ENDR("conjugate_monomial");
}




OP find_schur(a,b) OP a,b;
/* AK 161001 */
/* return OP pointer to SCHUR element
   with partition eq b, b is monom or partition */
{
    INT erg = OK;
    OP z,p;
    CTO(SCHUR,"find_schur(1)",a);
    CTTO(PARTITION,MONOM,"find_schur(2)",b);
    if (S_O_K(b) == MONOM)
        {
        CTO(PARTITION,"find_schur(2b)",S_MO_S(b));
        p = S_MO_S(b);
        }
    else p = b;

    z = a;
    while (z != NULL)
        {
        if (EQ(p,S_S_S(z))) return S_L_S(z);
        z = S_S_N(z);
        }
    return NULL;
    ENDO("find_schur");
}

OP find_monomial(a,b) OP a,b;
/* AK 161001 */
/* return OP pointer to MONOMIAL element
   with partition eq b, b is monom or partition */
{
    INT erg = OK;
    OP z,p;
    CTO(MONOMIAL,"find_monomial(1)",a);
    CTTO(PARTITION,MONOM,"find_monomial(2)",b);
    if (S_O_K(b) == MONOM)
        {
        CTO(PARTITION,"find_monomial(2b)",S_MO_S(b));
        p = S_MO_S(b);
        }
    else p = b;

    z = a;
    while (z != NULL)
        {
        if (EQ(p,S_S_S(z))) return S_L_S(z);
        z = S_S_N(z);
        }
    return NULL;
    ENDO("find_monomial");
}





#define FINDMAX_SF(a,cf)\
    {\
    OP z,res=NULL;\
    if (cf == NULL) cf = comp;\
    FORALL(z,a, {\
        if (res ==NULL) res = z;\
        else if ( (*cf)(S_MO_S(z),S_MO_S(res)) > 0 ) res = z;\
        } );\
    return res;\
    }

OP findmax_schur(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns maximum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(SCHUR,HASHTABLE,"findmax_schur(1)",a);
    FINDMAX_SF(a,cf);
    ENDO("findmax_schur");
}


OP findmax_monomial(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns maximum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,MONOMIAL,"find_monomial",a);
    FINDMAX_SF(a,cf);
    ENDO("findmax_monomial");
}

OP findmax_powsym(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns maximum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,POWSYM,"find_powsym",a);
    FINDMAX_SF(a,cf);
    ENDO("findmax_powsym");
}

OP findmax_elmsym(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns maximum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,ELMSYM,"find_elmsym",a);
    FINDMAX_SF(a,cf);
    ENDO("findmax_elmsym");
}

OP findmax_homsym(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns maximum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,HOMSYM,"find_homsym",a);
    FINDMAX_SF(a,cf);
    ENDO("findmax_homsym");
}

#define FINDMIN_SF(a,cf)\
 {\
    OP z,res=NULL;\
    if (cf == NULL) cf = comp;\
    FORALL(z,a,\
        {\
        if (res ==NULL) res = z;\
        else if ( (*cf)(S_MO_S(z),S_MO_S(res)) < 0 ) res = z;\
        });\
    return res;\
}

OP findmin_monomial(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns minimum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,MONOMIAL,"findmin_monomial(1)",a);
    FINDMIN_SF(a,cf);
    ENDO("findmin_monomial");
}

OP findmin_schur(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns minimum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,SCHUR,"findmin_schur(1)",a);
    FINDMIN_SF(a,cf);
    ENDO("findmin_schur");
}

OP findmin_elmsym(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns minimum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,ELMSYM,"findmin_elmsym(1)",a);
    FINDMIN_SF(a,cf);
    ENDO("findmin_elmsym");
}

OP findmin_homsym(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns minimum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,HOMSYM,"findmin_homsym(1)",a);
    FINDMIN_SF(a,cf);
    ENDO("findmin_homsym");
}

OP findmin_powsym(a,cf) OP a; INT (*cf)();
/* AK 161001 */
/* returns minimum according to comp function */
/* comp function operates on PARTITIONS */
/* if cf == NULL comp is used */
{
    INT erg = OK;
    CTTO(HASHTABLE,POWSYM,"findmin_powsym(1)",a);
    FINDMIN_SF(a,cf);
    ENDO("findmin_powsym");
}

INT m_forall_monomials_in_a(a,b,c,f,partf) OP a,b,c,f; INT (*partf)();
/* basic routine for multiplication */
{
    INT erg = OK;
    OP ff,z;
    ff = CALLOCOBJECT();

    FORALL (z,a, {
            if (EINSP(f))
                erg += (*partf)(S_MO_S(z),b,c,S_MO_K(z));
            else {
                FREESELF(ff);
                MULT(f,S_MO_K(z),ff);
                erg += (*partf)(S_MO_S(z),b,c,ff);
                }
            } );
    FREEALL(ff);

    ENDR("m_forall_monomials_in_a");
}


INT t_forall_monomials_in_a(a,b,f,partf) OP a,b,f; INT (*partf)();
/* basic routine for multiplication */
{
    INT erg = OK;
    OP ff,z;
    ff = CALLOCOBJECT();

    FORALL (z,a, {
            if (EINSP(f))
                erg += (*partf)(S_MO_S(z),b,S_MO_K(z));
            else {
                FREESELF(ff);
                MULT(f,S_MO_K(z),ff);
                erg += (*partf)(S_MO_S(z),b,ff);
                }
            } );
    FREEALL(ff);

    ENDR("t_forall_monomials_in_a");
}



INT m_forall_monomials_in_b(a,b,c,f,partf) OP a,b,c,f; INT (*partf)();
/* basic routine for multiplication */
{
    INT erg = OK;
    OP ff,z;
    ff = CALLOCOBJECT();

    FORALL (z,b, {
            CTO(ANYTYPE,"m_forall_monomials_in_b(i4)",f);
            CTO(MONOM,"m_forall_monomials_in_b(z)",z);
            if (EINSP(f))
                {
                erg += (*partf)(a,S_MO_S(z),c,S_MO_K(z));
                }
            else {
                FREESELF(ff);
                MULT(f,S_MO_K(z),ff);
                erg += (*partf)(a,S_MO_S(z),c,ff);
                }
            } );
    FREEALL(ff);

    ENDR("m_forall_monomials_in_b");
}

INT m_forall_monomials_in_ab(a,b,c,f,partf) OP a,b,c,f; INT (*partf)();
/* basic routine for multiplication */
{
    INT erg = OK;
    OP ff,z,y;
    ff = CALLOCOBJECT();

    FORALL (y,a, {
    FORALL (z,b, {
            FREESELF(ff);
            MULT(S_MO_K(z),S_MO_K(y),ff);
            if (not EINSP(f)) {
                MULT_APPLY(f,ff);
                }
            erg += (*partf)(S_MO_S(y),S_MO_S(z),c,ff);
            } );
            } );
    FREEALL(ff);

    ENDR("m_forall_monomials_in_b");
}

INT t_loop_partition(a,b,f,intf,multf,multapplyf) OP a,b,f;
      INT (*intf)(); INT (*multf)(); INT (*multapplyf)();
/* computes the decomposition of a partition by looping over all
   parts of the partition and multiplying */
{
    INT erg = OK;
    CTO(PARTITION,"t_loop_partition(1)",a);
    CTTTTTTO(HASHTABLE,POWSYM,SCHUR,HOMSYM,ELMSYM,MONOMIAL,"t_loop_partition(2)",b);
    if (S_PA_LI(a) == 0) {
        OP m;
        m = CALLOCOBJECT();
        b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),m);
        erg += first_partition(cons_null,S_MO_S(m));
        COPY(f,S_MO_K(m));
        INSERT_HOMSYMMONOM_(m,b); /* also for other types working */
        goto ende;
        }
    else if (S_PA_LI(a) == 1)
        {
        (*intf)(S_PA_I(a,0),b,f);
        goto ende;
        }
    else {
        OP l,d;
        INT i;
        l = CALLOCOBJECT();
        d = CALLOCOBJECT();
        init_hashtable(l);
        init_hashtable(d);

        erg += (*intf)(S_PA_I(a,0),d,f);

        for (i=1;i<S_PA_LI(a)-1;i++)
            {
            erg += (*intf)(S_PA_I(a,i),l,cons_eins);
            (*multapplyf)(l,d);
            CLEAR_HASHTABLE(l);
            }

        erg += (*intf)(S_PA_I(a,S_PA_LI(a)-1),l,cons_eins);
        erg += (*multf)(d,l,b);
        FREEALL(d);
        FREEALL(l);
        }


ende:
    ENDR("t_loop_partition");
}




INT t_splitpart(a,b,f,tf,mf) OP a,b,f; INT (*tf)(), (*mf)();
/* to trans using spliting of partition */
/* both parts are transferred to the new basis and
   multiplication afterwards */
{
        INT erg = OK;
        OP m1,m2,h1,h2;
        CTO(PARTITION,"t_splitpart(1)",a);
        SYMCHECK((S_PA_LI(a) <2),"t_splitpart:partition with < 2 parts ");

        m1 = CALLOCOBJECT();
        m2 = CALLOCOBJECT();
        erg += splitpart(a,m1,m2);
        NEW_HASHTABLE(h1);
        NEW_HASHTABLE(h2);
        erg += (*tf)(m1,h1,cons_eins);
        erg += (*tf)(m2,h2,cons_eins);
        erg += (*mf)(h1,h2,b,f);
        FREEALL(m1);
        FREEALL(m2);
        FREEALL(h1);
        FREEALL(h2);
        ENDR("t_splitpart");
}

INT p_splitpart(a,b,c,f,tf,mf) OP a,b,c,f; INT (*tf)(), (*mf)();
/* to compute plethysm using spliting of partition */
/* for both parts the plethysm is computed
   multiplication afterwards */
/* works for outer homsym/powsym or elmsym
   H_ab [C] = H_a [C] * H_b[C] */
/* AK 061201 */
{
        INT erg = OK;
        CTO(PARTITION,"p_splitpart(1)",a);
        SYMCHECK((S_PA_LI(a) <2),"p_splitpart:partition with < 2 parts ");

        if (S_PA_II(a,0) == S_PA_II(a,S_PA_LI(a)-1) )
            {
            INT i,l = S_PA_LI(a);
            OP h1,h2;
            h1 = CALLOCOBJECT();init_hashtable(h1);
            h2 = CALLOCOBJECT();
            M_I_I(1,S_PA_L(a));
            erg += (*tf)(a,b,h1,cons_eins);
            M_I_I(l,S_PA_L(a));
            COPY(h1,h2);

            for (i=2;i<l;i++)
                {
                OP h3;
                h3 = CALLOCOBJECT();init_hashtable(h3);
                SWAP(h3,h2);
                erg += (*mf)(h1,h3,h2,cons_eins);
                FREEALL(h3);
                }

            erg += (*mf)(h1,h2,c,f);
            FREEALL(h1);
            FREEALL(h2);
            goto ende;
            }
        else {
            OP m1,m2,h1,h2;
            m1 = CALLOCOBJECT();
            m2 = CALLOCOBJECT();
            splitpart(a,m1,m2);
            NEW_HASHTABLE(h1);
            NEW_HASHTABLE(h2);
            erg += (*tf)(m1,b,h1,cons_eins);
            erg += (*tf)(m2,b,h2,cons_eins);
            erg += (*mf)(h1,h2,c,f);
            FREEALL(m1);
            FREEALL(m2);
            FREEALL(h1);
            FREEALL(h2);
            goto ende;
            }
ende:
        ENDR("p_splitpart");
}


INT p_splitpart2(a,b,c,f,tf,mf) OP a,b,c,f; INT (*tf)(), (*mf)();
/* to compute plethysm using spliting of partition */
/* for both parts the plethysm is computed
   multiplication afterwards */
/* works for outer powsym
   P_I [h_ab] = P_I [h_a] * P_I[h_b] */
/* AK 061201 */
{
        INT erg = OK;
        OP m1,m2,h1,h2;
        CTO(PARTITION,"p_splitpart(2)",b);
        SYMCHECK((S_PA_LI(b) <2),"p_splitpart:partition with < 2 parts ");

        m1 = CALLOCOBJECT();
        m2 = CALLOCOBJECT();
        splitpart(b,m1,m2);
        h1 = CALLOCOBJECT();init_hashtable(h1);
        h2 = CALLOCOBJECT();init_hashtable(h2);
        erg += (*tf)(a,m1,h1,cons_eins);
        erg += (*tf)(a,m2,h2,cons_eins);
        erg += (*mf)(h1,h2,c,f);
        FREEALL(m1);
        FREEALL(m2);
        FREEALL(h1);
        FREEALL(h2);
        ENDR("p_splitpart2");
}

INT txx_null__faktor(b,f) OP b,f;
/* b = b + symfunc_0 * f */
{
    INT erg = OK;
    OP m;
    CTO(ANYTYPE,"txx_null__faktor(2)",f);
    CTTTTTTO(SCHUR,MONOMIAL,HASHTABLE,HOMSYM,ELMSYM,POWSYM,"txx_null__faktor(1)",b);
    m = CALLOCOBJECT();
    erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),m);
    COPY(f,S_MO_K(m));
    erg += first_partition(cons_null,S_MO_S(m));
    if (S_O_K(b) == HASHTABLE)
        insert_scalar_hashtable(m,b,add_koeff,eq_monomsymfunc,hash_monompartition);
    else
        INSERT_LIST(m,b,add_koeff,comp_monomschur);

    ENDR("txx_null__faktor");
}



INT mxx_null__(b,c,f) OP b,c,f;
/* c = c +  symfunc_b * f */
/* called from several m.._null__ routines */
{
    INT erg = OK;
    CTTTTTTTTO(PARTITION,INTEGER,SCHUR,MONOMIAL,HASHTABLE,HOMSYM,ELMSYM,POWSYM,"mxx_null__(1)",b);
    CTTTTTTO(SCHUR,MONOMIAL,HASHTABLE,HOMSYM,ELMSYM,POWSYM,"mxx_null__(2)",c);
    CTO(ANYTYPE,"mxx_null__(3)",f);


    if (S_O_K(b) == INTEGER) {
        OP m;
        m = CALLOCOBJECT();
        b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),m);
        COPY(f,S_MO_K(m));
        erg += first_partition(b,S_MO_S(m));
        if (S_O_K(c) == HASHTABLE)
            insert_scalar_hashtable(m,c,add_koeff,eq_monomsymfunc,hash_monompartition);
        else
            INSERT_LIST(m,c,add_koeff,comp_monomschur);
        }
    else if (S_O_K(b) == PARTITION)
        {
        _NULL_PARTITION_(b,c,f);
        }
    else
        {
        OP m;
        m = CALLOCOBJECT();
        COPY(b,m);
        if (not EINSP(f))
            {
            MULT_APPLY(f,m);
            }
        if (S_O_K(c) == HASHTABLE)
            INSERT_HASHTABLE(m,c,add_koeff,eq_monomsymfunc,hash_monompartition);
        else
            INSERT_LIST(m,c,add_koeff,comp_monomschur);
        }
    ENDR("mxx_null");
}


INT giambelli_matrix(a,b) OP a,b;
/* computes the  matrix of giambelli
   a determinant of hooks */
/* AK 260603 */
{
    INT erg = OK;
    INT i,j,k;
    OP c,l,r;
    CTO(PARTITION,"giambelli_matrix(1)",a);
    c=CALLOCOBJECT();
    erg += t_VECTOR_FROBENIUS(a,c);
    r=S_V_I(S_PA_S(c),0);
    l=S_V_I(S_PA_S(c),1);
    erg += m_lh_m(S_V_L(l) ,S_V_L(r), b);
    for (i=0;i<S_M_HI(b);i++)
    for (j=0;j<S_M_LI(b);j++)
        {
        OP d,e;
        d = CALLOCOBJECT();
        e = CALLOCOBJECT();
        erg += m_il_integervector(S_V_II(l,j)+1,d);
        for (k=0;k<S_V_II(l,j);k++) M_I_I(1,S_V_I(d,k));
        M_I_I(S_V_II(r,i)+1,S_V_I(d,k));
        erg += b_ks_pa(VECTOR,d,e);
        erg += b_pa_s(e,S_M_IJ(b,i,j));
        }
    FREEALL(c);
    ENDR("giambelli_matrix");
}



INT ribbon_matrix(a,b) OP a,b;
/* computes the  matrix of ribbons,
   whose determinant is the schur function
*/
/* AK 270603 */
{
    INT erg = OK;
    INT i,j;
    OP c;
    CTO(PARTITION,"ribbon_matrix(1)",a);
    erg += durfee_size_part(a,b);
    erg += m_lh_m(b,b,b);
    c=CALLOCOBJECT();
    for (i=0;i<S_M_HI(b);i++)
    for (j=0;j<S_M_LI(b);j++)
        {
        erg += ribbon_partition(a,i,j,c);
        erg += part_part_skewschur(S_SPA_G(c),S_SPA_K(c),S_M_IJ(b,i,j));
        }
    FREEALL(c);
    ENDR("ribbon_matrix");
}


INT jacobitrudimatrix(a,b) OP a,b;
/* build the matrix of jacobi trudi */
/* AK 010703 V2.0 */
{
    INT i,j,k;
    INT erg = OK;
    CTTTO(INTEGERVECTOR,PARTITION,SKEWPARTITION,"jacobitrudimatrix(1)",a);
    if (S_O_K(a) == PARTITION)
        {
        m_lh_nm(S_PA_L(a),S_PA_L(a),b);
        for (i=0;i<S_M_HI(b);i++)
        for (j=0;j<S_M_LI(b);j++) {
             k = S_PA_II(a,S_PA_LI(a)-1-i)+j-i;
             if (k<0) continue;
             m_int_pa(k,S_M_IJ(b,i,j));
             m_pa_s(S_M_IJ(b,i,j),S_M_IJ(b,i,j));
             }
        }
    else if (S_O_K(a) == INTEGERVECTOR) /* AK 010703 */
        {
        m_lh_nm(S_V_L(a),S_V_L(a),b);
        for (i=0;i<S_M_HI(b);i++)
        for (j=0;j<S_M_LI(b);j++) {
             k = S_V_II(a,i)+j-i;
             if (k<0) continue;
             m_int_pa(k,S_M_IJ(b,i,j));
             m_pa_s(S_M_IJ(b,i,j),S_M_IJ(b,i,j));
             }
        }
    else if (S_O_K(a) = SKEWPARTITION) /* AK 010703 */
        {
        OP g = S_SPA_G(a);
        OP kl = S_SPA_K(a);
        m_lh_nm(S_PA_L(g),S_PA_L(g),b);
        for (i=0;i<S_M_HI(b);i++)
        for (j=0;j<S_M_LI(b);j++) {
             k = S_PA_II(g,S_PA_LI(g)-1-i)+j-i;
             /* now subtract smaller part */
             if (S_PA_LI(kl)-j-1 >= 0) k = k - S_PA_II(kl,S_PA_LI(kl)-j-1);
             if (k<0) continue;
             m_int_pa(k,S_M_IJ(b,i,j));
             m_pa_s(S_M_IJ(b,i,j),S_M_IJ(b,i,j));
             }
        }
    ENDR("jacobitrudimatrix");
}




