/* file: macro.h symmetrica source code */
#ifndef MACRO_H

#ifdef SYMMAGMA
#define SYM_MALLOC(a) mem_malloc(a)
#else
#define SYM_MALLOC(a) SYM_malloc(a)
#endif

#ifdef SYMMAGMA
#define SYM_FREE(a) mem_free(a)
#else
#define SYM_FREE(a) SYM_free(a)
#endif


#define S_O_S(a) (((OP)(a))->ob_self)
#define S_O_K(a) (((OP)(a))->ob_kind)
#define C_O_K(a,b) ((a)->ob_kind = (OBJECTKIND)(b))
#define C_O_S(a,b) S_O_S(a)=(b)
#define B_KS_O(a,b,c) do { C_O_S(c,b); C_O_K(c,a); } while(0)

#define EMPTYP(a) ((a)->ob_kind == (OBJECTKIND)0)

#define EQ_INTEGER(a,b) \
    ( S_O_K(b) == INTEGER ?(S_I_I(a) == S_I_I(b)):(comp_integer(a,b) == 0))
#define EQ_PARTITION(a,b) \
    ( S_O_K(b) == PARTITION ? eq_partition_partition(a,b) :FALSE )
#define EQ_LONGINT(a,b) \
    ( S_O_K(b) == LONGINT ? eq_longint_longint(a,b) :(comp_longint(a,b) == 0))
#define EQ_FF(a,b) \
    ( comp_ff(a,b) == 0 )
#define EQ(a,b) \
    (\
        S_O_K(a) == INTEGER? EQ_INTEGER(a,b) : \
          (\
          S_O_K(a) == LONGINT? EQ_LONGINT(a,b) :\
            ( \
            S_O_K(a) == PARTITION?EQ_PARTITION(a,b): \
                (\
                S_O_K(a)==FF?EQ_FF(a,b): eq(a,b) \
                )\
            ) \
          )\
    )



#define S_I_I(a) (((a)->ob_self).ob_INT)
#define C_I_I(a,b) (a)->ob_self.ob_INT = (INT)(b)
#define M_I_I(a,b)  (\
            ((b)->ob_self.ob_INT = (INT)(a)),\
             ((b)->ob_kind = INTEGER)\
            )

#define DEC_INTEGER(a) (((a)->ob_self).ob_INT --)
#define INC_INTEGER(a) (((a)->ob_self).ob_INT ++)


#define COPY_INTEGER(a,b) M_I_I(S_I_I(a),b)
#define FREESELF_INTEGER(a) C_O_K((a),0)
#define NULLP_INTEGER(a) (S_I_I(a) == 0)
#define POSP_INTEGER(a) (S_I_I(a) > (INT)0)
#define NEGP_INTEGER(a) (S_I_I(a) < (INT)0)
#define EINSP_INTEGER(a) (S_I_I(a) == (INT)1)
#define NEGEINSP_INTEGER(a) (S_I_I(a) == (INT)-1)
#define COMP_INTEGER_INTEGER(a,b) ((S_I_I(a) == S_I_I(b))? 0L :\
    ((S_I_I(a) > S_I_I(b))? 1L : -1L) )
#define INTEGERP(a) (S_O_K(a) == INTEGER)


#define S_V_S(a) ((((a)->ob_self).ob_vector)->v_self)
#define C_V_S(a,b) (((((a)->ob_self).ob_vector)->v_self) = ((OP) b))
#define S_V_L(a) ((((a)->ob_self).ob_vector)->v_length)
#define C_V_L(a,b) (((((a)->ob_self).ob_vector)->v_length) = (b))
#define S_V_I(a,i) (((((a)->ob_self).ob_vector)->v_self)+(i))
#define C_V_I(a,i,b) ( *(((((a)->ob_self).ob_vector)->v_self)+(i)) = *(b))
#define S_V_II(a,i) (((((((a)->ob_self).ob_vector)->v_self)+(i))\
            ->ob_self).ob_INT)
#define S_V_LI(a) ((((((a)->ob_self).ob_vector)->v_length)->ob_self).ob_INT)

#define B_O_V(a,b) \
do { \
struct vector *callocvectorstruct();\
FREESELF(b);C_O_K(b,VECTOR);b->ob_self.ob_vector = callocvectorstruct();\
C_V_S(b,a);C_V_L(b,CALLOCOBJECT()); M_I_I(1,S_V_L(b)); } while(0)



#define INTEGRALP(a) (       (S_O_K(a) == LONGINT) \
            || (S_O_K(a) == INTEGER)  \
                      )
#define VECTORP_CO1(a) (       (S_O_K(a) == VECTOR) \
            || (S_O_K(a) == WORD) \
            || (S_O_K(a) == QUEUE) \
            || (S_O_K(a) == INTEGERVECTOR) \
            || (S_O_K(a) == COMPOSITION) \
            || (S_O_K(a) == HASHTABLE) \
            || (S_O_K(a) == LAURENT) \
            || (S_O_K(a) == KRANZ) \
            || (S_O_K(a) == SUBSET) \
            || (S_O_K(a) == FF) \
        )

#define VECTORP(a) ((a == NULL)? FALSE: VECTORP_CO1(a))

#define LISTP_CO1(a) (         (S_O_K(a) == LIST) \
            || (S_O_K(a) == SCHUR) \
            || (S_O_K(a) == GRAL) \
            || (S_O_K(a) == POLYNOM) \
            || (S_O_K(a) == SCHUBERT) \
            || (S_O_K(a) == MONOPOLY) \
            || (S_O_K(a) == ELM_SYM) \
            || (S_O_K(a) == POW_SYM) \
            || (S_O_K(a) == MONOMIAL) \
            || (S_O_K(a) == HOM_SYM) \
        )
#define LISTP(a) ((a == NULL)? FALSE: LISTP_CO1(a))

#define POLYP(a) (          \
            (S_O_K(a) == SCHUR) \
            || (S_O_K(a) == GRAL) \
            || (S_O_K(a) == POLYNOM) \
            || (S_O_K(a) == SCHUBERT) \
            || (S_O_K(a) == MONOPOLY) \
            || (S_O_K(a) == ELM_SYM) \
            || (S_O_K(a) == POW_SYM) \
            || (S_O_K(a) == MONOMIAL) \
            || (S_O_K(a) == HOM_SYM) \
        )



#define S_L_S(a) ((((a)->ob_self).ob_list)->l_self)
#define C_L_S(a,b) (((((a)->ob_self).ob_list)->l_self) = (OP)(b))
#define S_L_N(a) ((((a)->ob_self).ob_list)->l_next)
#define C_L_N(a,b) (((((a)->ob_self).ob_list)->l_next) = (OP)(b))

#define S_PA_C(a) ((unsigned char *)  (S_PA_S(a)))
#define S_PA_K(a) ((a)->ob_self.ob_partition->pa_kind)
#define C_PA_K(a,b) (((((a)->ob_self).ob_partition)->pa_kind) = b)
#define C_PA_HASH(a,b) (((((a)->ob_self).ob_partition)->pa_hash) = b)
#define S_PA_S(a) ((a)->ob_self.ob_partition->pa_self)
#define S_PA_HASH(a) ((a)->ob_self.ob_partition->pa_hash)
#define C_PA_S(a,b) ((a)->ob_self.ob_partition->pa_self = b)
#define S_PA_I(a,i) S_V_I((a)->ob_self.ob_partition->pa_self,i)
#define S_PA_II(a,i)  \
((S_O_K(a) == CHARPARTITION || S_O_K(a) == CHAR_AUG_PART) ?   \
        (INT)(S_PA_C(a)[i+1]) : \
        S_V_II((a)->ob_self.ob_partition->pa_self,i))
#define S_PA_CII(a,i) (S_PA_C(a)[i+1])
#define S_PA_CI(a,i) (S_PA_C(a)+i+1)
#define S_PA_CL(a) (S_PA_C(a)[0])
#define S_PA_L(a) S_V_L(S_PA_S(a))
#define S_PA_LI(a) \
((S_O_K(a) == CHARPARTITION || S_O_K(a) == CHAR_AUG_PART) ?  \
        (INT)(S_PA_C(a)[0]) : \
        S_V_LI(S_PA_S(a)))
#define INC_PARTITION(a) inc_vector(S_PA_S(a))
#define DEC_PARTITION(a) dec_integervector(S_PA_S(a))

#ifdef FAST
#define PART_CHECK_KIND(t,a,b)
#else
#define PART_CHECK_KIND(t,a,b)\
    CTO(PARTITION,t,a);\
    if (S_PA_K(a) != b)\
        wrong_kind_part(t,a,b);
#endif



extern long partition_speichersize,partition_speicherindex,mem_counter_part;
extern struct partition **partition_speicher;

#define FREEPARTITION(d)\
FREE_MEMMANAGER(struct partition *,\
                    partition_speicher,\
                    partition_speicherindex,\
                    partition_speichersize,\
                    mem_counter_part,\
                    d)

#define B_KS_PA(a,b,c) \
    do {\
    C_O_K(c,PARTITION);\
    CALLOC_MEMMANAGER(struct partition,\
                      partition_speicher,\
                      partition_speicherindex,\
                      mem_counter_part,\
                      (c ->ob_self).ob_partition);\
    /* (c ->ob_self).ob_partition=CALLOCPARTITION();*/\
    C_PA_K(c,a);\
    C_PA_S(c,b);\
    C_PA_HASH(c,-1L); \
    } while (0)

#define M_KL_PA(a,b,c) \
do { B_KS_PA(a,CALLOCOBJECT(),c);\
        erg += m_l_v(b,S_PA_S(c));\
        C_O_K(S_PA_S(c),INTEGERVECTOR) );\
while(0)

#define B_KL_PA(a,b,c) \
do { B_KS_PA(a,CALLOCOBJECT(),c);\
        erg += b_l_v(b,S_PA_S(c));\
        C_O_K(S_PA_S(c),INTEGERVECTOR);\
} while(0)




#define S_P_K(a) ((((a)->ob_self).ob_permutation)->p_kind)
#define C_P_K(a,b) (((((a)->ob_self).ob_permutation)->p_kind) = b)
#define S_P_S(a) ((((a)->ob_self).ob_permutation)->p_self)
#define C_P_S(a,b) (((((a)->ob_self).ob_permutation)->p_self) = b)
#define S_P_I(a,i) S_V_I(((((a)->ob_self).ob_permutation)->p_self),(i))
#define S_P_II(a,i) S_V_II(((((a)->ob_self).ob_permutation)->p_self),(i))
#define S_P_L(a) S_V_L(S_P_S(a))
#define S_P_LI(a) S_V_LI(S_P_S(a))

/* kranz produkt AK 120804 */
#define S_KR_G(a) S_V_I(a,0)
#define S_KR_GI(a,i) S_P_I(S_V_I(a,0),i)
#define S_KR_GL(a) S_P_L(S_V_I(a,0))
#define S_KR_GLI(a) S_P_LI(S_V_I(a,0))
#define S_KR_V(a) S_V_I(a,1)
#define S_KR_VL(a) S_V_L(S_V_I(a,1))
#define S_KR_VLI(a) S_V_LI(S_V_I(a,1))
#define S_KR_I(a,i) S_V_I(S_V_I(a,1),i)

#define S_M_L(a) ((((a)->ob_self).ob_matrix)->m_length)
#define C_M_L(a,b) (((((a)->ob_self).ob_matrix)->m_length)=(b))
#define S_M_LI(a) ((((((a)->ob_self).ob_matrix)->m_length)->ob_self).ob_INT)
#define S_M_H(a) ((((a)->ob_self).ob_matrix)->m_height)
#define S_M_HASH(a) ((((a)->ob_self).ob_matrix)->m_hash)
#define C_M_H(a,b) (((((a)->ob_self).ob_matrix)->m_height)=(b))
#define S_M_HI(a) ((((((a)->ob_self).ob_matrix)->m_height)->ob_self).ob_INT)
#define S_M_S(a) ((((a)->ob_self).ob_matrix)->m_self)
#define C_M_S(a,b) (((((a)->ob_self).ob_matrix)->m_self)=(b))
#define C_M_HASH(a,b) (((((a)->ob_self).ob_matrix)->m_hash)=(b))
#define S_M_IJ(a,i,j) ( ((((a)->ob_self).ob_matrix)->m_self)\
     + ((((((a)->ob_self).ob_matrix)->m_length)->ob_self).ob_INT)\
     * (i) + (j) )
#define S_M_IJI(a,i,j) S_I_I(S_M_IJ(a,i,j))
#define MATRIXP(a) ((S_O_K(a) == MATRIX) || (S_O_K(a) == KRANZTYPUS)\
    || (S_O_K(a) == KOSTKA) || (S_O_K(a) == INTEGERMATRIX) )

#define S_MO_S(a) (((a->ob_self).ob_monom)->mo_self)
#define C_MO_S(a,b) ((((a->ob_self).ob_monom)->mo_self)=(b))
#define S_MO_K(a) (((a->ob_self).ob_monom)->mo_koeff)
#define C_MO_K(a,b) ((((a->ob_self).ob_monom)->mo_koeff)=(b))
#define S_MO_KI(a) ((((a->ob_self).ob_monom)->mo_koeff)->ob_self.ob_INT)
#define S_MO_SI(a,i) S_V_I(S_MO_S(a),(i))
#define S_MO_SII(a,i) S_V_II(S_MO_S(a),(i))
#define S_MO_SL(a) S_V_L(S_MO_S(a))
#define S_MO_SLI(a) S_V_LI(S_MO_S(a))
#define COPY_MONOM(a,b) M_SK_MO(S_MO_S(a),S_MO_K(a),b)

#define B_SK_MO(a,b,c)\
do { \
      C_O_K(c,MONOM);\
      CALLOC_MEMMANAGER(struct monom, monom_speicher, monom_speicherindex,mem_counter_monom ,(c->ob_self).ob_monom);\
      C_MO_S(c,a) ; \
      C_MO_K(c,b); \
} while(0)

extern long monom_speicherindex,mem_counter_monom,monom_speichersize;
extern struct monom **monom_speicher;

#define FREEMONOM(v) \
FREE_MEMMANAGER(struct monom *,\
 monom_speicher,\
 monom_speicherindex,\
 monom_speichersize,\
 mem_counter_monom,\
 v)


#define FREESELF_MONOM(a)\
do {\
    if (S_O_K(S_MO_S(a)) == PARTITION) \
        { erg+= freeself_partition(S_MO_S(a)); }\
    else if (S_O_K(S_MO_S(a)) == INTEGERMATRIX) \
        { erg+= freeself_integermatrix(S_MO_S(a)); }\
    else erg += freeself(S_MO_S(a));\
    FREE_EMPTY_OBJECT(S_MO_S(a));\
    if (S_O_K(S_MO_K(a)) == INTEGER) C_O_K(S_MO_K(a),EMPTY);\
    else if (S_O_K(S_MO_K(a)) == LONGINT) erg += freeself_longint(S_MO_K(a));\
    else if (S_O_K(S_MO_K(a)) == BRUCH) erg += freeself_bruch(S_MO_K(a));\
    else if (S_O_K(S_MO_K(a)) == FF) erg += freeself_ff(S_MO_K(a));\
    else  erg += freeself(S_MO_K(a));\
    FREE_EMPTY_OBJECT(S_MO_K(a));\
    FREEMONOM(((a)->ob_self).ob_monom);\
    C_O_K(a,EMPTY);\
} while(0)

#define S_B_I(a) ((INT)((((a)->ob_self).ob_bruch)->b_info))
#define S_B_O(a) ((OP)((((a)->ob_self).ob_bruch)->b_oben))
#define S_B_OI(a) ((INT)((((((a)->ob_self).ob_bruch)->b_oben)->ob_self).ob_INT))
#define S_B_U(a) ((OP)((((a)->ob_self).ob_bruch)->b_unten))
#define S_B_UI(a) ((INT)((((((a)->ob_self).ob_bruch)->b_unten)->ob_self).ob_INT))
#define FREESELF_BRUCH(a) (freeall(S_B_O(a)),freeall(S_B_U(a)),\
        free(S_O_S(a).ob_bruch),C_O_K(a,0),OK)
#define C_B_I(a,b) (((((a)->ob_self).ob_bruch)->b_info)=(INT)(b))
#define C_B_O(a,b) (((((a)->ob_self).ob_bruch)->b_oben)=(b))
#define C_B_U(a,b) (((((a)->ob_self).ob_bruch)->b_unten)=(b))
#define EINSP_BRUCH(a) (EQ(S_B_O((a)),S_B_U((a))))
#define BRUCHP(a) (S_O_K(a) == BRUCH)

#define S_S_S(a)  ((((((((a)->ob_self).ob_list)->l_self))\
        ->ob_self).ob_monom)->mo_self)
#define S_S_SL(a) S_PA_L( ((((((((a)->ob_self).ob_list)->l_self))\
        ->ob_self).ob_monom)->mo_self)\
            )
#define S_S_SLI(a) S_PA_LI( ((((((((a)->ob_self).ob_list)->l_self))\
        ->ob_self).ob_monom)->mo_self)\
            )
#define S_S_SI(a,i) S_PA_I( ((((((((a)->ob_self).ob_list)->l_self))\
        ->ob_self).ob_monom)->mo_self)\
        ,(i))
#define S_S_SII(a,i) S_PA_II( ((((((((a)->ob_self).ob_list)->l_self))\
        ->ob_self).ob_monom)->mo_self)\
        ,(i))
#define S_S_K(a) (S_MO_K(((((a)->ob_self).ob_list)->l_self)))
#define C_S_K(a,b) (C_MO_K(((((a)->ob_self).ob_list)->l_self),(b)))
#define C_S_S(a,b) (C_MO_S(((((a)->ob_self).ob_list)->l_self),(b)))
#define S_S_KI(a) (S_MO_KI(((((a)->ob_self).ob_list)->l_self)))
#define S_S_N(a) (S_L_N(a))
#define C_S_N(a,b) (C_L_N((a),(b)))

#define S_SCH_S(a) (S_MO_S(S_L_S(a)))
#define S_SCH_SL(a) (S_P_L(S_MO_S(S_L_S(a))))
#define S_SCH_SLI(a) (S_P_LI(S_MO_S(S_L_S(a))))
#define S_SCH_SI(a,i) S_P_I(S_SCH_S(a),(i))
#define S_SCH_SII(a,i) S_P_II(S_SCH_S(a),(i))
#define S_SCH_K(a) (S_MO_K(S_L_S(a)))
#define C_SCH_K(a,b) (C_MO_K(S_L_S(a),(b)))
#define S_SCH_KI(a) (S_MO_KI(S_L_S(a)))
#define S_SCH_N(a) (S_L_N(a))
#define C_SCH_N(a,b) (C_L_N((a),(b)))


#define POLYNOMP(a) (S_O_K(a) == POLYNOM)
#define S_PO_S(a) (S_MO_S(S_L_S(a)))
#define S_PO_SI(a,i) (S_V_I(S_MO_S(S_L_S(a)),i))
#define S_PO_SIJ(a,i,j) (S_M_IJ(S_MO_S(S_L_S(a)),i,j))
#define S_PO_SIJI(a,i,j) (S_M_IJI(S_MO_S(S_L_S(a)),i,j))
#define S_PO_SII(a,i) (S_V_II(S_MO_S(S_L_S(a)),i))
#define S_PO_SL(a) (S_V_L(S_MO_S(S_L_S(a))))
#define S_PO_SLI(a) (S_V_LI(S_MO_S(S_L_S(a))))
#define S_PO_K(a) (S_MO_K(S_L_S(a)))
#define S_PO_KI(a) (S_MO_KI(S_L_S(a)))
#define S_PO_N(a) (S_L_N(a))
#define C_PO_N(a,b) (C_L_N((a),(b)))
#define C_PO_K(a,b) (C_MO_K(S_L_S(a),(b)))
#define COPY_POLYNOM(a,b) copy_list(a,b)
#define M_SK_MO(a,b,c) (b_sk_mo(callocobject(),callocobject(),c),\
            copy(b,S_MO_K(c)),\
            copy(a,S_MO_S(c)))

#define S_SPA_G(a) ((((a)->ob_self).ob_skewpartition)->spa_gross)
#define S_SPA_GL(a) S_PA_L((((a)->ob_self).ob_skewpartition)->spa_gross)
#define S_SPA_GLI(a) S_PA_LI((((a)->ob_self).ob_skewpartition)->spa_gross)
#define S_SPA_GI(a,i) S_PA_I(((((a)->ob_self).ob_skewpartition)->spa_gross),i)
#define S_SPA_GII(a,i) S_PA_II(((((a)->ob_self).ob_skewpartition)->spa_gross),i)
#define S_SPA_GS(a) S_PA_S((((a)->ob_self).ob_skewpartition)->spa_gross)
#define S_SPA_K(a) ((((a)->ob_self).ob_skewpartition)->spa_klein)
#define S_SPA_KL(a) S_PA_L((((a)->ob_self).ob_skewpartition)->spa_klein)
#define S_SPA_KLI(a) S_PA_LI((((a)->ob_self).ob_skewpartition)->spa_klein)
#define S_SPA_KI(a,i) S_PA_I(S_SPA_K(a),i)
#define S_SPA_KII(a,i) S_PA_II(S_SPA_K(a),i)
#define S_SPA_KS(a) S_PA_S((((a)->ob_self).ob_skewpartition)->spa_klein)

#define S_T_S(a) ((((a)->ob_self).ob_tableaux)->t_self)
#define S_T_U(a) ((((a)->ob_self).ob_tableaux)->t_umriss)
#define S_T_IJ(a,i,j) S_M_IJ(S_T_S(a),i,j)
#define S_T_H(a) S_M_H(S_T_S(a))
#define S_T_HI(a) S_M_HI(S_T_S(a))
#define S_T_L(a) S_M_L(S_T_S(a))
#define S_T_LI(a) S_M_LI(S_T_S(a))
#define S_T_IJI(a,i,j) S_M_IJI(S_T_S(a),i,j)
#define S_T_UG(a) S_SPA_G(S_T_U(a))
#define S_T_UGLI(a) S_SPA_GLI(S_T_U(a))
#define S_T_UGL(a) S_SPA_GL(S_T_U(a))
#define S_T_UGII(a,i) S_SPA_GII(S_T_U(a),i)
#define S_T_UGI(a,i) S_SPA_GI(S_T_U(a),i)
#define S_T_UK(a) S_SPA_K(S_T_U(a))
#define S_T_UKLI(a) S_SPA_KLI(S_T_U(a))
#define S_T_UKL(a) S_SPA_KL(S_T_U(a))
#define S_T_UKII(a,i) S_SPA_KII(S_T_U(a),i)
#define S_T_UKI(a,i) S_SPA_KI(S_T_U(a),i)
#define S_T_UL(a) S_PA_L(S_T_U(a))
#define S_T_ULI(a) S_PA_LI(S_T_U(a))
#define S_T_UII(a,i) S_PA_II(S_T_U(a),i)
#define S_T_UI(a,i) S_PA_I(S_T_U(a),i)
#define TABLEAUXP(a) (S_O_K(a) == TABLEAUX)

#define S_W_I(a,i) S_V_I(a,i)
#define S_W_II(a,i) S_V_II(a,i)
#define S_W_L(a) S_V_L(a)
#define S_W_LI(a) S_V_LI(a)
#define S_W_S(a) S_V_S(a)
#define S_LA_I(a,i) S_V_I(a,i)
#define S_LA_II(a,i) S_V_II(a,i)
#define S_LA_L(a) S_V_L(a)
#define S_LA_LI(a) S_V_LI(a)
#define S_LA_S(a) S_V_S(a)
#define b_l_w(a,b) (b_l_v(a,b),C_O_K(b,WORD),OK)
#define m_l_w(a,b) (m_l_v(a,b),C_O_K(b,WORD),OK)
#define m_il_w(a,b) (m_il_v(a,b),C_O_K(b,WORD),OK)
#define m_il_nw(a,b) (m_il_nv(a,b),C_O_K(b,WORD),OK)
#define m_l_nw(a,b) (m_l_nv(a,b),C_O_K(b,WORD),OK)
#define b_l_la(a,b) (b_l_v(a,b),C_O_K(b,LAURENT),OK)
#define m_l_la(a,b) (m_l_v(a,b),C_O_K(b,LAURENT),OK)
#define m_il_la(a,b) (m_il_v(a,b),C_O_K(b,LAURENT),OK)
#define m_il_nla(a,b) (m_il_nv(a,b),C_O_K(b,LAURENT),OK)
#define m_l_nla(a,b) (m_l_nv(a,b),C_O_K(b,LAURENT),OK)

#define S_SC_D(a) ((((a)->ob_self).ob_symchar)->sy_dimension)
#define S_SC_DI(a) S_I_I((((a)->ob_self).ob_symchar)->sy_dimension)
#define S_SC_W(a) ((((a)->ob_self).ob_symchar)->sy_werte)
#define S_SC_WI(a,i) S_V_I(((((a)->ob_self).ob_symchar)->sy_werte),i)
#define S_SC_WII(a,i) S_V_II(((((a)->ob_self).ob_symchar)->sy_werte),i)
#define S_SC_P(a) ((((a)->ob_self).ob_symchar)->sy_parlist)
#define S_SC_PI(a,i) S_V_I(((((a)->ob_self).ob_symchar)->sy_parlist),i)
#define S_SC_PLI(a) S_V_LI(((((a)->ob_self).ob_symchar)->sy_parlist))
#define S_SC_WLI(a) S_V_LI(((((a)->ob_self).ob_symchar)->sy_werte))

/* MD */
#define S_N_S(a) ((((a)->ob_self).ob_number)->n_self)
#define C_N_S(a,b) (((((a)->ob_self).ob_number)->n_self) = (b))
#define S_N_D(a) (((((a)->ob_self).ob_number)->n_data).o_data)
#define C_N_D(a,b) ((((((a)->ob_self).ob_number)->n_data).o_data) = (b))
#define S_N_DC(a) (((((a)->ob_self).ob_number)->n_data).c_data)
#define S_N_DCI(a) ((((((a)->ob_self).ob_number)->n_data).c_data)->index)
#define S_N_DCII(a) S_I_I(S_N_DCI(a))
#define S_N_DCD(a) ((((((a)->ob_self).ob_number)->n_data).c_data)->deg)
#define S_N_DCP(a) ((((((a)->ob_self).ob_number)->n_data).c_data)->poly)
#define    OBJECTREAD_CYCLO(f,a)    objectread_number((f),(a),CYCLOTOMIC)
#define    OBJECTREAD_SQRAD(f,a)    objectread_number((f),(a),SQ_RADICAL)
#define    EINSP_MONOPOLY(a)    eq_fieldobject_int(MONOPOLY,(a),1L)
#define    EINSP_CYCLO(a)        eq_fieldobject_int(CYCLOTOMIC,(a),1L)
#define    EINSP_SQRAD(a)        eq_fieldobject_int(SQ_RADICAL,(a),1L)
#define    NEGEINSP_MONOPOLY(a)    eq_fieldobject_int(MONOPOLY,(a),-1L)
#define    NEGEINSP_CYCLO(a)    eq_fieldobject_int(CYCLOTOMIC,(a),-1L)
#define    NEGEINSP_SQRAD(a)    eq_fieldobject_int(SQ_RADICAL,(a),-1L)




#define evenp(a) even(a)
#define oddp(a) odd(a)
#define first_ym(a,b) ym_min(a,b)

#define addinvers_schubert(a,b) addinvers_polynom(a,b)
#define addinvers_schur(a,b) addinvers_polynom(a,b)
#define einsp_schur(a) einsp_symfunc(a)
#define add_apply_schur(a,b) add_apply_symfunc(a,b)
#define add_apply_schur_schur(a,b) add_apply_symfunc_symfunc(a,b)
#define m_part_schur(a,b) m_skn_s(a,cons_eins,NULL,b)

/* for MONOMIAL */
#define m_v_mon(a,b) (m_v_s(a,b),C_O_K(b,MONOMIAL),OK)
#define b_skn_mon(a,b,c,d) (b_skn_s(a,b,c,d),C_O_K(d,MONOMIAL),OK)
#define m_skn_mon(a,b,c,d) (m_skn_s(a,b,c,d),C_O_K(d,MONOMIAL),OK)
#define addinvers_monomial(a,b) addinvers_polynom(a,b)
#define add_apply_monomial(a,b) add_apply_symfunc(a,b)
#define add_apply_monomial_monomial(a,b) add_apply_symfunc_symfunc(a,b)
#define m_part_monomial(a,b) m_skn_mon(a,cons_eins,NULL,b)
#define einsp_monomial(a) einsp_symfunc(a)

/* for ELM_SYM */
#define m_v_e(a,b) (m_v_s(a,b),C_O_K(b,ELM_SYM),OK)
#define add_monom_elmsym(a,b,c) add_monom_schur(a,b,c)
#define addinvers_elmsym(a,b) addinvers_polynom(a,b)
#define b_skn_e(a,b,c,d) (b_skn_s(a,b,c,d),C_O_K(d,ELM_SYM),OK)
#define m_skn_e(a,b,c,d) (m_skn_s(a,b,c,d),C_O_K(d,ELM_SYM),OK)
#define add_apply_elmsym(a,b) add_apply_symfunc(a,b)
#define add_apply_elmsym_elmsym(a,b) add_apply_symfunc_symfunc(a,b)
#define m_part_elmsym(a,b) m_skn_e(a,cons_eins,NULL,b)
#define einsp_elmsym(a) einsp_symfunc(a)

/* for POW_SYM */
#define m_v_ps(a,b) (m_v_s(a,b),C_O_K(b,POW_SYM),OK)
#define add_monom_powsym(a,b,c) add_monom_schur(a,b,c)
#define addinvers_powsym(a,b) addinvers_polynom(a,b)
#define b_skn_ps(a,b,c,d) (b_skn_s(a,b,c,d),C_O_K(d,POW_SYM),OK)
#define m_skn_ps(a,b,c,d) (m_skn_s(a,b,c,d),C_O_K(d,POW_SYM),OK)
#define compute_powsym_with_alphabet(a,b,c) \
    compute_power_with_alphabet(a,b,c)
#define add_apply_powsym(a,b) add_apply_symfunc(a,b)
#define add_apply_powsym_powsym(a,b) add_apply_symfunc_symfunc(a,b)
#define m_part_powsym(a,b) m_skn_ps(a,cons_eins,NULL,b)
#define einsp_powsym(a) einsp_symfunc(a)

/* for HOM_SYM */
#define m_v_h(a,b) (m_v_s(a,b),C_O_K(b,HOM_SYM),OK)
#define add_monom_homsym(a,b,c) add_monom_schur(a,b,c)
#define addinvers_homsym(a,b) addinvers_schur(a,b)
#define b_skn_h(a,b,c,d) (b_skn_s(a,b,c,d),C_O_K(d,HOM_SYM),OK)
#define m_skn_h(a,b,c,d) (m_skn_s(a,b,c,d),C_O_K(d,HOM_SYM),OK)
#define compute_homsym_with_alphabet(a,b,c) \
    compute_complete_with_alphabet(a,b,c)
#define add_apply_homsym(a,b) add_apply_symfunc(a,b)
#define add_apply_homsym_homsym(a,b) add_apply_symfunc_symfunc(a,b)
#define m_part_homsym(a,b) m_skn_h(a,cons_eins,NULL,b)
#define einsp_homsym(a) einsp_symfunc(a)

/* for nc.c */
#define S_NC_GL(a) S_V_I(a,0L)
#define S_NC_C(a) S_V_I(a,1L)
#define SYM_GL(a) (S_V_II(a,0L)==1L) /* true falls sym */
#define ALT_GL(a) (S_V_II(a,0L)==2L) /* true falls alt */
#define KRANZ_GL(a) (S_V_II(a,0L)==3L) /* true falls kranz */
#define CYCLIC_GL(a) (S_V_II(a,0L)==4L) /* true falls cyclic */
#define GLNQ_GL(a) (S_V_II(a,0L)==5L) /* true falls gl(n,q) */
#define S_GL_SYM_A(a) S_V_I(a,1L)
#define S_GL_ALT_A(a) S_V_I(a,1L)
#define S_GL_CYCLIC_A(a) S_V_I(a,1L)
#define S_GL_KRANZ_A(a) S_GL_SYM_A(S_V_I(S_V_I(a,1L),0L))
#define S_GL_KRANZ_GLA(a) (S_V_I(S_V_I(a,1L),0L))
#define S_GL_KRANZ_B(a) S_GL_SYM_A(S_V_I(S_V_I(a,1L),1L))
#define S_GL_KRANZ_GLB(a) (S_V_I(S_V_I(a,1L),1L))
#define S_GL_GLNQ_N(a) (S_V_I(S_V_I(a,1L),0L))
#define S_GL_GLNQ_Q(a) (S_V_I(S_V_I(a,1L),1L))



/* for ga.c */

#define m_skn_gral(a,b,c,d) ( m_skn_po(a,b,c,d), C_O_K(d,GRAL), OK )
#define b_skn_gral(a,b,c,d) ( b_skn_po(a,b,c,d), C_O_K(d,GRAL), OK )
#define hplus_hecke(a,b) hplus(a,b)

/* for ff.c */
#define S_FF_C(a) S_V_I(a,0)
#define S_FF_CI(a) S_V_II(a,0)
#define S_FF_IP(a) S_O_S(S_V_I(a,1)).ob_INTpointer
#define C_FF_IP(a,p) S_O_S(S_V_I(a,1)).ob_INTpointer=(INT*)p
#define S_FF_II(a,i) *((S_O_S(S_V_I(a,1)).ob_INTpointer) + i)
#define S_FF_DI(a) S_FF_II(a,0)
#define S_FF_ME(a) S_V_I(a,2)
#define S_FF_MEI(a) S_V_II(a,2)

/* for galois.c */
#define S_GR_CI(a) S_V_II(a,1)
#define S_GR_C(a) S_V_I(a,1)
#define S_GR_D(a) S_V_I(a,0)
#define S_GR_DI(a) S_V_II(a,0)



/* for longint */

extern long loc_index, loc_size,loc_counter;
extern struct loc **loc_speicher;

extern long longint_speicherindex,mem_counter_loc,longint_speichersize;
extern struct longint **longint_speicher;



#define FREE_LOC(a) \
    FREE_MEMMANAGER(struct loc *,loc_speicher,loc_index,loc_size,loc_counter,a)

#define LOCSGN(lx) ( ((lx)->w2 || (lx)->w1 || (lx)->w0 ) ? 1 : 0 )

#define LOCNULL(lx) ((lx)->w0=0,(lx)->w1=0,(lx)->w2=0)

#define LOCHOLE(aloc) \
do {\
CALLOC_MEMMANAGER(struct loc,loc_speicher,loc_index,loc_counter,*(aloc));\
LOCNULL(*(aloc));\
(*(aloc))->nloc = NULL;\
} while(0)

#define GANZDEFINIERE(x) \
do {\
 (x)->signum = (signed char)0;\
    (x)->laenge = (INT)1;\
     (x)->floc = NULL;\
      LOCHOLE(&((x)->floc)) ;\
} while (0)

#define INTGANZ(x) \
    ( x->signum < 0 ?\
      - x->floc->w0 - x->floc->w1 * LO_B - x->floc->w2 * LO_B * LO_B :\
      (x->floc->w0&BMINUSEINS)\
              +(x->floc->w1&BMINUSEINS) * LO_B\
              +(x->floc->w2&BMINUSEINS) * LO_B * LO_B )

#define T_LONGINT_INT(a)\
    if ((S_O_S(a).ob_longint) ->laenge == 1)\
        if (S_O_S(a).ob_longint ->floc ->w2 <= 1) /* AK 051101 */\
        {\
        INT wert = INTGANZ(S_O_S(a).ob_longint);\
        FREESELF(a);\
        M_I_I(wert,a);\
        }


#define INIT_LONGINT(l) \
 do { \
 C_O_K(l,LONGINT); \
 CALLOC_MEMMANAGER(struct longint, longint_speicher,\
                   longint_speicherindex,mem_counter_loc,(l->ob_self).ob_longint);\
 GANZDEFINIERE(S_O_S(l).ob_longint);\
 } while(0)



#define t_INTEGER_LONGINT(a,b) m_i_longint(a,b)

#ifdef FAST
#define SYMCHECK(a,b)
#define COP(text,object)
#define CTO(type,text,object)
#define CTTO(type,type2,text,object)
#define CTTTO(type,type2,type3,text,object)
#define CTTTTO(type,type2,type3,type4,text,object)
#define CTTTTTO(type,type2,type3,type4,type5,text,object)
#define CTTTTTTO(type,type2,type3,type4,type5,type6,text,object)
#define CTTTTTTTO(type,type2,type3,type4,type5,type6,type7,text,object)
#define CTTTTTTTTO(type,type2,type3,type4,type5,type6,type7,type8,text,object)
#define TCTO(type,text,object)
#define TCTTO(type,type2,text,object)
#define TCTTTO(type,type2,type3,text,object)
#define TCTTTTO(type,type2,type3,type4,text,object)
#define TCTTTTTO(type,type2,type3,type4,type5,text,object)
#define TCTTTTTTO(type,type2,type3,type4,type5,type6,text,object)
#define TCTTTTTTTO(type,type2,type3,type4,type5,type6,type7,text,object)
#define TCTTTTTTTTO(type,type2,type3,type4,type5,type6,type7,type8,text,object)
#else
#define SYMCHECK(a,b) if (a) { erg += error(b); goto endr_ende; }
#define COP(b,a) if (((void *)a) == NULL)  { erg += null_object(b); goto endr_ende;}
#define CTO(type,text,object) \
        if (type == INTTYPE);\
        else { \
            COP(text,object);\
            if (type == ANYTYPE);\
            else if(S_O_K(object) != type) {\
                erg += wrong_type_oneparameter(text,object);\
                goto endr_ende;\
               }\
            }

#define CTTO(type,type2,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&(S_O_K(object) != type2)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende;\
        }

#define CTTTO(type,type2,type3,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&(S_O_K(object) != type2)\
                &&(S_O_K(object) != type3)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende; }
#define CTTTTO(type,type2,type3,type4,text,object) \
        COP(text,object);\
        if(               (S_O_K(object) != type)\
                                &&(S_O_K(object) != type2)\
                &&(S_O_K(object) != type3) \
                                &&(S_O_K(object) != type4)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende; }
#define CTTTTTO(type,type2,type3,type4,type5,text,object) \
        COP(text,object);\
        if(       (S_O_K(object) != type)\
                &&(S_O_K(object) != type2)\
                &&(S_O_K(object) != type3)\
                &&(S_O_K(object) != type4)\
                &&(S_O_K(object) != type5)        ) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende; }
#define CTTTTTTO(type,type2,type3,type4,type5,type6,text,object) \
        COP(text,object);\
        if(       (S_O_K(object) != type)\
                &&(S_O_K(object) != type2)\
                &&(S_O_K(object) != type3)\
                &&(S_O_K(object) != type4)\
                &&(S_O_K(object) != type5)\
                &&(S_O_K(object) != type6)        ) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende; }

#define CTTTTTTTO(type,type2,type3,type4,type5,type6,type7,text,object) \
        COP(text,object);\
        if(       (S_O_K(object) != type)\
                &&(S_O_K(object) != type2)\
                &&(S_O_K(object) != type3)\
                &&(S_O_K(object) != type4)\
                &&(S_O_K(object) != type5)\
                &&(S_O_K(object) != type7)\
                &&(S_O_K(object) != type6)        ) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende; }

#define CTTTTTTTTO(type,type2,type3,type4,type5,type6,type7,type8,text,object) \
        COP(text,object);\
        if(       (S_O_K(object) != type)\
                &&(S_O_K(object) != type2)\
                &&(S_O_K(object) != type3)\
                &&(S_O_K(object) != type4)\
                &&(S_O_K(object) != type5)\
                &&(S_O_K(object) != type7)\
                &&(S_O_K(object) != type8)\
                &&(S_O_K(object) != type6)        ) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende; }



#define TCTO(type,text,object) \
       CTO(type,text,object);\
       if (type == INTTYPE)  printf("%s= %ld\n",text,(INT)object);\
       else { printf("%s=",text);println(object);}

#define TCTTO(type,type2,text,object) \
       CTTO(type,type2,text,object);\
       printf("%s=",text);println(object);

#define TCTTTO(type,type2,type3,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&\
            (S_O_K(object) != type3)&&(S_O_K(object) != type2)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende;\
       }\
       printf("%s=",text);println(object);

#define TCTTTTO(type,type2,type3,type4,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&(S_O_K(object) != type4)&&\
            (S_O_K(object) != type3)&&(S_O_K(object) != type2)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende;\
       }\
       printf("%s=",text);println(object);

#define TCTTTTTO(type,type2,type3,type4,type5,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&(S_O_K(object) != type4)&&\
            (S_O_K(object) != type5)&&\
            (S_O_K(object) != type3)&&(S_O_K(object) != type2)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende;\
       }\
       printf("%s=",text);println(object);

#define TCTTTTTTO(type,type2,type3,type4,type5,type6,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&(S_O_K(object) != type4)&&\
            (S_O_K(object) != type5)&&(S_O_K(object) != type6)&&\
            (S_O_K(object) != type3)&&(S_O_K(object) != type2)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende;\
       }\
       printf("%s=",text);println(object);

#define TCTTTTTTTO(type,type2,type3,type4,type5,type6,type7,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&(S_O_K(object) != type4)&&\
            (S_O_K(object) != type5)&&(S_O_K(object) != type6)&&\
            (S_O_K(object) != type7)&&\
            (S_O_K(object) != type3)&&(S_O_K(object) != type2)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende;\
       }\
       printf("%s=",text);println(object);

#define TCTTTTTTTTO(type,type2,type3,type4,type5,type6,type7,type8,text,object) \
        COP(text,object);\
        if((S_O_K(object) != type)&&(S_O_K(object) != type4)&&\
            (S_O_K(object) != type5)&&(S_O_K(object) != type6)&&\
            (S_O_K(object) != type7)&&(S_O_K(object) != type8)&&\
            (S_O_K(object) != type3)&&(S_O_K(object) != type2)) {\
            erg += wrong_type_oneparameter(text,object);\
            goto endr_ende;\
       }\
       printf("%s=",text);println(object);


#endif
#define C2R(a,b,t,c) if (check_result_2(a,b,t,c) \
    != NORESULT) goto strlabel
#define S2R(a,b,t,c) erg += store_result_2(a,b,t,c);strlabel:
#define C3R(a,b,d,t,c) if (check_result_3(a,b,d,t,c)\
     != NORESULT) goto strlabel
#define S3R(a,b,d,t,c)  erg += store_result_3(a,b,d,t,c);strlabel:
#define C5R(a,b,d,e,f,t,c) if (check_result_5(a,b,d,e,f,t,c)\
     != NORESULT) goto strlabel
#define S5R(a,b,d,e,f,t,c)  erg += store_result_5(a,b,d,e,f,t,c);strlabel:
#define C0R(t,c) if (check_result_0(t,c)\
     != NORESULT) goto strlabel
#define S0R(t,c) erg += store_result_0(t,c);strlabel:
#define C1R(a,t,c) if (check_result_1(a,t,c)\
     != NORESULT) goto strlabel
#define S1R(a,t,c) erg += store_result_1(a,t,c);strlabel:
#define WTO(text,b) erg += wrong_type_oneparameter(text,b)
#define WTT(text,b,c) erg += wrong_type_twoparameter(text,b,c)
#define EDC(b) error_during_computation_code(b,erg)
#define ENDR(a) endr_ende: if (erg != OK) EDC(a); return erg;
#define ENDO(a) endr_ende: if (erg != OK) EDC(a); return (OP) NULL;
#define ENDTYP(a,t) endr_ende: if (erg != OK) EDC(a); return (t) NULL;
#define NYI(b) not_yet_implemented(b)
#define NOP(a) null_object(a)
#define EOP(text,object) \
        COP(text,object);\
        if (EMPTYP(object)) {erg += empty_object(text); goto endr_ende;}
#define FATALERROR(text) fatal_error(text)
#define CH2D(a,b) if ((a) == (b) ) { erg += equal_2_error(); goto endr_ende;}

/*
#define CE2(a,b,f) if (check_equal_2(a,b,f,&erg) == EQUAL) goto endr_ende;
*/

#define CE2A(a,b,f) if (check_equal_2a(a,b,f,&erg) == EQUAL) goto endr_ende;

#define CE3(a,b,c,f) \
    if ((a==c) && (b == c))\
        {\
        OP checkequal3_d = callocobject();\
        *checkequal3_d = *c;\
        C_O_K(c,EMPTY);\
        erg += (*f)(checkequal3_d,checkequal3_d,c);\
        erg += freeall(checkequal3_d);\
        goto endr_ende;\
        }\
    else if (a==c)\
        {\
        OP checkequal3_d = callocobject();\
        *checkequal3_d = *c;\
        C_O_K(c,EMPTY);\
        erg += (*f)(checkequal3_d,b,c);\
        erg += freeall(checkequal3_d);\
        goto endr_ende;\
        }\
    else if (b==c)\
        {\
        OP checkequal3_d = callocobject();\
        *checkequal3_d = *c;\
        C_O_K(c,EMPTY);\
        erg += (*f)(a,checkequal3_d,c);\
        erg += freeall(checkequal3_d);\
        goto endr_ende;\
        }\
    else\
        {\
        if (c != NULL)\
            FREESELF(c);\
        }

#define CE4(a,b,c,d,f) \
    if (check_equal_4(a,b,c,d,\
    f,&erg) == EQUAL) goto endr_ende;
#define CE5(a,b,c,d,e,f) \
    if (check_equal_5(a,b,c,d,e,\
    f,&erg) == EQUAL) goto endr_ende;

#define GET_BIT_I(a,i) (((*a) >> (i))&1)
/*
#define GET_BV_I(a,i) ((*(((unsigned char *)S_V_S(a) ) + ((i)>>3))  >> ((i)&7))&1)
*/
#define GET_BV_I(a,i) ((*(((unsigned char *)S_V_S(a) ) + ((i)>>3))  >> (7-((i)&7))&1))
#define GET_BV_B(a,i) (*(((unsigned char *)S_V_S(a) ) + (i)>>3))
/*
#define UNSET_BV_I(a,i) (*((unsigned char *)S_V_S(a)  + ((i)>>3))  &= (~(1 << ((i)&7))))
*/
#define UNSET_BV_I(a,i) (\
    *((unsigned char *)S_V_S(a)  + ((i)>>3))  \
    &= (~(1 << (7-((i)&7))       )           )\
    )
#define UNSET_BIT_I(a,i) ((*(a))  &= (~(1 << (i))))
/*
#define SET_BV_I(a,i) ( *((unsigned char *)S_V_S(a)  + ((i)>>3))  |= (1 << ((i)&7)))
*/
#define SET_BV_I(a,i) ( *((unsigned char *)S_V_S(a)  + ((i)>>3))  |= (128 >> ((i)&7)))
#define SET_BIT_I(a,i) ( *(a)  |= (1 << (i)))
#define S_BV_LI(a) (S_V_LI(a) % 8 == 0 ? S_V_LI(a)>>3 : (S_V_LI(a)>>3) +1)


#define cons_two cons_zwei
#define ANFANG \
int main()\
{\
OP a,b,c,d,e,f,g,h;\
INT i,j,k;\
OP z=NULL;\
anfang();\
a=callocobject();\
b=callocobject();\
c=callocobject();\
d=callocobject();\
e=callocobject();\
f=callocobject();\
g=callocobject();\
h=callocobject();

#define ENDE \
freeall(a);\
freeall(b);\
freeall(c);\
freeall(d);\
freeall(e);\
freeall(f);\
freeall(g);\
freeall(h);\
ende();\
i=j=k=0;z=NULL;\
return 0;\
}
#define SYM_BEGIN ANFANG
#define SYM_END ENDE

#define inhalt(a,b) content(a,b)
#define inhalt_tableaux(a,b) content_tableaux(a,b)
#define inhalt_word(a,b) content_word(a,b)









/* FOR ALL macros */
#define FORALL_HASHTABLE_PRE060202(z,a,B)\
do { \
    OP FORALL_E;INT I,J; \
    for (FORALL_E=S_V_S(a), I=0; I<S_V_LI(a); I++,FORALL_E++) \
        {\
        if (S_O_K(FORALL_E)==VECTOR)\
           for (J=0;J<S_V_LI(FORALL_E) ;J++) \
               { \
               z = S_V_I(FORALL_E,J); \
               if (not EMPTYP(z)) {B;} \
               } \
        }\
} while(0)

#define FORALL_HASHTABLE(z,a,B)\
do { \
    OP FORALL_E;INT I,J; \
    for (FORALL_E=S_V_S(a), I=0; I<S_V_LI(a); I++,FORALL_E++) \
        {\
        if (S_O_K(FORALL_E)==VECTOR)\
           for (J=0,z=S_V_S(FORALL_E);J<S_V_LI(FORALL_E) ;J++,z++) \
               { if (not EMPTYP(z)) {B;} } \
        else if (S_I_I(FORALL_E) == -1) break;\
        else { I = S_I_I(FORALL_E)-1; FORALL_E=S_V_I(a,I); } \
        }\
} while(0)


#define FORALL(z,a,B)\
if (S_O_K(a) == HASHTABLE) \
do { \
    OP FORALL_E;INT I,J; \
    for (FORALL_E=S_V_S(a), I=0; I<S_V_LI(a); I++,FORALL_E++) \
        {\
        if (S_O_K(FORALL_E)==VECTOR)\
           for (J=0,z=S_V_S(FORALL_E);J<S_V_LI(FORALL_E) ;J++,z++) \
               { if (not EMPTYP(z)) {B;} }\
        else if (S_I_I(FORALL_E) == -1) break;\
        else { I = S_I_I(FORALL_E)-1; FORALL_E=S_V_I(a,I); } \
        }\
} while(0);\
else if (LISTP(a)) { \
            OP FORALL_E; \
            FORALL_E=a; \
            while (FORALL_E!=NULL) \
                { z = S_L_S(FORALL_E); \
                if (z != NULL) {B;}; \
                FORALL_E=S_L_N(FORALL_E); } \
            }\
else if (MATRIXP(a)) {\
    INT I=S_M_HI(a)*S_M_LI(a)-1;\
    for (z=S_M_S(a)+I;I>=0;I--,z--) {B;}\
}\
else if (VECTORP(a)) {\
    INT I=S_V_LI(a)-1;\
    for (z=S_V_S(a)+I;I>=0;I--,z--) {B;}\
}





/* new macros for main functions */

#define ABSOLUTE_INTEGER(a,b) \
    M_I_I((S_I_I(a) > 0 ? S_I_I(a) : -S_I_I(a)), b)

#define ADD_INTEGER(a,b,c) \
    if (S_O_K(b) == INTEGER) erg += add_integer_integer(a,b,c);\
    else if (S_O_K(b) == LONGINT) erg += add_longint_integer(b,a,c);\
    else if (S_O_K(b) == BRUCH) erg += add_bruch_integer(b,a,c);\
    else  erg += add_integer(a,b,c)

#define ADD(a,b,c) \
if (S_O_K(a) == INTEGER)    ADD_INTEGER(a,b,c); \
else if (S_O_K(a) == LONGINT)  \
    {\
    if (S_O_K(b) == INTEGER) erg += add_longint_integer(a,b,c);\
    else if (S_O_K(b) == LONGINT) erg += add_longint_longint(b,a,c);\
    else  erg += add_longint(a,b,c);\
    }\
else\
    erg += add(a,b,c)

#define ADD_APPLY_INTEGER(a,b)\
        if (S_O_K(b) == INTEGER) erg += add_apply_integer_integer(a,b);\
        else if (S_O_K(b) == LONGINT) erg += add_apply_integer_longint(a,b);\
        else if (S_O_K(b) == BRUCH) erg += add_apply_integer_bruch(a,b);\
        else  erg += add_apply_integer(a,b)

#define ADD_APPLY(a,b) \
    if (S_O_K(a) == INTEGER)  ADD_APPLY_INTEGER(a,b);\
    else if (S_O_K(a) == LONGINT)  \
        {\
        if (S_O_K(b) == INTEGER) erg += add_apply_longint_integer(a,b);\
        else if (S_O_K(b) == LONGINT) erg += add_apply_longint_longint(a,b);\
        else  erg += add_apply_longint(a,b);\
        }\
    else if (S_O_K(a) == BRUCH)  \
        {\
        if (S_O_K(b) == INTEGER) erg += add_apply_bruch_integer(a,b);\
        else if (S_O_K(b) == BRUCH) erg += add_apply_bruch_bruch(a,b);\
        else  erg += add_apply_bruch(a,b);\
        }\
    else if (S_O_K(a) == INTEGERVECTOR)  erg += add_apply_integervector(a,b);\
    else if (S_O_K(a) == POLYNOM)  erg += add_apply_polynom(a,b);\
    else if (S_O_K(a) == FF)  erg += add_apply_ff(a,b);\
    else\
        erg += add_apply(a,b)

#define ADDINVERS_INTEGER(a,b) M_I_I(-S_I_I(a),b)

#define ADDINVERS(a,b) \
if (S_O_K(a) == INTEGER) ADDINVERS_INTEGER(a,b);\
else if (S_O_K(a) == LONGINT) erg += addinvers_longint(a,b);\
else if (S_O_K(a) == BRUCH) erg += addinvers_bruch(a,b);\
else erg += addinvers(a,b)


#define ADDINVERS_APPLY_LONGINT(a) erg += (GANZNEG(S_O_S(a).ob_longint),OK)
#define ADDINVERS_APPLY_INTEGER(a) M_I_I( - S_I_I(a),a)

#define ADDINVERS_APPLY(a)\
if (S_O_K(a) == INTEGER) ADDINVERS_APPLY_INTEGER(a);\
else if (S_O_K(a) == LONGINT) ADDINVERS_APPLY_LONGINT(a);\
else if (S_O_K(a) == BRUCH) erg += addinvers_apply_bruch(a);\
else if (S_O_K(a) == MONOM) erg += addinvers_apply_monom(a);\
else erg += addinvers_apply(a)


#define BINOM_POSINTEGER_POSINTEGER(a,b,c)\
if (S_I_I(a) < BINOMLIMIT) {\
     M_I_I(  \
          (      (S_I_I(b)>S_I_I(a))     ? 0 :  binom_values [ S_I_I(a) ] [S_I_I(b)] )\
          ,c);\
    }\
else binom(a,b,c)


/*
#ifdef SYMMAGMA
#define CALLOCOBJECT() \
       ( (freeall_speicherposition >= 0L) ? \
         freeall_speicher[freeall_speicherposition--] : \
         mem_calloc(sizeof(struct object),1) )
#define CALLOCOBJECT() (OP)mem_calloc(sizeof(struct object),1)
#else
*/
/* freeall_speicherposition is next free object */
#define CALLOCOBJECT() \
       ( (freeall_speicherposition >= 0L) ? \
         freeall_speicher[freeall_speicherposition--] : \
         callocobject_fast() )
/*
#endif
*/

#define CALLOCOBJECT2(a,b) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT(); } while(0)

#define CALLOCOBJECT3(a,b,c) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT();c=CALLOCOBJECT(); } while(0)

#define CALLOCOBJECT4(a,b,c,d) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT();c=CALLOCOBJECT(); d=CALLOCOBJECT();} while(0)

#define CALLOCOBJECT5(a,b,c,d,e) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT();\
     c=CALLOCOBJECT(); d=CALLOCOBJECT();\
     e=CALLOCOBJECT();} while(0)

#define CALLOCOBJECT6(a,b,c,d,e,f) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT();\
     c=CALLOCOBJECT(); d=CALLOCOBJECT();\
     e=CALLOCOBJECT(); f=CALLOCOBJECT(); } while(0)

#define CALLOCOBJECT7(a,b,c,d,e,f,g) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT();\
     c=CALLOCOBJECT(); d=CALLOCOBJECT();\
     e=CALLOCOBJECT(); f=CALLOCOBJECT(); \
     g=CALLOCOBJECT();} while(0)

#define CALLOCOBJECT8(a,b,c,d,e,f,g,h) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT();\
     c=CALLOCOBJECT(); d=CALLOCOBJECT();\
     e=CALLOCOBJECT(); f=CALLOCOBJECT();h=CALLOCOBJECT(); \
     g=CALLOCOBJECT();} while(0)

#define CALLOCOBJECT9(a,b,c,d,e,f,g,h,i) \
do { a=CALLOCOBJECT(); b=CALLOCOBJECT(); h=CALLOCOBJECT(); i=CALLOCOBJECT();\
     c=CALLOCOBJECT(); d=CALLOCOBJECT();\
     e=CALLOCOBJECT(); f=CALLOCOBJECT(); \
     g=CALLOCOBJECT();} while(0)

#define COMP_INTEGER(a,b)\
    ( S_O_K(b) == INTEGER ? COMP_INTEGER_INTEGER(a,b) : \
       (S_O_K(b) == LONGINT ? - comp_longint_integer(b,a) : comp_integer(a,b) )\
    )

#define COMP(a,b)\
    ( S_O_K(a) == INTEGER ? COMP_INTEGER(a,b) :\
       (S_O_K(a) == LONGINT ? comp_longint(a,b) : \
          (S_O_K(a) == INTEGERMATRIX ? comp_integermatrix(a,b) : \
             comp(a,b) \
          )\
       )\
    )

#define COPY(a,b) \
    ( (S_O_K(a) == INTEGER) ? M_I_I(S_I_I(a),b):\
      ( (S_O_K(a) == LONGINT) ? copy_longint(a,b): \
        ( (S_O_K(a) == BRUCH) ? copy_bruch(a,b): \
          ( (S_O_K(a) == MONOM) ?  copy_monom(a,b): \
            ( (S_O_K(a) == PARTITION) ? copy_partition(a,b): \
              ( (S_O_K(a) == HASHTABLE)? copy_hashtable(a,b): \
                ( (S_O_K(a) == MATRIX)? copy_matrix(a,b): \
                  ( (S_O_K(a) == INTEGERMATRIX)? copy_integermatrix(a,b): \
                      copy(a,b)\
                  )\
                )\
              )\
            )\
          )\
        )\
      )\
    )

#define CLEVER_COPY_INTEGER(a,b) \
do { switch(S_O_K(b)) {\
        case INTEGER:\
        case EMPTY:\
             M_I_I(S_I_I(a),b);\
             break;\
        default:\
             FREESELF(b);\
             M_I_I(S_I_I(a),b);\
             break;\
        }\
} while (0)

#define CLEVER_COPY_FF(a,b) \
do { switch(S_O_K(b)) {\
        case INTEGER:\
             C_O_K(b,EMPTY);\
        case EMPTY:\
             erg += copy_ff(a,b);\
             break;\
        case FF:\
             {\
             INT *ap,*bp,i;\
             COPY(S_FF_C(a),S_FF_C(b));\
             COPY(S_FF_ME(a),S_FF_ME(b));\
             ap =S_FF_IP(a);bp=S_FF_IP(b);\
             if (*ap != *bp)\
                 bp = (INT*) SYM_realloc(bp,(S_FF_DI(a)+1)*sizeof(INT));\
             for(i=0;i<=ap[0];i++) bp[i]=ap[i];\
             C_FF_IP(b,bp);\
             };\
             break;   \
        default:\
             FREESELF(b);\
             erg += copy_ff(a,b);\
             break;\
        }\
} while (0)
#define CLEVER_COPY_BRUCH(a,b) \
do { switch(S_O_K(b)) {\
        case INTEGER:\
             C_O_K(b,EMPTY);\
        case EMPTY:\
             erg += copy_bruch(a,b);\
             break;\
        case BRUCH:\
             FREESELF(S_B_O(b));\
             FREESELF(S_B_U(b));\
             COPY(S_B_O(a),S_B_O(b));\
             COPY(S_B_U(a),S_B_U(b));\
             C_B_I(b,S_B_I(a));\
             break;\
        default:\
             FREESELF(b);\
             erg += copy_bruch(a,b);\
             break;\
        }\
} while (0)

#define CLEVER_COPY_LONGINT(a,b)\
do { switch(S_O_K(b)) {\
        case INTEGER:\
             C_O_K(b,EMPTY);\
        case EMPTY:\
             erg += copy_longint(a,b);\
             break;\
        default:\
             FREESELF(b);\
             erg += copy_longint(a,b);\
             break;\
        }\
} while(0)

#define CLEVER_COPY_PARTITION(a,b)\
do { switch(S_O_K(b)) {\
        case INTEGER:\
             C_O_K(b,EMPTY);\
        case EMPTY:\
             erg += copy_partition(a,b);\
             break;\
        case PARTITION:\
             FREESELF(S_PA_S(b));\
             if (S_O_K(S_PA_S(a))==INTEGERVECTOR) copy_integervector(S_PA_S(a),S_PA_S(b));\
             else COPY(S_PA_S(a),S_PA_S(b));\
             C_PA_K(b, S_PA_K(a));\
             C_PA_HASH(b, S_PA_HASH(a));\
             break;\
        default:\
             FREESELF(b);\
             erg += copy_partition(a,b);\
             break;\
        }\
} while(0)



#define CLEVER_COPY(a,b) \
do { if (S_O_K(a) == INTEGER) CLEVER_COPY_INTEGER(a,b);\
else if (S_O_K(a) == LONGINT) CLEVER_COPY_LONGINT(a,b);\
else if (S_O_K(a) == BRUCH) CLEVER_COPY_BRUCH(a,b);\
else if (S_O_K(a) == PARTITION) CLEVER_COPY_PARTITION(a,b);\
else if (S_O_K(a) == FF) CLEVER_COPY_FF(a,b);\
else { erg += copy(a,b); }\
} while(0)



#define DEC(a)\
if (S_O_K(a) == INTEGER) DEC_INTEGER(a);\
else if (S_O_K(a) == LONGINT) erg+= dec_longint(a);\
else dec(a)

#define EINSP_LONGINT(a)\
    (\
    ((S_O_S(a).ob_longint) ->floc ->w0 == 1) \
    &&\
    ((S_O_S(a).ob_longint) ->floc ->w1 == 0)\
    &&\
    ((S_O_S(a).ob_longint) ->floc ->w2 == 0)\
    &&\
    ((S_O_S(a).ob_longint) ->signum == 1)\
    &&\
    ((S_O_S(a).ob_longint) ->laenge == 1)\
    )


#define EINSP(a)\
    ( S_O_K(a) == INTEGER ? EINSP_INTEGER(a): \
       ( S_O_K(a) == LONGINT ? EINSP_LONGINT(a) : \
          ( S_O_K(a) == BRUCH ?  EINSP_BRUCH(a) : \
            ( einsp(a) ) \
          ) \
       ) \
    )

#define EVEN_INTEGER(a) (S_I_I(a) % 2 == 0)
#define EVEN_LONGINT(a) \
     (\
     ((S_O_S(a).ob_longint) ->signum == 0) \
     ||\
     ((S_O_S(a).ob_longint) ->floc ->w0 % 2 == 0)\
     )

#define EVEN(a) \
    ( S_O_K(a) == INTEGER ? (EVEN_INTEGER(a)) : \
       ( S_O_K(a) == LONGINT ? EVEN_LONGINT(a) : even(a) )\
    )

#define FREE_EMPTY_OBJECT(a)\
    do {\
    CTO(EMPTY,"FREE_EMPTY_OBJECT(1)",a);\
    if (freeall_speichersize+SPEICHERSIZE <freeall_speichersize_max)\
        {\
        if (freeall_speicherposition+1  == freeall_speichersize) \
            {\
            freeall_speicher = (OP *) \
                SYM_realloc(freeall_speicher,\
                            (freeall_speichersize+SPEICHERSIZE)*sizeof(OP));\
            if (freeall_speicher == NULL) {\
                erg += error("no more memory in freeall");\
                goto endr_ende;\
                }\
            freeall_speichersize = freeall_speichersize+SPEICHERSIZE;\
            }\
        freeall_speicher[++freeall_speicherposition] = a;\
        }\
    else SYM_FREE(a);\
    } while(0)

#define FREEALL_INTEGER(a) do { C_O_K(a,EMPTY); FREE_EMPTY_OBJECT(a); } while(0)
#define FREESELF_INTEGERVECTOR(a) \
    do {\
        extern INT freevectorstruct();\
        if (S_V_LI(a) == (INT)1) FREEALL_INTEGER(S_V_S(a));\
        else if (S_V_LI(a) > (INT)0) SYM_free(S_V_S(a));\
        FREEALL_INTEGER(S_V_L(a));\
        freevectorstruct(S_O_S(a).ob_vector);\
        C_O_K(a,EMPTY);\
    } while(0)
#define FREEALL_INTEGERVECTOR(a) do { \
    FREESELF_INTEGERVECTOR(a); FREE_EMPTY_OBJECT(a); } while(0)


#define FREESELF(a) \
    if (S_O_K(a) == EMPTY);\
    else if (S_O_K(a) == INTEGER) C_O_K(a,EMPTY);\
    else if (S_O_K(a) == LONGINT)  erg += freeself_longint(a);  \
    else if (S_O_K(a) == BRUCH)  erg += freeself_bruch(a);  \
    else if (S_O_K(a) == PARTITION)  erg += freeself_partition(a);  \
    else if (S_O_K(a) == MATRIX)  erg += freeself_matrix(a);  \
    else if (S_O_K(a) == INTEGERMATRIX)  erg += freeself_integermatrix(a);  \
    else if (S_O_K(a) == MONOM)  FREESELF_MONOM(a);  \
    else if (S_O_K(a) == INTEGERVECTOR)  FREESELF_INTEGERVECTOR(a);  \
    else if (S_O_K(a) == VECTOR)  erg += freeself_vector(a);  \
    else if (S_O_K(a) == HASHTABLE)  erg += freeself_hashtable(a);  \
    else if (LISTP(a))  erg += freeself_list(a);  \
    else if (S_O_K(a) == PERMUTATION)  erg += freeself_permutation(a);  \
    else if (S_O_K(a) == SKEWPARTITION)  erg += freeself_skewpartition(a);  \
    else if (S_O_K(a) == FF)  erg += freeself_ff(a);  \
    else erg += freeself(a)

#define FREESELF2(a,b) do { FREESELF(a); FREESELF(b); } while(0)
#define FREESELF2(a,b) do { FREESELF(a); FREESELF(b); } while(0)
#define FREESELF3(a,b,c) do { FREESELF(a); FREESELF(b); FREESELF(c); } while(0)
#define FREESELF4(a,b,c,d) do { FREESELF2(a,b); FREESELF2(c,d); } while(0)
#define FREESELF5(a,b,c,d,e) do { FREESELF2(a,b); FREESELF3(c,d,e); } while(0)
#define FREESELF6(a,b,c,d,e,f) do { FREESELF3(a,b,f); FREESELF3(c,d,e); } while(0)
#define FREESELF7(a,b,c,d,e,f,g) do { FREESELF4(a,b,f,g); FREESELF3(c,d,e); } while(0)

#define FREEALL(a) do { FREESELF(a); FREE_EMPTY_OBJECT(a); } while(0)
#define FREEALL2(a,b) do { FREEALL(a); FREEALL(b); } while(0)
#define FREEALL3(a,b,c) do { FREEALL(a); FREEALL(b); FREEALL(c); } while(0)
#define FREEALL4(a,b,c,d) do { FREEALL(a); FREEALL(b); FREEALL(c);FREEALL(d); } while(0)
#define FREEALL5(a,b,c,d,e) do { FREEALL2(a,b); FREEALL3(c,d,e); } while(0)
#define FREEALL6(a,b,c,d,e,f) do { FREEALL3(a,b,f); FREEALL3(c,d,e); } while(0)
#define FREEALL7(a,b,c,d,e,f,g) do { FREEALL4(a,b,f,g); FREEALL3(c,d,e); } while(0)
#define FREEALL8(a,b,c,d,e,f,g,h) do { FREEALL4(a,b,f,g); FREEALL4(c,d,e,h); } while(0)
#define FREEALL9(a,b,c,d,e,f,g,h,i) do { FREEALL4(a,b,f,g); FREEALL5(c,d,e,h,i); } while(0)

#define GANZDIV_INTEGER(a,b,c) \
    if (S_O_K(b) == INTEGER) M_I_I(S_I_I(a)/S_I_I(b),c);\
    else if (S_O_K(b) == LONGINT) erg += ganzdiv_integer_longint(a,b,c);\
    else  erg += ganzdiv_integer(a,b,c)

#define GANZDIV_LONGINT(a,b,c) \
    if (S_O_K(b) == INTEGER) erg += ganzdiv_longint_integer(a,b,c);\
    else if (S_O_K(b) == LONGINT) erg += ganzdiv_longint_longint(a,b,c);\
    else  erg += ganzdiv_longint(a,b,c)

#define GANZDIV(a,b,c) \
if (S_O_K(a) == INTEGER)    GANZDIV_INTEGER(a,b,c); \
else if (S_O_K(a) == LONGINT)  GANZDIV_LONGINT(a,b,c); \
else\
    erg += ganzdiv(a,b,c)


#define GANZDIV_APPLY_INTEGER(a,b) \
    if (S_O_K(b) == INTEGER) M_I_I(S_I_I(a)/S_I_I(b),a);\
    else ganzdiv_apply_integer(a,b)

#define GANZDIV_APPLY_LONGINT(a,b) \
    if (S_O_K(b) == INTEGER) erg += ganzdiv_apply_longint_integer(a,b);\
    else if (S_O_K(b) == LONGINT) erg += ganzdiv_apply_longint_longint(a,b);\
    else ganzdiv_apply_longint(a,b)

#define GANZDIV_APPLY(a,b) \
if (S_O_K(a) == INTEGER) GANZDIV_APPLY_INTEGER(a,b);\
else if (S_O_K(a) == LONGINT) GANZDIV_APPLY_LONGINT(a,b);\
else\
    erg += ganzdiv_apply(a,b);


#define GGT_INTEGER(a,b,c) \
    if (S_O_K(b) == INTEGER) erg += ggt_integer_integer(a,b,c);\
    else if (S_O_K(b) == LONGINT) erg += ggt_integer_longint(a,b,c);\
    else  erg += ggt_integer(a,b,c)

#define GGT_LONGINT(a,b,c) \
    if (S_O_K(b) == INTEGER) erg += ggt_integer_longint(b,a,c);\
    else if (S_O_K(b) == LONGINT) erg += ggt_longint_longint(a,b,c);\
    else  erg += ggt_longint(a,b,c)

#define GGT(a,b,c) \
if (S_O_K(a) == INTEGER)    GGT_INTEGER(a,b,c); \
else if (S_O_K(a) == LONGINT)  GGT_LONGINT(a,b,c); \
else\
    erg += ggt(a,b,c)

#define HALF_APPLY_INTEGER(a) M_I_I((S_I_I(a) >> 1),a)

#define HALF_APPLY(a)\
if (S_O_K(a) == INTEGER)\
    HALF_APPLY_INTEGER(a);\
else if (S_O_K(a) == LONGINT)\
    half_apply_longint(a);\
else\
    erg += half_apply(a)\


#define HASH_INTEGERVECTOR(a,res)\
    if (S_V_LI(a) == 0) res=4711;\
    else {\
        INT hash_integer_vector_i;\
        res = S_V_II(a,0);\
        for (hash_integer_vector_i=1;hash_integer_vector_i<S_V_LI(a);hash_integer_vector_i++)\
            {\
            res *= 4711;\
            res += S_V_II(a,hash_integer_vector_i);\
            }\
        }

#define HASH_MONOMPARTITION(a) \
    ( S_PA_HASH(S_MO_S(a)) == -1 ? hash_partition(S_MO_S(a)) : S_PA_HASH(S_MO_S(a))  )
#define HASH_PARTITION(a) \
    ( S_PA_HASH(a) == -1 ? hash_partition(a) : S_PA_HASH(a)  )
#define HASH_MONOM(a) \
    ( S_O_K(S_MO_S(a)) == PARTITION ? HASH_PARTITION(S_MO_S(a)) : hash(S_MO_S(a)) )
#define HASH(a) \
    ( S_O_K(a) == MONOM ?  HASH_MONOM(a): \
      ( S_O_K(a) == INTEGER ? S_I_I(a) : hash(a) )\
    )

#define INC(a)\
if (S_O_K(a) == INTEGER) INC_INTEGER(a);\
else if (S_O_K(a) == LONGINT) erg+= inc_longint(a);\
else inc(a)

#define INSERT_BINTREE(a,b,eh,cf) insert_bintree(a,b,eh,cf)
#define INSERT_HASHTABLE(a,b,eh,cf,hf) \
    if (S_O_K(a) == HASHTABLE) insert_hashtable_hashtable(a,b,eh,cf,hf);\
    else if (S_O_K(a) == SCHUR) insert_schur_hashtable(a,b,eh,cf,hf);\
    else if (S_O_K(a) == MONOMIAL) insert_monomial_hashtable(a,b,eh,cf,hf);\
    else if (S_O_K(a) == ELMSYM) insert_elmsym_hashtable(a,b,eh,cf,hf);\
    else if (S_O_K(a) == POWSYM) insert_powsym_hashtable(a,b,eh,cf,hf);\
    else if (S_O_K(a) == HOMSYM) insert_homsym_hashtable(a,b,eh,cf,hf);\
    else insert_scalar_hashtable(a,b,eh,cf,hf)
#define INSERT_LIST(a,b,eh,cf) \
    if (LISTP(a)) insert_list_list(a,b,eh,cf); \
    else insert_list(a,b,eh,cf)

#define INSERT(a,b,eh,cf) \
    if (S_O_K(b) == HASHTABLE) INSERT_HASHTABLE(a,b,eh,cf,hash); \
    else if (S_O_K(b) == BINTREE) INSERT_BINTREE(a,b,eh,cf); \
    else INSERT_LIST(a,b,eh,cf)

#define INSERT_SCHURMONOM_(m,c)\
    if (S_O_K(c) == HASHTABLE)\
       INSERT_HASHTABLE(m,c,add_koeff,eq_monomsymfunc,hash_monompartition);\
    else if (S_O_K(c) == SCHUR)\
       INSERT_LIST(m,c,add_koeff,comp_monomschur);\
    else if (S_O_K(c) == BINTREE)\
       INSERT_BINTREE(m,c,add_koeff,comp_monomschur);\
    else\
       WTO("INSERT_SCHURMONOM_(2)",c)

#define INSERT_POWSYMMONOM_(m,c)\
    if (S_O_K(c) == HASHTABLE)\
       INSERT_HASHTABLE(m,c,add_koeff,eq_monomsymfunc,hash_monompartition);\
    else if (S_O_K(c) == POWSYM)\
       INSERT_LIST(m,c,add_koeff,comp_monompowsym);\
    else if (S_O_K(c) == BINTREE)\
       INSERT_BINTREE(m,c,add_koeff,comp_monompowsym);\
    else\
       WTO("INSERT_POWSYMMONOM_(2)",c)

#define INSERT_ELMSYMMONOM_(m,c)\
    if (S_O_K(c) == HASHTABLE)\
       INSERT_HASHTABLE(m,c,add_koeff,eq_monomsymfunc,hash_monompartition);\
    else if (S_O_K(c) == ELMSYM)\
       INSERT_LIST(m,c,add_koeff,comp_monomelmsym);\
    else if (S_O_K(c) == BINTREE)\
       INSERT_BINTREE(m,c,add_koeff,comp_monomelmsym);\
    else\
       WTO("INSERT_ELMSYMMONOM_(2)",c)

#define INSERT_HOMSYMMONOM_(m,c)\
    if (S_O_K(c) == HASHTABLE)\
       INSERT_HASHTABLE(m,c,add_koeff,eq_monomsymfunc,hash_monompartition);\
    else if (S_O_K(c) == HOMSYM)\
       INSERT_LIST(m,c,add_koeff,comp_monomhomsym);\
    else if (S_O_K(c) == BINTREE)\
       INSERT_BINTREE(m,c,add_koeff,comp_monomhomsym);\
    else\
       WTO("INSERT_HOMSYMMONOM_(2)",c)





#define INTLOGPOS5_4(ai)   (( ai >= 10000L)     ? 5 : 4)
#define INTLOGPOS5_3(ai)   (( ai >= 1000L)      ? INTLOGPOS5_4(ai) : 3)
#define INTLOGPOS2_1(ai)   (( ai >= 10L)        ? 2 : 1)
#define INTLOGPOS5_1(ai)   (( ai >= 100L)       ? INTLOGPOS5_3(ai) : INTLOGPOS2_1(ai))
#define INTLOGPOS7_6(ai)   (( ai >= 1000000L)   ? 7 : 6)
#define INTLOGPOS10_9(ai)  (( ai >= 1000000000L)? 10 : 9)
#define INTLOGPOS10_8(ai)  (( ai >= 100000000L) ? INTLOGPOS10_9(ai) : 8)
#define INTLOGPOS10_6(ai)  (( ai >= 10000000L ) ? INTLOGPOS10_8(ai) : INTLOGPOS7_6(ai) )
#define INTLOGPOS(ai)      (( ai >= 100000L )   ? INTLOGPOS10_6(ai) : INTLOGPOS5_1(ai) )
#define INTLOG(a) ( (S_I_I(a)) >= 0 ? INTLOGPOS(S_I_I(a)) : INTLOGPOS(-S_I_I(a)) )


#define INVERS_INTEGER(a,b) do {\
          if (S_I_I(a)==1)  M_I_I(1,b); else \
          if (S_I_I(a)==-1)  M_I_I(-1,b); else \
          { b_ou_b(CALLOCOBJECT(),CALLOCOBJECT(),b); \
            M_I_I(1,S_B_O(b)); \
            M_I_I(S_I_I(a),S_B_U(b)); \
            C_B_I(b,GEKUERZT); \
          }\
        } while(0)

#define INVERS_LONGINT(a,b) \
      do {  b_ou_b(CALLOCOBJECT(),CALLOCOBJECT(),b); M_I_I(1,S_B_O(b)); \
       copy_longint(a,S_B_U(b)); C_B_I(b,GEKUERZT) ; } while(0)

#define INVERS_BRUCH(a,b) \
      do {  b_ou_b(CALLOCOBJECT(),CALLOCOBJECT(),b); COPY(S_B_O(a),S_B_U(b)); \
      COPY(S_B_U(a), S_B_O(b)); C_B_I(b,S_B_I(a)) ;} while(0)

#define INVERS(a,b) \
if (S_O_K(a) == INTEGER)  INVERS_INTEGER(a,b);\
else if (S_O_K(a) == LONGINT)  INVERS_LONGINT(a,b);\
else if (S_O_K(a) == BRUCH)   INVERS_BRUCH(a,b);\
else\
    erg += invers(a,b)

#define KUERZEN(a)\
if (S_O_K(S_B_O(a)) == INTEGER) {\
    if (S_O_K(S_B_U(b)) == INTEGER) erg += kuerzen_integer_integer(a);\
    else if (S_O_K(S_B_U(b)) == LONGINT) erg += kuerzen_integer_longint(a);\
    else erg += krz(a);\
}\
else if (S_O_K(S_B_O(a)) == LONGINT) {\
    if (S_O_K(S_B_U(b)) == INTEGER) erg += kuerzen_longint_integer(a);\
    else if (S_O_K(S_B_U(b)) == LONGINT) erg += kuerzen_longint_longint(a);\
    else erg += krz(a);\
}\
else krz(a)


#define MOD_APPLY_INTEGER(a,b) \
    if (S_O_K(b) == INTEGER) M_I_I(S_I_I(a) % S_I_I(b), a);\
    else if (S_O_K(b) == LONGINT) erg += mod_apply_integer_longint(a,b);\
    else mod_apply(a,b)

#define MOD_APPLY(a,b) \
if (S_O_K(a) == INTEGER)  MOD_APPLY_INTEGER(a,b);\
else if (S_O_K(a) == LONGINT)  erg += mod_apply_longint(a,b);\
else erg += mod_apply(a,b)


#define MULT_APPLY_INTEGER_INTEGER(a,b) \
    if ( NULLP_INTEGER(a) || NULLP_INTEGER(b) ) \
            { \
            M_I_I(0,b); \
            } \
    else if ( (INTLOG(a) + INTLOG(b)) > 9L )\
            {\
            erg += t_int_longint(b,b);\
            erg += mult_apply_integer_longint(a,b);\
            }\
    else M_I_I(S_I_I(a)*S_I_I(b),b)

#define MULT_APPLY_INTEGER(a,b) \
    if (S_O_K(b) == INTEGER) MULT_APPLY_INTEGER_INTEGER(a,b);\
    else if (S_O_K(b) == LONGINT) erg += mult_apply_integer_longint(a,b);\
    else if (S_O_K(b) == BRUCH) erg += mult_apply_integer_bruch(a,b);\
    else if (S_O_K(b) == MONOM) erg += mult_apply_integer_monom(a,b);\
    else if (POLYP(b)) erg += mult_apply_integer_polynom(a,b);\
    else if (S_O_K(b) == HASHTABLE) erg += mult_apply_integer_hashtable(a,b);\
    else  erg += mult_apply_integer(a,b)

#define MULT_APPLY_LONGINT(a,b)\
    if (S_O_K(b) == INTEGER) erg += mult_apply_longint_integer(a,b);\
    else if (S_O_K(b) == LONGINT) erg += mult_apply_longint_longint(a,b);\
    else if (S_O_K(b) == BRUCH) erg += mult_apply_longint_bruch(a,b);\
    else if (POLYP(b)) erg += mult_apply_longint_polynom(a,b);\
    else  erg += mult_apply_longint(a,b)

#define MULT_APPLY_BRUCH(a,b)\
    if (S_O_K(b) == INTEGER) erg += mult_apply_bruch_integer(a,b);\
    else if (S_O_K(b) == LONGINT) erg += mult_apply_bruch_longint(a,b);\
    else if (S_O_K(b) == BRUCH) erg += mult_apply_bruch_bruch(a,b);\
    else if (POLYP(b)) erg += mult_apply_bruch_polynom(a,b);\
    else if (S_O_K(b) == HASHTABLE) erg += mult_apply_bruch_hashtable(a,b);\
    else  erg += mult_apply_bruch(a,b)

#define MULT_APPLY(a,b) \
if (S_O_K(a) == INTEGER)  MULT_APPLY_INTEGER(a,b);\
else if (S_O_K(a) == LONGINT)  MULT_APPLY_LONGINT(a,b);\
else if (S_O_K(a) == BRUCH)   MULT_APPLY_BRUCH(a,b);\
else if (S_O_K(a) == POLYNOM)   erg += mult_apply_polynom(a,b);\
else if (S_O_K(a) == FF)   erg += mult_apply_ff(a,b);\
else\
    erg += mult_apply(a,b)

#define MULT_INTEGER_INTEGER(a,b,c) \
    if (INTLOG(a) + INTLOG(b) > 9) {\
        OP mii_c= CALLOCOBJECT();\
        erg += t_int_longint(a,mii_c);\
        erg += mult_longint_integer(mii_c,b,c);\
        FREEALL(mii_c);\
        }\
    else\
        M_I_I(S_I_I(a)*S_I_I(b),c)

#define MULT_INTEGER(a,b,c) \
    if (S_O_K(b) == INTEGER) MULT_INTEGER_INTEGER(a,b,c);\
    else if (S_O_K(b) == LONGINT) erg += mult_longint_integer(b,a,c);\
    else if (S_O_K(b) == BRUCH) erg += mult_bruch_integer(b,a,c);\
    else if (S_O_K(b) == CYCLOTOMIC) erg += mult_scalar_cyclo(a,b,c);\
    else  erg += mult_integer(a,b,c)

#define MULT_LONGINT(a,b,c) \
    if (S_O_K(b) == INTEGER) erg += mult_longint_integer(a,b,c);\
    else if (S_O_K(b) == LONGINT) erg += mult_longint_longint(b,a,c);\
    else if (S_O_K(b) == CYCLOTOMIC) erg += mult_scalar_cyclo(a,b,c);\
    else  erg += mult_longint(a,b,c)

#define MULT_BRUCH(a,b,c) \
    if (S_O_K(b) == INTEGER) erg += mult_bruch_integer(a,b,c);\
    else if (S_O_K(b) == LONGINT) erg += mult_bruch_longint(a,b,c);\
    else if (S_O_K(b) == BRUCH) erg += mult_bruch_bruch(a,b,c);\
    else if (S_O_K(b) == CYCLOTOMIC) erg += mult_scalar_cyclo(a,b,c);\
    else  erg += mult_bruch(a,b,c)

#define MULT(a,b,c) \
if (S_O_K(a) == INTEGER)    MULT_INTEGER(a,b,c); \
else if (S_O_K(a) == LONGINT)  MULT_LONGINT(a,b,c); \
else if (S_O_K(a) == BRUCH)  MULT_BRUCH(a,b,c);\
else if (S_O_K(a) == CYCLOTOMIC)  mult_cyclo(a,b,c);\
else if (S_O_K(a) == FF)  mult_ff(a,b,c);\
else if (S_O_K(a) == SQ_RADICAL)  mult_sqrad(a,b,c);\
else\
    erg += mult(a,b,c)

#define MULT_SCALAR_MONOMLIST(a,b,c)\
if ((NULLP(a))|| (NULLP(b))) erg += init(S_O_K(b),c);\
else erg += trans2formlist(a,b,c,mult)

#define CLEVER_MULT_INTEGER(a,b,c)\
do {\
FREESELF(c); MULT_INTEGER(a,b,c);\
} while(0)

#define CLEVER_MULT_LONGINT(a,b,c)\
do {\
FREESELF(c); MULT_LONGINT(a,b,c);\
} while(0)

#define CLEVER_MULT_BRUCH(a,b,c)\
if (S_O_K(b) == BRUCH) {\
    switch (S_O_K(c)) {\
        case INTEGER: C_O_K(c,EMPTY);\
        case EMPTY:\
        case BRUCH: erg += mult_bruch_bruch(a,b,c); break;\
        default: FREESELF(c); mult_bruch_bruch(a,b,c); break;\
        }\
    }\
else do { FREESELF(c); MULT_BRUCH(a,b,c); } while(0)

#define CLEVER_MULT_FF(a,b,c)\
if (S_O_K(b) == FF) {\
    switch (S_O_K(c)) {\
        case INTEGER: C_O_K(c,EMPTY);\
        case EMPTY:\
        case BRUCH: erg += mult_ff_ff(a,b,c); break;\
        default: FREESELF(c); mult_ff_ff(a,b,c); break;\
        }\
    }\
else do { FREESELF(c); erg += mult_ff(a,b,c); } while(0)

#define CLEVER_MULT(a,b,c) \
if (S_O_K(a) == INTEGER)    CLEVER_MULT_INTEGER(a,b,c); \
else if (S_O_K(a) == LONGINT)  CLEVER_MULT_LONGINT(a,b,c); \
else if (S_O_K(a) == BRUCH)  CLEVER_MULT_BRUCH(a,b,c);\
else if (S_O_K(a) == FF)  CLEVER_MULT_FF(a,b,c);\
else\
    do { FREESELF(c); MULT(a,b,c); } while(0)


#define NEGEINSP_LONGINT(a)\
    (\
    ((S_O_S(a).ob_longint) ->floc ->w0 == 1) \
    &&\
    ((S_O_S(a).ob_longint) ->floc ->w1 == 0)\
    &&\
    ((S_O_S(a).ob_longint) ->floc ->w2 == 0)\
    &&\
    ((S_O_S(a).ob_longint) ->signum == -1)\
    &&\
    ((S_O_S(a).ob_longint) ->laenge == 1)\
    )

#define NEGEINSP(a)\
    ( S_O_K(a) == INTEGER ? NEGEINSP_INTEGER(a): \
       ( S_O_K(a) == LONGINT ? NEGEINSP_LONGINT(a) : \
            ( negeinsp(a) ) \
       ) \
    )

#define NEGP_LONGINT(a)  (GANZSIGNUM(S_O_S(a).ob_longint) == (signed char)-1)

#define NEGP(a) \
    ( S_O_K(a) == INTEGER ? (NEGP_INTEGER(a)) : \
       ( S_O_K(a) == LONGINT ? NEGP_LONGINT(a) : negp(a) )\
    )

#define NEW_HASHTABLE(c)\
    do { c = CALLOCOBJECT(); erg += init_hashtable(c); } while(0)

#define CLEAR_HASHTABLE(c) /* removes all entries in a hashtable */ \
    do { OP z,zz;INT i,j;\
         for (i=0,z=S_V_S(c);i<S_V_LI(c);i++,z++)\
             { \
             if (not EMPTYP(z)) {\
                 for (j=0,zz=S_V_S(z);j<S_V_LI(z);j++,zz++) FREESELF(zz);\
                 /*FREESELF_INTEGERVECTOR(z);*/ C_I_I(S_V_L(z),1);\
                 }\
             /* C_I_I(z,-1);*/ \
             else if (S_I_I(z) == -1) break;\
             else { i = S_I_I(z)-1; z = S_V_I(c,i); }\
             }\
         M_I_I(0,S_V_I(c,S_V_LI(c)));} while(0)



#define NEW_INTEGER(m,i) do { m=CALLOCOBJECT(); M_I_I(i,m); } while(0)
#define NEW_TABLEAU(t,um) do { t=CALLOCOBJECT(); m_u_t(a,t); } while(0)
#define NEW_VECTOR(v,i)\
    do { v = CALLOCOBJECT(); erg += m_il_v(i,v); } while(0)
#define NEW_INTEGERVECTOR(v,i)\
    do { v = CALLOCOBJECT(); erg += m_il_integervector(i,v); } while(0)
#define NEW_HOMSYM(a)\
    do { a = CALLOCOBJECT(); erg += b_sn_l(NULL,NULL,a); \
         C_O_K(a,HOMSYM); } while(0)


#define NULLP_LONGINT(a)  (GANZSIGNUM(S_O_S(a).ob_longint) == (signed char) 0)
#define NULLP_HASHTABLE(a)  (S_V_II(a,S_V_LI(a)) == 0)
#define NULLP_BRUCH(a) \
   (S_O_K(S_B_O(a)) == INTEGER ? NULLP_INTEGER(S_B_O(a)) : \
       ( S_O_K(S_B_O(a)) == LONGINT  ?  NULLP_LONGINT(S_B_O(a)) : nullp(S_B_O(a))\
       )\
   )


#define NULLP(a)\
    ( S_O_K(a) == INTEGER ? NULLP_INTEGER(a) : \
       ( S_O_K(a) == LONGINT ? NULLP_LONGINT(a) : \
         ( S_O_K(a) == HASHTABLE ? NULLP_HASHTABLE(a) : \
          ( S_O_K(a) == BRUCH ?  nullp_bruch(a) : \
            ( S_O_K(a) == FF ? nullp_ff(a):  \
               ( POLYP(a) ? nullp_polynom(a): nullp(a) ) \
            )\
          ) \
         )\
       ) \
    )

#define ODD_INTEGER(a) (S_I_I(a) % 2 == 1)
#define ODD_LONGINT(a) \
     ((S_O_S(a).ob_longint) ->floc ->w0 % 2 == 1)

#define ODD(a) \
    ( S_O_K(a) == INTEGER ? (ODD_INTEGER(a)) : \
       ( S_O_K(a) == LONGINT ? ODD_LONGINT(a) : odd(a) )\
    )

#define PARTITION_WEIGHT(a,i) \
do { \
    OP z; \
    INT j; \
    for(j=S_PA_LI(a),i=0,z=S_V_S(S_PA_S(a));j>0;j--,z++) \
        i+=S_I_I(z); \
} while(0)

#define MAXPARTI(a) ((S_PA_LI(a) == 0) ? 0 : S_PA_II(a,S_PA_LI(a)-1) )

#define POSP_LONGINT(a)  (GANZSIGNUM(S_O_S(a).ob_longint) == (signed char)1)

#define POSP(a) \
    ( S_O_K(a) == INTEGER ? (POSP_INTEGER(a)) : \
       ( S_O_K(a) == LONGINT ? POSP_LONGINT(a) : posp(a) )\
    )

#define SWAP(a,b) do { \
        struct object swap_object;  \
        swap_object = *a; \
        *a = *b; \
        *b = swap_object; \
        } while(0)

#define CE2(a,b,f) \
    if (a==b)  {\
        OP checkequal2_c = CALLOCOBJECT();\
        *checkequal2_c = *b;\
        C_O_K(b,EMPTY);\
        erg += (*f)(checkequal2_c,b);\
        FREEALL(checkequal2_c);\
        goto endr_ende;\
        }\
    else   FREESELF(b)
/* used for transfunctions */
#define TCE2(a,b,f,typ) \
    if (a==b)  {\
        OP checkequal2_c = CALLOCOBJECT();\
        *checkequal2_c = *b;\
        C_O_K(b,EMPTY);\
        erg += (*f)(checkequal2_c,b);\
        FREEALL(checkequal2_c);\
        goto endr_ende;\
        }\
    else   if ( (S_O_K(b) != HASHTABLE) && (S_O_K(b) != typ) ) \
        FREESELF(b)

#define ADD_KOEFF(a,b) \
    ADD_APPLY(S_MO_K(a), S_MO_K(b));\
    if (NULLP(S_MO_K(b)))\
        FREESELF_MONOM(b)



#define M_FORALL_MONOMIALS_IN_AB(a,b,c,f,partf)\
{\
OP ff,z,y;\
ff = CALLOCOBJECT();\
    FORALL (y,a, {\
    FORALL (z,b, {\
            FREESELF(ff);\
            MULT(S_MO_K(z),S_MO_K(y),ff);\
            if (not EINSP(f))\
                {\
                MULT_APPLY(f,ff);\
                }\
            erg += (*partf)(S_MO_S(y),S_MO_S(z),c,ff);\
            } );\
            } );\
    FREEALL(ff);\
}

#define M2_FORALL_MONOMIALS_IN_AB(a,b,c,f,m,partf)\
{\
OP ff,z,y;\
ff = CALLOCOBJECT();\
    FORALL (y,a, {\
    FORALL (z,b, {\
            FREESELF(ff);\
            MULT(S_MO_K(z),S_MO_K(y),ff);\
            if (not EINSP(f))\
                {\
                MULT_APPLY(f,ff);\
                }\
            erg += (*partf)(S_MO_S(y),S_MO_S(z),c,ff,m);\
            } );\
            } );\
    FREEALL(ff);\
}

#define M3_FORALL_MONOMIALS_IN_AB(a,b,c,f,m,l,partf)\
{\
OP ff,z,y;\
ff = CALLOCOBJECT();\
    FORALL (y,a, {\
    FORALL (z,b, {\
            FREESELF(ff);\
            MULT(S_MO_K(z),S_MO_K(y),ff);\
            if (not EINSP(f))\
                {\
                MULT_APPLY(f,ff);\
                }\
            erg += (*partf)(S_MO_S(y),S_MO_S(z),c,ff,m,l);\
            } );\
            } );\
    FREEALL(ff);\
}



#define M_FORALL_MONOMIALS_IN_B(a,b,c,f,partf)\
{\
OP ff,z;\
if (not EINSP(f)) {\
    ff = CALLOCOBJECT();\
    FORALL (z,b, {\
            FREESELF(ff);\
            MULT(S_MO_K(z),f,ff);\
            erg += (*partf)(a,S_MO_S(z),c,ff);\
            } );\
    FREEALL(ff);\
    }\
else {\
    FORALL (z,b, {\
       erg += (*partf)(a,S_MO_S(z),c,S_MO_K(z));\
            } );\
    }\
}
#define M2_FORALL_MONOMIALS_IN_B(a,b,c,f,m,partf)\
{\
OP ff,z;\
if (not EINSP(f)) {\
    ff = CALLOCOBJECT();\
    FORALL (z,b, {\
            MULT(S_MO_K(z),f,ff);\
            erg += (*partf)(a,S_MO_S(z),c,ff,m);\
            FREESELF(ff);\
            } );\
    FREEALL(ff);\
    }\
else {\
    FORALL (z,b, {\
       erg += (*partf)(a,S_MO_S(z),c,S_MO_K(z),m);\
            } );\
    }\
}

#define M3_FORALL_MONOMIALS_IN_B(a,b,c,f,m,l,partf)\
{\
OP ff,z;\
if (not EINSP(f)) {\
    ff = CALLOCOBJECT();\
    FORALL (z,b, {\
            MULT(S_MO_K(z),f,ff);\
            erg += (*partf)(a,S_MO_S(z),c,ff,m,l);\
            FREESELF(ff);\
            } );\
    FREEALL(ff);\
    }\
else {\
    FORALL (z,b, {\
       erg += (*partf)(a,S_MO_S(z),c,S_MO_K(z),m,l);\
            } );\
    }\
}



#define M_FORALL_MONOMIALS_IN_A(a,b,c,f,partf)\
{\
OP ff,z;\
if (not EINSP(f)) {\
    ff = CALLOCOBJECT();\
    FORALL (z,a, {\
            MULT(S_MO_K(z),f,ff);\
            erg += (*partf)(S_MO_S(z),b,c,ff);\
            FREESELF(ff);\
            } );\
    FREEALL(ff);\
    }\
else {\
    FORALL (z,a, {\
       erg += (*partf)(S_MO_S(z),b,c,S_MO_K(z));\
            } );\
    }\
}
#define M2_FORALL_MONOMIALS_IN_A(a,b,c,f,m,partf)\
{\
OP ff,z;\
if (not EINSP(f)) {\
    ff = CALLOCOBJECT();\
    FORALL (z,a, {\
            MULT(S_MO_K(z),f,ff);\
            erg += (*partf)(S_MO_S(z),b,c,ff,m);\
            FREESELF(ff);\
            } );\
    FREEALL(ff);\
    }\
else {\
    FORALL (z,a, {\
       erg += (*partf)(S_MO_S(z),b,c,S_MO_K(z),m);\
            } );\
    }\
}

#define M3_FORALL_MONOMIALS_IN_A(a,b,c,f,m,l,partf)\
{\
OP ff,z;\
if (not EINSP(f)) {\
    ff = CALLOCOBJECT();\
    FORALL (z,a, {\
            MULT(S_MO_K(z),f,ff);\
            erg += (*partf)(S_MO_S(z),b,c,ff,m,l);\
            FREESELF(ff);\
            } );\
    FREEALL(ff);\
    }\
else {\
    FORALL (z,a, {\
       erg += (*partf)(S_MO_S(z),b,c,S_MO_K(z),m,l);\
            } );\
    }\
}



#define T_FORALL_MONOMIALS_IN_A(a,b,f,partf)\
{\
OP ff,z;\
if (not EINSP(f)) {\
    ff = CALLOCOBJECT();\
    FORALL (z,a, {\
            MULT(S_MO_K(z),f,ff);\
            erg += (*partf)(S_MO_S(z),b,ff);\
            FREESELF(ff);\
            } );\
    FREEALL(ff);\
    }\
else {\
    FORALL (z,a, {\
       erg += (*partf)(S_MO_S(z),b,S_MO_K(z));\
            } );\
    }\
}

#define _NULL_PARTITION_(b,c,f) \
do { OP m;\
CTO(PARTITION,"_NULL_PARTITION_(1)",b);\
m=CALLOCOBJECT(); \
erg += b_sk_mo(CALLOCOBJECT(),CALLOCOBJECT(),m);\
erg += copy_partition(b,S_MO_S(m));\
COPY(f,S_MO_K(m));\
if (S_O_K(c)==HASHTABLE)\
insert_scalar_hashtable(m,c,add_koeff,eq_monomsymfunc,hash_monompartition);\
else \
INSERT_LIST(m,c,add_koeff,comp_monommonomial);\
} while(0)

#define WEIGHT_HASHTABLE(a) S_V_II(a,S_V_LI(a))



#define NEQ(a,b) (! EQ(a,b))
#define GR(a,b) (COMP(a,b) > (INT)0)
#define GT(a,b) (COMP(a,b) > (INT)0)
#define GE(a,b) (COMP(a,b) >= (INT)0)
#define LE(a,b) (COMP(a,b) <= (INT)0)
#define LT(a,b) (COMP(a,b) < (INT)0)


/* speicher managment */

#define FREE_MEMMANAGER(t,s,i,g,c,v)\
do {\
    c--;\
    if ((i+1) == g) {\
       if (g+SPEICHERSIZE <freeall_speichersize_max)\
       {\
       if (g == 0) {\
           s = (t*) SYM_MALLOC(SPEICHERSIZE * sizeof(t));\
           SYMCHECK(s == NULL,"no memory");\
           g = SPEICHERSIZE;\
           }\
       else {\
           s = (t*) SYM_realloc (s, (g+SPEICHERSIZE) * sizeof(t));\
           SYMCHECK(s == NULL,"no memory");\
           g += SPEICHERSIZE;\
           }\
       s[++i] = (v);\
       }\
       else SYM_FREE(v);\
       }\
    else s[++i] = (v);\
} while(0)

#define CALLOC_MEMMANAGER(t,s,i,c,v)\
do {\
    c++;\
    if (i>=0) v = s[i--];\
    else v = (t*) SYM_MALLOC(sizeof(t));\
} while(0)

#define ANFANG_MEMMANAGER(s,i,g,c) \
s = NULL; i = -1; c = 0; g = 0

#define ENDE_MEMMANAGER(s,i,g,c,text) \
if (no_banner != TRUE) { SYMCHECK(c, text);} \
if (s != NULL) { \
    INT jj;\
    for (jj = 0; jj<=i;jj++) SYM_FREE(s[jj]);\
    SYM_FREE(s);\
    s = NULL;\
    }\
i = -1;\
g = 0


#define MACRO_H
#endif /* MACRO_H */
