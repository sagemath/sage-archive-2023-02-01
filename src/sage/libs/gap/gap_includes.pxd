# distutils: libraries = gap gmp m
###############################################################################
#       Copyright (C) 2009, William Stein <wstein@gmail.com>
#       Copyright (C) 2012, Volker Braun <vbraun.name@gmail.com>
#
#   Distributed under the terms of the GNU General Public License (GPL)
#   as published by the Free Software Foundation; either version 2 of
#   the License, or (at your option) any later version.
#                   http://www.gnu.org/licenses/
###############################################################################


cdef extern from "<gap/system.h>":
    ctypedef char Char "libGAP_Char"
    ctypedef int Int "libGAP_Int"
    ctypedef unsigned char UChar "libGAP_UChar"

cdef extern from "<gap/libgap.h>":
    void libgap_initialize(int argc, char** argv)
    ctypedef void(*libgap_gasman_callback_ptr)()
    void libgap_set_gasman_callback(libgap_gasman_callback_ptr callback)
    ctypedef void(*libgap_error_func_ptr)(char* msg)
    void libgap_set_error_handler(libgap_error_func_ptr error_handler)
    void libgap_call_error_handler()
    void libgap_finalize()
    void libgap_start_interaction(char* inputline)
    char* libgap_get_output()
    char* libgap_get_error()
    void libgap_finish_interaction()
    void libgap_mark_stack_bottom()
    void libgap_enter()
    void libgap_exit()

cdef extern from "<gap/code.h>":
    ctypedef unsigned int Stat "libGAP_Stat"
    ctypedef Stat* PtrBody "libGAP_PtrBody"

cdef extern from "<gap/gap.h>":
    ctypedef unsigned int UInt "libGAP_UInt"
    ctypedef void* ExecStatus "libGAP_ExecStatus"
    void ViewObjHandler "libGAP_ViewObjHandler"(void*)
    void InitializeGap "libGAP_InitializeGap"(int*, char** argv)
    void set_system_variables "libGAP_set_system_variables"(char**, char**)
    cdef UInt Last "libGAP_Last"
    cdef UInt Last2 "libGAP_Last2"
    cdef UInt Last3 "libGAP_Last3"
    cdef ExecStatus STATUS_END "libGAP_STATUS_END"
    cdef ExecStatus STATUS_RETURN_VAL "libGAP_STATUS_RETURN_VAL"
    cdef ExecStatus STATUS_RETURN_VOID "libGAP_STATUS_RETURN_VOID"
    cdef ExecStatus STATUS_TNM "libGAP_STATUS_TNM"
    cdef ExecStatus STATUS_QUIT "libGAP_STATUS_QUIT"
    cdef ExecStatus STATUS_EOF "libGAP_STATUS_EOF"
    cdef ExecStatus STATUS_ERROR "libGAP_STATUS_ERROR"
    cdef ExecStatus STATUS_QQUIT "libGAP_STATUS_QQUIT"

cdef extern from "<gap/objects.h>":
    ctypedef void* Obj "libGAP_Obj"
    bint IS_MUTABLE_OBJ "libGAP_IS_MUTABLE_OBJ"(Obj obj)
    bint IS_COPYABLE_OBJ "libGAP_IS_COPYABLE_OBJ"(Obj obj)
    Obj SHALLOW_COPY_OBJ "libGAP_SHALLOW_COPY_OBJ"(Obj obj)
    Obj CopyObj "libGAP_CopyObj"(Obj obj, int mut)

    bint IS_INTOBJ "libGAP_IS_INTOBJ"(Obj obj)
    Obj INTOBJ_INT "libGAP_INTOBJ_INT"(Int)
    Int INT_INTOBJ "libGAP_INT_INTOBJ"(Obj)
    UInt TNUM_OBJ "libGAP_TNUM_OBJ"(Obj obj)
    char* TNAM_OBJ "libGAP_TNAM_OBJ"(Obj obj)
    cdef int FIRST_REAL_TNUM "libGAP_FIRST_REAL_TNUM"
    cdef int FIRST_CONSTANT_TNUM "libGAP_FIRST_CONSTANT_TNUM"
    cdef int T_INT "libGAP_T_INT"
    cdef int T_INTPOS "libGAP_T_INTPOS"
    cdef int T_INTNEG "libGAP_T_INTNEG"
    cdef int T_RAT "libGAP_T_RAT"
    cdef int T_CYC "libGAP_T_CYC"
    cdef int T_FFE "libGAP_T_FFE"
    cdef int T_PERM2 "libGAP_T_PERM2"
    cdef int T_PERM4 "libGAP_T_PERM4"
    cdef int T_BOOL "libGAP_T_BOOL"
    cdef int T_CHAR "libGAP_T_CHAR"
    cdef int T_FUNCTION "libGAP_T_FUNCTION"
    cdef int T_FLAGS "libGAP_T_FLAGS"
    cdef int T_MACFLOAT "libGAP_T_MACFLOAT"
    cdef int T_RESERVED_BY_GAP "libGAP_T_RESERVED_BY_GAP"
    cdef int LAST_CONSTANT_TNUM "libGAP_LAST_CONSTANT_TNUM"
    cdef int IMMUTABLE "libGAP_IMMUTABLE"
    cdef int FIRST_IMM_MUT_TNUM "libGAP_FIRST_IMM_MUT_TNUM"
    cdef int FIRST_RECORD_TNUM "libGAP_FIRST_RECORD_TNUM"
    cdef int T_PREC "libGAP_T_PREC"
    cdef int LAST_RECORD_TNUM "libGAP_LAST_RECORD_TNUM"
    cdef int FIRST_LIST_TNUM "libGAP_FIRST_LIST_TNUM"
    cdef int FIRST_PLIST_TNUM "libGAP_FIRST_PLIST_TNUM"
    cdef int T_PLIST "libGAP_T_PLIST"
    cdef int T_PLIST_NDENSE "libGAP_T_PLIST_NDENSE"
    cdef int T_PLIST_DENSE "libGAP_T_PLIST_DENSE"
    cdef int T_PLIST_DENSE_NHOM "libGAP_T_PLIST_DENSE_NHOM"
    cdef int T_PLIST_DENSE_NHOM_SSORT "libGAP_T_PLIST_DENSE_NHOM_SSORT"
    cdef int T_PLIST_DENSE_NHOM_NSORT "libGAP_T_PLIST_DENSE_NHOM_NSORT"
    cdef int T_PLIST_EMPTY "libGAP_T_PLIST_EMPTY"
    cdef int T_PLIST_HOM "libGAP_T_PLIST_HOM"
    cdef int T_PLIST_HOM_NSORT "libGAP_T_PLIST_HOM_NSORT"
    cdef int T_PLIST_HOM_SSORT "libGAP_T_PLIST_HOM_SSORT"
    cdef int T_PLIST_TAB "libGAP_T_PLIST_TAB"
    cdef int T_PLIST_TAB_NSORT "libGAP_T_PLIST_TAB_NSORT"
    cdef int T_PLIST_TAB_SSORT "libGAP_T_PLIST_TAB_SSORT"
    cdef int T_PLIST_TAB_RECT "libGAP_T_PLIST_TAB_RECT"
    cdef int T_PLIST_TAB_RECT_NSORT "libGAP_T_PLIST_TAB_RECT_NSORT"
    cdef int T_PLIST_TAB_RECT_SSORT "libGAP_T_PLIST_TAB_RECT_SSORT"
    cdef int T_PLIST_CYC "libGAP_T_PLIST_CYC"
    cdef int T_PLIST_CYC_NSORT "libGAP_T_PLIST_CYC_NSORT"
    cdef int T_PLIST_CYC_SSORT "libGAP_T_PLIST_CYC_SSORT"
    cdef int T_PLIST_FFE "libGAP_T_PLIST_FFE"
    cdef int LAST_PLIST_TNUM "libGAP_LAST_PLIST_TNUM"
    cdef int T_RANGE_NSORT "libGAP_T_RANGE_NSORT"
    cdef int T_RANGE_SSORT "libGAP_T_RANGE_SSORT"
    cdef int T_BLIST "libGAP_T_BLIST"
    cdef int T_BLIST_NSORT "libGAP_T_BLIST_NSORT"
    cdef int T_BLIST_SSORT "libGAP_T_BLIST_SSORT"
    cdef int T_STRING "libGAP_T_STRING"
    cdef int T_STRING_NSORT "libGAP_T_STRING_NSORT"
    cdef int T_STRING_SSORT "libGAP_T_STRING_SSORT"
    cdef int LAST_LIST_TNUM "libGAP_LAST_LIST_TNUM"
    cdef int LAST_IMM_MUT_TNUM "libGAP_LAST_IMM_MUT_TNUM"
    cdef int FIRST_EXTERNAL_TNUM "libGAP_FIRST_EXTERNAL_TNUM"
    cdef int T_COMOBJ "libGAP_T_COMOBJ"
    cdef int T_POSOBJ "libGAP_T_POSOBJ"
    cdef int T_DATOBJ "libGAP_T_DATOBJ"
    cdef int T_WPOBJ "libGAP_T_WPOBJ"
    cdef int LAST_EXTERNAL_TNUM "libGAP_LAST_EXTERNAL_TNUM"
    cdef int LAST_REAL_TNUM "libGAP_LAST_REAL_TNUM"
    cdef int LAST_VIRTUAL_TNUM "libGAP_LAST_VIRTUAL_TNUM"
    cdef int FIRST_COPYING_TNUM "libGAP_FIRST_COPYING_TNUM"
    cdef int COPYING "libGAP_COPYING"
    cdef int LAST_COPYING_TNUM "libGAP_LAST_COPYING_TNUM"
    cdef int FIRST_TESTING_TNUM "libGAP_FIRST_TESTING_TNUM"
    cdef int TESTING "libGAP_TESTING"
    cdef int LAST_TESTING_TNUM "libGAP_LAST_TESTING_TNUM"

cdef extern from "<gap/read.h>":
    void* ReadEvalCommand "libGAP_ReadEvalCommand"(Obj context, UInt *dualSemicolon)
    void* ReadEvalFile "libGAP_ReadEvalFile"()
    void* ReadEvalResult "libGAP_ReadEvalResult"
    bint READ_ERROR "libGAP_READ_ERROR"()

cdef extern from "<gap/scanner.h>":
    void ClearError "libGAP_ClearError"()
    UInt NrError "libGAP_NrError"
    UInt Symbol "libGAP_Symbol"
    void GetSymbol "libGAP_GetSymbol"()
    void Match "libGAP_Match" (UInt symbol, char* msg, UInt skipto)
    int S_ILLEGAL "libGAP_S_ILLEGAL"
    int S_IDENT "libGAP_S_IDENT"
    int S_UNBIND "libGAP_S_UNBIND"
    int S_ISBOUND "libGAP_S_ISBOUND"
    int S_TRYNEXT "libGAP_S_TRYNEXT"
    int S_INFO "libGAP_S_INFO"
    int S_ASSERT "libGAP_S_ASSERT"
    int S_SAVEWS "libGAP_S_SAVEWS"
    int S_LOADWS "libGAP_S_LOADWS"
    int S_LBRACK "libGAP_S_LBRACK"
    int S_LBRACE "libGAP_S_LBRACE"
    int S_BLBRACK "libGAP_S_BLBRACK"
    int S_BLBRACE "libGAP_S_BLBRACE"
    int S_RBRACK "libGAP_S_RBRACK"
    int S_RBRACE "libGAP_S_RBRACE"
    int S_DOT "libGAP_S_DOT"
    int S_BDOT "libGAP_S_BDOT"
    int S_LPAREN "libGAP_S_LPAREN"
    int S_RPAREN "libGAP_S_RPAREN"
    int S_COMMA "libGAP_S_COMMA"
    int S_DOTDOT "libGAP_S_DOTDOT"
    int S_COLON "libGAP_S_COLON"
    int S_PARTIALINT "libGAP_S_PARTIALINT"
    int S_INT "libGAP_S_INT"
    int S_TRUE "libGAP_S_TRUE"
    int S_FALSE "libGAP_S_FALSE"
    int S_CHAR "libGAP_S_CHAR"
    int S_STRING "libGAP_S_STRING"
    int S_PARTIALSTRING "libGAP_S_PARTIALSTRING"
    int S_REC "libGAP_S_REC"
    int S_FUNCTION "libGAP_S_FUNCTION"
    int S_LOCAL "libGAP_S_LOCAL"
    int S_END "libGAP_S_END"
    int S_MAPTO "libGAP_S_MAPTO"
    int S_MULT "libGAP_S_MULT"
    int S_DIV "libGAP_S_DIV"
    int S_MOD "libGAP_S_MOD"
    int S_POW "libGAP_S_POW"
    int S_PLUS "libGAP_S_PLUS"
    int S_MINUS "libGAP_S_MINUS"
    int S_EQ "libGAP_S_EQ"
    int S_LT "libGAP_S_LT"
    int S_GT "libGAP_S_GT"
    int S_NE "libGAP_S_NE"
    int S_LE "libGAP_S_LE"
    int S_GE "libGAP_S_GE"
    int S_IN "libGAP_S_IN"
    int S_NOT "libGAP_S_NOT"
    int S_AND "libGAP_S_AND"
    int S_OR "libGAP_S_OR"
    int S_ASSIGN "libGAP_S_ASSIGN"
    int S_IF "libGAP_S_IF"
    int S_FOR "libGAP_S_FOR"
    int S_WHILE "libGAP_S_WHILE"
    int S_REPEAT "libGAP_S_REPEAT"
    int S_THEN "libGAP_S_THEN"
    int S_ELIF "libGAP_S_ELIF"
    int S_ELSE "libGAP_S_ELSE"
    int S_FI "libGAP_S_FI"
    int S_DO "libGAP_S_DO"
    int S_OD "libGAP_S_OD"
    int S_UNTIL "libGAP_S_UNTIL"
    int S_BREAK "libGAP_S_BREAK"
    int S_RETURN "libGAP_S_RETURN"
    int S_QUIT "libGAP_S_QUIT"
    int S_QQUIT "libGAP_S_QQUIT"
    int S_CONTINUE "libGAP_S_CONTINUE"
    int S_SEMICOLON "libGAP_S_SEMICOLON"
    int S_EOF "libGAP_S_EOF"

cdef extern from "<gap/gvars.h>":
    UInt GVarName "libGAP_GVarName"(char* name)
    void AssGVar "libGAP_AssGVar"(UInt gvar, Obj val)
    Obj VAL_GVAR "libGAP_VAL_GVAR"(UInt gvar)

cdef extern from "<gap/string.h>":
    char* CSTR_STRING "libGAP_CSTR_STRING"(Obj list)
    int GET_LEN_STRING "libGAP_GET_LEN_STRING"(Obj list)
    bint IS_STRING "libGAP_IS_STRING"(Obj obj)
    bint IsStringConv "libGAP_IsStringConv"(Obj obj)
    bint ConvString "libGAP_ConvString"(Obj obj)
    void C_NEW_STRING "libGAP_C_NEW_STRING"(Obj new_gap_string, int length, char* c_string)

cdef extern from "<gap/gasman.h>":
    void InitGlobalBag "libGAP_InitGlobalBag"(Obj* addr, char* cookie)
    Obj NewBag "libGAP_NewBag"(UInt type, UInt size)
    void CHANGED_BAG "libGAP_CHANGED_BAG"(Obj bag)
    void MARK_BAG "libGAP_MARK_BAG"(Obj bag)
    bint IS_MARKED_ALIVE "libGAP_IS_MARKED_ALIVE"(Obj bag)
    bint IS_MARKED_DEAD "libGAP_IS_MARKED_DEAD"(Obj bag)
    bint IS_MARKED_HALFDEAD "libGAP_IS_MARKED_HALFDEAD"(Obj bag)
    cdef UInt NrAllBags "libGAP_NrAllBags"
    cdef UInt SizeAllBags "libGAP_SizeAllBags"
    cdef UInt NrLiveBags "libGAP_NrLiveBags"
    cdef UInt SizeLiveBags "libGAP_SizeLiveBags"
    cdef UInt NrDeadBags "libGAP_NrDeadBags"
    cdef UInt SizeDeadBags "libGAP_SizeDeadBags"
    cdef UInt NrHalfDeadBags "libGAP_NrHalfDeadBags"
    UInt CollectBags "libGAP_CollectBags"(UInt size, UInt full)
    void CallbackForAllBags "libGAP_CallbackForAllBags"(void (*func)(Obj))
    char* TNAM_BAG "libGAP_TNAM_BAG"(Obj obj)
    UInt TNUM_BAG "libGAP_TNUM_BAG"(Obj)
    UInt SIZE_BAG "libGAP_SIZE_BAG"(Obj)
    void CheckMasterPointers "libGAP_CheckMasterPointers"()
    Obj* MptrBags "libGAP_MptrBags"
    Obj* YoungBags "libGAP_YoungBags"
    Obj* OldBags "libGAP_OldBags"
    Obj* AllocBags "libGAP_AllocBags"
    Obj* MarkedBags "libGAP_MarkedBags"
    Obj* ChangedBags "libGAP_ChangedBags"

# in gasman.c but not declared in gasman.h
cdef extern Obj* StopBags "libGAP_StopBags"
cdef extern Obj* EndBags "libGAP_EndBags"

cdef extern from "<gap/ariths.h>":
    Obj SUM "libGAP_SUM" (Obj, Obj)
    Obj DIFF "libGAP_DIFF"(Obj, Obj)
    Obj PROD "libGAP_PROD"(Obj, Obj)
    Obj QUO "libGAP_QUO"(Obj, Obj)
    Obj POW "libGAP_POW"(Obj, Obj)
    Obj MOD "libGAP_MOD"(Obj, Obj)
    Obj CALL_0ARGS "libGAP_CALL_0ARGS"(Obj f)              # 0 arguments
    Obj CALL_1ARGS "libGAP_CALL_1ARGS"(Obj f, Obj a1)      # 1 argument
    Obj CALL_2ARGS "libGAP_CALL_2ARGS"(Obj f, Obj a1, Obj a2)
    Obj CALL_3ARGS "libGAP_CALL_3ARGS"(Obj f, Obj a1, Obj a2, Obj a3)
    Obj CALL_4ARGS "libGAP_CALL_4ARGS"(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4)
    Obj CALL_5ARGS "libGAP_CALL_5ARGS"(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4, Obj a5)
    Obj CALL_6ARGS "libGAP_CALL_6ARGS"(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4, Obj a5, Obj a6)
    Obj CALL_XARGS "libGAP_CALL_XARGS"(Obj f, Obj args)   # more than 6 arguments
    bint EQ "libGAP_EQ"(Obj opL, Obj opR)
    bint LT "libGAP_LT"(Obj opL, Obj opR)

cdef extern from "<gap/calls.h>":
    bint IS_FUNC "libGAP_IS_FUNC"(Obj)

cdef extern from "<gap/plist.h>":
    Obj NEW_PLIST "libGAP_NEW_PLIST"(int type, int len)
    bint IS_PLIST "libGAP_IS_PLIST"(Obj lst)
    int LEN_PLIST "libGAP_LEN_PLIST"(Obj lst)
    Obj ELM_PLIST "libGAP_ELM_PLIST"(Obj lst, int pos)

cdef extern from "<gap/lists.h>":
    void UNB_LIST "libGAP_UNB_LIST"(Obj lst, int pos)
    bint IS_LIST "libGAP_IS_LIST"(Obj lst)
    int LEN_LIST "libGAP_LEN_LIST"(Obj lst)
    Obj ELM_LIST "libGAP_ELM_LIST"(Obj lst, int pos)
    void ASS_LIST "libGAP_ASS_LIST"(Obj lst, int pos, Obj elt)

cdef extern from "<gap/listfunc.h>":
    void AddList "libGAP_AddList"(Obj list, Obj obj)
    void AddPlist "libGAP_AddPlist"(Obj list, Obj obj)

cdef extern from "<gap/macfloat.h>":
    bint IS_MACFLOAT "libGAP_IS_MACFLOAT"(Obj obj)
    double VAL_MACFLOAT "libGAP_VAL_MACFLOAT"(Obj obj)

cdef extern from "<gap/records.h>":
    char* NAME_RNAM "libGAP_NAME_RNAM"(UInt rnam)
    UInt RNamIntg "libGAP_RNamIntg"(int i)
    bint IS_REC "libGAP_IS_REC"(Obj obj)
    Obj ELM_REC "libGAP_ELM_REC"(Obj rec, UInt rnam)
    UInt RNamName "libGAP_RNamName"(Char* name)

cdef extern from "<gap/precord.h>":
    Obj NEW_PREC "libGAP_NEW_PREC"(int len)
    int LEN_PREC "libGAP_LEN_PREC"(Obj rec)
    int GET_RNAM_PREC "libGAP_GET_RNAM_PREC"(Obj rec, int i)
    Obj GET_ELM_PREC "libGAP_GET_ELM_PREC"(Obj rec, int i)
    void AssPRec "libGAP_AssPRec"(Obj rec, UInt rnam, Obj val)
    void UnbPRec "libGAP_UnbPRec"(Obj rec, UInt rnam)
    bint IsbPRec "libGAP_IsbPRec"(Obj rec, UInt rnam)
    Obj ElmPRec "libGAP_ElmPRec"(Obj rec, UInt rnam)

cdef extern from "<gap/cyclotom.h>":
    pass

cdef extern from "<gap/bool.h>":
    cdef Obj libGAP_True
    cdef Obj libGAP_False

cdef extern from "<gap/vars.h>":
     cdef int T_LVARS "libGAP_T_LVARS"
     Obj BottomLVars "libGAP_BottomLVars"
