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
    ctypedef char Char "Char"
    ctypedef int Int "Int"
    ctypedef unsigned char UChar "UChar"

cdef extern from "<gap/libgap-api.h>":
#    void libgap_initialize(int argc, char** argv)
#    void libgap_set_gasman_callback(libgap_gasman_callback_ptr callback)
#    ctypedef void(*libgap_error_func_ptr)(char* msg)
#    void libgap_set_error_handler(libgap_error_func_ptr error_handler)
#    void libgap_call_error_handler()
#    void libgap_finalize()
    ctypedef void (*CallbackFunc)()
    void GAP_Initialize(int argc, char ** argv, char ** env, 
        CallbackFunc, CallbackFunc)
    void libgap_start_interaction(char* inputline)
    char* libgap_get_output()
    char* libgap_get_error()
    void libgap_mark_stack_bottom()

cdef extern from "<gap/code.h>":
    ctypedef unsigned int Stat "Stat"
    ctypedef Stat* PtrBody "PtrBody"

cdef extern from "<gap/gap.h>":
    ctypedef unsigned int UInt "UInt"
    ctypedef void* ExecStatus "ExecStatus"
    void ViewObjHandler "ViewObjHandler"(void*)
    void InitializeGap "InitializeGap"(int*, char** argv)
    void set_system_variables "set_system_variables"(char**, char**)
    cdef UInt Last
    cdef UInt Last2
    cdef UInt Last3
    cdef ExecStatus STATUS_END
    cdef ExecStatus STATUS_RETURN_VAL
    cdef ExecStatus STATUS_RETURN_VOID
    cdef ExecStatus STATUS_TNM
    cdef ExecStatus STATUS_QUIT
    cdef ExecStatus STATUS_EOF
    cdef ExecStatus STATUS_ERROR
    cdef ExecStatus STATUS_QQUIT

cdef extern from "<gap/objects.h>":
    ctypedef void* Obj "Obj"
    bint IS_MUTABLE_OBJ(Obj obj)
    bint IS_COPYABLE_OBJ(Obj obj)
    Obj SHALLOW_COPY_OBJ(Obj obj)
    Obj CopyObj(Obj obj, int mut)

    bint IS_INTOBJ(Obj obj)
    Obj INTOBJ_INT(Int)
    Int INT_INTOBJ(Obj)
    UInt TNUM_OBJ(Obj obj)
    char* TNAM_OBJ(Obj obj)
    cdef int FIRST_REAL_TNUM
    cdef int FIRST_CONSTANT_TNUM
    cdef int T_INTPOS
    cdef int T_INTNEG
    cdef int T_RAT
    cdef int T_CYC
    cdef int T_FFE
    cdef int T_PERM2
    cdef int T_PERM4
    cdef int T_BOOL
    cdef int T_CHAR
    cdef int T_FUNCTION
    cdef int T_FLAGS
    cdef int T_MACFLOAT
    cdef int T_RESERVED_BY_GAP
    cdef int LAST_CONSTANT_TNUM
    cdef int IMMUTABLE
    cdef int FIRST_IMM_MUT_TNUM
    cdef int FIRST_RECORD_TNUM
    cdef int T_PREC
    cdef int LAST_RECORD_TNUM
    cdef int FIRST_LIST_TNUM
    cdef int FIRST_PLIST_TNUM
    cdef int T_PLIST
    cdef int T_PLIST_NDENSE
    cdef int T_PLIST_DENSE
    cdef int T_PLIST_DENSE_NHOM
    cdef int T_PLIST_DENSE_NHOM_SSORT
    cdef int T_PLIST_DENSE_NHOM_NSORT
    cdef int T_PLIST_EMPTY
    cdef int T_PLIST_HOM
    cdef int T_PLIST_HOM_NSORT
    cdef int T_PLIST_HOM_SSORT "T_PLIST_HOM_SSORT"
    cdef int T_PLIST_TAB "T_PLIST_TAB"
    cdef int T_PLIST_TAB_NSORT "T_PLIST_TAB_NSORT"
    cdef int T_PLIST_TAB_SSORT "T_PLIST_TAB_SSORT"
    cdef int T_PLIST_TAB_RECT "T_PLIST_TAB_RECT"
    cdef int T_PLIST_TAB_RECT_NSORT "T_PLIST_TAB_RECT_NSORT"
    cdef int T_PLIST_TAB_RECT_SSORT "T_PLIST_TAB_RECT_SSORT"
    cdef int T_PLIST_CYC "T_PLIST_CYC"
    cdef int T_PLIST_CYC_NSORT "T_PLIST_CYC_NSORT"
    cdef int T_PLIST_CYC_SSORT "T_PLIST_CYC_SSORT"
    cdef int T_PLIST_FFE "T_PLIST_FFE"
    cdef int LAST_PLIST_TNUM "LAST_PLIST_TNUM"
    cdef int T_RANGE_NSORT "T_RANGE_NSORT"
    cdef int T_RANGE_SSORT "T_RANGE_SSORT"
    cdef int T_BLIST "T_BLIST"
    cdef int T_BLIST_NSORT "T_BLIST_NSORT"
    cdef int T_BLIST_SSORT "T_BLIST_SSORT"
    cdef int T_STRING "T_STRING"
    cdef int T_STRING_NSORT "T_STRING_NSORT"
    cdef int T_STRING_SSORT "T_STRING_SSORT"
    cdef int LAST_LIST_TNUM "LAST_LIST_TNUM"
    cdef int LAST_IMM_MUT_TNUM "LAST_IMM_MUT_TNUM"
    cdef int FIRST_EXTERNAL_TNUM "FIRST_EXTERNAL_TNUM"
    cdef int T_COMOBJ "T_COMOBJ"
    cdef int T_POSOBJ "T_POSOBJ"
    cdef int T_DATOBJ "T_DATOBJ"
    cdef int T_WPOBJ "T_WPOBJ"
    cdef int LAST_EXTERNAL_TNUM "LAST_EXTERNAL_TNUM"
    cdef int LAST_REAL_TNUM "LAST_REAL_TNUM"
    cdef int LAST_VIRTUAL_TNUM "LAST_VIRTUAL_TNUM"
    cdef int FIRST_COPYING_TNUM "FIRST_COPYING_TNUM"
    cdef int COPYING "COPYING"
    cdef int LAST_COPYING_TNUM "LAST_COPYING_TNUM"
    cdef int FIRST_TESTING_TNUM "FIRST_TESTING_TNUM"
    cdef int TESTING "TESTING"
    cdef int LAST_TESTING_TNUM "LAST_TESTING_TNUM"

cdef extern from "<gap/integer.h>":
    Int IS_INT(Obj)

cdef extern from "<gap/read.h>":
    void* ReadEvalCommand "ReadEvalCommand"(Obj context, UInt *dualSemicolon)
    void* ReadEvalFile "ReadEvalFile"()
    void* ReadEvalResult "ReadEvalResult"
    bint READ_ERROR "READ_ERROR"()

cdef extern from "<gap/scanner.h>":
    void ClearError "ClearError"()
    UInt NrError "NrError"
    UInt Symbol "Symbol"
    void GetSymbol "GetSymbol"()
    void Match "Match" (UInt symbol, char* msg, UInt skipto)
    int S_ILLEGAL "S_ILLEGAL"
    int S_IDENT "S_IDENT"
    int S_UNBIND "S_UNBIND"
    int S_ISBOUND "S_ISBOUND"
    int S_TRYNEXT "S_TRYNEXT"
    int S_INFO "S_INFO"
    int S_ASSERT "S_ASSERT"
    int S_SAVEWS "S_SAVEWS"
    int S_LOADWS "S_LOADWS"
    int S_LBRACK "S_LBRACK"
    int S_LBRACE "S_LBRACE"
    int S_BLBRACK "S_BLBRACK"
    int S_BLBRACE "S_BLBRACE"
    int S_RBRACK "S_RBRACK"
    int S_RBRACE "S_RBRACE"
    int S_DOT "S_DOT"
    int S_BDOT "S_BDOT"
    int S_LPAREN "S_LPAREN"
    int S_RPAREN "S_RPAREN"
    int S_COMMA "S_COMMA"
    int S_DOTDOT "S_DOTDOT"
    int S_COLON "S_COLON"
    int S_PARTIALINT "S_PARTIALINT"
    int S_INT "S_INT"
    int S_TRUE "S_TRUE"
    int S_FALSE "S_FALSE"
    int S_CHAR "S_CHAR"
    int S_STRING "S_STRING"
    int S_PARTIALSTRING "S_PARTIALSTRING"
    int S_REC "S_REC"
    int S_FUNCTION "S_FUNCTION"
    int S_LOCAL "S_LOCAL"
    int S_END "S_END"
    int S_MAPTO "S_MAPTO"
    int S_MULT "S_MULT"
    int S_DIV "S_DIV"
    int S_MOD "S_MOD"
    int S_POW "S_POW"
    int S_PLUS "S_PLUS"
    int S_MINUS "S_MINUS"
    int S_EQ "S_EQ"
    int S_LT "S_LT"
    int S_GT "S_GT"
    int S_NE "S_NE"
    int S_LE "S_LE"
    int S_GE "S_GE"
    int S_IN "S_IN"
    int S_NOT "S_NOT"
    int S_AND "S_AND"
    int S_OR "S_OR"
    int S_ASSIGN "S_ASSIGN"
    int S_IF "S_IF"
    int S_FOR "S_FOR"
    int S_WHILE "S_WHILE"
    int S_REPEAT "S_REPEAT"
    int S_THEN "S_THEN"
    int S_ELIF "S_ELIF"
    int S_ELSE "S_ELSE"
    int S_FI "S_FI"
    int S_DO "S_DO"
    int S_OD "S_OD"
    int S_UNTIL "S_UNTIL"
    int S_BREAK "S_BREAK"
    int S_RETURN "S_RETURN"
    int S_QUIT "S_QUIT"
    int S_QQUIT "S_QQUIT"
    int S_CONTINUE "S_CONTINUE"
    int S_SEMICOLON "S_SEMICOLON"
    int S_EOF "S_EOF"

cdef extern from "<gap/gvars.h>":
    UInt GVarName "GVarName"(char* name)
    void AssGVar "AssGVar"(UInt gvar, Obj val)
    Obj VAL_GVAR "VAL_GVAR"(UInt gvar)

cdef extern from "<gap/stringobj.h>":
    char* CSTR_STRING "CSTR_STRING"(Obj list)
    int GET_LEN_STRING "GET_LEN_STRING"(Obj list)
    bint IS_STRING "IS_STRING"(Obj obj)
    bint IsStringConv "IsStringConv"(Obj obj)
    bint ConvString "ConvString"(Obj obj)
    void C_NEW_STRING "C_NEW_STRING"(Obj new_gap_string, int length, char* c_string)

cdef extern from "<gap/gasman.h>":
    void InitGlobalBag "InitGlobalBag"(Obj* addr, char* cookie)
    Obj NewBag "NewBag"(UInt type, UInt size)
    void CHANGED_BAG "CHANGED_BAG"(Obj bag)
    void MARK_BAG "MARK_BAG"(Obj bag)
    bint IS_MARKED_ALIVE "IS_MARKED_ALIVE"(Obj bag)
    bint IS_MARKED_DEAD "IS_MARKED_DEAD"(Obj bag)
    bint IS_MARKED_HALFDEAD "IS_MARKED_HALFDEAD"(Obj bag)
    cdef UInt NrAllBags "NrAllBags"
    cdef UInt SizeAllBags "SizeAllBags"
    cdef UInt NrLiveBags "NrLiveBags"
    cdef UInt SizeLiveBags "SizeLiveBags"
    cdef UInt NrDeadBags "NrDeadBags"
    cdef UInt SizeDeadBags "SizeDeadBags"
    cdef UInt NrHalfDeadBags "NrHalfDeadBags"
    UInt CollectBags "CollectBags"(UInt size, UInt full)
    void CallbackForAllBags "CallbackForAllBags"(void (*func)(Obj))
    char* TNAM_BAG "TNAM_BAG"(Obj obj)
    UInt TNUM_BAG "TNUM_BAG"(Obj)
    UInt SIZE_BAG "SIZE_BAG"(Obj)
    void CheckMasterPointers "CheckMasterPointers"()
    Obj* MptrBags "MptrBags"
    Obj* YoungBags "YoungBags"
    Obj* OldBags "OldBags"
    Obj* AllocBags "AllocBags"
    Obj* MarkedBags "MarkedBags"
    Obj* ChangedBags "ChangedBags"

# in gasman.c but not declared in gasman.h
cdef extern Obj* StopBags "StopBags"
cdef extern Obj* EndBags "EndBags"

cdef extern from "<gap/ariths.h>":
    Obj SUM "SUM" (Obj, Obj)
    Obj DIFF "DIFF"(Obj, Obj)
    Obj PROD "PROD"(Obj, Obj)
    Obj QUO "QUO"(Obj, Obj)
    Obj POW "POW"(Obj, Obj)
    Obj MOD "MOD"(Obj, Obj)
    Obj CALL_0ARGS "CALL_0ARGS"(Obj f)              # 0 arguments
    Obj CALL_1ARGS "CALL_1ARGS"(Obj f, Obj a1)      # 1 argument
    Obj CALL_2ARGS "CALL_2ARGS"(Obj f, Obj a1, Obj a2)
    Obj CALL_3ARGS "CALL_3ARGS"(Obj f, Obj a1, Obj a2, Obj a3)
    Obj CALL_4ARGS "CALL_4ARGS"(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4)
    Obj CALL_5ARGS "CALL_5ARGS"(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4, Obj a5)
    Obj CALL_6ARGS "CALL_6ARGS"(Obj f, Obj a1, Obj a2, Obj a3,
                                 Obj a4, Obj a5, Obj a6)
    Obj CALL_XARGS "CALL_XARGS"(Obj f, Obj args)   # more than 6 arguments
    bint EQ "EQ"(Obj opL, Obj opR)
    bint LT "LT"(Obj opL, Obj opR)

cdef extern from "<gap/calls.h>":
    bint IS_FUNC "IS_FUNC"(Obj)

cdef extern from "<gap/plist.h>":
    Obj NEW_PLIST "NEW_PLIST"(int type, int len)
    bint IS_PLIST "IS_PLIST"(Obj lst)
    int LEN_PLIST "LEN_PLIST"(Obj lst)
    Obj ELM_PLIST "ELM_PLIST"(Obj lst, int pos)

cdef extern from "<gap/lists.h>":
    void UNB_LIST "UNB_LIST"(Obj lst, int pos)
    bint IS_LIST "IS_LIST"(Obj lst)
    int LEN_LIST "LEN_LIST"(Obj lst)
    Obj ELM_LIST "ELM_LIST"(Obj lst, int pos)
    void ASS_LIST "ASS_LIST"(Obj lst, int pos, Obj elt)

cdef extern from "<gap/listfunc.h>":
    void AddList "AddList"(Obj list, Obj obj)
    void AddPlist "AddPlist"(Obj list, Obj obj)

cdef extern from "<gap/macfloat.h>":
    bint IS_MACFLOAT "IS_MACFLOAT"(Obj obj)
    double VAL_MACFLOAT "VAL_MACFLOAT"(Obj obj)

cdef extern from "<gap/records.h>":
    char* NAME_RNAM "NAME_RNAM"(UInt rnam)
    UInt RNamIntg "RNamIntg"(int i)
    bint IS_REC "IS_REC"(Obj obj)
    Obj ELM_REC "ELM_REC"(Obj rec, UInt rnam)
    UInt RNamName "RNamName"(Char* name)

cdef extern from "<gap/precord.h>":
    Obj NEW_PREC "NEW_PREC"(int len)
    int LEN_PREC "LEN_PREC"(Obj rec)
    int GET_RNAM_PREC "GET_RNAM_PREC"(Obj rec, int i)
    Obj GET_ELM_PREC "GET_ELM_PREC"(Obj rec, int i)
    void AssPRec "AssPRec"(Obj rec, UInt rnam, Obj val)
    void UnbPRec "UnbPRec"(Obj rec, UInt rnam)
    bint IsbPRec "IsbPRec"(Obj rec, UInt rnam)
    Obj ElmPRec "ElmPRec"(Obj rec, UInt rnam)

cdef extern from "<gap/cyclotom.h>":
    pass

cdef extern from "<gap/bool.h>":
    cdef Obj GAP_True "True"
    cdef Obj GAP_False "False"


cdef extern from "<gap/gapstate.h>":
    Obj BottomLVars "BottomLVars"

cdef extern from "<gap/vars.h>":
    cdef int T_LVARS "T_LVARS"
