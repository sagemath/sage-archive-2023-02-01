/****************************************************************************
**
*W  gap.c                       GAP source                       Frank Celler
*W                                                         & Martin Schoenert
**
*H  @(#)$Id: gap.c,v 4.185.2.7 2007/10/05 14:05:17 gap Exp $
**
*Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
*Y  (C) 1998 School Math and Comp. Sci., University of St.  Andrews, Scotland
*Y  Copyright (C) 2002 The GAP Group
**
**  This file contains the various read-eval-print loops and  related  stuff.
*/
#include        <stdio.h>
#include        <setjmp.h>              /* jmp_buf, setjmp, longjmp        */
#include        <string.h>              /* memcpy */

#include        "system.h"              /* system dependent part           */

const char * Revision_gap_c =
"@(#)$Id: gap.c,v 4.185.2.7 2007/10/05 14:05:17 gap Exp $";

extern char * In;

#include        "gasman.h"              /* garbage collector               */
#include        "objects.h"             /* objects                         */
#include        "scanner.h"             /* scanner                         */

#define INCLUDE_DECLARATION_PART
#include        "gap.h"                 /* error handling, initialisation  */
#undef  INCLUDE_DECLARATION_PART

#include        "read.h"                /* reader                          */

#include        "gvars.h"               /* global variables                */
#include        "calls.h"               /* generic call mechanism          */
#include        "opers.h"               /* generic operations              */

#include        "ariths.h"              /* basic arithmetic                */

#include        "integer.h"             /* integers                        */
#include        "rational.h"            /* rationals                       */
#include        "cyclotom.h"            /* cyclotomics                     */
#include        "finfield.h"            /* finite fields and ff elements   */

#include        "bool.h"                /* booleans                        */
#include        "float.h"               /* machine doubles                 */
#include        "permutat.h"            /* permutations                    */

#include        "records.h"             /* generic records                 */
#include        "precord.h"             /* plain records                   */

#include        "lists.h"               /* generic lists                   */
#include        "listoper.h"            /* operations for generic lists    */
#include        "listfunc.h"            /* functions for generic lists     */
#include        "plist.h"               /* plain lists                     */
#include        "set.h"                 /* plain sets                      */
#include        "vector.h"              /* functions for plain vectors     */
#include        "vecffe.h"              /* functions for fin field vectors */
#include        "blister.h"             /* boolean lists                   */
#include        "range.h"               /* ranges                          */
#include        "string.h"              /* strings                         */
#include        "vecgf2.h"              /* functions for GF2 vectors       */
#include        "vec8bit.h"             /* functions for other compressed
					   GF(q) vectors                   */

#include        "objfgelm.h"            /* objects of free groups          */
#include        "objpcgel.h"            /* objects of polycyclic groups    */
#include        "objscoll.h"            /* single collector                */
#include        "objccoll.h"            /* combinatorial collector         */
#include        "objcftl.h"             /* from the left collect           */

#include        "dt.h"                  /* deep thought                    */
#include        "dteval.h"              /* deep though evaluation          */

#include        "sctable.h"             /* structure constant table        */
#include        "costab.h"              /* coset table                     */
#include        "tietze.h"              /* tietze helper functions         */

#include        "code.h"                /* coder                           */

#include        "vars.h"                /* variables                       */
#include        "exprs.h"               /* expressions                     */
#include        "stats.h"               /* statements                      */
#include        "funcs.h"               /* functions                       */


#include        "intrprtr.h"            /* interpreter                     */

#include        "compiler.h"            /* compiler                        */

#include        "compstat.h"            /* statically linked modules       */

#include        "saveload.h"            /* saving and loading              */

#include        "streams.h"             /* streams package                 */
#include        "sysfiles.h"            /* file input/output               */
#include        "weakptr.h"             /* weak pointers                   */

#ifdef GAPMPI
#include        "gapmpi.h"              /* ParGAP/MPI			   */
#endif

#ifdef SYS_IS_MAC_MWC
#include        "macintr.h"              /* Mac interrupt handlers	      */
#endif

#include        "iostream.h"

/****************************************************************************
**

*V  Last  . . . . . . . . . . . . . . . . . . . . . . global variable  'last'
**
**  'Last',  'Last2', and 'Last3'  are the  global variables 'last', 'last2',
**  and  'last3', which are automatically  assigned  the result values in the
**  main read-eval-print loop.
*/
UInt Last;


/****************************************************************************
**
*V  Last2 . . . . . . . . . . . . . . . . . . . . . . global variable 'last2'
*/
UInt Last2;


/****************************************************************************
**
*V  Last3 . . . . . . . . . . . . . . . . . . . . . . global variable 'last3'
*/
UInt Last3;


/****************************************************************************
**
*V  Time  . . . . . . . . . . . . . . . . . . . . . . global variable  'time'
**
**  'Time' is the global variable 'time', which is automatically assigned the
**  time the last command took.
*/
UInt Time;


/****************************************************************************
**
*V  BreakOnError  . . . . . . . . . . . . . . . . . . . . . . enter breakloop
*/
UInt BreakOnError = 1;

/****************************************************************************
**
*V  ErrorCount  . . . . . . . .  . . . . .how many times have we had an error
**                             note that this includes cases where the break
**                             loop was skipped.
*/
UInt ErrorCount = 0;


/****************************************************************************
**
*F  ViewObjHandler  . . . . . . . . . handler to view object and catch errors
**
**  This is the function actually called in Read-Eval-View loops.
**  We might be in trouble if the library has not (yet) loaded and so ViewObj
**  is not yet defined, or the fallback methods not yet installed. To avoid
**  this problem, we check, and use PrintObj if there is a problem
**
**  This function also supplies the \n after viewing.
*/
UInt ViewObjGVar;

void ViewObjHandler ( Obj obj )
{
  volatile Obj        func;
  jmp_buf             readJmpError;

  /* get the function                                                    */
  func = ValAutoGVar(ViewObjGVar);

  /* if non-zero use this function, otherwise use `PrintObj'             */
  memcpy( readJmpError, ReadJmpError, sizeof(jmp_buf) );
  if ( ! READ_ERROR() ) {
    if ( func == 0 || TNUM_OBJ(func) != T_FUNCTION ) {
      PrintObj(obj);
    }
    else {
      ViewObj( obj );
    }
    Pr( "\n", 0L, 0L );
    memcpy( ReadJmpError, readJmpError, sizeof(jmp_buf) );
  }
  else {
    memcpy( ReadJmpError, readJmpError, sizeof(jmp_buf) );
  }
}


/****************************************************************************
**
*F  main( <argc>, <argv> )  . . . . . . .  main program, read-eval-print loop
*/
Obj AtExitFunctions;

Obj AlternativeMainLoop = 0;

UInt SaveOnExitFileGVar;

UInt QUITTINGGVar;

Obj OnGapPromptHook = 0;

static char **sysargv;
static char **sysenviron;

#ifdef COMPILECYGWINDLL
int realmain (
#else
int main (
#endif
	  int                 argc,
          char *              argv [],
          char *              environ [] )

{
  ExecStatus          status;                 /* result of ReadEvalCommand*/
  UInt                type;                   /* result of compile       */
  UInt                time;                   /* start time              */
  Obj                 func;                   /* function (compiler)     */
  Int4                crc;                    /* crc of file to compile  */
  volatile UInt       i;                      /* loop variable           */
  Obj                 SaveOnExitFile;         /* contents of the GVar    */

  sysargv = argv;
  sysenviron = environ;

  /* initialize everything                                               */
  InitializeGap( &argc, argv );
  if (UserHasQUIT)		/* maybe the user QUIT from the initial
				   read of init.g */
    goto finalize;

  /* maybe compile                                                       */
  if ( SyCompilePlease ) {
    if ( ! OpenInput(SyCompileInput) ) {
      SyExit(1);
    }
    func = READ_AS_FUNC();
    crc  = SyGAPCRC(SyCompileInput);
    if (SyStrlen(SyCompileOptions) != 0)
      SetCompileOpts(SyCompileOptions);
    type = CompileFunc(
		       SyCompileOutput,
		       func,
		       SyCompileName,
		       crc,
		       SyCompileMagic1 );
    if ( type == 0 )
      SyExit( 1 );
    SyExit( 0 );
  }

  if (AlternativeMainLoop != (Obj) 0)
    {
      if (!IS_FUNC(AlternativeMainLoop))
	{
	  Pr("#E AlternativeMainLoop set to non-function, ignoring\n",0L,0L);
	}
      else
	{
	  ClearError();
	  if (READ_ERROR())
	    {
	      if (UserHasQUIT)
		goto finalize;
	      ClearError();
	    }
	  ExecBegin( BottomLVars );
	  CALL_0ARGS(AlternativeMainLoop);
	  /* If we ever get to here then the AlternativeMainLoop function actual
	     exited. In this case, we are done */
	  ExecEnd(0);
	  goto finalize;
	}
    }

  /* read-eval-print loop                                                */
  while ( 1 ) {

    /* start the stopwatch                                             */
    time = SyTime();

    /* read and evaluate one command                                   */
    Prompt = "gap> ";
    ClearError();

    /* here is a hook: */
    if (OnGapPromptHook) {
      if (!IS_FUNC(OnGapPromptHook))
	{
	  Pr("#E OnGapPromptHook set to non-function, ignoring\n",0L,0L);
	}
      else
        {
          Call0ArgsInNewReader(OnGapPromptHook);
          /* Recover from a potential break loop: */
          Prompt = "gap> ";
          ClearError();
        }
    }

    /* now enter the read-eval-command loop: */
    status = ReadEvalCommand();
    if (UserHasQUIT)
      break;

    /* stop the stopwatch                                              */
    AssGVar( Time, INTOBJ_INT( SyTime() - time ) );

    /* handle ordinary command                                         */
    if ( status == STATUS_END && ReadEvalResult != 0 ) {

      /* remember the value in 'last' and the time in 'time'         */
      AssGVar( Last3, VAL_GVAR( Last2 ) );
      AssGVar( Last2, VAL_GVAR( Last  ) );
      AssGVar( Last,  ReadEvalResult   );

      /* print the result                                            */
      if ( ! DualSemicolon ) {
	ViewObjHandler( ReadEvalResult );
      }

    }

    /* handle return-value or return-void command                      */
    else if ( status & (STATUS_RETURN_VAL | STATUS_RETURN_VOID) ) {
      Pr( "'return' must not be used in main read-eval-print loop\n",
	  0L, 0L );
    }

    /* handle quit command or <end-of-file>                            */
    else if ( status & (STATUS_EOF | STATUS_QUIT ) ) {
      break;
    }

    /* handle QUIT */
    else if (status & (STATUS_QQUIT)) {
      UserHasQUIT = 1;
      break;
    }

    /* stop the stopwatch                                          */
    AssGVar( Time, INTOBJ_INT( SyTime() - time ) );

    if (UserHasQuit)
      {
	FlushRestOfInputLine();
	UserHasQuit = 0;	/* quit has done its job if we are here */
      }

  }

 finalize:
  /* The QUITTING variable is made available to the AtExitFunction */
  MakeReadWriteGVar(QUITTINGGVar);
  AssGVar(QUITTINGGVar, UserHasQUIT ? True : False );
  MakeReadOnlyGVar(QUITTINGGVar);

  /* call the exit functions                                             */
  BreakOnError = 0;

  for ( i = 1;  i <= LEN_PLIST(AtExitFunctions);  i++ ) {
    if ( setjmp(ReadJmpError) == 0 ) {
      func = ELM_PLIST( AtExitFunctions, i );
      CALL_0ARGS(func);
    }
  }

  /* Possibly save the workspace */
  if ( !UserHasQUIT
       && (SaveOnExitFile = VAL_GVAR(SaveOnExitFileGVar))
       && SaveOnExitFile != False)
    {
      if (IsStringConv(SaveOnExitFile))
	{
	  if (SaveWorkspace(SaveOnExitFile) == True)
	    Pr("Workspace saved in %s\n",
	       (Int)CSTR_STRING(SaveOnExitFile), 0);
	  else
	    Pr("Attempt to save in %s failed\n",
	       (Int)CSTR_STRING(SaveOnExitFile), 0);
	}
      else
	Pr("SaveOnExitFile is a %s not a filename\n",
	   (Int)TNAM_OBJ(SaveOnExitFile), 0);
    }

  /* exit to the operating system, the return is there to please lint    */
  SyExit(0);
  return 0;
}


/****************************************************************************
**
*F  FuncID_FUNC( <self>, <val1> ) . . . . . . . . . . . . . . . return <val1>
*/
Obj FuncID_FUNC (
		 Obj                 self,
		 Obj                 val1 )
{
  return val1;
}


/****************************************************************************
**
*F  FuncRuntime( <self> ) . . . . . . . . . . . . internal function 'Runtime'
**
**  'FuncRuntime' implements the internal function 'Runtime'.
**
**  'Runtime()'
**
**  'Runtime' returns the time spent since the start of GAP in  milliseconds.
**  How much time execution of statements take is of course system dependent.
**  The accuracy of this number is also system dependent.
*/
Obj FuncRuntime (
		 Obj                 self )
{
  return INTOBJ_INT( SyTime() );
}



Obj FuncRUNTIMES( Obj     self)
{
  Obj    res;
#if HAVE_GETRUSAGE
  res = NEW_PLIST(T_PLIST, 4);
  SET_LEN_PLIST(res, 4);
  SET_ELM_PLIST(res, 1, INTOBJ_INT( SyTime() ));
  SET_ELM_PLIST(res, 2, INTOBJ_INT( SyTimeSys() ));
  SET_ELM_PLIST(res, 3, INTOBJ_INT( SyTimeChildren() ));
  SET_ELM_PLIST(res, 4, INTOBJ_INT( SyTimeChildrenSys() ));
#else
  res = INTOBJ_INT( SyTime() );
#endif
  return res;

}



/****************************************************************************
**
*F  FuncSizeScreen( <self>, <args> )  . . . .  internal function 'SizeScreen'
**
**  'FuncSizeScreen'  implements  the  internal  function 'SizeScreen' to get
**  or set the actual screen size.
**
**  'SizeScreen()'
**
**  In this form 'SizeScreen'  returns the size of the  screen as a list with
**  two  entries.  The first is  the length of each line,  the  second is the
**  number of lines.
**
**  'SizeScreen( [ <x>, <y> ] )'
**
**  In this form 'SizeScreen' sets the size of the screen.  <x> is the length
**  of each line, <y> is the number of lines.  Either value may  be  missing,
**  to leave this value unaffected.  Note that those parameters can  also  be
**  set with the command line options '-x <x>' and '-y <y>'.
*/
Obj FuncSizeScreen (
		    Obj                 self,
		    Obj                 args )
{
  Obj                 size;           /* argument and result list        */
  Obj                 elm;            /* one entry from size             */
  UInt                len;            /* length of lines on the screen   */
  UInt                nr;             /* number of lines on the screen   */

  /* check the arguments                                                 */
  while ( ! IS_SMALL_LIST(args) || 1 < LEN_LIST(args) ) {
    args = ErrorReturnObj(
			  "Function: number of arguments must be 0 or 1 (not %d)",
			  LEN_LIST(args), 0L,
			  "you can replace the argument list <args> via 'return <args>;'" );
  }

  /* get the arguments                                                   */
  if ( LEN_LIST(args) == 0 ) {
    size = NEW_PLIST( T_PLIST, 0 );
    SET_LEN_PLIST( size, 0 );
  }

  /* otherwise check the argument                                        */
  else {
    size = ELM_LIST( args, 1 );
    while ( ! IS_SMALL_LIST(size) || 2 < LEN_LIST(size) ) {
      size = ErrorReturnObj(
			    "SizeScreen: <size> must be a list of length 2",
			    0L, 0L,
			    "you can replace <size> via 'return <size>;'" );
    }
  }

  /* extract the length                                                  */
  if ( LEN_LIST(size) < 1 || ELM0_LIST(size,1) == 0 ) {
    len = SyNrCols;
  }
  else {
    elm = ELMW_LIST(size,1);
    while ( TNUM_OBJ(elm) != T_INT ) {
      elm = ErrorReturnObj(
			   "SizeScreen: <x> must be an integer",
			   0L, 0L,
			   "you can replace <x> via 'return <x>;'" );
    }
    len = INT_INTOBJ( elm );
    if ( len < 20  )  len = 20;
    if ( 256 < len )  len = 256;
  }

  /* extract the number                                                  */
  if ( LEN_LIST(size) < 2 || ELM0_LIST(size,2) == 0 ) {
    nr = SyNrRows;
  }
  else {
    elm = ELMW_LIST(size,2);
    while ( TNUM_OBJ(elm) != T_INT ) {
      elm = ErrorReturnObj(
			   "SizeScreen: <y> must be an integer",
			   0L, 0L,
			   "you can replace <y> via 'return <y>;'" );
    }
    nr = INT_INTOBJ( elm );
    if ( nr < 10 )  nr = 10;
  }

  /* set length and number                                               */
  SyNrCols = len;
  SyNrRows = nr;

  /* make and return the size of the screen                              */
  size = NEW_PLIST( T_PLIST, 2 );
  SET_LEN_PLIST( size, 2 );
  SET_ELM_PLIST( size, 1, INTOBJ_INT(len) );
  SET_ELM_PLIST( size, 2, INTOBJ_INT(nr)  );
  return size;

}


/****************************************************************************
**
*F  FuncWindowCmd( <self>, <args> ) . . . . . . . .  execute a window command
*/
static Obj WindowCmdString;

Obj FuncWindowCmd (
		   Obj	      	    self,
		   Obj             args )
{
  Obj             tmp;
  Obj       	    list;
  Int             len;
  Int             n,  m;
  Int             i;
  Char *          ptr;
  Char *          qtr;

  /* check arguments                                                     */
  while ( ! IS_SMALL_LIST(args) ) {
    args = ErrorReturnObj( "argument list must be a list (not a %s)",
			   (Int)TNAM_OBJ(args), 0L,
			   "you can replace the argument list <args> via 'return <args>;'" );

  }
  tmp = ELM_LIST(args,1);
  while ( ! IsStringConv(tmp) || 3 != LEN_LIST(tmp) ) {
    while ( ! IsStringConv(tmp) ) {
      tmp = ErrorReturnObj( "<cmd> must be a string (not a %s)",
			    (Int)TNAM_OBJ(tmp), 0L,
			    "you can replace <cmd> via 'return <cmd>;'" );
    }
    if ( 3 != LEN_LIST(tmp) ) {
      tmp = ErrorReturnObj( "<cmd> must be a string of length 3",
			    0L, 0L,
			    "you can replace <cmd> via 'return <cmd>;'" );
    }
  }

  /* compute size needed to store argument string                        */
  len = 13;
  for ( i = 2;  i <= LEN_LIST(args);  i++ )
    {
      tmp = ELM_LIST( args, i );
      while ( TNUM_OBJ(tmp) != T_INT && ! IsStringConv(tmp) ) {
	tmp = ErrorReturnObj(
			     "%d. argument must be a string or integer (not a %s)",
			     i, (Int)TNAM_OBJ(tmp),
			     "you can replace the argument <arg> via 'return <arg>;'" );
	SET_ELM_PLIST( args, i, tmp );
      }
      if ( TNUM_OBJ(tmp) == T_INT )
	len += 12;
      else
	len += 12 + LEN_LIST(tmp);
    }
  if ( SIZE_OBJ(WindowCmdString) <= len ) {
    ResizeBag( WindowCmdString, 2*len+1 );
  }

  /* convert <args> into an argument string                              */
  ptr  = (Char*) CSTR_STRING(WindowCmdString);
  *ptr = '\0';

  /* first the command name                                              */
  SyStrncat( ptr, CSTR_STRING( ELM_LIST(args,1) ), 3 );
  ptr += 3;

  /* and now the arguments                                               */
  for ( i = 2;  i <= LEN_LIST(args);  i++ )
    {
      tmp = ELM_LIST(args,i);

      if ( TNUM_OBJ(tmp) == T_INT ) {
	*ptr++ = 'I';
	m = INT_INTOBJ(tmp);
	for ( m = (m<0)?-m:m;  0 < m;  m /= 10 )
	  *ptr++ = (m%10) + '0';
	if ( INT_INTOBJ(tmp) < 0 )
	  *ptr++ = '-';
	else
	  *ptr++ = '+';
      }
      else {
	*ptr++ = 'S';
	m = LEN_LIST(tmp);
	for ( ; 0 < m;  m/= 10 )
	  *ptr++ = (m%10) + '0';
	*ptr++ = '+';
	qtr = CSTR_STRING(tmp);
	for ( m = LEN_LIST(tmp);  0 < m;  m-- )
	  *ptr++ = *qtr++;
      }
    }
  *ptr = 0;

  /* now call the window front end with the argument string              */
  qtr = CSTR_STRING(WindowCmdString);
  ptr = SyWinCmd( qtr, SyStrlen(qtr) );
  len = SyStrlen(ptr);

  /* now convert result back into a list                                 */
  list = NEW_PLIST( T_PLIST, 11 );
  SET_LEN_PLIST( list, 0 );
  i = 1;
  while ( 0 < len ) {
    if ( *ptr == 'I' ) {
      ptr++;
      for ( n=0,m=1; '0' <= *ptr && *ptr <= '9'; ptr++,m *= 10,len-- )
	n += (*ptr-'0') * m;
      if ( *ptr++ == '-' )
	n *= -1;
      len -= 2;
      AssPlist( list, i, INTOBJ_INT(n) );
    }
    else if ( *ptr == 'S' ) {
      ptr++;
      for ( n=0,m=1;  '0' <= *ptr && *ptr <= '9';  ptr++,m *= 10,len-- )
	n += (*ptr-'0') * m;
      ptr++; /* ignore the '+' */
      /*CCC tmp = NEW_STRING(n);
      *CSTR_STRING(tmp) = '\0';
      SyStrncat( CSTR_STRING(tmp), ptr, n ); CCC*/
      C_NEW_STRING(tmp, n, ptr);
      ptr += n;
      len -= n+2;
      AssPlist( list, i, tmp );
    }
    else {
      ErrorQuit( "unknown return value '%s'", (Int)ptr, 0 );
      return 0;
    }
    i++;
  }

  /* if the first entry is one signal an error */
  if ( ELM_LIST(list,1) == INTOBJ_INT(1) ) {
    /*CCCtmp = NEW_STRING(15);
      SyStrncat( CSTR_STRING(tmp), "window system: ", 15 );CCC*/
    C_NEW_STRING(tmp, 15, "window system: ");
    SET_ELM_PLIST( list, 1, tmp );
    SET_LEN_PLIST( list, i-1 );
    return FuncError( 0, list );
  }
  else {
    for ( m = 1;  m <= i-2;  m++ )
      SET_ELM_PLIST( list, m, ELM_PLIST(list,m+1) );
    SET_LEN_PLIST( list, i-2 );
    return list;
  }
}


/****************************************************************************
**

*F * * * * * * * * * * * * * * error functions * * * * * * * * * * * * * * *
*/



/****************************************************************************
**

*F  FuncDownEnv( <self>, <level> )  . . . . . . . . .  change the environment
*/
UInt ErrorLevel;

Obj  ErrorLVars0;
Obj  ErrorLVars;
Int  ErrorLLevel;

extern Obj BottomLVars;


void DownEnvInner( Int depth )
{
  /* if we really want to go up                                          */
  if ( depth < 0 && -ErrorLLevel <= -depth ) {
    depth = 0;
    ErrorLVars = ErrorLVars0;
    ErrorLLevel = 0;
  }
  else if ( depth < 0 ) {
    depth = -ErrorLLevel + depth;
    ErrorLVars = ErrorLVars0;
    ErrorLLevel = 0;
  }

  /* now go down                                                         */
  while ( 0 < depth
	  && ErrorLVars != BottomLVars
	  && PTR_BAG(ErrorLVars)[2] != BottomLVars ) {
    ErrorLVars = PTR_BAG(ErrorLVars)[2];
    ErrorLLevel--;
    depth--;
  }
}

Obj FuncDownEnv (
		 Obj                 self,
		 Obj                 args )
{
  Int                 depth;

  if ( LEN_LIST(args) == 0 ) {
    depth = 1;
  }
  else if ( LEN_LIST(args) == 1 && IS_INTOBJ( ELM_PLIST(args,1) ) ) {
    depth = INT_INTOBJ( ELM_PLIST( args, 1 ) );
  }
  else {
    ErrorQuit( "usage: DownEnv( [ <depth> ] )", 0L, 0L );
    return 0;
  }
  if ( ErrorLVars == 0 ) {
    Pr( "not in any function\n", 0L, 0L );
    return 0;
  }

  DownEnvInner( depth);

  /* return nothing                                                      */
  return 0;
}

Obj FuncUpEnv (
	       Obj                 self,
	       Obj                 args )
{
  Int                 depth;
  if ( LEN_LIST(args) == 0 ) {
    depth = 1;
  }
  else if ( LEN_LIST(args) == 1 && IS_INTOBJ( ELM_PLIST(args,1) ) ) {
    depth = INT_INTOBJ( ELM_PLIST( args, 1 ) );
  }
  else {
    ErrorQuit( "usage: UpEnv( [ <depth> ] )", 0L, 0L );
    return 0;
  }
  if ( ErrorLVars == 0 ) {
    Pr( "not in any function\n", 0L, 0L );
    return 0;
  }

  DownEnvInner(-depth);
  return 0;
}

/****************************************************************************
**
*F  FuncWhere( <self>, <depth> )  . . . . . . . . . . . .  print stack frames
*/
Obj FuncWhere (
	       Obj                 self,
	       Obj                 args )
{
  Obj                 currLVars;
  Int                 depth;
  Expr                call;

#ifndef NO_BRK_CALLS

  /* evaluate the argument                                               */
  if ( LEN_LIST(args) == 0 ) {
    depth = 5;
  }
  else if ( LEN_LIST(args) == 1 && IS_INTOBJ( ELM_PLIST(args,1) ) ) {
    depth = INT_INTOBJ( ELM_PLIST( args, 1 ) );
    if ( depth == 0 ) {
      /* We do this to avoid the enclosing " called from" and "..."
         when the depth is zero                                         */
      depth = -1;
    }
  }
  else {
    ErrorQuit( "usage: Where( [ <depth> ] )", 0L, 0L );
    return 0;
  }

  currLVars = CurrLVars;

  if ( ErrorLVars != 0  && ErrorLVars != BottomLVars ) {
    SWITCH_TO_OLD_LVARS( ErrorLVars );
    SWITCH_TO_OLD_LVARS( BRK_CALL_FROM() );
    if ( 0 < depth ) {
      Pr( " called from\n", 0L, 0L );
    }
    while ( CurrLVars != BottomLVars && 0 < depth ) {
      call = BRK_CALL_TO();
      if ( call == 0 ) {
	Pr( "<compiled or corrupted call value> ", 0L, 0L );
      }
#if T_PROCCALL_0ARGS
      else if ( T_PROCCALL_0ARGS <= TNUM_STAT(call)
		&& TNUM_STAT(call)  <= T_PROCCALL_XARGS ) {
#else
      else if ( TNUM_STAT(call)  <= T_PROCCALL_XARGS ) {
#endif
	PrintStat( call );
      }
      else if ( T_FUNCCALL_0ARGS <= TNUM_EXPR(call)
                && TNUM_EXPR(call)  <= T_POW ) {
        PrintExpr( call );
      }
      Pr( " called from\n", 0L, 0L );
      SWITCH_TO_OLD_LVARS( BRK_CALL_FROM() );
      depth--;
    }
    if ( 0 < depth ) {
      Pr( "<function>( <arguments> ) called from read-eval-loop\n",
          0L, 0L );
    }
    else if ( depth == 0 ) {
      Pr( "...\n", 0L, 0L );
    }
  }
  else {
    Pr( "not in any function\n", 0L, 0L );
  }

  SWITCH_TO_OLD_LVARS( currLVars );

#endif

  return 0;
}

/****************************************************************************
**
*F  FuncCallFuncTrapError( <self>, <func> )
**
*/


Obj FuncCallFuncTrapError( Obj self, Obj func)
  {
    jmp_buf readJmpError;
    ClearError();

    /* Save the old error long-jump */
    memcpy( readJmpError, ReadJmpError, sizeof(jmp_buf) );
    if (READ_ERROR())
      {
	/* Get out any old how */
	if (UserHasQUIT)
	  longjmp( readJmpError, 1);
	ExecEnd(1);
	ClearError();
	return True;
      }

    ExecBegin( BottomLVars );
    CALL_0ARGS( func );
    ExecEnd(0);
    return False;

  }

/****************************************************************************
**
*F  ErrorMode( <msg>, <arg1>, <arg2>, <args>, <msg2>, <mode> )
*/

    Obj OnBreak;			/* a Fopy of the global OnBreak,
				         which by default is set to Where. */
    Obj OnBreakMessage;			/* a Fopy of the global
                                         OnBreakMessage.                   */
    Obj OnQuit;		                /* a Copy of the global OnQuit     */

Obj ErrorHandler = (Obj) 0;		/* not yet settable from GAP level */

 Obj FuncSetErrorHandler( Obj self, Obj Handler)
   {
     Obj handler;
     if (ErrorHandler != (Obj) 0)
       handler = ErrorHandler;
     else
       handler = Fail;
     ErrorHandler = Handler;
     return handler;
   }

UInt UserHasQuit = 0;
UInt UserHasQUIT = 0;

Obj ErrorMode (
    const Char *        msg,
    Int                 arg1,
    Int                 arg2,
    Obj                 args,
    const Char *        msg2,
    Char                mode )
{
    Obj                 errorLVars0;
    Obj                 errorLVars;
    UInt                errorLLevel;
    ExecStatus          status;
    char                prompt [16];
    Obj                 errorHandler;

    ErrorCount++;
    if (ErrorCount >= ((UInt)1)<<NR_SMALL_INT_BITS)
      ErrorCount = 0;

    /* open the standard error output file                                 */
    OpenOutput( "*errout*" );
    ErrorLevel += 1;
    errorLVars0 = ErrorLVars0;
    ErrorLVars0 = CurrLVars;
    errorLVars  = ErrorLVars;
    ErrorLVars  = CurrLVars;
    errorLLevel = ErrorLLevel;
    ErrorLLevel = 0;

    /* ignore all errors when testing or quitting                          */
    if ( ( TestInput != 0 && TestOutput == Output ) || ! BreakOnError ) {
        if ( msg != (Char*)0 ) {
            Pr( msg, arg1, arg2 );
        }
        else if ( args != (Obj)0 ) {
            Pr( "Error, ", 0L, 0L );
            FuncPrint( (Obj)0, args );
        }
        Pr( "\n", 0L, 0L );

	ErrorLevel -= 1;
	ErrorLVars0 = errorLVars0;
	ErrorLVars = errorLVars;
	ErrorLLevel = errorLLevel;
	ClearError();
	CloseOutput();
        ReadEvalError();
    }

    /* See if we have an Error Handler */
    if (ErrorHandler != 0 && IS_FUNC(ErrorHandler))
      {
	Obj mess;
	Obj mess2;
	Obj ret;
	Int len;
	errorHandler = ErrorHandler;
	ErrorHandler = (Obj)0;

	/* Transform the messages into GAP strings */
	if (msg)
	  {
	    len = SyStrlen(msg);
	    C_NEW_STRING(mess, len, msg);
	  }
	else
	  {
	    mess = NEW_STRING(0);
	  }
	if (msg2)
	  {
	    len = SyStrlen(msg2);
	    C_NEW_STRING(mess2, len, msg2);
	  }
	else
	  {
	    mess2 = NEW_STRING(0);
	  }

	/* Now call the handler */
	if (args != 0)
	  ret=CALL_4ARGS(errorHandler,mess, args, mess2, ObjsChar[(Int)mode]);
	else
	  ret=CALL_4ARGS(errorHandler,mess, Fail, mess2, ObjsChar[(Int)mode]);

	/* Now handle the return, allowing for the mode */
	if (ret == True)
	  {
	    ReadEvalError();
	  }
	else if (IS_PLIST(ret))
	  {
	    if (LEN_PLIST(ret) == 0)
	      {
		if (mode == 'x')
		  return 0;
		else
		  {
		    Pr("%%E Error handler tried to return null in mode %c\n",mode,0);
		    if (mode == 'v')
		      return Fail;
		    if (mode == 'q' || mode == 'm')
		      {
			ReadEvalError();
		      }
		    Pr("Panic: impossible error mode %c\n",mode,0);
		  }
	      }
	    else if (LEN_PLIST(ret) > 1)
	      {
		Pr("%%E Error handler returned %d objects, all but the first are ignored\n",
		   LEN_PLIST(ret),0);
	      }

	    if (mode == 'v')
	      return ELM_PLIST(ret,1);

	    Pr("%%E Error handler tried to return an object in mode %c\n",mode,0);
	    if (mode == 'x')
	      {
		return 0;
	      }
	    if (mode == 'q' || mode == 'm')
	      {
		ReadEvalError();
	      }
	    Pr("Panic: impossible error mode %c\n",mode,0);
	  }
	else
	  {
	    Pr("%%E Error handler returned bad value or nothing, treating as true\n", 0L, 0L);
	  }
	ReadEvalError();
      }

    /* print the error message                                             */
    if ( msg != (Char*)0 ) {
        Pr( msg, arg1, arg2 );
    }
    else if ( args != (Obj)0 ) {
        Pr( "Error, ", 0L, 0L );
        FuncPrint( (Obj)0, args );
    }

    /* print the location                                                  */
    if ( CurrStat != 0 && msg != (Char*)0) {
	/* only print the current command if the `msg' variable is not 0.
	 * This will avoid printing the `Error("blabla")' line again */
        Pr( " at\n", 0L, 0L );
	PrintStat( CurrStat );
        Pr( "\n", 0L, 0L );
    }
    else if ( mode != 'f' ) {
        Pr( "\n", 0L, 0L );
    }

    /* try to open input for a break loop                                  */
    if ( mode == 'q' || ! OpenInput( "*errin*") ) {
        ErrorLevel -= 1;
        ErrorLVars0 = errorLVars0;
        ErrorLVars = errorLVars;
        ErrorLLevel = errorLLevel;
        CloseOutput();
        ReadEvalError();
    }
    ClearError();

    /* Call the OnBreak function. This can't be done earlier, because
       there seems to be no safe way to handle errors in it, unless we
       are about to enter a break r-e-v loop*/
    CALL_0ARGS(OnBreak);

    /* print the second message                                            */
    Pr( "Entering break read-eval-print loop ...\n", 0L, 0L );

    if ( mode == 'f' ) {
        /* If called via (Func)Error call the OnBreakMessage which the
           user may have customised                                        */
        CALL_0ARGS(OnBreakMessage);
    }
    else {
        Pr( "you can 'quit;' to quit to outer loop", 0L, 0L );
	if (mode != 'm')
	  Pr( ", or\n%s to continue", (Int)msg2, 0L );
	Pr( "\n",0L,0L);
    }

    /* read-eval-print loop                                                */
    while ( 1 ) {

        /* read and evaluate one command                                   */
        if ( ErrorLevel == 1 ) {
            Prompt = "brk> ";
        }
        else {
            prompt[0] = 'b';
            prompt[1] = 'r';
            prompt[2] = 'k';
            prompt[3] = '_';
            prompt[4] = ErrorLevel / 10 + '0';
            prompt[5] = ErrorLevel % 10 + '0';
            prompt[6] = '>';
            prompt[7] = ' ';
            prompt[8] = '\0';
            Prompt = prompt;
        }

        /* read and evaluate one command                                   */
        ClearError();
        DualSemicolon = 0;
        status = ReadEvalDebug();
	UserHasQuit = 0;	/* it is enough for quit
				 to have got us here */

        /* handle ordinary command                                         */
        if ( status == STATUS_END && ReadEvalResult != 0 ) {

            /* remember the value in 'last'                                */
            AssGVar( Last,  ReadEvalResult   );

            /* print the result                                            */
            if ( ! DualSemicolon ) {
                ViewObjHandler( ReadEvalResult );
            }

        }

        /* handle return-value                                             */
        else if ( status == STATUS_RETURN_VAL ) {
            if ( mode == 'v' ) {
                ErrorLevel -= 1;
                ErrorLVars0 = errorLVars0;
                ErrorLVars = errorLVars;
                ErrorLLevel = errorLLevel;
                CloseInput();
                ClearError();
                CloseOutput();
                return ReadEvalResult;
            }
            else {
                Pr( "'return <value>;' cannot be used in this break-loop\n",
                    0L, 0L );
            }
        }

        /* handle return-void                                             */
        else if ( status == STATUS_RETURN_VOID ) {
            if ( mode == 'x' || mode == 'f' ) {
                ErrorLevel -= 1;
                ErrorLVars0 = errorLVars0;
                ErrorLVars = errorLVars;
                ErrorLLevel = errorLLevel;
                CloseInput();
                ClearError();
                CloseOutput();
                return (Obj)0;
            }
            else {
                Pr( "'return;' cannot be used in this break-loop\n",
                    0L, 0L );
            }
        }

        /* handle quit command or <end-of-file>                            */
        else if ( status == STATUS_EOF || status == STATUS_QUIT ) {
          CALL_0ARGS(OnQuit);
	  UserHasQuit = 1;
	  break;
        }
	else if ( status == STATUS_QQUIT ) {
	  UserHasQUIT = 1;
	  break;
	}
	if (UserHasQUIT)
	  break;

    }

    /* return to the outer read-eval-print loop                            */
    ErrorLevel -= 1;
    ErrorLVars0 = errorLVars0;
    ErrorLVars = errorLVars;
    ErrorLLevel = errorLLevel;
    CloseInput();
    ClearError();
    CloseOutput();
    ReadEvalError();

    /* this is just to please GNU cc, 'ReadEvalError' never returns        */
    return 0;
}


/****************************************************************************
**
*F  ErrorQuit( <msg>, <arg1>, <arg2> )  . . . . . . . . . . .  print and quit
*/
void ErrorQuit (
    const Char *        msg,
    Int                 arg1,
    Int                 arg2 )
{
    ErrorMode( msg, arg1, arg2, (Obj)0, (Char*)0, 'q' );
}


/****************************************************************************
**
*F  ErrorQuitBound( <name> )  . . . . . . . . . . . . . . .  unbound variable
*/
void ErrorQuitBound (
    Char *              name )
{
    ErrorQuit(
        "variable '%s' must have an assigned value",
        (Int)name, 0L );
}


/****************************************************************************
**
*F  ErrorQuitFuncResult() . . . . . . . . . . . . . . . . must return a value
*/
void ErrorQuitFuncResult ( void )
{
    ErrorQuit(
        "function must return a value",
        0L, 0L );
}


/****************************************************************************
**
*F  ErrorQuitIntSmall( <obj> )  . . . . . . . . . . . . . not a small integer
*/
void ErrorQuitIntSmall (
    Obj                 obj )
{
    ErrorQuit(
        "<obj> must be a small integer (not a %s)",
        (Int)TNAM_OBJ(obj), 0L );
}


/****************************************************************************
**
*F  ErrorQuitIntSmallPos( <obj> ) . . . . . . .  not a positive small integer
*/
void ErrorQuitIntSmallPos (
    Obj                 obj )
{
    ErrorQuit(
        "<obj> must be a positive small integer (not a %s)",
        (Int)TNAM_OBJ(obj), 0L );
}

/****************************************************************************
**
*F  ErrorQuitIntPos( <obj> ) . . . . . . .  not a positive small integer
*/
void ErrorQuitIntPos (
    Obj                 obj )
{
    ErrorQuit(
        "<obj> must be a positive integer (not a %s)",
        (Int)TNAM_OBJ(obj), 0L );
}


/****************************************************************************
**
*F  ErrorQuitBool( <obj> )  . . . . . . . . . . . . . . . . . . not a boolean
*/
void ErrorQuitBool (
    Obj                 obj )
{
    ErrorQuit(
        "<obj> must be 'true' or 'false' (not a %s)",
        (Int)TNAM_OBJ(obj), 0L );
}


/****************************************************************************
**
*F  ErrorQuitFunc( <obj> )  . . . . . . . . . . . . . . . . .  not a function
*/
void ErrorQuitFunc (
    Obj                 obj )
{
    ErrorQuit(
        "<obj> must be a function (not a %s)",
        (Int)TNAM_OBJ(obj), 0L );
}


/****************************************************************************
**
*F  ErrorQuitNrArgs( <narg>, <args> ) . . . . . . . wrong number of arguments
*/
void ErrorQuitNrArgs (
    Int                 narg,
    Obj                 args )
{
    ErrorQuit(
        "Function Calls: number of arguments must be %d (not %d)",
        narg, LEN_PLIST( args ) );
}

/****************************************************************************
**
*F  ErrorQuitRange3( <first>, <second>, <last> ) . . divisibility
*/
void ErrorQuitRange3 (
		      Obj                 first,
		      Obj                 second,
		      Obj                 last)
{
    ErrorQuit(
        "Range expression <last>-<first> must be divisible by <second>-<first>, not %d %d",
        INT_INTOBJ(last)-INT_INTOBJ(first), INT_INTOBJ(second)-INT_INTOBJ(first) );
}


/****************************************************************************
**
*F  ErrorReturnObj( <msg>, <arg1>, <arg2>, <msg2> ) . .  print and return obj
*/
Obj ErrorReturnObj (
    const Char *        msg,
    Int                 arg1,
    Int                 arg2,
    const Char *        msg2 )
{
    return ErrorMode( msg, arg1, arg2, (Obj)0, msg2, 'v' );
}


/****************************************************************************
**
*F  ErrorReturnVoid( <msg>, <arg1>, <arg2>, <msg2> )  . . .  print and return
*/
void ErrorReturnVoid (
    const Char *        msg,
    Int                 arg1,
    Int                 arg2,
    const Char *        msg2 )
{
    ErrorMode( msg, arg1, arg2, (Obj)0, msg2, 'x' );
}

/****************************************************************************
**
*F  ErrorMayQuit( <msg>, <arg1>, <arg2> )  . . .  print and return
*/
void ErrorMayQuit (
    const Char *        msg,
    Int                 arg1,
    Int                 arg2)
{
    ErrorMode( msg, arg1, arg2, (Obj)0, (Char *)0, 'm' );
}


/****************************************************************************
**
*F  FuncError( <self>, <args> ) . . . . . . . . . . . . . . . signal an error
**
*/
Obj FuncError (
    Obj                 self,
    Obj                 args )
{
    return ErrorMode( (Char*)0, 0L, 0L, args, (Char*)0, 'f' );
}

/****************************************************************************
**
*F  FuncErrorCount( <self> ) . . . . . . . . . . . . .return the error count
**
*/

Obj FuncErrorCount( Obj self )
{
  return INTOBJ_INT(ErrorCount);
}

/****************************************************************************
**

*F * * * * * * * * * functions for creating the init file * * * * * * * * * *
*/



/****************************************************************************
**

*F  Complete( <list> )  . . . . . . . . . . . . . . . . . . . complete a file
*/
Obj  CompNowFuncs;
UInt CompNowCount;
Obj  CompLists;
Obj  CompThenFuncs;

#define COMP_THEN_OFFSET        2

void Complete (
    Obj                 list )
{
    Obj                 filename;
    UInt                type;
    Int4                crc;
    Int4                crc1;

    /* get the filename                                                    */
    filename = ELM_PLIST( list, 1 );

    /* and the crc value                                                   */
    crc = INT_INTOBJ( ELM_PLIST( list, 2 ) );

    /* check the crc value                                                 */
    if ( SyCheckCompletionCrcRead ) {
        crc1 = SyGAPCRC( CSTR_STRING(filename) );
        if ( crc != crc1 ) {
            ErrorQuit(
 "Error, Rebuild completion files! (Crc value of\n\"%s\" does not match.)",
                (Int)CSTR_STRING(filename), 0L );
            return;
        }
    }

    /* try to open the file                                                */
    if ( ! OpenInput( CSTR_STRING(filename) ) ) {
        return;
    }
    ClearError();

    /* switch on the buffer for faster reading */
    SySetBuffering(Input->file);

    /* we are now completing                                               */
    if ( SyDebugLoading ) {
        Pr( "#I  completing '%s'\n", (Int)CSTR_STRING(filename), 0L );
    }
    CompNowFuncs = list;
    CompNowCount = COMP_THEN_OFFSET;

    /* now do the reading                                                  */
    while ( 1 ) {
        type = ReadEvalCommand();
        if ( type == STATUS_RETURN_VAL || type == STATUS_RETURN_VOID ) {
            Pr( "'return' must not be used in file read-eval loop",
                0L, 0L );
        }
        else if ( type == STATUS_QUIT || type == STATUS_EOF ) {
            break;
        }
    }

    /* thats it for completing                                             */
    CompNowFuncs = 0;
    CompNowCount = 0;

    /* close the input file again, and return 'true'                       */
    if ( ! CloseInput() ) {
        ErrorQuit(
            "Panic: COMPLETE cannot close input, this should not happen",
            0L, 0L );
    }
    ClearError();
}


/****************************************************************************
**
*F  DoComplete<i>args( ... )  . . . . . . . . . .  handler to complete a file
*/
Obj DoComplete0args (
    Obj                 self )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_0ARGS( self );
}

Obj DoComplete1args (
    Obj                 self,
    Obj                 arg1 )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_1ARGS( self, arg1 );
}

Obj DoComplete2args (
    Obj                 self,
    Obj                 arg1,
    Obj                 arg2 )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_2ARGS( self, arg1, arg2 );
}

Obj DoComplete3args (
    Obj                 self,
    Obj                 arg1,
    Obj                 arg2,
    Obj                 arg3 )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_3ARGS( self, arg1, arg2, arg3 );
}

Obj DoComplete4args (
    Obj                 self,
    Obj                 arg1,
    Obj                 arg2,
    Obj                 arg3,
    Obj                 arg4 )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_4ARGS( self, arg1, arg2, arg3, arg4 );
}

Obj DoComplete5args (
    Obj                 self,
    Obj                 arg1,
    Obj                 arg2,
    Obj                 arg3,
    Obj                 arg4,
    Obj                 arg5 )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_5ARGS( self, arg1, arg2, arg3, arg4, arg5 );
}

Obj DoComplete6args (
    Obj                 self,
    Obj                 arg1,
    Obj                 arg2,
    Obj                 arg3,
    Obj                 arg4,
    Obj                 arg5,
    Obj                 arg6 )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_6ARGS( self, arg1, arg2, arg3, arg4, arg5, arg6 );
}

Obj DoCompleteXargs (
    Obj                 self,
    Obj                 args )
{
    COMPLETE_FUNC( self );
    if ( IS_UNCOMPLETED_FUNC(self) ) {
        ErrorQuit( "panic: completion did not define function",
                   0, 0 );
        return 0;
    }
    return CALL_XARGS( self, args );
}


/****************************************************************************
**
*F  FuncCOM_FILE( <self>, <filename>, <crc> ) . . . . . . . . .  set filename
*/
Obj FuncCOM_FILE (
    Obj                 self,
    Obj                 filename,
    Obj                 crc )
{
    Int                 len;
    StructInitInfo *    info;
    Int4                crc1;
    Int4                crc2;
    Char                result[256];
    Int                 res;


    /* check the argument                                                  */
    while ( ! IsStringConv(filename) ) {
        filename = ErrorReturnObj(
            "<filename> must be a string (not a %s)",
            (Int)TNAM_OBJ(filename), 0L,
            "you can replace <filename> via 'return <filename>;'" );
    }
    while ( ! IS_INTOBJ(crc) ) {
        crc = ErrorReturnObj(
            "<crc> must be a small integer (not a %s)",
            (Int)TNAM_OBJ(crc), 0L,
            "you can replace <crc> via 'return <crc>;'" );
    }

    /* check if have a statically or dynamically loadable module           */
    crc1 = INT_INTOBJ(crc);
    res  = SyFindOrLinkGapRootFile(CSTR_STRING(filename), crc1, result, 256, &info);

    /* not found                                                           */
    if ( res == 0 ) {
        ErrorQuit( "cannot find module or file '%s'",
                   (Int)CSTR_STRING(filename), 0L );
        return Fail;
    }

    /* dynamically linked                                                  */
    else if ( res == 1 ) {
        if ( SyDebugLoading ) {
            Pr( "#I  READ_GAP_ROOT: loading '%s' dynamically\n",
                (Int)CSTR_STRING(filename), 0L );
        }
        res  = info->initKernel(info);
	UpdateCopyFopyInfo();
        res  = res || info->initLibrary(info);
        if ( res ) {
            Pr( "#W  init functions returned non-zero exit code\n", 0L, 0L );
        }
	info->isGapRootRelative = 1;
	RecordLoadedModule(info, CSTR_STRING(filename));
        return INTOBJ_INT(1);
    }

    /* statically linked                                                   */
    else if ( res == 2 ) {
        if ( SyDebugLoading ) {
            Pr( "#I  READ_GAP_ROOT: loading '%s' statically\n",
                (Int)CSTR_STRING(filename), 0L );
        }
        res  = info->initKernel(info);
	UpdateCopyFopyInfo();
        res  = res || info->initLibrary(info);
        if ( res ) {
            Pr( "#W  init functions returned non-zero exit code\n", 0L, 0L );
        }
	info->isGapRootRelative = 1;
	RecordLoadedModule(info, CSTR_STRING(filename));
        return INTOBJ_INT(2);
    }


    /* we have to read the GAP file                                        */
    else if ( res == 3 ) {

        /* compute the crc value of the original and compare               */
        if ( SyCheckCompletionCrcComp ) {
            crc2 = SyGAPCRC(result);
            if ( crc1 != crc2 ) {
                return INTOBJ_INT(4);
            }
        }
        /*CCC filename = NEW_STRING( SyStrlen(result) );
	  SyStrncat( CSTR_STRING(filename), result, SyStrlen(result) );CCC*/
	len = SyStrlen(result);
	C_NEW_STRING(filename, len, result);

        CompThenFuncs = NEW_PLIST( T_PLIST, COMP_THEN_OFFSET );
        SET_LEN_PLIST( CompThenFuncs, COMP_THEN_OFFSET );
        SET_ELM_PLIST( CompThenFuncs, 1, filename );
        SET_ELM_PLIST( CompThenFuncs, 2, INTOBJ_INT(crc1) );

        len = LEN_PLIST( CompLists );
        GROW_PLIST(    CompLists, len+1 );
        SET_LEN_PLIST( CompLists, len+1 );
        SET_ELM_PLIST( CompLists, len+1, CompThenFuncs );
        CHANGED_BAG(   CompLists );

        return INTOBJ_INT(3);
    }

    /* we have to read the GAP file, crc mismatch                          */
    else if ( res == 4 ) {
        return INTOBJ_INT(4);
    }

    /* don't know                                                          */
    else {
        ErrorQuit( "unknown result code %d from 'SyFindGapRoot'", res, 0L );
        return Fail;
    }
}


/****************************************************************************
**
*F  FuncCOM_FUN( <self>, <num> )  . . . . . . . . make a completable function
*/
static Obj StringUncompleted;
static Obj EmptyList;

Obj FuncCOM_FUN (
    Obj                 self,
    Obj                 num )
{
    Obj                 func;
    Int                 n;

    /* if the file is not yet completed then make a new function           */
    n = INT_INTOBJ(num) + COMP_THEN_OFFSET;
    if ( LEN_PLIST( CompThenFuncs ) < n ) {

        /* make the function                                               */
        func = NewFunctionT( T_FUNCTION, SIZE_FUNC, EmptyList, -1,
                             StringUncompleted, 0 );
        HDLR_FUNC( func, 0 ) = DoComplete0args;
        HDLR_FUNC( func, 1 ) = DoComplete1args;
        HDLR_FUNC( func, 2 ) = DoComplete2args;
        HDLR_FUNC( func, 3 ) = DoComplete3args;
        HDLR_FUNC( func, 4 ) = DoComplete4args;
        HDLR_FUNC( func, 5 ) = DoComplete5args;
        HDLR_FUNC( func, 6 ) = DoComplete6args;
        HDLR_FUNC( func, 7 ) = DoCompleteXargs;
        BODY_FUNC( func )    = CompThenFuncs;

        /* add the function to the list of functions to complete           */
        GROW_PLIST(    CompThenFuncs, n );
        SET_LEN_PLIST( CompThenFuncs, n );
        SET_ELM_PLIST( CompThenFuncs, n, func );
        CHANGED_BAG(   CompThenFuncs );

    }

    /* return the function                                                 */
    return ELM_PLIST( CompThenFuncs, n );
}


/****************************************************************************
**
*F  FuncMAKE_INIT( <out>, <in>, ... ) . . . . . . . . . .  generate init file
**  XXX  This is not correct with long integers or strings which are long
**  or contain zero characters ! (FL) XXX
*/
#define MAKE_INIT_GET_SYMBOL                    \
    do {                                        \
        symbol = Symbol;                        \
        value[0] = '\0';                        \
        SyStrncat( value, Value, 1023 );        \
        if ( Symbol != S_EOF )  GetSymbol();    \
    } while (0)


Obj FuncMAKE_INIT (
    Obj                 self,
    Obj                 output,
    Obj                 filename )
{
    volatile UInt       level;
    volatile UInt       symbol;
    Char                value [1024];
    volatile UInt       funcNum;
    jmp_buf             readJmpError;

    /* check the argument                                                  */
    if ( ! IsStringConv( filename ) ) {
        ErrorQuit( "%d.th argument must be a string (not a %s)",
                   (Int)TNAM_OBJ(filename), 0L );
    }

    /* try to open the output                                              */
    if ( ! OpenAppend(CSTR_STRING(output)) ) {
        ErrorQuit( "cannot open '%s' for output",
                   (Int)CSTR_STRING(output), 0L );
    }

    /* try to open the file                                                */
    if ( ! OpenInput( CSTR_STRING(filename) ) ) {
        CloseOutput();
        ErrorQuit( "'%s' must exist and be readable",
                   (Int)CSTR_STRING(filename), 0L );
    }
    ClearError();

    /* where is this stuff                                                 */
    funcNum = 1;

    /* read the file                                                       */
    GetSymbol();
    MAKE_INIT_GET_SYMBOL;
    while ( symbol != S_EOF ) {

        memcpy( readJmpError, ReadJmpError, sizeof(jmp_buf) );
        if ( READ_ERROR() ) {
            memcpy( ReadJmpError, readJmpError, sizeof(jmp_buf) );
            CloseInput();
            CloseOutput();
            ReadEvalError();
        }
        memcpy( ReadJmpError, readJmpError, sizeof(jmp_buf) );

        /* handle function beginning and ending                            */
        if ( symbol == S_FUNCTION ) {
            Pr( "COM_FUN(%d)", funcNum++, 0L );
            MAKE_INIT_GET_SYMBOL;
            level = 0;
            while ( level != 0 || symbol != S_END ) {
                if ( symbol == S_FUNCTION )
                    level++;
                if ( symbol == S_END )
                    level--;
                MAKE_INIT_GET_SYMBOL;
            }
            MAKE_INIT_GET_SYMBOL;
        }

        /* handle -> expressions                                           */
        else if ( symbol == S_IDENT && Symbol == S_MAPTO ) {
            Pr( "COM_FUN(%d)", funcNum++, 0L );
            symbol = Symbol;  if ( Symbol != S_EOF )  GetSymbol();
            MAKE_INIT_GET_SYMBOL;
            level = 0;
            while ( level != 0
                 || (symbol != S_RBRACK  && symbol != S_RBRACE
                 && symbol != S_RPAREN  && symbol != S_COMMA
                 && symbol != S_DOTDOT  && symbol != S_SEMICOLON) )
            {
                 if ( symbol == S_LBRACK  || symbol == S_LBRACE
                   || symbol == S_LPAREN  || symbol == S_FUNCTION
                   || symbol == S_BLBRACK || symbol == S_BLBRACE )
                     level++;
                 if ( symbol == S_RBRACK  || symbol == S_RBRACE
                   || symbol == S_RPAREN  || symbol == S_END )
                     level--;
                 MAKE_INIT_GET_SYMBOL;
            }
        }

        /* handle the other symbols                                        */
        else {

            switch ( symbol ) {
            case S_IDENT:    Pr( "%I",      (Int)value, 0L );  break;
            case S_UNBIND:   Pr( "Unbind",  0L, 0L );  break;
            case S_ISBOUND:  Pr( "IsBound", 0L, 0L );  break;

            case S_LBRACK:   Pr( "[",       0L, 0L );  break;
            case S_RBRACK:   Pr( "]",       0L, 0L );  break;
            case S_LBRACE:   Pr( "{",       0L, 0L );  break;
            case S_RBRACE:   Pr( "}",       0L, 0L );  break;
            case S_DOT:      Pr( ".",       0L, 0L );  break;
            case S_LPAREN:   Pr( "(",       0L, 0L );  break;
            case S_RPAREN:   Pr( ")",       0L, 0L );  break;
            case S_COMMA:    Pr( ",",       0L, 0L );  break;
            case S_DOTDOT:   Pr( "%>..%<",  0L, 0L );  break;

            case S_BDOT:     Pr( "!.",      0L, 0L );  break;
            case S_BLBRACK:  Pr( "![",      0L, 0L );  break;
            case S_BLBRACE:  Pr( "!{",      0L, 0L );  break;

            case S_INT:      Pr( "%s",      (Int)value, 0L );  break;
            case S_TRUE:     Pr( "true",    0L, 0L );  break;
            case S_FALSE:    Pr( "false",   0L, 0L );  break;
            case S_CHAR:     Pr( "'%c'",    (Int)value[0], 0L );  break;
            case S_STRING:   Pr( "\"%S\"",  (Int)value, 0L );  break;

            case S_REC:      Pr( "rec",     0L, 0L );  break;

            case S_FUNCTION: /* handled above */       break;
            case S_LOCAL:    /* shouldn't happen */    break;
            case S_END:      /* handled above */       break;
            case S_MAPTO:    /* handled above */       break;

            case S_MULT:     Pr( "*",       0L, 0L );  break;
            case S_DIV:      Pr( "/",       0L, 0L );  break;
            case S_MOD:      Pr( " mod ",   0L, 0L );  break;
            case S_POW:      Pr( "^",       0L, 0L );  break;

            case S_PLUS:     Pr( "+",       0L, 0L );  break;
            case S_MINUS:    Pr( "-",       0L, 0L );  break;

            case S_EQ:       Pr( "=",       0L, 0L );  break;
            case S_LT:       Pr( "<",       0L, 0L );  break;
            case S_GT:       Pr( ">",       0L, 0L );  break;
            case S_NE:       Pr( "<>",      0L, 0L );  break;
            case S_LE:       Pr( "<=",      0L, 0L );  break;
            case S_GE:       Pr( ">=",      0L, 0L );  break;
            case S_IN:       Pr( " in ",    0L, 0L );  break;

            case S_NOT:      Pr( "not ",    0L, 0L );  break;
            case S_AND:      Pr( " and ",   0L, 0L );  break;
            case S_OR:       Pr( " or ",    0L, 0L );  break;

            case S_ASSIGN:   Pr( ":=",      0L, 0L );  break;

            case S_IF:       Pr( "if ",     0L, 0L );  break;
            case S_FOR:      Pr( "for ",    0L, 0L );  break;
            case S_WHILE:    Pr( "while ",  0L, 0L );  break;
            case S_REPEAT:   Pr( "repeat ", 0L, 0L );  break;

            case S_THEN:     Pr( " then\n", 0L, 0L );  break;
            case S_ELIF:     Pr( "elif ",   0L, 0L );  break;
            case S_ELSE:     Pr( "else\n",  0L, 0L );  break;
            case S_FI:       Pr( "fi",      0L, 0L );  break;
            case S_DO:       Pr( " do\n",   0L, 0L );  break;
            case S_OD:       Pr( "od",      0L, 0L );  break;
            case S_UNTIL:    Pr( "until ",  0L, 0L );  break;

            case S_BREAK:    Pr( "break",   0L, 0L );  break;
            case S_RETURN:   Pr( "return ", 0L, 0L );  break;
            case S_QUIT:     Pr( "quit",    0L, 0L );  break;

            case S_SEMICOLON: Pr( ";\n",    0L, 0L );  break;

            default: CloseInput();
                     CloseOutput();
                     ClearError();
                     ErrorQuit( "unknown symbol %d", (Int)symbol, 0L );

            }

            /* get the next symbol                                         */
            MAKE_INIT_GET_SYMBOL;
        }
    }

    /* close the input file again                                          */
    if ( ! CloseInput() ) {
        ErrorQuit(
            "Panic: MAKE_INIT cannot close input, this should not happen",
            0L, 0L );
    }
    ClearError();

    /* close the output file                                               */
    CloseOutput();

    return 0;
}


/****************************************************************************
**

*F * * * * * * * * * functions for dynamical/static modules * * * * * * * * *
*/



/****************************************************************************
**

*F  FuncGAP_CRC( <self>, <name> ) . . . . . . . create a crc value for a file
*/
Obj FuncGAP_CRC (
    Obj                 self,
    Obj                 filename )
{
    /* check the argument                                                  */
    while ( ! IsStringConv( filename ) ) {
        filename = ErrorReturnObj(
            "<filename> must be a string (not a %s)",
            (Int)TNAM_OBJ(filename), 0L,
            "you can replace <filename> via 'return <filename>;'" );
    }

    /* compute the crc value                                               */
    return INTOBJ_INT( SyGAPCRC( CSTR_STRING(filename) ) );
}


/****************************************************************************
**
*F  FuncLOAD_DYN( <self>, <name>, <crc> ) . . .  try to load a dynamic module
*/
Obj FuncLOAD_DYN (
    Obj                 self,
    Obj                 filename,
    Obj                 crc )
{
    InitInfoFunc        init;
    StructInitInfo *    info;
    Obj                 crc1;
    Int                 res;

    /* check the argument                                                  */
    while ( ! IsStringConv( filename ) ) {
        filename = ErrorReturnObj(
            "<filename> must be a string (not a %s)",
            (Int)TNAM_OBJ(filename), 0L,
            "you can replace <filename> via 'return <filename>;'" );
    }
    while ( ! IS_INTOBJ(crc) && crc!=False ) {
        crc = ErrorReturnObj(
            "<crc> must be a small integer or 'false' (not a %s)",
            (Int)TNAM_OBJ(crc), 0L,
            "you can replace <crc> via 'return <crc>;'" );
    }

    /* try to read the module                                              */
    init = SyLoadModule( CSTR_STRING(filename) );
    if ( (Int)init == 1 )
        ErrorQuit( "module '%s' not found", (Int)CSTR_STRING(filename), 0L );
    else if ( (Int) init == 3 )
        ErrorQuit( "symbol 'Init_Dynamic' not found", 0L, 0L );
    else if ( (Int) init == 5 )
        ErrorQuit( "forget symbol failed", 0L, 0L );

    /* no dynamic library support                                          */
    else if ( (Int) init == 7 ) {
        if ( SyDebugLoading ) {
            Pr( "#I  LOAD_DYN: no support for dynamical loading\n", 0L, 0L );
        }
        return False;
    }

    /* get the description structure                                       */
    info = (*init)();
    if ( info == 0 )
        ErrorQuit( "call to init function failed", 0L, 0L );

    /* check the crc value                                                 */
    if ( crc != False ) {
        crc1 = INTOBJ_INT( info->crc );
        if ( ! EQ( crc, crc1 ) ) {
            if ( SyDebugLoading ) {
                Pr( "#I  LOAD_DYN: crc values do not match, gap ", 0L, 0L );
                PrintInt( crc );
                Pr( ", dyn ", 0L, 0L );
                PrintInt( crc1 );
                Pr( "\n", 0L, 0L );
            }
            return False;
        }
    }

    /* link and init me                                                    */
    info->isGapRootRelative = 0;
    res = (info->initKernel)(info);
    UpdateCopyFopyInfo();

    /* Start a new executor to run the outer function of the module
       in global context */
    ExecBegin( BottomLVars );
    res = res || (info->initLibrary)(info);
    ExecEnd(res ? STATUS_ERROR : STATUS_END);
    if ( res ) {
        Pr( "#W  init functions returned non-zero exit code\n", 0L, 0L );
    }
    RecordLoadedModule(info, CSTR_STRING(filename));

    return True;
}


/****************************************************************************
**
*F  FuncLOAD_STAT( <self>, <name>, <crc> )  . . . . try to load static module
*/
Obj FuncLOAD_STAT (
    Obj                 self,
    Obj                 filename,
    Obj                 crc )
{
    StructInitInfo *    info = 0;
    Obj                 crc1;
    Int                 k;
    Int                 res;

    /* check the argument                                                  */
    while ( ! IsStringConv( filename ) ) {
        filename = ErrorReturnObj(
            "<filename> must be a string (not a %s)",
            (Int)TNAM_OBJ(filename), 0L,
            "you can replace <filename> via 'return <filename>;'" );
    }
    while ( !IS_INTOBJ(crc) && crc!=False ) {
        crc = ErrorReturnObj(
            "<crc> must be a small integer or 'false' (not a %s)",
            (Int)TNAM_OBJ(crc), 0L,
            "you can replace <crc> via 'return <crc>;'" );
    }

    /* try to find the module                                              */
    for ( k = 0;  CompInitFuncs[k];  k++ ) {
        info = (*(CompInitFuncs[k]))();
        if ( info == 0 ) {
            continue;
        }
        if ( ! SyStrcmp( CSTR_STRING(filename), info->name ) ) {
            break;
        }
    }
    if ( CompInitFuncs[k] == 0 ) {
        if ( SyDebugLoading ) {
            Pr( "#I  LOAD_STAT: no module named '%s' found\n",
                (Int)CSTR_STRING(filename), 0L );
        }
        return False;
    }

    /* check the crc value                                                 */
    if ( crc != False ) {
        crc1 = INTOBJ_INT( info->crc );
        if ( ! EQ( crc, crc1 ) ) {
            if ( SyDebugLoading ) {
                Pr( "#I  LOAD_STAT: crc values do not match, gap ", 0L, 0L );
                PrintInt( crc );
                Pr( ", stat ", 0L, 0L );
                PrintInt( crc1 );
                Pr( "\n", 0L, 0L );
            }
            return False;
        }
    }

    /* link and init me                                                    */
    info->isGapRootRelative = 0;
    res = (info->initKernel)(info);
    UpdateCopyFopyInfo();
    /* Start a new executor to run the outer function of the module
       in global context */
    ExecBegin( BottomLVars );
    res = res || (info->initLibrary)(info);
    ExecEnd(res ? STATUS_ERROR : STATUS_END);
    if ( res ) {
        Pr( "#W  init functions returned non-zero exit code\n", 0L, 0L );
    }
    RecordLoadedModule(info, CSTR_STRING(filename));

    return True;
}


/****************************************************************************
**
*F  FuncSHOW_STAT() . . . . . . . . . . . . . . . . . . . show static modules
*/
Obj FuncSHOW_STAT (
    Obj                 self )
{
    Obj                 modules;
    Obj                 name;
    StructInitInfo *    info;
    Int                 k;
    Int                 im;
    Int                 len;

    /* count the number of install modules                                 */
    for ( k = 0,  im = 0;  CompInitFuncs[k];  k++ ) {
        info = (*(CompInitFuncs[k]))();
        if ( info == 0 ) {
            continue;
        }
        im++;
    }

    /* make a list of modules with crc values                              */
    modules = NEW_PLIST( T_PLIST, 2*im );
    SET_LEN_PLIST( modules, 2*im );

    for ( k = 0,  im = 1;  CompInitFuncs[k];  k++ ) {
        info = (*(CompInitFuncs[k]))();
        if ( info == 0 ) {
            continue;
        }
        /*CCC name = NEW_STRING( SyStrlen(info->name) );
	  SyStrncat( CSTR_STRING(name), info->name, SyStrlen(info->name) );CCC*/
	len = SyStrlen(info->name);
	C_NEW_STRING(name, len, info->name);

        SET_ELM_PLIST( modules, im, name );

        /* compute the crc value                                           */
        SET_ELM_PLIST( modules, im+1, INTOBJ_INT( info->crc ) );
        im += 2;
    }

    return modules;
}


/****************************************************************************
**
*F  FuncLoadedModules( <self> ) . . . . . . . . . . . list all loaded modules
*/
Obj FuncLoadedModules (
    Obj                 self )
{
    Int                 i;
    StructInitInfo *    m;
    Obj                 str;
    Obj                 list;

    /* create a list                                                       */
    list = NEW_PLIST( T_PLIST, NrModules * 3 );
    SET_LEN_PLIST( list, NrModules * 3 );
    for ( i = 0;  i < NrModules;  i++ ) {
        m = Modules[i];
        if ( m->type == MODULE_BUILTIN ) {
            SET_ELM_PLIST( list, 3*i+1, ObjsChar[(Int)'b'] );
	    CHANGED_BAG(list);
            C_NEW_STRING( str, SyStrlen(m->name), m->name );
            SET_ELM_PLIST( list, 3*i+2, str );
            SET_ELM_PLIST( list, 3*i+3, INTOBJ_INT(m->version) );
        }
        else if ( m->type == MODULE_DYNAMIC ) {
            SET_ELM_PLIST( list, 3*i+1, ObjsChar[(Int)'d'] );
	    CHANGED_BAG(list);
            C_NEW_STRING( str, SyStrlen(m->name), m->name );
            SET_ELM_PLIST( list, 3*i+2, str );
	    CHANGED_BAG(list);
            C_NEW_STRING( str, SyStrlen(m->filename), m->filename );
            SET_ELM_PLIST( list, 3*i+3, str );
        }
        else if ( m->type == MODULE_STATIC ) {
            SET_ELM_PLIST( list, 3*i+1, ObjsChar[(Int)'s'] );
	    CHANGED_BAG(list);
            C_NEW_STRING( str, SyStrlen(m->name), m->name );
            SET_ELM_PLIST( list, 3*i+2, str );
	    CHANGED_BAG(list);
            C_NEW_STRING( str, SyStrlen(m->filename), m->filename );
            SET_ELM_PLIST( list, 3*i+3, str );
        }
    }
    return CopyObj( list, 0 );
}


/****************************************************************************
**


*F * * * * * * * * * * * * * * debug functions  * * * * * * * * * * * * * * *
*/

/****************************************************************************
**

*F  FuncGASMAN( <self>, <args> )  . . . . . . . . .  expert function 'GASMAN'
**
**  'FuncGASMAN' implements the internal function 'GASMAN'
**
**  'GASMAN( "display" | "clear" | "collect" | "message" | "partial" )'
*/
Obj FuncGASMAN (
    Obj                 self,
    Obj                 args )
{
    Obj                 cmd;            /* argument                        */
    UInt                i,  k;          /* loop variables                  */
    Char                buf[100];

    /* check the argument                                                  */
    while ( ! IS_SMALL_LIST(args) || LEN_LIST(args) == 0 ) {
        args = ErrorReturnObj(
            "usage: GASMAN( \"display\"|\"displayshort\"|\"clear\"|\"collect\"|\"message\"|\"partial\" )",
            0L, 0L,
            "you can replace the argument list <args> via 'return <args>;'" );
    }

    /* loop over the arguments                                             */
    for ( i = 1; i <= LEN_LIST(args); i++ ) {

        /* evaluate and check the command                                  */
        cmd = ELM_PLIST( args, i );
again:
        while ( ! IsStringConv(cmd) ) {
           cmd = ErrorReturnObj(
               "GASMAN: <cmd> must be a string (not a %s)",
               (Int)TNAM_OBJ(cmd), 0L,
               "you can replace <cmd> via 'return <cmd>;'" );
       }

        /* if request display the statistics                               */
        if ( SyStrcmp( CSTR_STRING(cmd), "display" ) == 0 ) {
            Pr( "%40s ", (Int)"type",  0L          );
            Pr( "%8s %8s ",  (Int)"alive", (Int)"kbyte" );
            Pr( "%8s %8s\n",  (Int)"total", (Int)"kbyte" );
            for ( k = 0; k < 256; k++ ) {
                if ( InfoBags[k].name != 0 ) {
                    buf[0] = '\0';
                    SyStrncat( buf, InfoBags[k].name, 40 );
                    Pr("%40s ",    (Int)buf, 0L );
                    Pr("%8d %8d ", (Int)InfoBags[k].nrLive,
                                   (Int)(InfoBags[k].sizeLive/1024));
                    Pr("%8d %8d\n",(Int)InfoBags[k].nrAll,
                                   (Int)(InfoBags[k].sizeAll/1024));
                }
            }
        }

        /* if request give a short display of the statistics                */
        if ( SyStrcmp( CSTR_STRING(cmd), "displayshort" ) == 0 ) {
            Pr( "%40s ", (Int)"type",  0L          );
            Pr( "%8s %8s ",  (Int)"alive", (Int)"kbyte" );
            Pr( "%8s %8s\n",  (Int)"total", (Int)"kbyte" );
            for ( k = 0; k < 256; k++ ) {
                if ( InfoBags[k].name != 0 &&
                     (InfoBags[k].nrLive != 0 ||
                      InfoBags[k].sizeLive != 0 ||
                      InfoBags[k].nrAll != 0 ||
                      InfoBags[k].sizeAll != 0) ) {
                    buf[0] = '\0';
                    SyStrncat( buf, InfoBags[k].name, 40 );
                    Pr("%40s ",    (Int)buf, 0L );
                    Pr("%8d %8d ", (Int)InfoBags[k].nrLive,
                                   (Int)(InfoBags[k].sizeLive/1024));
                    Pr("%8d %8d\n",(Int)InfoBags[k].nrAll,
                                   (Int)(InfoBags[k].sizeAll/1024));
                }
            }
        }

        /* if request display the statistics                               */
        else if ( SyStrcmp( CSTR_STRING(cmd), "clear" ) == 0 ) {
            for ( k = 0; k < 256; k++ ) {
#ifdef GASMAN_CLEAR_TO_LIVE
                InfoBags[k].nrAll    = InfoBags[k].nrLive;
                InfoBags[k].sizeAll  = InfoBags[k].sizeLive;
#else
                InfoBags[k].nrAll    = 0;
                InfoBags[k].sizeAll  = 0;
#endif
            }
        }

        /* or collect the garbage                                          */
        else if ( SyStrcmp( CSTR_STRING(cmd), "collect" ) == 0 ) {
            CollectBags(0,1);
        }

        /* or collect the garbage                                          */
        else if ( SyStrcmp( CSTR_STRING(cmd), "partial" ) == 0 ) {
            CollectBags(0,0);
        }

        /* or display information about global bags                        */
        else if ( SyStrcmp( CSTR_STRING(cmd), "global" ) == 0 ) {
            for ( i = 0;  i < GlobalBags.nr;  i++ ) {
                if ( *(GlobalBags.addr[i]) != 0 ) {
                    Pr( "%50s: %12d bytes\n", (Int)GlobalBags.cookie[i],
                        (Int)SIZE_BAG(*(GlobalBags.addr[i])) );
                }
            }
        }

        /* or finally toggle Gasman messages                               */
        else if ( SyStrcmp( CSTR_STRING(cmd), "message" ) == 0 ) {
            SyMsgsFlagBags = (SyMsgsFlagBags + 1) % 3;
        }

        /* otherwise complain                                              */
        else {
            cmd = ErrorReturnObj(
                "GASMAN: <cmd> must be %s or %s",
                (Int)"\"display\" or \"clear\" or \"global\" or ",
                (Int)"\"collect\" or \"partial\" or \"message\"",
                "you can replace <cmd> via 'return <cmd>;'" );
            goto again;
        }
    }

    /* return nothing, this function is a procedure                        */
    return 0;
}

Obj FuncGASMAN_STATS(Obj self)
{
  Obj res;
  Obj row;
  Obj entry;
  UInt i,j;
  Int x;
  res = NEW_PLIST(T_PLIST_TAB_RECT + IMMUTABLE, 2);
  SET_LEN_PLIST(res, 2);
  for (i = 1; i <= 2; i++)
    {
      row = NEW_PLIST(T_PLIST_CYC + IMMUTABLE, 6);
      SET_ELM_PLIST(res, i, row);
      CHANGED_BAG(res);
      SET_LEN_PLIST(row, 6);
      for (j = 1; j <= 6; j++)
	{
	  x = SyGasmanNumbers[i-1][j];

	  /* convert x to GAP integer. x may be too big to be a small int */
	  if (x < (1L << NR_SMALL_INT_BITS))
	    entry = INTOBJ_INT(x);
	  else
	    entry = SUM( PROD(INTOBJ_INT(x >> (NR_SMALL_INT_BITS/2)),
			      INTOBJ_INT(1 << (NR_SMALL_INT_BITS/2))),
			 INTOBJ_INT( x % ( 1 << (NR_SMALL_INT_BITS/2))));
	  SET_ELM_PLIST(row, j, entry);
	}
    }
  return res;
}

Obj FuncGASMAN_MESSAGE_STATUS( Obj self )
{
  return INTOBJ_INT(SyMsgsFlagBags);
}

Obj FuncGASMAN_LIMITS( Obj self )
{
  Obj list;
  list = NEW_PLIST(T_PLIST_CYC+IMMUTABLE, 3);
  SET_LEN_PLIST(list,3);
  SET_ELM_PLIST(list, 1, INTOBJ_INT(SyStorMin));
  SET_ELM_PLIST(list, 2, INTOBJ_INT(SyStorMax));
  SET_ELM_PLIST(list, 3, INTOBJ_INT(SyStorKill));
  return list;
}

/****************************************************************************
**
*F  FuncSHALLOW_SIZE( <self>, <obj> ) . . . .  expert function 'SHALLOW_SIZE'
*/
Obj FuncSHALLOW_SIZE (
    Obj                 self,
    Obj                 obj )
{
  if (IS_INTOBJ(obj) || IS_FFE(obj))
    return INTOBJ_INT(0);
  else
    return INTOBJ_INT( SIZE_BAG( obj ) );
}


/****************************************************************************
**
*F  FuncTNUM_OBJ( <self>, <obj> ) . . . . . . . .  expert function 'TNUM_OBJ'
*/

Obj FuncTNUM_OBJ (
    Obj                 self,
    Obj                 obj )
{
    Obj                 res;
    Obj                 str;
    Int                 len;
    const Char *        cst;

    res = NEW_PLIST( T_PLIST, 2 );
    SET_LEN_PLIST( res, 2 );

    /* set the type                                                        */
    SET_ELM_PLIST( res, 1, INTOBJ_INT( TNUM_OBJ(obj) ) );
    cst = TNAM_OBJ(obj);
    /*CCC    str = NEW_STRING( SyStrlen(cst) );
      SyStrncat( CSTR_STRING(str), cst, SyStrlen(cst) );CCC*/
    len = SyStrlen(cst);
    C_NEW_STRING(str, len, cst);
    SET_ELM_PLIST( res, 2, str );

    /* and return                                                          */
    return res;
}

Obj FuncTNUM_OBJ_INT (
    Obj                 self,
    Obj                 obj )
{


    return  INTOBJ_INT( TNUM_OBJ(obj) ) ;
}

/****************************************************************************
**
*F  FuncXTNUM_OBJ( <self>, <obj> )  . . . . . . . expert function 'XTNUM_OBJ'
*/
Obj FuncXTNUM_OBJ (
    Obj                 self,
    Obj                 obj )
{
    Obj                 res;
    Obj                 str;

    res = NEW_PLIST( T_PLIST, 2 );
    SET_LEN_PLIST( res, 2 );
    SET_ELM_PLIST( res, 1, Fail );
    C_NEW_STRING(str, 16, "xtnums abolished");
    SET_ELM_PLIST(res, 2,str);
    /* and return                                                          */
    return res;
}


/****************************************************************************
**
*F  FuncOBJ_HANDLE( <self>, <obj> ) . . . . . .  expert function 'OBJ_HANDLE'
*/
Obj FuncOBJ_HANDLE (
    Obj                 self,
    Obj                 obj )
{
    UInt                hand;
    UInt                prod;
    Obj                 rem;

    if ( IS_INTOBJ(obj) ) {
        return (Obj)INT_INTOBJ(obj);
    }
    else if ( TNUM_OBJ(obj) == T_INTPOS ) {
        hand = 0;
        prod = 1;
        while ( EQ( obj, INTOBJ_INT(0) ) == 0 ) {
            rem  = RemInt( obj, INTOBJ_INT( 1 << 16 ) );
            obj  = QuoInt( obj, INTOBJ_INT( 1 << 16 ) );
            hand = hand + prod * INT_INTOBJ(rem);
            prod = prod * ( 1 << 16 );
        }
        return (Obj) hand;
    }
    else {
        ErrorQuit( "<handle> must be a positive integer", 0L, 0L );
        return 0;
    }
}


/****************************************************************************
**
*F  FuncHANDLE_OBJ( <self>, <obj> ) . . . . . .  expert function 'HANDLE_OBJ'
*/
Obj FuncHANDLE_OBJ (
    Obj                 self,
    Obj                 obj )
{
    Obj                 hnum;
    Obj                 prod;
    Obj                 tmp;
    UInt                hand;

    hand = (UInt) obj;
    hnum = INTOBJ_INT(0);
    prod = INTOBJ_INT(1);
    while ( 0 < hand ) {
        tmp  = PROD( prod, INTOBJ_INT( hand & 0xffff ) );
        prod = PROD( prod, INTOBJ_INT( 1 << 16 ) );
        hnum = SUM(  hnum, tmp );
        hand = hand >> 16;
    }
    return hnum;
}

Obj FuncMASTER_POINTER_NUMBER(Obj self, Obj o)
{
    if ((void **) o >= (void **) MptrBags && (void **) o < (void **) OldBags) {
        return INTOBJ_INT( ((void **) o - (void **) MptrBags) + 1 );
    } else {
        return INTOBJ_INT( 0 );
    }
}

Obj FuncFUNC_BODY_SIZE(Obj self, Obj f)
{
    Obj body;
    if (TNUM_OBJ(f) != T_FUNCTION) return Fail;
    body = BODY_FUNC(f);
    if (body == 0) return INTOBJ_INT(0);
    else return INTOBJ_INT( SIZE_BAG( body ) );
}

/****************************************************************************
**
*F  FuncSWAP_MPTR( <self>, <obj1>, <obj2> ) . . . . . . . swap master pointer
**
**  Never use this function unless you are debugging.
*/
Obj FuncSWAP_MPTR (
    Obj                 self,
    Obj                 obj1,
    Obj                 obj2 )
{
    if ( TNUM_OBJ(obj1) == T_INT || TNUM_OBJ(obj1) == T_FFE ) {
        ErrorQuit("SWAP_MPTR: <obj1> must not be an integer or ffe", 0L, 0L);
        return 0;
    }
    if ( TNUM_OBJ(obj2) == T_INT || TNUM_OBJ(obj2) == T_FFE ) {
        ErrorQuit("SWAP_MPTR: <obj2> must not be an integer or ffe", 0L, 0L);
        return 0;
    }

    SwapMasterPoint( obj1, obj2 );
    return 0;
}


/****************************************************************************
**

*F * * * * * * * * * * * * * initialize package * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*F  FillInVersion( <module>, <rev_c>, <rev_h> ) . . .  fill in version number
*/
static UInt ExtractRevision (
    const Char *                rev,
    const Char * *              name )
{
    const Char *                p;
    const Char *                major;
    const Char *                minor;
    UInt                        ver1;
    UInt                        ver2;

    /* store the revision strings                                          */

    /* the revision string is "@(#)Id: filename.x,v major.minor ..."       */
    p = rev;
    while ( *p && *p != ':' )  p++;
    if ( *p )  p++;
    while ( *p && *p == ' ' )  p++;
    *name = p;
    while ( *p && *p != ' ' )  p++;
    while ( *p && *p == ' ' )  p++;
    major = p;
    while ( *p && *p != '.' )  p++;
    if ( *p )  p++;
    while ( *p && *p == '.' )  p++;
    minor = p;

    /* the version is MMmmm, that is 2 digits major, 3 digits minor        */
    ver1 = 0;
    while ( '0' <= *major && *major <= '9' ) {
        ver1 = ver1 * 10 + (UInt)( *major - '0' );
        major++;
    }
    ver2 = 0;
    while ( '0' <= *minor && *minor <= '9' ) {
        ver2 = ver2 * 10 + (UInt)( *minor - '0' );
        minor++;
    }

    return ver1 * 1000 + ver2;
}


void FillInVersion (
    StructInitInfo *            module )
{
    const Char *                p;
    const Char *                q;
    const Char *                name;
    const Char *                rev_c;
    const Char *                rev_h;
    UInt                        c_ver;
    UInt                        h_ver;

    /* store revision entries                                              */
    rev_c = module->revision_c;
    rev_h = module->revision_h;

    /* extract the filename and version entry from <rev_c>                 */
    c_ver = ExtractRevision( rev_c, &name );
    if ( module->name ) {
        p = name;
        q = module->name;
        while ( *p && *q && *p == *q ) { p++; q++; }
        if ( *q || *p != '.' ) {
            FPUTS_TO_STDERR( "#W  corrupt version info '" );
            FPUTS_TO_STDERR( rev_c );
            FPUTS_TO_STDERR( "'\n" );
        }
    }
    h_ver = ExtractRevision( rev_h, &name );
    if ( module->name ) {
        p = name;
        q = module->name;
        while ( *p && *q && *p == *q ) { p++; q++; }
        if ( *q || *p != '.' ) {
            FPUTS_TO_STDERR( "#W  corrupt version info '" );
            FPUTS_TO_STDERR( rev_h );
            FPUTS_TO_STDERR( "'\n" );
        }
    }
    module->version = c_ver*100000+h_ver;
}


/****************************************************************************
**
*F  RequireModule( <calling>, <required>, <version> ) . . . .  require module
*/
void RequireModule (
    StructInitInfo *            module,
    const Char *                required,
    UInt                        version )
{
}


/****************************************************************************
**
*F  InitBagNamesFromTable( <table> )  . . . . . . . . .  initialise bag names
*/
void InitBagNamesFromTable (
    StructBagNames *            tab )
{
    Int                         i;

    for ( i = 0;  tab[i].tnum != -1;  i++ ) {
        InfoBags[tab[i].tnum].name = tab[i].name;
    }
}


/****************************************************************************
**
*F  InitClearFiltsTNumsFromTable( <tab> ) . . .  initialise clear filts tnums
*/
void InitClearFiltsTNumsFromTable (
    Int *               tab )
{
    Int                 i;

    for ( i = 0;  tab[i] != -1;  i += 2 ) {
        ClearFiltsTNums[tab[i]] = tab[i+1];
    }
}


/****************************************************************************
**
*F  InitHasFiltListTNumsFromTable( <tab> )  . . initialise tester filts tnums
*/
void InitHasFiltListTNumsFromTable (
    Int *               tab )
{
    Int                 i;

    for ( i = 0;  tab[i] != -1;  i += 3 ) {
        HasFiltListTNums[tab[i]][tab[i+1]] = tab[i+2];
    }
}


/****************************************************************************
**
*F  InitSetFiltListTNumsFromTable( <tab> )  . . initialise setter filts tnums
*/
void InitSetFiltListTNumsFromTable (
    Int *               tab )
{
    Int                 i;

    for ( i = 0;  tab[i] != -1;  i += 3 ) {
        SetFiltListTNums[tab[i]][tab[i+1]] = tab[i+2];
    }
}


/****************************************************************************
**
*F  InitResetFiltListTNumsFromTable( <tab> )  initialise unsetter filts tnums
*/
void InitResetFiltListTNumsFromTable (
    Int *               tab )
{
    Int                 i;

    for ( i = 0;  tab[i] != -1;  i += 3 ) {
        ResetFiltListTNums[tab[i]][tab[i+1]] = tab[i+2];
    }
}


/****************************************************************************
**
*F  InitGVarFiltsFromTable( <tab> ) . . . . . . . . . . . . . . . new filters
*/
void InitGVarFiltsFromTable (
    StructGVarFilt *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        AssGVar( GVarName( tab[i].name ),
            NewFilterC( tab[i].name, 1, tab[i].argument, tab[i].handler ) );
        MakeReadOnlyGVar( GVarName( tab[i].name ) );
    }
}


/****************************************************************************
**
*F  InitGVarAttrsFromTable( <tab> ) . . . . . . . . . . . . .  new attributes
*/
void InitGVarAttrsFromTable (
    StructGVarAttr *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
       AssGVar( GVarName( tab[i].name ),
         NewAttributeC( tab[i].name, 1, tab[i].argument, tab[i].handler ) );
       MakeReadOnlyGVar( GVarName( tab[i].name ) );
    }
}


/****************************************************************************
**
*F  InitGVarPropsFromTable( <tab> ) . . . . . . . . . . . . .  new properties
*/
void InitGVarPropsFromTable (
    StructGVarProp *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
       AssGVar( GVarName( tab[i].name ),
         NewPropertyC( tab[i].name, 1, tab[i].argument, tab[i].handler ) );
       MakeReadOnlyGVar( GVarName( tab[i].name ) );
    }
}


/****************************************************************************
**
*F  InitGVarOpersFromTable( <tab> ) . . . . . . . . . . . . .  new operations
*/
void InitGVarOpersFromTable (
    StructGVarOper *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        AssGVar( GVarName( tab[i].name ), NewOperationC( tab[i].name,
            tab[i].nargs, tab[i].args, tab[i].handler ) );
        MakeReadOnlyGVar( GVarName( tab[i].name ) );
    }
}


/****************************************************************************
**
*F  InitGVarFuncsFromTable( <tab> ) . . . . . . . . . . . . . . new functions
*/
void InitGVarFuncsFromTable (
    StructGVarFunc *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        AssGVar( GVarName( tab[i].name ), NewFunctionC( tab[i].name,
            tab[i].nargs, tab[i].args, tab[i].handler ) );
        MakeReadOnlyGVar( GVarName( tab[i].name ) );
    }
}


/****************************************************************************
**
*F  InitHdlrFiltsFromTable( <tab> ) . . . . . . . . . . . . . . . new filters
*/
void InitHdlrFiltsFromTable (
    StructGVarFilt *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        InitHandlerFunc( tab[i].handler, tab[i].cookie );
        InitFopyGVar( tab[i].name, tab[i].filter );
    }
}


/****************************************************************************
**
*F  InitHdlrAttrsFromTable( <tab> ) . . . . . . . . . . . . .  new attributes
*/
void InitHdlrAttrsFromTable (
    StructGVarAttr *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        InitHandlerFunc( tab[i].handler, tab[i].cookie );
        InitFopyGVar( tab[i].name, tab[i].attribute );
    }
}


/****************************************************************************
**
*F  InitHdlrPropsFromTable( <tab> ) . . . . . . . . . . . . .  new properties
*/
void InitHdlrPropsFromTable (
    StructGVarProp *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        InitHandlerFunc( tab[i].handler, tab[i].cookie );
        InitFopyGVar( tab[i].name, tab[i].property );
    }
}


/****************************************************************************
**
*F  InitHdlrOpersFromTable( <tab> ) . . . . . . . . . . . . .  new operations
*/
void InitHdlrOpersFromTable (
    StructGVarOper *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        InitHandlerFunc( tab[i].handler, tab[i].cookie );
        InitFopyGVar( tab[i].name, tab[i].operation );
    }
}


/****************************************************************************
**
*F  InitHdlrFuncsFromTable( <tab> ) . . . . . . . . . . . . . . new functions
*/
void InitHdlrFuncsFromTable (
    StructGVarFunc *    tab )
{
    Int                 i;

    for ( i = 0;  tab[i].name != 0;  i++ ) {
        InitHandlerFunc( tab[i].handler, tab[i].cookie );
    }
}


/****************************************************************************
**
*F  ImportGVarFromLibrary( <name>, <address> )  . . .  import global variable
*/
typedef struct {
    const Char *                name;
    Obj *                       address;
} StructImportedGVars;

#ifndef MAX_IMPORTED_GVARS
#define MAX_IMPORTED_GVARS      1024
#endif

static StructImportedGVars ImportedGVars[MAX_IMPORTED_GVARS];
static Int NrImportedGVars = 0;

void ImportGVarFromLibrary(
    const Char *        name,
    Obj *               address )
{
    if ( NrImportedGVars == 1024 ) {
        Pr( "#W  warning: too many imported GVars\n", 0L, 0L );
    }
    else {
        ImportedGVars[NrImportedGVars].name    = name;
        ImportedGVars[NrImportedGVars].address = address;
        NrImportedGVars++;
    }
    if ( address != 0 ) {
        InitCopyGVar( name, address );
    }
}


/****************************************************************************
**
*F  ImportFuncFromLibrary( <name>, <address> )  . . .  import global function
*/
static StructImportedGVars ImportedFuncs[MAX_IMPORTED_GVARS];
static Int NrImportedFuncs = 0;


void ImportFuncFromLibrary(
    const Char *        name,
    Obj *               address )
{
    if ( NrImportedFuncs == 1024 ) {
        Pr( "#W  warning: too many imported Funcs\n", 0L, 0L );
    }
    else {
        ImportedFuncs[NrImportedFuncs].name    = name;
        ImportedFuncs[NrImportedFuncs].address = address;
        NrImportedFuncs++;
    }
    if ( address != 0 ) {
        InitFopyGVar( name, address );
    }
}


/****************************************************************************
**
*F  FuncExportToKernelFinished( <self> )  . . . . . . . . . . check functions
*/
Obj FuncExportToKernelFinished (
    Obj             self )
{
    UInt            i;
    Int             errs = 0;
    Obj             val;

    for ( i = 0;  i < NrImportedGVars;  i++ ) {
        if ( ImportedGVars[i].address == 0 ) {
            val = ValAutoGVar(GVarName(ImportedGVars[i].name));
            if ( val == 0 ) {
                errs++;
                if ( ! SyQuiet ) {
                    Pr( "#W  global variable '%s' has not been defined\n",
                        (Int)ImportedFuncs[i].name, 0L );
                }
            }
        }
        else if ( *ImportedGVars[i].address == 0 ) {
            errs++;
            if ( ! SyQuiet ) {
                Pr( "#W  global variable '%s' has not been defined\n",
                    (Int)ImportedGVars[i].name, 0L );
            }
        }
        else {
            MakeReadOnlyGVar(GVarName(ImportedGVars[i].name));
        }
    }

    for ( i = 0;  i < NrImportedFuncs;  i++ ) {
        if (  ImportedFuncs[i].address == 0 ) {
            val = ValAutoGVar(GVarName(ImportedFuncs[i].name));
            if ( val == 0 || ! IS_FUNC(val) ) {
                errs++;
                if ( ! SyQuiet ) {
                    Pr( "#W  global function '%s' has not been defined\n",
                        (Int)ImportedFuncs[i].name, 0L );
                }
            }
        }
        else if ( *ImportedFuncs[i].address == ErrorMustEvalToFuncFunc
          || *ImportedFuncs[i].address == ErrorMustHaveAssObjFunc )
        {
            errs++;
            if ( ! SyQuiet ) {
                Pr( "#W  global function '%s' has not been defined\n",
                    (Int)ImportedFuncs[i].name, 0L );
            }
        }
        else {
            MakeReadOnlyGVar(GVarName(ImportedFuncs[i].name));
        }
    }

    return errs == 0 ? True : False;
}


/****************************************************************************
**
*F  FuncSleep( <self>, <secs> )
**
*/

Obj FuncSleep( Obj self, Obj secs )
{
  Int  s;

  while( ! IS_INTOBJ(secs) )
    secs = ErrorReturnObj( "<secs> must be a small integer", 0L, 0L,
                           "you can replace <secs> via 'return <secs>;'" );


  if ( (s = INT_INTOBJ(secs)) > 0)
    SySleep((UInt)s);

  /* either we used up the time, or we were interrupted. */
  if (SyIsIntr())
    ErrorReturnVoid("user interrupt in sleep", 0L, 0L,
		    "you can 'return;' as if the sleep was finished");

  return (Obj) 0;
}

/****************************************************************************
**
*F  FuncQUIT_GAP()
**
*/

Obj FuncQUIT_GAP( Obj self )
{
  UserHasQUIT = 1;
  ReadEvalError();
  return (Obj)0;
}

/****************************************************************************
**
*F  MakeOptionsRecord()
**
**  Assemble a GAP readable record of all the command-line options
**
*/

Obj MakeOptionsRecord( void )
{
  Obj optrec, list;
  UInt i,j,optcount,len;
  Char name[2], opt;
  Char *intarg;
  Obj thearg;

  /* export the command line arguments                                   */
  optrec = NEW_PREC(255);
  j = 1;
  name[1] = '\0';
  for (opt = '\1'; opt != '\177'; opt++)
    {
      optcount = getOptionCount( opt );
      if (optcount != 0)
	{
	  name[0] = opt;
	  SET_RNAM_PREC( optrec, j, RNamName(name));
	  list = NEW_PLIST(T_PLIST+IMMUTABLE, optcount);

	  for (i = 0; i < optcount; i++)
	    {
	      intarg = getOptionArg(opt, i);
	      if (intarg == NULL)
		SET_ELM_PLIST(list, i+1, False);
	      else
		{
		  len = SyStrlen(intarg);
		  thearg = NEW_STRING(len);
		  SyStrncat(CSTR_STRING(thearg), intarg, len);
		  SET_ELM_PLIST(list,i+1, thearg);
		  CHANGED_BAG(list);
		}
	    }
	  SET_LEN_PLIST(list, optcount);
	  SET_ELM_PREC( optrec, j++, list);
	  CHANGED_BAG(optrec);
	}
    }
  ResizeBag(optrec, 2*sizeof(Obj)*j);
  RetypeBag(optrec, T_PREC+IMMUTABLE);
  return optrec;
}

/****************************************************************************
**
*V  Revisions . . . . . . . . . . . . . . . . . .  record of revision numbers
*/
Obj Revisions;


/****************************************************************************
**
*V  GVarFuncs . . . . . . . . . . . . . . . . . . list of functions to export
*/
static StructGVarFunc GVarFuncs [] = {

    { "Runtime", 0, "",
      FuncRuntime, "src/gap.c:Runtime" },

    { "RUNTIMES", 0, "",
      FuncRUNTIMES, "src/gap.c:RUNTIMES" },

    { "SizeScreen", -1, "args",
      FuncSizeScreen, "src/gap.c:SizeScreen" },

    { "ID_FUNC", 1, "object",
      FuncID_FUNC, "src/gap.c:ID_FUNC" },

    { "ExportToKernelFinished", 0, "",
      FuncExportToKernelFinished, "src/gap.c:ExportToKernelFinished" },

    { "DownEnv", -1, "args",
      FuncDownEnv, "src/gap.c:DownEnv" },

    { "UpEnv", -1, "args",
      FuncUpEnv, "src/gap.c:UpEnv" },

    { "Where", -1, "args",
      FuncWhere, "src/gap.c:Where" },

    { "Error", -1, "args",
      FuncError, "src/gap.c:Error" },

    { "COM_FILE", 2, "filename, crc",
      FuncCOM_FILE, "src/gap.c:COM_FILE" },

    { "COM_FUN", 1, "number",
      FuncCOM_FUN, "src/gap.c:COM_FUN" },

    { "MAKE_INIT", 2, "output, input",
      FuncMAKE_INIT, "src/gap.c:MAKE_INIT" },

    { "GAP_CRC", 1, "filename",
      FuncGAP_CRC, "src/gap.c:GAP_CRC" },

    { "LOAD_DYN", 2, "filename, crc",
      FuncLOAD_DYN, "src/gap.c:LOAD_DYN" },

    { "LOAD_STAT", 2, "filename, crc",
      FuncLOAD_STAT, "src/gap.c:LOAD_STAT" },

    { "SHOW_STAT", 0, "",
      FuncSHOW_STAT, "src/gap.c:SHOW_STAT" },

    { "GASMAN", -1, "args",
      FuncGASMAN, "src/gap.c:GASMAN" },

    { "GASMAN_STATS", 0, "",
      FuncGASMAN_STATS, "src/gap.c:GASMAN_STATS" },

    { "GASMAN_MESSAGE_STATUS", 0, "",
      FuncGASMAN_MESSAGE_STATUS, "src/gap.c:GASMAN_MESSAGE_STATUS" },

    { "GASMAN_LIMITS", 0, "",
      FuncGASMAN_LIMITS, "src/gap.c:GASMAN_LIMITS" },

    { "SHALLOW_SIZE", 1, "object",
      FuncSHALLOW_SIZE, "src/gap.c:SHALLOW_SIZE" },

    { "TNUM_OBJ", 1, "object",
      FuncTNUM_OBJ, "src/gap.c:TNUM_OBJ" },

    { "TNUM_OBJ_INT", 1, "object",
      FuncTNUM_OBJ_INT, "src/gap.c:TNUM_OBJ_INT" },

    { "XTNUM_OBJ", 1, "object",
      FuncXTNUM_OBJ, "src/gap.c:XTNUM_OBJ" },

    { "OBJ_HANDLE", 1, "object",
      FuncOBJ_HANDLE, "src/gap.c:OBJ_HANDLE" },

    { "HANDLE_OBJ", 1, "object",
      FuncHANDLE_OBJ, "src/gap.c:HANDLE_OBJ" },

    { "SWAP_MPTR", 2, "obj1, obj2",
      FuncSWAP_MPTR, "src/gap.c:SWAP_MPTR" },

    { "LoadedModules", 0, "",
      FuncLoadedModules, "src/gap.c:LoadedModules" },

    { "WindowCmd", 1, "arg-list",
      FuncWindowCmd, "src/gap.c:WindowCmd" },

    { "ErrorCount", 0, "",
      FuncErrorCount, "src/gap.c:ErrorCount" },

    { "Sleep", 1, "secs",
      FuncSleep, "src/gap.c:Sleep" },

    { "QUIT_GAP", 0, "",
      FuncQUIT_GAP, "src/gap.c:QUIT_GAP" },

    { "SetErrorHandler", 1, "handler",
      FuncSetErrorHandler, "src/gap.c:SetErrorHandler" },

    { "CallFuncTrapError", 1, "func",
      FuncCallFuncTrapError, "src/gap.c:FuncCallFuncTrapError" },

    { "MASTER_POINTER_NUMBER", 1, "ob",
      FuncMASTER_POINTER_NUMBER, "src/gap.c:MASTER_POINTER_NUMBER" },

    { "FUNC_BODY_SIZE", 1, "f",
      FuncFUNC_BODY_SIZE, "src/gap.c:FUNC_BODY_SIZE" },

    { 0 }

};


/****************************************************************************
**

*F  InitKernel( <module> )  . . . . . . . . initialise kernel data structures
*/
static Int InitKernel (
    StructInitInfo *    module )
{
    /* init the completion function                                        */
    InitGlobalBag( &CompNowFuncs,      "src/gap.c:CompNowFuncs"      );
    InitGlobalBag( &CompThenFuncs,     "src/gap.c:CompThenFuncs"     );
    InitGlobalBag( &CompLists,         "src/gap.c:CompLists"         );
    InitGlobalBag( &StringUncompleted, "src/gap.c:StringUncompleted" );
    InitGlobalBag( &EmptyList,         "src/gap.c:EmptyList"         );

    InitGlobalBag( &Revisions,         "src/gap.c:Revisions"         );

    /* list of exit functions                                              */
    InitGlobalBag( &AtExitFunctions, "src/gap.c:AtExitFunctions" );
    InitGlobalBag( &WindowCmdString, "src/gap.c:WindowCmdString" );

    /* init filters and functions                                          */
    InitHdlrFuncsFromTable( GVarFuncs );

    /* use short cookies to save space in saved workspace                  */
    InitHandlerFunc( DoComplete0args, "c0" );
    InitHandlerFunc( DoComplete1args, "c1" );
    InitHandlerFunc( DoComplete2args, "c2" );
    InitHandlerFunc( DoComplete3args, "c3" );
    InitHandlerFunc( DoComplete4args, "c4" );
    InitHandlerFunc( DoComplete5args, "c5" );
    InitHandlerFunc( DoComplete6args, "c6" );
    InitHandlerFunc( DoCompleteXargs, "cX" );


    /* establish Fopy of ViewObj                                           */
    ImportFuncFromLibrary(  "ViewObj", 0L );

    /* Also of OnBreak and OnBreakMessage, but we don't want them made
       ReadOnly                                                            */
    InitFopyGVar(  "OnBreak", &OnBreak );
    InitFopyGVar(  "OnBreakMessage", &OnBreakMessage );

    /* Initialize some hooks: */
    InitCopyGVar("OnGapPromptHook",&OnGapPromptHook);
    InitCopyGVar("OnQuit",&OnQuit);
#if !SYS_MAC_MWC
#if HAVE_SELECT
    InitCopyGVar("OnCharReadHookActive",&OnCharReadHookActive);
    InitCopyGVar("OnCharReadHookInFds",&OnCharReadHookInFds);
    InitCopyGVar("OnCharReadHookInFuncs",&OnCharReadHookInFuncs);
    InitCopyGVar("OnCharReadHookOutFds",&OnCharReadHookOutFds);
    InitCopyGVar("OnCharReadHookOutFuncs",&OnCharReadHookOutFuncs);
    InitCopyGVar("OnCharReadHookExcFds",&OnCharReadHookExcFds);
    InitCopyGVar("OnCharReadHookExcFuncs",&OnCharReadHookExcFuncs);
#endif
#endif

    /* If a package or .gaprc or file read from the command line
       sets this to a function, then we want to know                       */
    InitCopyGVar(  "AlternativeMainLoop", &AlternativeMainLoop );

    InitGlobalBag(&ErrorHandler, "gap.c: ErrorHandler");

    /* return success                                                      */
    return 0;
}


/****************************************************************************
**
*F  PostRestore( <module> ) . . . . . . . . . . . . . after restore workspace
*/
static Int PostRestore (
    StructInitInfo *    module )
{
    UInt var;
    Obj optrec;

    optrec = MakeOptionsRecord();
    var = GVarName("SY_RESTORE_OPTIONS");
    MakeReadWriteGVar(var);
    AssGVar(var, optrec);
    MakeReadOnlyGVar(var);

    /* create a revision record                                            */
    Revisions = NEW_PREC(0);
    var = GVarName( "Revision" );
    MakeReadWriteGVar(var);
    AssGVar( var, Revisions );
    MakeReadOnlyGVar(var);

    /* library name and other stuff                                        */
    var = GVarName( "DEBUG_LOADING" );
    MakeReadWriteGVar(var);
    AssGVar( var, (SyDebugLoading ? True : False) );
    MakeReadOnlyGVar(var);

    /* construct the `ViewObj' variable                                    */
    ViewObjGVar = GVarName( "ViewObj" );

    /* construct the last and time variables                               */
    Last              = GVarName( "last"  );
    Last2             = GVarName( "last2" );
    Last3             = GVarName( "last3" );
    Time              = GVarName( "time"  );
    SaveOnExitFileGVar= GVarName( "SaveOnExitFile" );
    QUITTINGGVar      = GVarName( "QUITTING" );

    /* return success                                                      */
    return 0;
}


/****************************************************************************
**
*F  InitLibrary( <module> ) . . . . . . .  initialise library data structures
*/
static Int InitLibrary (
    StructInitInfo *    module )
{
    UInt                var, lenvec, lenstr, i;
    Obj                 optrec, tmp, str;


    optrec = MakeOptionsRecord();
    var = GVarName("SY_COMMAND_LINE_OPTIONS");
    AssGVar(var, optrec);
    MakeReadOnlyGVar(var);

    /* make command line and environment available to GAP level       */
    for (lenvec=0; sysargv[lenvec]; lenvec++);
    tmp = NEW_PLIST( T_PLIST+IMMUTABLE, lenvec );
    SET_LEN_PLIST( tmp, lenvec );
    for (i = 0; i<lenvec; i++) {
      lenstr = SyStrlen(sysargv[i]);
      str = NEW_STRING(lenstr);
      SyStrncat(CSTR_STRING(str), sysargv[i], lenstr);
      SET_LEN_STRING(str, lenstr);
      SET_ELM_PLIST(tmp, i+1, str);
      CHANGED_BAG(tmp);
    }
    var = GVarName("SYSTEM_COMMAND_LINE");
    MakeReadWriteGVar(var);
    AssGVar(var, tmp);
    MakeReadOnlyGVar(var);

    for (lenvec=0; sysenviron[lenvec]; lenvec++);
    tmp = NEW_PLIST( T_PLIST+IMMUTABLE, lenvec );
    SET_LEN_PLIST( tmp, lenvec );
    for (i = 0; i<lenvec; i++) {
      lenstr = SyStrlen(sysenviron[i]);
      str = NEW_STRING(lenstr);
      SyStrncat(CSTR_STRING(str), sysenviron[i], lenstr);
      SET_LEN_STRING(str, lenstr);
      SET_ELM_PLIST(tmp, i+1, str);
      CHANGED_BAG(tmp);
    }
    var = GVarName("SYSTEM_ENVIRONMENT");
    MakeReadWriteGVar(var);
    AssGVar(var, tmp);
    MakeReadOnlyGVar(var);

    /* init the completion function                                        */
    CompLists = NEW_PLIST( T_PLIST, 0 );
    SET_LEN_PLIST( CompLists, 0 );

    /* list of exit functions                                              */
    AtExitFunctions = NEW_PLIST( T_PLIST, 0 );
    SET_LEN_PLIST( AtExitFunctions, 0 );
    var = GVarName( "AT_EXIT_FUNCS" );
    AssGVar( var, AtExitFunctions );
    MakeReadOnlyGVar(var);

    /* we are not (yet) bailing out after QUIT */
    QUITTINGGVar      = GVarName( "QUITTING" );
    AssGVar(QUITTINGGVar, False);
    MakeReadOnlyGVar( QUITTINGGVar );

    /* share between uncompleted functions                                 */
    C_NEW_STRING( StringUncompleted, 11, "uncompleted" );
    RESET_FILT_LIST( StringUncompleted, FN_IS_MUTABLE );
    EmptyList = NEW_PLIST( T_PLIST+IMMUTABLE, 0 );
    SET_LEN_PLIST( EmptyList, 0 );

    /* init filters and functions                                          */
    InitGVarFuncsFromTable( GVarFuncs );

    /* create windows command buffer                                       */
    WindowCmdString = NEW_STRING( 1000 );



    /* return success                                                      */
    return PostRestore( module );
}


/****************************************************************************
**
*F  InitInfoGap() . . . . . . . . . . . . . . . . . . table of init functions
*/
static StructInitInfo module = {
    MODULE_BUILTIN,                     /* type                           */
    "gap",                              /* name                           */
    0,                                  /* revision entry of c file       */
    0,                                  /* revision entry of h file       */
    0,                                  /* version                        */
    0,                                  /* crc                            */
    InitKernel,                         /* initKernel                     */
    InitLibrary,                        /* initLibrary                    */
    0,                                  /* checkInit                      */
    0,                                  /* preSave                        */
    0,                                  /* postSave                       */
    PostRestore                         /* postRestore                    */
};

StructInitInfo * InitInfoGap ( void )
{
    module.revision_c = Revision_gap_c;
    module.revision_h = Revision_gap_h;
    FillInVersion( &module );
    return &module;
}


/****************************************************************************
**

*V  InitFuncsBuiltinModules . . . . .  list of builtin modules init functions
*/
static InitInfoFunc InitFuncsBuiltinModules[] = {

    /* global variables                                                    */
    InitInfoGVars,

    /* objects                                                             */
    InitInfoObjects,

    /* scanner, reader, interpreter, coder, caller, compiler               */
    InitInfoScanner,
    InitInfoRead,
    InitInfoCalls,
    InitInfoExprs,
    InitInfoStats,
    InitInfoCode,
    InitInfoVars,       /* must come after InitExpr and InitStats */
    InitInfoFuncs,
    InitInfoOpers,
    InitInfoIntrprtr,
    InitInfoCompiler,

    /* arithmetic operations                                               */
    InitInfoAriths,
    InitInfoInt,
    InitInfoRat,
    InitInfoCyc,
    InitInfoFinfield,
    InitInfoPermutat,
    InitInfoBool,
    InitInfoFloat,

    /* record packages                                                     */
    InitInfoRecords,
    InitInfoPRecord,

    /* list packages                                                       */
    InitInfoLists,
    InitInfoListOper,
    InitInfoListFunc,
    InitInfoPlist,
    InitInfoSet,
    InitInfoVector,
    InitInfoVecFFE,
    InitInfoBlist,
    InitInfoRange,
    InitInfoString,
    InitInfoGF2Vec,
    InitInfoVec8bit,

    /* free and presented groups                                           */
    InitInfoFreeGroupElements,
    InitInfoCosetTable,
    InitInfoTietze,
    InitInfoPcElements,
    InitInfoSingleCollector,
    InitInfoCombiCollector,
    InitInfoPcc,
    InitInfoDeepThought,
    InitInfoDTEvaluation,

    /* algebras                                                            */
    InitInfoSCTable,

    /* save and load workspace, weak pointers                              */
    InitInfoWeakPtr,
    InitInfoSaveLoad,

    /* input and output                                                    */
    InitInfoStreams,
    InitInfoSysFiles,
    InitInfoIOStream,

    /* main module                                                         */
    InitInfoGap,

#ifdef GAPMPI
    /* ParGAP/MPI module						   */
    InitInfoGapmpi,
#endif

    0
};


/****************************************************************************
**
*F  Modules . . . . . . . . . . . . . . . . . . . . . . . . . list of modules
*/
#ifndef MAX_MODULES
#define MAX_MODULES     1000
#endif


#ifndef MAX_MODULE_FILENAMES
#define MAX_MODULE_FILENAMES (MAX_MODULES*50)
#endif

Char LoadedModuleFilenames[MAX_MODULE_FILENAMES];
Char *NextLoadedModuleFilename = LoadedModuleFilenames;


StructInitInfo * Modules [ MAX_MODULES ];
UInt NrModules = 0;
UInt NrBuiltinModules = 0;


/****************************************************************************
**
*F  RecordLoadedModule( <module> )  . . . . . . . . store module in <Modules>
*/

void RecordLoadedModule (
    StructInitInfo *        info,
    Char *filename )
{
  UInt len;
    if ( NrModules == MAX_MODULES ) {
        Pr( "panic: no room to record module\n", 0L, 0L );
    }
    len = SyStrlen(filename);
    if (NextLoadedModuleFilename + len + 1
	> LoadedModuleFilenames+MAX_MODULE_FILENAMES) {
      Pr( "panic: no room for module filename\n", 0L, 0L );
    }
    *NextLoadedModuleFilename = '\0';
    SyStrncat(NextLoadedModuleFilename,filename, len);
    info->filename = NextLoadedModuleFilename;
    NextLoadedModuleFilename += len +1;
    Modules[NrModules++] = info;
}


/****************************************************************************
**

*F  SET_REVISION( <file>, <revision> )  . . . . . . . . . enter revision info
*/
#define SET_REVISION( file, revision ) \
  do { \
      UInt                    rev_rnam; \
      Obj                     rev_str; \
      rev_rnam = RNamName(file); \
      C_NEW_STRING( rev_str, SyStrlen(revision), (revision) ); \
      RESET_FILT_LIST( rev_str, FN_IS_MUTABLE ); \
      AssPRec( Revisions, rev_rnam, rev_str ); \
  } while (0)


/****************************************************************************
**
*F  InitializeGap() . . . . . . . . . . . . . . . . . . . . . . intialize GAP
**
**  Each module  (builtin  or compiled) exports  a sturctures  which contains
**  information about the name, version, crc, init function, save and restore
**  functions.
**
**  The init process is split into three different functions:
**
**  `InitKernel':   This function setups the   internal  data structures  and
**  tables,   registers the global bags  and   functions handlers, copies and
**  fopies.  It is not allowed to create objects, gvar or rnam numbers.  This
**  function is used for both starting and restoring.
**
**  `InitLibrary': This function creates objects,  gvar and rnam number,  and
**  does  assignments of auxillary C   variables (for example, pointers  from
**  objects, length of hash lists).  This function is only used for starting.
**
**  `PostRestore': Everything in  `InitLibrary' execpt  creating objects.  In
**  general    `InitLibrary'  will  create    all objects    and  then  calls
**  `PostRestore'.  This function is only used when restoring.
*/
extern TNumMarkFuncBags TabMarkFuncBags [ 256 ];

static Obj POST_RESTORE_FUNCS;

void InitializeGap (
    int *               pargc,
    char *              argv [] )
{
    UInt                type;
    UInt                i;
    Int                 ret;


    /* initialize the basic system and gasman                              */
#ifdef GAPMPI
    /* ParGAP/MPI needs to call MPI_Init() first to remove command line args */
    InitGapmpi( pargc, &argv, &BreakOnError );
#endif

    InitSystem( *pargc, argv );

    InitBags( SyAllocBags, SyStorMin,
              0, (Bag*)(((UInt)pargc/SyStackAlign)*SyStackAlign), SyStackAlign,
              SyCacheSize, 0, SyAbortBags );
    InitMsgsFuncBags( SyMsgsBags );


    /* get info structures for the build in modules                        */
    for ( i = 0;  InitFuncsBuiltinModules[i];  i++ ) {
        if ( NrModules == MAX_MODULES ) {
            FPUTS_TO_STDERR( "panic: too many builtin modules\n" );
            SyExit(1);
        }
        Modules[NrModules++] = InitFuncsBuiltinModules[i]();
#       ifdef DEBUG_LOADING
            FPUTS_TO_STDERR( "#I  InitInfo(builtin " );
            FPUTS_TO_STDERR( Modules[NrModules-1]->name );
            FPUTS_TO_STDERR( ")\n" );
#       endif
    }
    NrBuiltinModules = NrModules;

    /* call kernel initialisation                                          */
    for ( i = 0;  i < NrBuiltinModules;  i++ ) {
        if ( Modules[i]->initKernel ) {
#           ifdef DEBUG_LOADING
                FPUTS_TO_STDERR( "#I  InitKernel(builtin " );
                FPUTS_TO_STDERR( Modules[i]->name );
                FPUTS_TO_STDERR( ")\n" );
#           endif
            ret =Modules[i]->initKernel( Modules[i] );
            if ( ret ) {
                FPUTS_TO_STDERR( "#I  InitKernel(builtin " );
                FPUTS_TO_STDERR( Modules[i]->name );
                FPUTS_TO_STDERR( ") returned non-zero value\n" );
            }
        }
    }

    InitGlobalBag(&POST_RESTORE_FUNCS, "gap.c: POST_RESTORE_FUNCS");
    InitCopyGVar( "POST_RESTORE_FUNCS", &POST_RESTORE_FUNCS);

    /* you should set 'COUNT_BAGS' as well                                 */
#   ifdef DEBUG_LOADING
        if ( SyRestoring ) {
            Pr( "#W  after setup\n", 0L, 0L );
            Pr( "#W  %36s ", (Int)"type",  0L          );
            Pr( "%8s %8s ",  (Int)"alive", (Int)"kbyte" );
            Pr( "%8s %8s\n",  (Int)"total", (Int)"kbyte" );
            for ( i = 0;  i < 256;  i++ ) {
                if ( InfoBags[i].name != 0 && InfoBags[i].nrAll != 0 ) {
                    char    buf[41];

                    buf[0] = '\0';
                    SyStrncat( buf, InfoBags[i].name, 40 );
                    Pr("#W  %36s ",    (Int)buf, 0L );
                    Pr("%8d %8d ", (Int)InfoBags[i].nrLive,
                       (Int)(InfoBags[i].sizeLive/1024));
                    Pr("%8d %8d\n",(Int)InfoBags[i].nrAll,
                       (Int)(InfoBags[i].sizeAll/1024));
                }
            }
        }
#   endif

#ifdef SYS_IS_MAC_MWC
	ActivateIntr ();
#endif

    /* and now for a special hack                                          */
    for ( i = LAST_CONSTANT_TNUM+1; i <= LAST_REAL_TNUM; i++ ) {
        TabMarkFuncBags[ i+COPYING ] = TabMarkFuncBags[ i ];
    }

    /* if we are restoring, load the workspace and call the post restore   */
    if ( SyRestoring ) {
        LoadWorkspace(SyRestoring);
        for ( i = 0;  i < NrModules;  i++ ) {
            if ( Modules[i]->postRestore ) {
#               ifdef DEBUG_LOADING
                    FPUTS_TO_STDERR( "#I  PostRestore(builtin " );
                    FPUTS_TO_STDERR( Modules[i]->name );
                    FPUTS_TO_STDERR( ")\n" );
#               endif
                ret = Modules[i]->postRestore( Modules[i] );
                if ( ret ) {
                    FPUTS_TO_STDERR( "#I  PostRestore(builtin " );
                    FPUTS_TO_STDERR( Modules[i]->name );
                    FPUTS_TO_STDERR( ") returned non-zero value\n" );
                }
            }
        }
	SyRestoring = NULL;

	/* call the post restore functions */
	if (POST_RESTORE_FUNCS != (Obj) 0
	    && IS_SMALL_LIST(POST_RESTORE_FUNCS))
        {
	    UInt l;
	    UInt j;
	    Obj func;
	    Obj res;
	    l = LEN_LIST(POST_RESTORE_FUNCS);
	    for (j = 1; j <= l; j++)
	      {
		func = ELM0_LIST(POST_RESTORE_FUNCS, j);
		if (func != (Obj) 0 && IS_FUNC(func))
		  {
		    res = CALL_0ARGS(func);
		    if (res == Fail)
		      {
			FPUTS_TO_STDERR("panic -- post restore function returned fail");
			SyExit(j);
		      }
		  }
	      }
        }
    }

    /* otherwise call library initialisation                               */
    else {
        WarnInitGlobalBag = 1;
#       ifdef DEBUG_HANDLER_REGISTRATION
            CheckAllHandlers();
#       endif

	SyInitializing = 1;
        for ( i = 0;  i < NrBuiltinModules;  i++ ) {
            if ( Modules[i]->initLibrary ) {
#               ifdef DEBUG_LOADING
                    FPUTS_TO_STDERR( "#I  InitLibrary(builtin " );
                    FPUTS_TO_STDERR( Modules[i]->name );
                    FPUTS_TO_STDERR( ")\n" );
#               endif
                ret = Modules[i]->initLibrary( Modules[i] );
                if ( ret ) {
                    FPUTS_TO_STDERR( "#I  InitLibrary(builtin " );
                    FPUTS_TO_STDERR( Modules[i]->name );
                    FPUTS_TO_STDERR( ") returned non-zero value\n" );
                }
            }
        }
        WarnInitGlobalBag = 0;
    }

    /* check initialisation                                                */
    for ( i = 0;  i < NrModules;  i++ ) {
        if ( Modules[i]->checkInit ) {
#           ifdef DEBUG_LOADING
                FPUTS_TO_STDERR( "#I  CheckInit(builtin " );
                FPUTS_TO_STDERR( Modules[i]->name );
                FPUTS_TO_STDERR( ")\n" );
#           endif
            ret = Modules[i]->checkInit( Modules[i] );
            if ( ret ) {
                FPUTS_TO_STDERR( "#I  CheckInit(builtin " );
                FPUTS_TO_STDERR( Modules[i]->name );
                FPUTS_TO_STDERR( ") returned non-zero value\n" );
            }
        }
    }

    /* create a revision record (overwrite a restored one)                 */
    for ( i = 0;  i < NrBuiltinModules;  i++ ) {
        Char buf[30];

        buf[0] = 0;
        SyStrncat( buf, Modules[i]->name, 27 );
        SyStrncat( buf, "_c", 2 );
        SET_REVISION( buf, Modules[i]->revision_c );
        buf[0] = 0;
        SyStrncat( buf, Modules[i]->name, 27 );
        SyStrncat( buf, "_h", 2 );
        SET_REVISION( buf, Modules[i]->revision_h );
    }

    /* add revisions for files which are not modules                       */
    {
        SET_REVISION( "system_c", Revision_system_c );
        SET_REVISION( "system_h", Revision_system_h );
        SET_REVISION( "gasman_c", Revision_gasman_c );
        SET_REVISION( "gasman_h", Revision_gasman_h );
    }

    /* read the init files                                                 */
    if ( SySystemInitFile[0] ) {
        if ( READ_GAP_ROOT(SySystemInitFile) == 0 ) {
            if ( ! SyQuiet ) {
                Pr( "gap: hmm, I cannot find '%s' maybe",
                    (Int)SySystemInitFile, 0L );
                Pr( " use option '-l <gaproot>'?\n If you ran the GAP\
 binary directly, try running the 'gap.sh' or 'gap.bat' script instead.", 0L, 0L );
            }
        }
    }
    SyInitializing = 0;
    for ( i = 0; i < sizeof(SyInitfiles)/sizeof(SyInitfiles[0]); i++ ) {
        if ( SyInitfiles[i][0] != '\0' ) {
            if ( OpenInput( SyInitfiles[i] ) ) {
                ClearError();
                while ( 1 ) {
                    type = ReadEvalCommand();
                    if ( type == STATUS_RETURN_VAL || type == STATUS_RETURN_VOID ) {
                        Pr("'return' must not be used in file",0L,0L);
                    }
                    else if ( type == STATUS_QUIT || type == STATUS_EOF ) {
                        break;
                    }
		    else if (type == STATUS_QQUIT)
		      {
			UserHasQUIT = 1;
			break;
		      }
		    else if (type == STATUS_ERROR)
		      {
			Pr("Error reading initial file \"%s\" abandoning remaining files\n",
			   (Int)SyInitfiles[i], 0);
			UserHasQuit = 0; /* enough to have got to here */
			break;
		      }
		    else if (type != STATUS_END)
		      {
			Pr("Unexpected status %d from reading initial file \"%s\"\n",type, (Int)SyInitfiles[i]);
			break;
		      }
                }
                CloseInput();
                ClearError();
		if (UserHasQUIT)
		  break;
            }
            else {
                Pr( "Error, file \"%s\" must exist and be readable\n",
                    (Int)SyInitfiles[i], 0L );
            }
        }
    }

    if (SyBreakSuppress)
      BreakOnError = 0;

}


/****************************************************************************
**

*E  gap.c . . . . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
*/






