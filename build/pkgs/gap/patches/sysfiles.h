/****************************************************************************
**
*W  sysfiles.h                  GAP source                       Frank Celler
*W                                                         & Martin Schoenert
*W                                                  & Burkhard Hoefling (MAC)
**
*H  @(#)$Id: sysfiles.h,v 4.34.6.1 2008/09/03 15:50:39 sal Exp $
**
*Y  Copyright (C)  1996,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
*Y  (C) 1998 School Math and Comp. Sci., University of St.  Andrews, Scotland
*Y  Copyright (C) 2002 The GAP Group
**
**  The  file 'system.c'  declares  all operating system  dependent functions
**  except file/stream handling which is done in "sysfiles.h".
*/
#ifdef  INCLUDE_DECLARATION_PART
const char * Revision_sysfiles_h =
   "@(#)$Id: sysfiles.h,v 4.34.6.1 2008/09/03 15:50:39 sal Exp $";
#endif


#ifndef SYS_STDIO_H                     /* standard input/output functions */
# include <stdio.h>
# define SYS_STDIO_H
#endif


/****************************************************************************
**


*F * * * * * * * * * * * * * * dynamic loading  * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*F  SyFindOrLinkGapRootFile( <filename>, <res>, <len> ) . . . .  load or link
**
**  'SyFindOrLinkGapRootFile'  tries to find a GAP  file in the root area and
**  check  if   there is a corresponding    statically  or dynamically linked
**  module.  If the CRC matches this module  is loaded otherwise the filename
**  is returned.
*/
extern Int SyFindOrLinkGapRootFile (
            Char *          filename,
            Int4            crc_gap,
            Char *          result,
            Int             len,
        StructInitInfo **   info_result );


/****************************************************************************
**
*F  SyGAPCRC( <name> )  . . . . . . . . . . . . . . . . . . crc of a GAP file
**
**  This function should  be clever and handle  white spaces and comments but
**  one has to certain that such characters are not ignored in strings.
*/
extern Int4 SyGAPCRC(
            Char *          name );


/****************************************************************************
**
*F  SyLoadModule( <name> )  . . . . . . . . . . . . .  load a compiled module
*/
extern InitInfoFunc SyLoadModule(
            Char *          name );


/****************************************************************************
**
*V  SyCanLoadDynamicModules . . . true if system supports dynamical libraries
*/
#if SYS_MAC_MWC
# include <Types.h>
extern Boolean SyCanLoadDynamicModules;
#endif


/****************************************************************************
**

*F * * * * * * * * * * * * * * * window handler * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*F  syWinPut( <fid>, <cmd>, <str> ) . . . . send a line to the window handler
**
**  'syWinPut'  send the command   <cmd> and the  string  <str> to the window
**  handler associated with the  file identifier <fid>.   In the string <str>
**  '@'  characters are duplicated, and   control characters are converted to
**  '@<chr>', e.g., <newline> is converted to '@J'.
*/
extern void syWinPut (
            Int                 fid,
            const Char *    cmd,
            const Char *    str );


/****************************************************************************
**
*F  SyWinCmd( <str>, <len> )  . . . . . . . . . . . .  . execute a window cmd
**
**  'SyWinCmd' send   the  command <str> to  the   window  handler (<len>  is
**  ignored).  In the string <str> '@' characters are duplicated, and control
**  characters  are converted to  '@<chr>', e.g.,  <newline> is converted  to
**  '@J'.  Then  'SyWinCmd' waits for  the window handlers answer and returns
**  that string.
*/
extern Char * SyWinCmd (
            const Char *    str,
            UInt                len );


/****************************************************************************
**


*F * * * * * * * * * * * * * * * * open/close * * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*V  syBuf . . . . . . . . . . . . . .  buffer and other info for files, local
**
**  'syBuf' is  a array used as  buffers for  file I/O to   prevent the C I/O
**  routines  from   allocating their  buffers  using  'malloc',  which would
**  otherwise confuse Gasman.
*/
#if SYS_MAC_MWC
/* on the Mac, we need some additional info about the file
** fromDoc should really be of type DocumentRecord*, but but this would require "macedit.h",
** which gives trouble with some of the other  GAP source files under CW Pro 2,
** because it refuses to include "List.h" if it has already read <List.h>
*/
#include <Files.h>    /* for FSSpec, SIgnedByte and Boolean (from MacTypes.h) */

#define MIN_BUFSIZ 16
typedef struct {
    short	     	fp;                     /* reference number for this file      */
#if 0
    FILE *      	echo;                   /* file pointer for the echo       */
#endif
    void *		 	fromDoc;				/* the document window from which it reads, if any */
    Int              bufno;                 /* if non-negative then this file has a buffer in
					                           syBuffers[bufno]; If negative, this file may not
					                           be buffered */
	FSSpec 			fsspec;					/* the file specs for this file */
	SignedByte		permission;
	Boolean			binary;                 /* binary file? */
    UInt       		isTTY;		      /* set in Fopen when this fid is a *stdin* or *errin*
				       and really is a tty*/
} SYS_SY_BUF;

#else
typedef struct {
  int         fp;                     /* file descriptor for this file      */
  int         echo;                   /* file descriptor for the echo       */
  UInt        pipe;                   /* file is really a pipe           */
  FILE       *pipehandle;             /* for pipes we need to remember the file handle */
  UInt       ateof;                   /* set to 1 by any read operation that hits eof
					 reset to 0 by a subsequent successful read */
  UInt        crlast;                 /* records that last character read was \r for
					 cygwin and othger systems that need end-of-line
					 hackery */
  Int        bufno;                   /* if non-negative then this file has a buffer in
					 syBuffers[bufno]; If negative, this file may not
					 be buffered */
  UInt       isTTY;		      /* set in Fopen when this fid is a *stdin* or *errin*
				       and really is a tty*/
} SYS_SY_BUF;

#endif

#define SYS_FILE_BUF_SIZE 20000

typedef struct {
  Char buf[SYS_FILE_BUF_SIZE];
  UInt inuse;
  UInt bufstart;
  UInt buflen;
} SYS_SY_BUFFER;

extern SYS_SY_BUF syBuf [256];

extern SYS_SY_BUFFER syBuffers[32];

extern UInt SySetBuffering( UInt fid );

/****************************************************************************
**
*F  SyFileno( <fid> ) . . . . . . . . . . . . . . get operating system fileno
*/
#define SyFileno(fid)   (fid==-1?-1:syBuf[fid].fp)


/****************************************************************************
**

*F  SyFopen( <name>, <mode> ) . . . . . . . .  open the file with name <name>
**
**  The function 'SyFopen'  is called to open the file with the name  <name>.
**  If <mode> is "r" it is opened for reading, in this case  it  must  exist.
**  If <mode> is "w" it is opened for writing, it is created  if  neccessary.
**  If <mode> is "a" it is opened for appending, i.e., it is  not  truncated.
**
**  'SyFopen' returns an integer used by the scanner to  identify  the  file.
**  'SyFopen' returns -1 if it cannot open the file.
**
**  The following standard files names and file identifiers  are  guaranteed:
**  'SyFopen( "*stdin*", "r")' returns 0 identifying the standard input file.
**  'SyFopen( "*stdout*","w")' returns 1 identifying the standard outpt file.
**  'SyFopen( "*errin*", "r")' returns 2 identifying the brk loop input file.
**  'SyFopen( "*errout*","w")' returns 3 identifying the error messages file.
**
**  If it is necessary to adjust the  filename  this  should  be  done  here.
**  Right now GAP does not read nonascii files, but if this changes sometimes
**  'SyFopen' must adjust the mode argument to open the file in binary  mode.
*/
extern Int SyFopen (
            Char *              name,
            Char *              mode );


/****************************************************************************
**
*F  SyFclose( <fid> ) . . . . . . . . . . . . . . . . .  close the file <fid>
**
**  'SyFclose' closes the file with the identifier <fid>  which  is  obtained
**  from 'SyFopen'.
*/
extern Int SyFclose (
            Int                 fid );


/****************************************************************************
**
*F  SyIsEndOfFile( <fid> )  . . . . . . . . . . . . . . . end of file reached
*/
extern Int SyIsEndOfFile (
    Int                 fid );


/****************************************************************************
**
*F  SyFgets( <line>, <lenght>, <fid> )  . . . . .  get a line from file <fid>
**
**  'SyFgets' is called to read a line from the file  with  identifier <fid>.
**  'SyFgets' (like 'fgets') reads characters until either  <length>-1  chars
**  have been read or until a <newline> or an  <eof> character is encoutered.
**  It retains the '\n' (unlike 'gets'), if any, and appends '\0' to  <line>.
**  'SyFgets' returns <line> if any char has been read, otherwise '(char*)0'.
**
**  'SyFgets'  allows to edit  the input line if the  file  <fid> refers to a
**  terminal with the following commands:
**
**      <ctr>-A move the cursor to the beginning of the line.
**      <esc>-B move the cursor to the beginning of the previous word.
**      <ctr>-B move the cursor backward one character.
**      <ctr>-F move the cursor forward  one character.
**      <esc>-F move the cursor to the end of the next word.
**      <ctr>-E move the cursor to the end of the line.
**
**      <ctr>-H, <del> delete the character left of the cursor.
**      <ctr>-D delete the character under the cursor.
**      <ctr>-K delete up to the end of the line.
**      <esc>-D delete forward to the end of the next word.
**      <esc>-<del> delete backward to the beginning of the last word.
**      <ctr>-X delete entire input line, and discard all pending input.
**      <ctr>-Y insert (yank) a just killed text.
**
**      <ctr>-T exchange (twiddle) current and previous character.
**      <esc>-U uppercase next word.
**      <esc>-L lowercase next word.
**      <esc>-C capitalize next word.
**
**      <tab>   complete the identifier before the cursor.
**      <ctr>-L insert last input line before current character.
**      <ctr>-P redisplay the last input line, another <ctr>-P will redisplay
**              the line before that, etc.  If the cursor is not in the first
**              column only the lines starting with the string to the left of
**              the cursor are taken. The history is limitied to ~8000 chars.
**      <ctr>-N Like <ctr>-P but goes the other way round through the history
**      <esc>-< goes to the beginning of the history.
**      <esc>-> goes to the end of the history.
**      <ctr>-O accept this line and perform a <ctr>-N.
**
**      <ctr>-V enter next character literally.
**      <ctr>-U execute the next command 4 times.
**      <esc>-<num> execute the next command <num> times.
**      <esc>-<ctr>-L repaint input line.
**
**  Not yet implemented commands:
**
**      <ctr>-S search interactive for a string forward.
**      <ctr>-R search interactive for a string backward.
**      <esc>-Y replace yanked string with previously killed text.
**      <ctr>-_ undo a command.
**      <esc>-T exchange two words.
*/
#if !SYS_MAC_MWC
#if HAVE_SELECT
extern Obj OnCharReadHookActive;  /* if bound the hook is active */
extern Obj OnCharReadHookInFds;   /* a list of UNIX file descriptors */
extern Obj OnCharReadHookInFuncs; /* a list of GAP functions */
extern Obj OnCharReadHookOutFds;  /* a list of UNIX file descriptors */
extern Obj OnCharReadHookOutFuncs;/* a list of GAP functions with 0 args */
extern Obj OnCharReadHookExcFds;  /* a list of UNIX file descriptors */
extern Obj OnCharReadHookExcFuncs;/* a list of GAP functions with 0 args */
#endif
#endif

extern Char * SyFgets (
            Char *              line,
            UInt                length,
            Int                 fid );


/****************************************************************************
**
*F  SyFputs( <line>, <fid> )  . . . . . . . .  write a line to the file <fid>
**
**  'SyFputs' is called to put the  <line>  to the file identified  by <fid>.
*/
extern void SyFputs (
            Char *              line,
            Int                 fid );


/****************************************************************************
**
*F  SyIsIntr()  . . . . . . . . . . . . . . . . check wether user hit <ctr>-C
**
**  'SyIsIntr' is called from the evaluator at  regular  intervals  to  check
**  wether the user hit '<ctr>-C' to interrupt a computation.
**
**  'SyIsIntr' returns 1 if the user typed '<ctr>-C' and 0 otherwise.
*/
extern void SyInstallAnswerIntr ( void );

extern UInt SyIsIntr ( void );


/****************************************************************************
**
*V  SyInFid . . .
**
*V  SyOutFid
*/
#if SYS_MAC_MWC
extern Int 			SyInFid, SyOutFid; /* for i/o redirection */
#endif

/****************************************************************************
**

*F * * * * * * * * * * * * * * * * * output * * * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*F  SyEchoch( <ch>, <fid> ) . . . . . . . . . . . echo a char to <fid>, local
*/
extern Int SyEchoch (
    Int                 ch,
    Int                 fid );


/****************************************************************************
**

*F * * * * * * * * * * * * * * * * * input  * * * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*F  SyFtell( <fid> )  . . . . . . . . . . . . . . . . . .  position of stream
*/
extern Int SyFtell (
    Int                 fid );


/****************************************************************************
**
*F  SyFseek( <fid>, <pos> )   . . . . . . . . . . . seek a position of stream
*/
extern Int SyFseek (
    Int                 fid,
    Int                 pos );


/****************************************************************************
**
*F  SyGetch( <fid> )  . . . . . . . . . . . . . . . . . get a char from <fid>
**
**  'SyGetch' reads a character from <fid>, which is switch to raw mode if it
**  is *stdin* or *errin*.
*/
extern Int SyGetch (
    Int                 fid );


/****************************************************************************
**
*F  SyGetc( <fid> ).  . . . . . . . . . . . . . . . . . get a char from <fid>
**
**  'SyGetc' reads a character from <fid>, without any translation or
**   interference
*/

extern Int SyGetc
(
    Int                 fid );

/****************************************************************************
**
*F  SyPutc( <fid>, <char> ).. . . . . . . . . . . . . . . put a char to <fid>
**
**  'SyPutc' writes a character to <fid>, without any translation or
**   interference
*/

extern Int SyPutc
(
    Int                 fid,
    Char                c );


/****************************************************************************
**

*F * * * * * * * * * * * * system error messages  * * * * * * * * * * * * * *
*/


/****************************************************************************
**
*V  SyLastMacErrorNo . . . . . . . . . . . . . .last error number, Macintosh
*/
#ifdef SYS_IS_MAC_MWC
extern OSErr SyLastMacErrorCode;
#endif


/****************************************************************************
**
*V  SyLastErrorNo . . . . . . . . . . . . . . . . . . . . . last error number
*/
extern Int SyLastErrorNo;


/****************************************************************************
**
*V  SyLastErrorMessage  . . . . . . . . . . . . . . . . .  last error message
*/
extern Char SyLastErrorMessage [ 1024 ];


/****************************************************************************
**
*F  SyClearErrorNo()  . . . . . . . . . . . . . . . . .  clear error messages
*/
extern void SyClearErrorNo ( void );


/****************************************************************************
**
*F  SySetErrorNo()  . . . . . . . . . . . . . . . . . . . . set error message
*/
extern void SySetErrorNo ( void );


/****************************************************************************
**

*F * * * * * * * * * * * * * file and execution * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*F  SyExec( <cmd> ) . . . . . . . . . . . execute command in operating system
**
**  'SyExec' executes the command <cmd> (a string) in the operating system.
**
**  'SyExec'  should call a command  interpreter  to execute the command,  so
**  that file name expansion and other common  actions take place.  If the OS
**  does not support this 'SyExec' should print a message and return.
**
**  For UNIX we can use 'system', which does exactly what we want.
*/
extern void SyExec (
    Char *              cmd );


/****************************************************************************
**
*F  SyExecuteProcess( <dir>, <prg>, <in>, <out>, <args> ) . . . . new process
**
**  Start  <prg> in  directory <dir>  with  standard input connected to <in>,
**  standard  output  connected to <out>   and arguments.  No  path search is
**  performed, the return  value of the process  is returned if the operation
**  system supports such a concept.
*/
extern UInt SyExecuteProcess (
    Char *                  dir,
    Char *                  prg,
    Int                     in,
    Int                     out,
    Char *                  args[] );


/****************************************************************************
**
*V  SyCanExec . . . . . . true if operating system supports launching another
**                        application
*/
#if SYS_MAC_MWC
extern Boolean		SyCanExec;
#endif


/****************************************************************************
**

*F  SyIsExistingFile( <name> )  . . . . . . . . . . . does file <name> exists
**
**  'SyIsExistingFile' returns 1 if the  file <name> exists and 0  otherwise.
**  It does not check if the file is readable, writable or excuteable. <name>
**  is a system dependent description of the file.
*/
extern Int SyIsExistingFile(
            Char * name );


/****************************************************************************
**
*F  SyIsReadableFile( <name> )  . . . . . . . . . . . is file <name> readable
**
**  'SyIsReadableFile'   returns 1  if the   file  <name> is   readable and 0
**  otherwise. <name> is a system dependent description of the file.
*/
extern Int SyIsReadableFile(
            Char * name );


/****************************************************************************
**
*F  SyIsWritable( <name> )  . . . . . . . . . . . is the file <name> writable
**
**  'SyIsWriteableFile'   returns 1  if the  file  <name>  is  writable and 0
**  otherwise. <name> is a system dependent description of the file.
*/
extern Int SyIsWritableFile(
            Char * name );


/****************************************************************************
**
*F  SyIsExecutableFile( <name> )  . . . . . . . . . is file <name> executable
**
**  'SyIsExecutableFile' returns 1 if the  file <name>  is  executable and  0
**  otherwise. <name> is a system dependent description of the file.
*/
extern Int SyIsExecutableFile(
            Char * name );


/****************************************************************************
**
*F  SyIsDirectoryPath( <name> ) . . . . . . . . .  is file <name> a directory
**
**  'SyIsDirectoryPath' returns 1 if the  file <name>  is a directory  and  0
**  otherwise. <name> is a system dependent description of the file.
*/
extern Int SyIsDirectoryPath ( Char * name );


/****************************************************************************
**
*F  SyRemoveFile( <name> )  . . . . . . . . . . . . . . .  remove file <name>
*/
extern Int SyRemoveFile ( Char * name );


/****************************************************************************
**
*F  SyFindGapRootFile( <filename> ) . . . . . . . .  find file in system area
*/
extern Char * SyFindGapRootFile (
            Char *          filename );


/****************************************************************************
**

*F * * * * * * * * * * * * * * * directories  * * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*F  SyTmpname() . . . . . . . . . . . . . . . . . return a temporary filename
**
**  'SyTmpname' creates and returns a new temporary name.
*/
extern Char * SyTmpname ( void );


/****************************************************************************
**
*F  SyTmpdir( <hint> )  . . . . . . . . . . . .  return a temporary directory
**
**  'SyTmpdir'  returns the directory   for  a temporary directory.  This  is
**  guaranteed  to be newly  created and empty  immediately after the call to
**  'SyTmpdir'. <hint> should be used by 'SyTmpdir' to  construct the name of
**  the directory (but 'SyTmpdir' is free to use only  a part of <hint>), and
**  must be a string of at most 8 alphanumerical characters.  Under UNIX this
**  would   usually   represent   '/usr/tmp/<hint>_<proc_id>_<cnt>/',   e.g.,
**  '/usr/tmp/guava_17188_1/'.
*/
extern Char * SyTmpdir ( Char * hint );

/****************************************************************************
**
*V  SyTmpVref . . . . . . . . . .. . . volume ref num for temporary directory
**
*V  SyTmpDirId . . . . . . . . . . . . . . . . dir id for temporary directory
*/
#if SYS_MAC_MWC
extern short 		SyTmpVref;
extern long 		SyTmpDirId;
#endif

/****************************************************************************
**
*F  SyFSMakeNewFSSpec ( <vol>, <dir>, <name>, <newFSSpec>) . .  get an unused
**                                                              Mac file spec
**  This returns a Mac file specification for a file in directory <dir> on
**  volume <vol>, whose name is based upon the string <name>
*/
#if SYS_MAC_MWC
OSErr SyFSMakeNewFSSpec (short vol, long dir, Str31 name, FSSpecPtr newFSSpec);
#endif

/****************************************************************************
**
*F  void getwindowsize( void )  . probe the OS for the window size and
**                               set SyNrRows and SyNrCols accordingly
*/

extern void getwindowsize( void );

extern void     InterruptExecStat ( void );

/***************************************************************************
**
*F HasAvailableBytes( <fid> ) returns positive if  a subsequent read to <fid>
**                            will read at least one byte without blocking
*/
extern Int HasAvailableBytes( UInt fid );

extern Char *SyFgetsSemiBlock (
    Char *              line,
    UInt                length,
    Int                 fid);

/***************************************************************************
**
*F SyWinPut( <fid>, <cmd>, <str> ) . . . . send a line to the window handler
**
**  'syWinPut'  send the command   <cmd> and the  string  <str> to the window
**  handler associated with the  file identifier <fid>.   In the string <str>
**  '@'  characters are duplicated, and   control characters are converted to
**  '@<chr>', e.g., <newline> is converted to '@J'.
*/

extern void syWinPut (
    Int                 fid,
    const Char *        cmd,
    const Char *        str );



/****************************************************************************
**

*F * * * * * * * * * * * * * initialize package * * * * * * * * * * * * * * *
*/

/****************************************************************************
**

*F  InitInfoSysFiles()  . . . . . . . . . . . . . . . table of init functions
*/
StructInitInfo * InitInfoSysFiles ( void );


/****************************************************************************
**

*E  sysfiles.h  . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
*/
