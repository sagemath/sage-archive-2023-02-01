/****************************************************************************
**
*W  saveload.c                  GAP source                       Steve Linton
**
*H  @(#)$Id: saveload.c,v 4.50 2002/06/17 11:37:21 sal Exp $
**
*Y  Copyright (C)  1997,  Lehrstuhl D fuer Mathematik,  RWTH Aachen,  Germany
*Y  (C) 1998 School Math and Comp. Sci., University of St.  Andrews, Scotland
*Y  Copyright (C) 2002 The GAP Group
**
**  This file contains the functions concerned with saving and loading
**  the workspace. There are support functions in gasman.c and elsewhere
**  throughout the kernel
*/
#include        "system.h"              /* system dependent part           */

const char * Revision_saveload_c =
   "@(#)$Id: saveload.c,v 4.50 2002/06/17 11:37:21 sal Exp $";

#include        <unistd.h>              /* write, read                     */

#include        "gasman.h"              /* garbage collector               */
#include        "objects.h"             /* objects                         */
#include        "bool.h"                /* booleans                        */
#include        "calls.h"               /* generic call mechanism          */
#include        "gap.h"                 /* error handling, initialisation  */
#include        "gvars.h"               /* global variables                */
#include        "streams.h"             /* streams                         */
#include        "string.h"              /* strings                         */
#include        "scanner.h"             /* scanner                         */
#include        "sysfiles.h"            /* file input/output               */
#include        "plist.h"               /* plain lists                     */
#include        "float.h"               /* floating points */
#include        "compstat.h"            /* statically compiled modules     */

#define INCLUDE_DECLARATION_PART
#include        "saveload.h"            /* saving and loading              */
#undef  INCLUDE_DECLARATION_PART


/***************************************************************************
**
** Temporary Stuff which will probably be revised to tie in with sysfiles
*/


static Int SaveFile = -1;
static UInt1 LoadBuffer[100000];
static UInt1* LBPointer = LoadBuffer;
static UInt1* LBEnd = LoadBuffer;

static Int OpenForSave( Obj fname )
{
  if (SaveFile != -1)
    {
      Pr("Already saving",0L,0L);
      return 1;
    }
  SaveFile = SyFopen(CSTR_STRING(fname), "wb");
  if (SaveFile == -1)
    {
      Pr("Couldn't open file %s to save workspace",
	 (UInt)CSTR_STRING(fname),0L);
      return 1;
    }
  LBPointer = LoadBuffer;
  LBEnd = LBPointer+sizeof(LoadBuffer);
  return 0;
}

#ifdef SYS_IS_MAC_MWC

static void CloseAfterSave( void )
{
  long count;

  if (SaveFile == -1)
    {
      Pr("Internal error -- this should never happen",0L,0L);
      SyExit(2);
    }
  count = LBPointer-LoadBuffer;
  FSWrite (syBuf[SaveFile].fp, &count, LoadBuffer);
  SyFclose(SaveFile);
  SaveFile = -1;
}

#else

static void CloseAfterSave( void )
{
  if (SaveFile == -1)
    {
      Pr("Internal error -- this should never happen",0L,0L);
      SyExit(2);
    }

  write(syBuf[SaveFile].fp, LoadBuffer, LBPointer-LoadBuffer);
  SyFclose(SaveFile);
  SaveFile = -1;
}

#endif
Int LoadFile = -1;


static void OpenForLoad( Char *fname )
{
  if (LoadFile != -1)
    {
      Pr("Internal error -- this should never happen\n",0L,0L);
      SyExit(2);
    }
  LoadFile = SyFopen(fname, "rb");
  if (LoadFile == -1)
    {
      Pr("Couldn't open saved workspace %s\n",(Int)fname,0L);
      SyExit(1);
    }
}


static void CloseAfterLoad( void )
{
  if (!LoadFile)
    {
      Pr("Internal error -- this should never happen\n",0L,0L);
      SyExit(2);
    }
  SyFclose(LoadFile);
  LoadFile = -1;
}

#ifdef SYS_IS_MAC_MWC

void SAVE_BYTE_BUF( void )
{
  long count;
  count = LBEnd-LoadBuffer;
  FSWrite (syBuf[SaveFile].fp, &count, LoadBuffer);
  LBPointer = LoadBuffer;
  return;
}

#else

void SAVE_BYTE_BUF( void )
{
  write(syBuf[SaveFile].fp, LoadBuffer, LBEnd-LoadBuffer);
  LBPointer = LoadBuffer;
  return;
}

#endif

#define SAVE_BYTE(byte) {if (LBPointer >= LBEnd) {SAVE_BYTE_BUF();} \
                          *LBPointer++ = (UInt1)(byte);}

Char * LoadByteErrorMessage = "Unexpected End of File in Load\n";

#ifdef SYS_IS_MAC_MWC

UInt1 LOAD_BYTE_BUF( void )
{
  OSErr ret;
  long count = sizeof(LoadBuffer);
  ret = FSRead (syBuf[LoadFile].fp, &count, LoadBuffer);
  if ((ret != noErr) && ((ret != eofErr || count == 0)))
    {
      Pr(LoadByteErrorMessage, 0L, 0L );
      SyExit(2);
    }
  LBEnd = LoadBuffer + count;
  LBPointer = LoadBuffer;
  return *LBPointer++;
}

#else

UInt1 LOAD_BYTE_BUF( void )
{
  Int ret;
  ret = read(syBuf[LoadFile].fp, LoadBuffer, 100000);
  if (ret <= 0)
    {
      Pr(LoadByteErrorMessage, 0L, 0L );
      SyExit(2);
    }
  LBEnd = LoadBuffer + ret;
  LBPointer = LoadBuffer;
  return *LBPointer++;
}

#endif

#define LOAD_BYTE()    (UInt1)((LBPointer >= LBEnd) ?\
                                  (LOAD_BYTE_BUF()) : (*LBPointer++))

/***************************************************************************
**
**  Low level saving routines
*/

void SaveUInt1( UInt1 data )
{
  SAVE_BYTE( data );
}

UInt1 LoadUInt1( void )
{
  return LOAD_BYTE( );
}

void SaveUInt2( UInt2 data )
{
  SAVE_BYTE( (UInt1) (data & 0xFF) );
  SAVE_BYTE( (UInt1) (data >> 8) );
}

UInt2 LoadUInt2 ( void )
{
  UInt2 res;
  res = (UInt2)LOAD_BYTE();
  res |= (UInt2)LOAD_BYTE()<<8;
  return res;
}

void SaveUInt4( UInt4 data )
{
  SAVE_BYTE( (UInt1) (data & 0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 8) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 16) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 24) &0xFF) );
}

UInt4 LoadUInt4 ( void )
{
  UInt4 res;
  res = (UInt)LOAD_BYTE();
  res |= (UInt)LOAD_BYTE() << 8;
  res |= (UInt)LOAD_BYTE() << 16;
  res |= (UInt)LOAD_BYTE() << 24;
  return res;
}


#ifdef SYS_IS_64_BIT

void SaveUInt8( UInt8 data )
{
  SAVE_BYTE( (UInt1) (data & 0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 8) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 16) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 24) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 32) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 40) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 48) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 56) &0xFF) );
}

UInt8 LoadUInt8 ( void )
{
  UInt8 res;
  res = (UInt)LOAD_BYTE();
  res |= (UInt)LOAD_BYTE() << 8;
  res |= (UInt)LOAD_BYTE() << 16;
  res |= (UInt)LOAD_BYTE() << 24;
  res |= (UInt)LOAD_BYTE() << 32;
  res |= (UInt)LOAD_BYTE() << 40;
  res |= (UInt)LOAD_BYTE() << 48;
  res |= (UInt)LOAD_BYTE() << 56;

  return res;
}


void SaveUInt( UInt data )
{
  SAVE_BYTE( (UInt1) (data & 0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 8) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 16) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 24) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 32) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 40) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 48) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 56) &0xFF) );
}

UInt8 LoadUInt ( void )
{
  UInt res;
  res = (UInt)LOAD_BYTE();
  res |= (UInt)LOAD_BYTE() << 8;
  res |= (UInt)LOAD_BYTE() << 16;
  res |= (UInt)LOAD_BYTE() << 24;
  res |= (UInt)LOAD_BYTE() << 32;
  res |= (UInt)LOAD_BYTE() << 40;
  res |= (UInt)LOAD_BYTE() << 48;
  res |= (UInt)LOAD_BYTE() << 56;

  return res;
}

#else

void SaveUInt( UInt data )
{
  SAVE_BYTE( (UInt1) (data & 0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 8) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 16) &0xFF) );
  SAVE_BYTE( (UInt1) ((data >> 24) &0xFF) );
}

UInt LoadUInt ( void )
{
  UInt res;
  res = (UInt)LOAD_BYTE();
  res |= (UInt)LOAD_BYTE() << 8;
  res |= (UInt)LOAD_BYTE() << 16;
  res |= (UInt)LOAD_BYTE() << 24;
  return res;
}

#endif

void SaveCStr( const Char * str)
{
  do {
    SAVE_BYTE( (UInt1) *str);
  } while (*(str++));
}

#include <assert.h>

void LoadCStr( Char *buf, UInt maxsize)
{
  UInt nread = 0;
  UInt1 c = 1;
  assert(maxsize > 0);
  while (c != '\0' && nread < maxsize )
    {
      c = LOAD_BYTE();
      *buf++ = (Char) c;
      nread++;
    }
  if (c != '\0')
    {
      Pr("Buffer overflow reading workspace\n",0L,0L);
      SyExit(1);
    }
}


/****************************************************************************
**
*F  SaveString( <string> )  . . . . . . . . . . . . . . . . . . save a string
**
*/
void SaveString ( Obj string )
{
  UInt i, len = GET_LEN_STRING(string);
  UInt1 *p = (UInt1*)CHARS_STRING(string);
  SaveUInt(len);
  for (i=0; i<len; i++)
    SAVE_BYTE(p[i]);
}

/****************************************************************************
**
*F  LoadString( <string> )
**
*/
void LoadString ( Obj string )
{
  UInt i, len;
  UInt1 c;
  UInt1 *p = (UInt1*)CHARS_STRING(string);
  len = LoadUInt();
  SET_LEN_STRING(string, len);
  for (i=0; i<len; i++) {
    c = LOAD_BYTE();
    p[i] = c;
  }
}

void SaveSubObj( Obj subobj )
{
  if (!subobj)
    SaveUInt(0);
  else if (IS_INTOBJ(subobj))
    SaveUInt((UInt) subobj);
  else if (IS_FFE(subobj))
    SaveUInt((UInt) subobj);
  else if ((((UInt)subobj & 3) != 0) ||
           subobj < (Bag)MptrBags ||
           subobj > (Bag)OldBags ||
           (Bag *)PTR_BAG(subobj) < OldBags)
    {
      Pr("#W bad bag id %d found, 0 saved\n", (Int)subobj, 0L);
      SaveUInt(0);
    }
  else
    SaveUInt(((UInt)((PTR_BAG(subobj))[-1])) << 2);

}

Obj LoadSubObj()
{
  UInt word = LoadUInt();
  if (word == 0)
    return (Obj) 0;
  if ((word & 0x3) == 1 || (word & 0x3) == 2)
    return (Obj) word;
  else
    return (Obj)(MptrBags + (word >> 2));
}

void SaveHandler( ObjFunc hdlr )
{
  if (hdlr == (ObjFunc)0)
    SaveCStr("");
  else
    SaveCStr((const Char *)CookieOfHandler(hdlr));
}


ObjFunc LoadHandler( )
{
  Char buf[256];
  LoadCStr(buf, 256);
  if (buf[0] == '\0')
    return (ObjFunc) 0;
  else
    return HandlerOfCookie(buf);
}


void SaveDouble( Double d)
{
  union { UInt1 buf[sizeof(Double)]; double d; } v;
  UInt i;
  v.d = d;
  for (i = 0; i < sizeof(Double); i++)
    SAVE_BYTE(v.buf[i]);
}

Double LoadDouble( void)
{
  union { UInt1 buf[sizeof(Double)]; double d; } v;
  UInt i;
  for (i = 0; i < sizeof(Double); i++)
    v.buf[i] = LOAD_BYTE();
  return v.d;
}

/***************************************************************************
**
**  Bag level saving routines
*/



static void SaveBagData (Bag bag)
{
  SaveUInt((UInt)PTR_BAG(bag)[-3]);
  SaveUInt((UInt)PTR_BAG(bag)[-2]);


  /* dispatch */
  (*(SaveObjFuncs[ TNUM_BAG( bag )]))(bag);

}



static void LoadBagData ( )
{
  Bag bag;
#ifndef NEWSHAPE
  UInt sizetype;
#else
  UInt size;
  UInt type;
#endif

  /* Recover the sizetype word */
#ifndef NEWSHAPE
  sizetype=LoadUInt();
#else
  type = LoadUInt();
  size = LoadUInt();
#endif

#ifdef DEBUG_LOADING
  {
    if (InfoBags[type].name == NULL)
      {
        Pr("Bad type %d, size %d\n",type,size);
        exit(1);
      }
  }

#endif

  /* Get GASMAN to set up the bag for me */
#ifndef NEWSHAPE
  bag = NextBagRestoring( sizetype );
#else
  bag = NextBagRestoring( size, type );
#endif

  /* despatch */
#ifndef NEWSHAPE
  (*(LoadObjFuncs[ sizetype & 0xFFL ]))(bag);
#else
  (*(LoadObjFuncs[ type ]))(bag);
#endif

  return;
}


/***************************************************************************
**

*F  WriteSaveHeader() . . . . .  and utility functions, and loading functions
**
*/

static void WriteEndiannessMarker( void )
{
  UInt x;
#ifdef SYS_IS_64_BIT
  x = 0x0102030405060708L;
#else
  x = 0x01020304L;
#endif
  SAVE_BYTE(((UInt1 *)&x)[0]);
  SAVE_BYTE(((UInt1 *)&x)[1]);
  SAVE_BYTE(((UInt1 *)&x)[2]);
  SAVE_BYTE(((UInt1 *)&x)[3]);
#if SYS_IS_64_BIT
  SAVE_BYTE(((UInt1 *)&x)[4]);
  SAVE_BYTE(((UInt1 *)&x)[5]);
  SAVE_BYTE(((UInt1 *)&x)[6]);
  SAVE_BYTE(((UInt1 *)&x)[7]);
#endif
}

static void CheckEndiannessMarker( void )
{
  UInt x;
  ((UInt1 *)&x)[0] = LOAD_BYTE();
  ((UInt1 *)&x)[1] = LOAD_BYTE();
  ((UInt1 *)&x)[2] = LOAD_BYTE();
  ((UInt1 *)&x)[3] = LOAD_BYTE();
#ifdef SYS_IS_64_BIT
  ((UInt1 *)&x)[4] = LOAD_BYTE();
  ((UInt1 *)&x)[5] = LOAD_BYTE();
  ((UInt1 *)&x)[6] = LOAD_BYTE();
  ((UInt1 *)&x)[7] = LOAD_BYTE();
  if (x != 0x0102030405060708L)
#else
  if (x != 0x01020304L)
#endif
    {
      Pr("Saved workspace with incompatible byte order\n",0L,0L);
      SyExit(1);
    }
}


/***************************************************************************
**
**  BagStats
*/

static FILE *file;

static void report( Bag bag)
{
  fprintf(file,"%li %li\n", (Int) TNUM_BAG(bag), (Int) SIZE_BAG(bag));
}

Obj BagStats(Obj self, Obj filename)
{
  file = fopen((Char *)CHARS_STRING(filename),"w");
  CallbackForAllBags(report);
  fclose(file);
  return (Obj) 0;
}

/***************************************************************************
**
**  Find Bags -- a useful debugging tool -- scan for a bag of specified
**   type and size and return it to the GAP level. Could be a problem
**  if the bag is not a valid GAP object -- eg a local variables bag or
**  a functions body.
*/


static UInt fb_minsize, fb_maxsize, fb_tnum;
static Bag hit;

static void ScanBag( Bag bag)
{
  if (hit == (Bag)0 &&
      SIZE_BAG(bag) >= fb_minsize &&
      SIZE_BAG(bag) <= fb_maxsize &&
      TNUM_BAG(bag) == fb_tnum)
    hit = bag;
  return;
}

Obj FuncFindBag( Obj self, Obj minsize, Obj maxsize, Obj tnum )
{
  hit = (Bag) 0;
  fb_minsize = INT_INTOBJ(minsize);
  fb_maxsize = INT_INTOBJ(maxsize);
  fb_tnum = INT_INTOBJ(tnum);
  CallbackForAllBags(ScanBag);
  return (hit != (Bag) 0) ? hit : Fail;
}


/***************************************************************************
**
*F  SaveWorkspace( <fname> )  . . . . .  save the workspace to the named file
**
**  'SaveWorkspace' is the entry point to the workspace saving. It is not
**  installed as a GAP function, but instead as a keyword, so that we can be
**  sure it is only being called from the top-most prompt level
**  The file saveload.tex in the dev directory describes the saved format
**  in more detail. Most of the work will be done from inside GASMAN, because
**  we need to fiddle with Bag internals somewhat
**
**  The return value is either True or Fail
*/

static UInt NextSaveIndex;

static void AddSaveIndex( Bag bag)
{
  PTR_BAG(bag)[-1] = (Obj)NextSaveIndex++;
}

/*  is this obsolete ???
static void CheckPlist( Bag bag)
{
  if (TNUM_BAG(bag) == 14 && sizeof(UInt)*((UInt)(PTR_BAG(bag)[0])) > SIZE_BAG(bag))
    Pr("Panic %d %d\n",sizeof(UInt)*((UInt)(PTR_BAG(bag)[0])), SIZE_BAG(bag));
  return;
}
*/

static void RemoveSaveIndex( Bag bag)
{
  PTR_BAG(bag)[-1] = bag;
}
static void WriteSaveHeader( void )
{
  UInt i;
  UInt globalcount = 0;

  SaveCStr("GAP workspace");
  SaveCStr(SyKernelVersion);

#ifdef SYS_IS_64_BIT
  SaveCStr("64 bit");
#else
  SaveCStr("32 bit");
#endif

  WriteEndiannessMarker();

  SaveCStr("Counts and Sizes");
  SaveUInt(NrModules - NrBuiltinModules);
  for (i = 0; i < GlobalBags.nr; i++)
    if (GlobalBags.cookie[i] != NULL)
      globalcount++;
  SaveUInt(globalcount);
  SaveUInt(NextSaveIndex);
  SaveUInt(AllocBags - OldBags);

  SaveCStr("Loaded Modules");

  for ( i = NrBuiltinModules; i < NrModules; i++)
    {
      SaveUInt(Modules[i]->type);
      SaveUInt(Modules[i]->isGapRootRelative);
      SaveCStr(Modules[i]->filename);
    }

  SaveCStr("Kernel to WS refs");
  for (i = 0; i < GlobalBags.nr; i++)
    {
      if (GlobalBags.cookie[i] != NULL)
	{
	  SaveCStr((const Char *)GlobalBags.cookie[i]);
	  SaveSubObj(*(GlobalBags.addr[i]));
	}
    }
}

static Obj ProtectFname;

Obj SaveWorkspace( Obj fname )
{

  Int i;
  if (!IsStringConv(fname))
    ErrorQuit("usage: SaveWorkspace( <filename> )",0,0);

  for (i = 0; i < NrModules; i++)
    if (Modules[i]->preSave != NULL &&
        (*(Modules[i]->preSave))(Modules[i]))
      {
        Pr("Failed to save workspace -- problem reported in %s\n",
           (Int)Modules[i]->name, 0L);
        for ( i--; i >= 0; i--)
          (*(Modules[i]->postSave))(Modules[i]);
        return Fail;
      }

  /* For some reason itanium GC seems unable to spot fname */
  ProtectFname = fname;
  /* Do a full garbage collection */
  CollectBags( 0, 1);

  ProtectFname = (Obj)0L;

  /* Add indices in link words of all bags, for saving inter-bag references */
  NextSaveIndex = 0;
  CallbackForAllBags( AddSaveIndex );

  /* Now do the work */
  if (!OpenForSave( fname ))
    {
      WriteSaveHeader();
      SaveCStr("Bag data");
      SortHandlers( 1 ); /* Sort by address to speed up CookieOfHandler */
      CallbackForAllBags( SaveBagData );
      CloseAfterSave();
    }

  /* Finally, reset all the link words */
  CallbackForAllBags( RemoveSaveIndex );

  /* Restore situation by calling all post-save methods */
  for (i = 0; i < NrModules; i++)
    if (Modules[i]->postSave != NULL)
      (*(Modules[i]->postSave))(Modules[i]);

  return True;
}




/***************************************************************************
**
*F  LoadWorkspace( <fname> )  . . . . .  load the workspace to the named file
**
**  'LoadWorkspace' is the entry point to the workspace saving. It is not
**  installed as a GAP function, but instead called from InitGap when the
**  -L commad-line flag is given
**
**  The file saveload.tex in the dev directory describes the saved format
**  in more detail. Most of the work will be done from inside GASMAN, because
**  we need to fiddle with Bag internals somewhat
**
*/


void LoadWorkspace( Char * fname )
{
  UInt nMods, nGlobs, nBags, i, maxSize, isGapRootRelative;
  UInt globalcount = 0;
  Char buf[256];
  Obj * glob;

  /* Open saved workspace  */
  OpenForLoad( fname );

  /* Check file header */

  LoadCStr(buf,256);
  if (SyStrncmp (buf, "GAP ", 4) != 0) {
     Pr("File %s does not appear to be a GAP workspae.\n", (long) fname, 0L);
     SyExit(1);
  }

  if (SyStrcmp (buf, "GAP workspace") == 0) {

     LoadCStr(buf,256);
     if (SyStrcmp (buf, SyKernelVersion) != 0) {
        Pr("This workspace is not compatible with GAP kernel (%s, present: %s).\n",
           (long)buf, (long)SyKernelVersion);
        SyExit(1);
     }

     LoadCStr(buf,256);
#ifdef SYS_IS_64_BIT
     if (SyStrcmp(buf,"64 bit") != 0)
#else
     if (SyStrcmp(buf,"32 bit") != 0)
#endif
        {
           Pr("This workspace was created by a %s version of GAP.\n", (long)buf, 0L);
           SyExit(1);
        }
  } else {

     /* try if it is an old workspace */

#ifdef SYS_IS_64_BIT
     if (SyStrcmp(buf,"GAP 4.0 beta 64 bit") != 0)
#else
     if (SyStrcmp(buf,"GAP 4.0 beta 32 bit") != 0)
#endif
        Pr("File %s probably isn't a GAP workspace.\n", (long)fname, 0L);
     else
        Pr("This workspace was created by an old version of GAP.\n", 0L, 0L);
     SyExit(1);
  }

  CheckEndiannessMarker();

  LoadCStr(buf,256);
  if (SyStrcmp(buf,"Counts and Sizes") != 0)
    {
      Pr("Bad divider\n",0L,0L);
      SyExit(1);
    }

  nMods = LoadUInt();
  nGlobs = LoadUInt();
  nBags = LoadUInt();
  maxSize = LoadUInt();

  /* Make sure there is enough room, and signal GASMAN that
     we are starting a restore */
  StartRestoringBags(nBags, maxSize);

  /* The restoring kernel must have at least as many compiled modules
     as the saving one. */
  LoadCStr(buf,256);
  if (SyStrcmp(buf,"Loaded Modules") != 0)
    {
      Pr("Bad divider\n",0L,0L);
      SyExit(1);
    }

  for (i = 0; i < nMods; i++)
    {
      UInt type = LoadUInt();
      isGapRootRelative = LoadUInt();
      LoadCStr(buf,256);
      if (isGapRootRelative)
        READ_GAP_ROOT( buf);
      else
	{
	  StructInitInfo *info = NULL;
 	  /* Search for user module static case first */
	  if (type == MODULE_STATIC) {
	    UInt k;
	    for ( k = 0;  CompInitFuncs[k];  k++ ) {
	      info = (*(CompInitFuncs[k]))();
	      if ( info == 0 ) {
		continue;
	      }
	      if ( ! SyStrcmp( buf, info->name ) ) {
		break;
	      }
	    }
	    if ( CompInitFuncs[k] == 0 ) {
	      Pr( "Static module %s not found in loading kernel\n",
		  (Int)buf, 0L );
	      SyExit(1);
	    }

	  } else {
	    /* and dynamic case */
	    InitInfoFunc init;

	    init = SyLoadModule(buf);

	    if ((Int)init == 1 || (Int) init == 3 || (Int) init == 5 || (Int) init == 7)
	      {
		Pr("Failed to load needed dynamic module %s, error code %d\n",
		   (Int)buf, (Int) init);
		SyExit(1);
	      }
	    info = (*init)();
	     if (info == 0 )
	       {
		Pr("Failed to init needed dynamic module %s, error code %d\n",
		   (Int)buf, (Int) info);
		SyExit(1);
	       }
	  }
	/* link and init me                                                    */
	info->isGapRootRelative = 0;
	(info->initKernel)(info);
	RecordLoadedModule(info, buf);
      }

    }

  /* Now the kernel variables that point into the workspace */
  LoadCStr(buf,256);
  if (SyStrcmp(buf,"Kernel to WS refs") != 0)
    {
      Pr("Bad divider\n",0L,0L);
      SyExit(1);
    }
  SortGlobals(2);               /* globals by cookie for quick
                                 lookup */
  for (i = 0; i < GlobalBags.nr; i++)
    {
      if (GlobalBags.cookie[i] != NULL)
	globalcount++;
      else
	*(GlobalBags.addr[i]) = (Bag) 0;
    }
  if (nGlobs != globalcount)
    {
      Pr("Wrong number of global bags in saved workspace %d %d\n",
         nGlobs, globalcount);
      SyExit(1);
    }
  for (i = 0; i < globalcount; i++)
    {
      LoadCStr(buf,256);
      glob = GlobalByCookie(buf);
      if (glob == (Obj *)0)
        {
          Pr("Global object cookie from workspace not found in kernel %s\n",
             (Int)buf,0L);
          SyExit(1);
        }
      *glob = LoadSubObj();
#ifdef DEBUG_LOADING
      Pr("Restored global %s\n",(Int)buf,0L);
#endif
    }

  LoadCStr(buf,256);
  if (SyStrcmp(buf,"Bag data") != 0)
    {
      Pr("Bad divider\n",0L,0L);
      SyExit(1);
    }

  SortHandlers(2);
  for (i = 0; i < nBags; i++)
    LoadBagData();

  FinishedRestoringBags();

  /* Post restore methods are called elsewhere */

  CloseAfterLoad();
}

#include        "finfield.h"            /* finite fields and ff elements   */

static void PrSavedObj( UInt x)
{
  if ((x & 3) == 1)
    Pr("Immediate  integer %d\n", INT_INTOBJ((Obj)x),0L);
  else if ((x & 3) == 2)
    Pr("Immedate FFE %d %d\n", VAL_FFE(x), SIZE_FF(FLD_FFE(x)));
  else
    Pr("Reference to bag number %d\n",x>>2,0L);
}

Obj FuncDumpWorkspace( Obj self, Obj fname )
{
  UInt nMods, nGlobs, nBags, i, relative;
  Char buf[256];
  OpenForLoad( CSTR_STRING(fname) );
  LoadCStr(buf,256);
  Pr("Header string: %s\n",(Int) buf, 0L);
#ifdef SYS_IS_64_BIT
  if (SyStrcmp(buf,"GAP 4.0 beta 64 bit") != 0)
#else
  if (SyStrcmp(buf,"GAP 4.0 beta 32 bit") != 0)
#endif
    ErrorQuit("Header is bad",0L,0L);
  CheckEndiannessMarker();
  LoadCStr(buf,256);
  Pr("Divider string: %s\n",(Int)buf,0L);
  if (SyStrcmp(buf,"Counts and Sizes") != 0)
    ErrorQuit("Bad divider",0L,0L);
  Pr("Loaded modules: %d\n",nMods = LoadUInt(), 0L);
  Pr("Global Bags   : %d\n",nGlobs = LoadUInt(), 0L);
  Pr("Total Bags    : %d\n",nBags = LoadUInt(), 0L);
  Pr("Maximum Size  : %d\n",sizeof(Bag)*LoadUInt(), 0L);
  LoadCStr(buf,256);
  Pr("Divider string: %s\n",(Int)buf, 0L);
  if (SyStrcmp(buf,"Loaded Modules") != 0)
    ErrorQuit("Bad divider",0L,0L);
  for (i = 0; i < nMods; i++)
    {
      relative = LoadUInt();
      if (relative)
	Pr("GAP root relative ", 0L, 0L);
      else
	Pr("absolute ", 0L, 0L);
      LoadCStr(buf,256);
      Pr("  %s\n",(Int)buf,0L);
    }
  LoadCStr(buf,256);
  Pr("Divider string: %s\n",(Int)buf,0L);
  if (SyStrcmp(buf,"Kernel to WS refs") != 0)
    ErrorQuit("Bad divider",0L,0L);
  for (i = 0; i < nGlobs; i++)
    {
      LoadCStr(buf,256);
      Pr("  %s ",(Int)buf,0L);
      PrSavedObj(LoadUInt());
    }
  LoadCStr(buf,256);
  Pr("Divider string: %s\n",(Int)buf,0L);
  if (SyStrcmp(buf,"Bag data") != 0)
    ErrorQuit("Bad divider",0L,0L);
  CloseAfterLoad();
  return (Obj) 0;
}


/****************************************************************************
**

*F * * * * * * * * * * * * * initialize package * * * * * * * * * * * * * * *
*/


/****************************************************************************
**

*V  GVarFuncs . . . . . . . . . . . . . . . . . . list of functions to export
*/
static StructGVarFunc GVarFuncs [] = {

    { "DumpWorkspace", 1, "fname",
      FuncDumpWorkspace, "src/saveload.c:DumpWorkspace" },

    { "FindBag", 3, "minsize, maxsize, tnum",
      FuncFindBag, "src/saveload.c:FindBag" },

    { "BagStats", 1, "filename",
      BagStats, "src/saveload.c:BagStats" },

    { 0 }

};


/****************************************************************************
**

*F  InitKernel( <module> )  . . . . . . . . initialise kernel data structures
*/
static Int InitKernel (
    StructInitInfo *    module )
{
    InitGlobalBag(&ProtectFname, "Protected Filename for SaveWorkspace");

    /* init filters and functions                                          */
    InitHdlrFuncsFromTable( GVarFuncs );

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
    /* Create dummy variable, to support tab-completion */
    (void)GVarName("SaveWorkspace");

    /* init filters and functions                                          */
    InitGVarFuncsFromTable( GVarFuncs );

    /* return success                                                      */
    return 0;
}


/****************************************************************************
**
*F  InitInfoSaveLoad()  . . . . . . . . . . . . . . . table of init functions
*/
static StructInitInfo module = {
    MODULE_BUILTIN,                     /* type                           */
    "saveload",                         /* name                           */
    0,                                  /* revision entry of c file       */
    0,                                  /* revision entry of h file       */
    0,                                  /* version                        */
    0,                                  /* crc                            */
    InitKernel,                         /* initKernel                     */
    InitLibrary,                        /* initLibrary                    */
    0,                                  /* checkInit                      */
    0,                                  /* preSave                        */
    0,                                  /* postSave                       */
    0                                   /* postRestore                    */
};

StructInitInfo * InitInfoSaveLoad ( void )
{
    module.revision_c = Revision_saveload_c;
    module.revision_h = Revision_saveload_h;
    FillInVersion( &module );
    return &module;
}


/****************************************************************************
**

*E  saveload.c  . . . . . . . . . . . . . . . . . . . . . . . . . . ends here
*/
