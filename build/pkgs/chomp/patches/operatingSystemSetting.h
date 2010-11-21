/// @addtogroup capd
/// @{

/////////////////////////////////////////////////////////////////////////////
/// @file operatingSystemSetting.h
///
/// @author The CAPD Group
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 2000-2005 by the CAPD Group.
//
// This file constitutes a part of the CAPD library,
// distributed under the terms of the GNU General Public License.
// Consult  http://capd.wsb-nlu.edu.pl/ for details.

// This file controls the operating system setting in which the krak and/or inerval
// packages are compiled.
// Note that you don't have to modify this file - most settings can be
// selected with the -D argument of your command-line compiler
// or in your project. However, if the auto-detection does not work
// in your case, you may want to modify some specific options.


#ifndef _CAPD_CAPD_OPERATINGSYSTEMSETTING_H_
#define _CAPD_CAPD_OPERATINGSYSTEMSETTING_H_

// ################ Manual operating system selection #######################
// Uncomment exactly one of the following defines if the automatic system
// selection doesn't work for you.
//
//#define DOS     // if working in DOS without Windows
//#define WIN31   // if working in Windows 3.1
//#define WIN95   // if working in Windows 95 or later
//#define SUN_OS  // if working in X Windows on SUN or SGI workstations
//#define LINUX   // if working on a PC under LINUX
//
// Some other old choices are: MAC, DOS_CH, SUN_CH, SUN_V but they may cause
// serious problems.



// ################## Automatic operating system selection ##################
// Determine the target platform. Note that the code in the entire package
// relies on the compiler settings like WIN95, LINUX, etc.

// auto-detect if not set before
#if !defined (DOS) && !defined (WIN31) && !defined (WIN95) && \
 !defined (SUN_OS) && !defined (LINUX) && !defined (MAC) && \
 !defined (DOS_CH) && !defined (SUN_CH) && !defined (SUN_V)

// Is the target operating system Windows 95 (or a newer Windows)?
#if defined (__WIN32__) || defined (_WIN32) || defined (WIN32)
#  define WIN95

// Is the target operating system Linux? Note: Linux on Sparc is not supported!
#elif defined (__linux__) || defined (__APPLE__)
#  define LINUX
#  ifdef __sparc__
#    error The Interval library does not work on the Linux/Sparc combination.
#  endif

// If 'sun' is defined, which it is by both Sun Studio and gcc, then assume Solaris.
#elif defined (sun)
#  define SUN_OS
#endif

#endif // auto detect operating system

// if the system has not been determined, you can't compile the package
#if !defined (DOS) && !defined (WIN31) && !defined (WIN95) && \
 !defined (SUN_OS) && !defined (LINUX) && !defined (MAC) && \
 !defined (DOS_CH) && !defined (SUN_CH) && !defined (SUN_V)
#error Your system cannot be determined automatically by the KRAK package.\
Please edit the configuration file appropriately.
#endif

// krak and interval packages use this in some low level routines
#if defined(DOS) || defined(WIN95)
  #define INTEL
  #define IBM
#endif


#endif // _CAPD_CAPD_OPERATINGSYSTEMSETTING_H_

/// @}
