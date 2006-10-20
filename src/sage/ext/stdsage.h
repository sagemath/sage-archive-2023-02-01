/*
General C (.h) code this is useful to include in any pyrex module.

Put

include 'relative/path/to/stdsage.pxi'

at the top of your Pyrex file.

*/

/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/


/*****************************************
          Memory management

 *****************************************/

#define sage_malloc  malloc
#define sage_free    free
#define sage_realloc realloc

