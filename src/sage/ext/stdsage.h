/*
General C (.h) code this is useful to include in any pyrex module.

These are mostly things that can't be done in Pyrex.
*/

/******************************************************************************
       Copyright (C) 2006 William Stein <wstein@gmail.com>

  Distributed under the terms of the GNU General Public License (GPL), Version 2.

  The full text of the GPL is available at:
                  http://www.gnu.org/licenses/

******************************************************************************/


/*****************************************
          For PARI
 *****************************************/

#define set_gel(x, n, z)  gel(x,n)=z;
