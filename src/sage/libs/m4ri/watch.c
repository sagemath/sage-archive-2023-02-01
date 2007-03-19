/******************************************************************************
*
*            M4RI: Method of the Four Russians Inversion
*
*       Copyright (C) 2007 Gregory Bard <gregory.bard@ieee.org>
*
*  Distributed under the terms of the GNU General Public License (GPL)
*
*    This code is distributed in the hope that it will be useful,
*    but WITHOUT ANY WARRANTY; without even the implied warranty of
*    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
*    General Public License for more details.
*
*  The full text of the GPL is available at:
*
*                  http://www.gnu.org/licenses/
******************************************************************************/

#include "watch.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>

clock_t watchstart;
clock_t watchstop;
long watchrunning;

double logtotal, logtotalsquare, logcount;

void seedWatch() {
  srand(clock()/10000);
}

double mysqrt(double n) {
  double guess=n/2;
  double miss, der;
  int i;

  if (n<0) {
    printf("\a Tried to take square root of %f.\n",n);
    exit(1);
  }

  for (i=0; i<40 ;i++) {
    miss=n-guess*guess;
    der=2.0*guess;
    guess+=miss/der;
  }

  return guess;
}

void startWatch() {
  watchstart=clock();
  watchrunning=YES;
}

void stopWatch() {
  watchstop=clock();
  watchrunning=NO;
}

clock_t getWatch() {
  if (watchrunning==YES)
    return clock()-watchstart;
  else
    return watchstop-watchstart;
}

void clearLogs() {
  logtotal=0;
  logtotalsquare=0;
  logcount=0;
}

void store(clock_t watch) {
  double value=(double)watch;

  logtotal+=value;
  logtotalsquare+=value*value;
  logcount++;
}

double getAverage() {
  return logtotal/(logcount*1.0);
}

double getSigma() {
  double mean = getAverage();
  double ex_sq = logtotalsquare/(logcount*1.0);
  double var = ex_sq-mean*mean;
  double sigma = mysqrt(var);
  return sigma;
}

long getCount() {
  return (long)logcount;
}




