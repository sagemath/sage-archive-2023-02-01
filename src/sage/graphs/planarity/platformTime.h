#ifndef PLATFORM_TIME
#define PLATFORM_TIME

/*
Planarity-Related Graph Algorithms Project
Copyright (c) 1997-2010, John M. Boyer
All rights reserved. Includes a reference implementation of the following:

* John M. Boyer. "Simplified O(n) Algorithms for Planar Graph Embedding,
  Kuratowski Subgraph Isolation, and Related Problems". Ph.D. Dissertation,
  University of Victoria, 2001.

* John M. Boyer and Wendy J. Myrvold. "On the Cutting Edge: Simplified O(n)
  Planarity by Edge Addition". Journal of Graph Algorithms and Applications,
  Vol. 8, No. 3, pp. 241-273, 2004.

* John M. Boyer. "A New Method for Efficiently Generating Planar Graph
  Visibility Representations". In P. Eades and P. Healy, editors,
  Proceedings of the 13th International Conference on Graph Drawing 2005,
  Lecture Notes Comput. Sci., Volume 3843, pp. 508-511, Springer-Verlag, 2006.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

* Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright notice, this
  list of conditions and the following disclaimer in the documentation and/or
  other materials provided with the distribution.

* Neither the name of the Planarity-Related Graph Algorithms Project nor the names
  of its contributors may be used to endorse or promote products derived from this
  software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#ifdef WIN32

#include <windows.h>
#include <winbase.h>

#define platform_time DWORD
#define platform_GetTime(timeVar) (timeVar = GetTickCount())
#define platform_GetDuration(startTime, endTime) ((double) (endTime-startTime) / 1000.0)

#else

#include <time.h>

typedef struct {
	clock_t hiresTime;
	time_t lowresTime;
} platform_time;

#define platform_GetTime(timeVar) (timeVar.hiresTime = clock(), timeVar.lowresTime = time(NULL))

// Many flavors of Unix have CLOCKS_PER_SEC at 1 million, and clock_t as a 4 byte long integer
// which means that the clock() construct has a resolution of only about 2000 seconds
// If we're getting a duration longer than that, then we fall back to the coarser time() measure

#define platform_GetDuration(startTime, endTime) ( \
		( (double) (endTime.lowresTime - startTime.lowresTime) ) > 2000 ? \
		( (double) (endTime.lowresTime - startTime.lowresTime) ) : \
		( (double) (endTime.hiresTime - startTime.hiresTime)) / CLOCKS_PER_SEC)

/*
#define platform_time clock_t
#define platform_GetTime() clock()
#define platform_GetDuration(startTime, endTime) (((double) (endTime - startTime)) / CLOCKS_PER_SEC)
*/

/*
#define platform_time time_t
#define platform_GetTime() time((time_t *)NULL)
#define platform_GetDuration(startTime, endTime) ((double) (endTime - startTime))
*/

#endif

#endif
