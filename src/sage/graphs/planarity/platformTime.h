/********************************************************************
Copyright 2005 John M. Boyer

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
 ********************************************************************/

#ifndef PLATFORM_TIME
#define PLATFORM_TIME

#ifdef WIN32

#include <windows.h>
#include <winbase.h>
#define platform_time DWORD
#define platform_GetTime GetTickCount
#define platform_GetDuration(startTime, endTime) ((double) (endTime-startTime) / 1000.0)

#else

#include <time.h>
#define platform_time time_t
#define platform_GetTime time(NULL)
#define platform_GetDuration(startTime, endTime) ((double) (endTime - startTime))

#endif

#endif
