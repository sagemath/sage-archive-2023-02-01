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

#ifndef STACK_H
#define STACK_H

typedef struct
{
        int *S;
        int Top, Size;
} stack;

typedef stack * stackP;

stackP sp_New(int);
void sp_Free(stackP *);

int  sp_Copy(stackP, stackP);

#ifndef SPEED_MACROS

int  sp_ClearStack(stackP);

int  sp_IsEmpty(stackP);
int  sp_NonEmpty(stackP);

int  sp_Push(stackP, int);
int  sp_Push2(stackP, int, int);

#define sp_Pop(theStack, a) sp__Pop(theStack, &(a))
#define sp_Pop2(theStack, a, b) sp__Pop2(theStack, &(a), &(b))

int  sp__Pop(stackP, int *);
int  sp__Pop2(stackP, int *, int *);

#else

#define sp_ClearStack(theStack) theStack->Top=0

#define sp_IsEmpty(theStack) !theStack->Top
#define sp_NonEmpty(theStack) theStack->Top

#define sp_Push(theStack, a) theStack->S[theStack->Top++] = a
#define sp_Push2(theStack, a, b) sp_Push(theStack, a); sp_Push(theStack, b)

#define sp_Pop(theStack, a) a=theStack->S[--theStack->Top]
#define sp_Pop2(theStack, a, b) sp_Pop(theStack, b);sp_Pop(theStack, a)

#endif

#endif
