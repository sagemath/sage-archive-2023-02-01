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

#include "appconst.h"
#include "stack.h"
#include <stdlib.h>

stackP sp_New(int Size)
{
stackP theStack;

     theStack = (stackP) malloc(sizeof(stack));

     if (theStack != NULL)
     {
         theStack->S = (int *) malloc(Size*sizeof(int));
         if (theStack->S == NULL)
         {
             free(theStack);
             theStack = NULL;
         }
     }

     if (theStack != NULL)
     {
         theStack->Size = Size;
         sp_ClearStack(theStack);
     }

     return theStack;
}

void sp_Free(stackP *pStack)
{
     if (pStack == NULL || *pStack == NULL) return;

     (*pStack)->Size = (*pStack)->Top = 0;

     if ((*pStack)->S != NULL)
          free((*pStack)->S);
     (*pStack)->S = NULL;
     free(*pStack);

     *pStack = NULL;
}

int  sp_Copy(stackP stackDst, stackP stackSrc)
{
stackP newStack = NULL;
int  I, *p;

     if (stackDst->Size == stackSrc->Size)
     {
         for (I=0; I < stackSrc->Top; I++)
              stackDst->S[I] = stackSrc->S[I];
     }

     else
     {
         newStack = sp_New(stackSrc->Size);
         if (newStack == NULL) return NOTOK;

         for (I=0; I < stackSrc->Top; I++)
              newStack->S[I] = stackSrc->S[I];

         p = stackDst->S;
         stackDst->S = newStack->S;
         newStack->S = p;
         newStack->Size = stackDst->Size;
         sp_Free(&newStack);
     }

     stackDst->Top = stackSrc->Top;
     stackDst->Size = stackSrc->Size;

     return OK;
}

#ifndef SPEED_MACROS

int  sp_ClearStack(stackP theStack)
{
     theStack->Top = 0;
     return OK;
}

int  sp_IsEmpty(stackP theStack)
{
     return !theStack->Top;
}

int  sp_NonEmpty(stackP theStack)
{
     return theStack->Top;
}

int  sp_Push(stackP theStack, int a)
{
//     if (theStack->Top >= theStack->Size)
//         return NOTOK;

     theStack->S[theStack->Top++] = a;
     return OK;
}

int  sp_Push2(stackP theStack, int a, int b)
{
//     if (theStack->Top + 1 >= theStack->Size)
//         return NOTOK;

     theStack->S[theStack->Top++] = a;
     theStack->S[theStack->Top++] = b;
     return OK;
}

int  sp__Pop(stackP theStack, int *pA)
{
//     if (theStack->Top <= 0)
//         return NOTOK;

     *pA = theStack->S[--theStack->Top];
     return OK;
}

int  sp__Pop2(stackP theStack, int *pA, int *pB)
{
//     if (theStack->Top <= 1)
//         return NOTOK;

     *pB = theStack->S[--theStack->Top];
     *pA = theStack->S[--theStack->Top];

     return OK;
}

#endif // SPEED_MACROS

