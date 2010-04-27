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

#include <stdlib.h>

#include "graphK23Search.private.h"
#include "graphK23Search.h"

extern int  _SearchForK23(graphP theGraph, int I);

extern int  _TestForK23GraphObstruction(graphP theGraph, int *degrees, int *imageVerts);
extern int  _getImageVertices(graphP theGraph, int *degrees, int maxDegree,
                              int *imageVerts, int maxNumImageVerts);
extern int  _TestSubgraph(graphP theSubgraph, graphP theGraph);

/* Forward declarations of overloading functions */

int  _K23Search_HandleBlockedEmbedIteration(graphP theGraph, int I);
int  _K23Search_EmbedPostprocess(graphP theGraph, int I, int edgeEmbeddingResult);
int  _K23Search_CheckEmbeddingIntegrity(graphP theGraph, graphP origGraph);
int  _K23Search_CheckObstructionIntegrity(graphP theGraph, graphP origGraph);

/* Forward declarations of functions used by the extension system */

void *_K23Search_DupContext(void *pContext, void *theGraph);
void _K23Search_FreeContext(void *);

/****************************************************************************
 * K23SEARCH_ID - the variable used to hold the integer identifier for this
 * extension, enabling this feature's extension context to be distinguished
 * from other features' extension contexts that may be attached to a graph.
 ****************************************************************************/

int K23SEARCH_ID = 0;

/****************************************************************************
 gp_AttachK23Search()

 This function adjusts the graph data structure to attach the K2,3 search
 feature.
 ****************************************************************************/

int  gp_AttachK23Search(graphP theGraph)
{
     K23SearchContext *context = NULL;

     // If the K2,3 search feature has already been attached to the graph
     // then there is no need to attach it again
     gp_FindExtension(theGraph, K23SEARCH_ID, (void *)&context);
     if (context != NULL)
     {
         return OK;
     }

     // Allocate a new extension context
     context = (K23SearchContext *) malloc(sizeof(K23SearchContext));
     if (context == NULL)
     {
         return NOTOK;
     }

     // Put the overload functions into the context function table.
     // gp_AddExtension will overload the graph's functions with these, and
     // return the base function pointers in the context function table
     memset(&context->functions, 0, sizeof(graphFunctionTable));

     context->functions.fpHandleBlockedEmbedIteration = _K23Search_HandleBlockedEmbedIteration;
     context->functions.fpEmbedPostprocess = _K23Search_EmbedPostprocess;
     context->functions.fpCheckEmbeddingIntegrity = _K23Search_CheckEmbeddingIntegrity;
     context->functions.fpCheckObstructionIntegrity = _K23Search_CheckObstructionIntegrity;

     // Store the K23 search context, including the data structure and the
     // function pointers, as an extension of the graph
     if (gp_AddExtension(theGraph, &K23SEARCH_ID, (void *) context,
                         _K23Search_DupContext, _K23Search_FreeContext,
                         &context->functions) != OK)
     {
         _K23Search_FreeContext(context);
         return NOTOK;
     }

     return OK;
}

/********************************************************************
 gp_DetachK23Search()
 ********************************************************************/

int gp_DetachK23Search(graphP theGraph)
{
    return gp_RemoveExtension(theGraph, K23SEARCH_ID);
}

/********************************************************************
 _K23Search_DupContext()
 ********************************************************************/

void *_K23Search_DupContext(void *pContext, void *theGraph)
{
     K23SearchContext *context = (K23SearchContext *) pContext;
     K23SearchContext *newContext = (K23SearchContext *) malloc(sizeof(K23SearchContext));

     if (newContext != NULL)
     {
         *newContext = *context;
     }

     return newContext;
}

/********************************************************************
 _K23Search_FreeContext()
 ********************************************************************/

void _K23Search_FreeContext(void *pContext)
{
     free(pContext);
}

/********************************************************************
 ********************************************************************/

int  _K23Search_HandleBlockedEmbedIteration(graphP theGraph, int I)
{
    if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK23)
        return _SearchForK23(theGraph, I);

    else
    {
        K23SearchContext *context = NULL;
        gp_FindExtension(theGraph, K23SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpHandleBlockedEmbedIteration(theGraph, I);
        }
    }

    return NOTOK;
}

/********************************************************************
 ********************************************************************/

int  _K23Search_EmbedPostprocess(graphP theGraph, int I, int edgeEmbeddingResult)
{
     // For K2,3 search, we just return the edge embedding result because the
     // search result has been obtained already.
     if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK23)
     {
         return edgeEmbeddingResult;
     }

     // When not searching for K2,3, we let the superclass do the work
     else
     {
        K23SearchContext *context = NULL;
        gp_FindExtension(theGraph, K23SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpEmbedPostprocess(theGraph, I, edgeEmbeddingResult);
        }
     }

     return NOTOK;
}

/********************************************************************
 ********************************************************************/

int  _K23Search_CheckEmbeddingIntegrity(graphP theGraph, graphP origGraph)
{
     if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK23)
     {
         return OK;
     }

     // When not searching for K2,3, we let the superclass do the work
     else
     {
        K23SearchContext *context = NULL;
        gp_FindExtension(theGraph, K23SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpCheckEmbeddingIntegrity(theGraph, origGraph);
        }
     }

     return NOTOK;
}

/********************************************************************
 ********************************************************************/

int  _K23Search_CheckObstructionIntegrity(graphP theGraph, graphP origGraph)
{
     // When searching for K2,3, we ensure that theGraph is a subgraph of
     // the original graph and that it contains a K2,3 homeomorph
     if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK23)
     {
         int  degrees[4], imageVerts[5];

         if (_TestSubgraph(theGraph, origGraph) != TRUE)
             return NOTOK;

         if (_getImageVertices(theGraph, degrees, 3, imageVerts, 5) != OK)
             return NOTOK;

         if (_TestForK23GraphObstruction(theGraph, degrees, imageVerts) == TRUE)
         {
             return OK;
         }

         return NOTOK;
     }

     // When not searching for K2,3, we let the superclass do the work
     else
     {
        K23SearchContext *context = NULL;
        gp_FindExtension(theGraph, K23SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpCheckObstructionIntegrity(theGraph, origGraph);
        }
     }

     return NOTOK;
}
