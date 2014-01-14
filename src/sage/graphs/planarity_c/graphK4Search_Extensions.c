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

#include "graphK4Search.private.h"
#include "graphK4Search.h"

extern int  _GetNextVertexOnExternalFace(graphP theGraph, int curVertex, int *pPrevLink);

extern int  _SearchForK4InBicomps(graphP theGraph, int I);
extern int  _SearchForK4InBicomp(graphP theGraph, K4SearchContext *context, int I, int R);

extern int _TestForCompleteGraphObstruction(graphP theGraph, int numVerts,
                                            int *degrees, int *imageVerts);

extern int  _getImageVertices(graphP theGraph, int *degrees, int maxDegree,
                              int *imageVerts, int maxNumImageVerts);

extern int  _TestSubgraph(graphP theSubgraph, graphP theGraph);

/* Forward declarations of local functions */

void _K4Search_ClearStructures(K4SearchContext *context);
int  _K4Search_CreateStructures(K4SearchContext *context);
int  _K4Search_InitStructures(K4SearchContext *context);

/* Forward declarations of overloading functions */

int  _K4Search_CreateFwdArcLists(graphP theGraph);
void _K4Search_CreateDFSTreeEmbedding(graphP theGraph);
void _K4Search_EmbedBackEdgeToDescendant(graphP theGraph, int RootSide, int RootVertex, int W, int WPrevLink);
int  _K4Search_MarkDFSPath(graphP theGraph, int ancestor, int descendant);
int  _K4Search_HandleBlockedEmbedIteration(graphP theGraph, int I);
int  _K4Search_HandleBlockedDescendantBicomp(graphP theGraph, int I, int RootVertex, int R, int *pRout, int *pW, int *pWPrevLink);
int  _K4Search_EmbedPostprocess(graphP theGraph, int I, int edgeEmbeddingResult);
int  _K4Search_CheckEmbeddingIntegrity(graphP theGraph, graphP origGraph);
int  _K4Search_CheckObstructionIntegrity(graphP theGraph, graphP origGraph);

void _K4Search_InitGraphNode(graphP theGraph, int I);
void _InitK4SearchGraphNode(K4SearchContext *context, int I);
void _K4Search_InitVertexRec(graphP theGraph, int I);
void _InitK4SearchVertexRec(K4SearchContext *context, int I);

int  _K4Search_InitGraph(graphP theGraph, int N);
void _K4Search_ReinitializeGraph(graphP theGraph);
int  _K4Search_EnsureArcCapacity(graphP theGraph, int requiredArcCapacity);

/* Forward declarations of functions used by the extension system */

void *_K4Search_DupContext(void *pContext, void *theGraph);
void _K4Search_FreeContext(void *);

/****************************************************************************
 * K4SEARCH_ID - the variable used to hold the integer identifier for this
 * extension, enabling this feature's extension context to be distinguished
 * from other features' extension contexts that may be attached to a graph.
 ****************************************************************************/

int K4SEARCH_ID = 0;

/****************************************************************************
 gp_AttachK4Search()

 This function adjusts the graph data structure to attach the K4 search
 feature.
 ****************************************************************************/

int  gp_AttachK4Search(graphP theGraph)
{
     K4SearchContext *context = NULL;

     // If the K4 search feature has already been attached to the graph,
     // then there is no need to attach it again
     gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);
     if (context != NULL)
     {
         return OK;
     }

     // Allocate a new extension context
     context = (K4SearchContext *) malloc(sizeof(K4SearchContext));
     if (context == NULL)
     {
         return NOTOK;
     }

     // First, tell the context that it is not initialized
     context->initialized = 0;

     // Save a pointer to theGraph in the context
     context->theGraph = theGraph;

     // Put the overload functions into the context function table.
     // gp_AddExtension will overload the graph's functions with these, and
     // return the base function pointers in the context function table
     memset(&context->functions, 0, sizeof(graphFunctionTable));

     context->functions.fpCreateFwdArcLists = _K4Search_CreateFwdArcLists;
     context->functions.fpCreateDFSTreeEmbedding = _K4Search_CreateDFSTreeEmbedding;
     context->functions.fpEmbedBackEdgeToDescendant = _K4Search_EmbedBackEdgeToDescendant;
     context->functions.fpMarkDFSPath = _K4Search_MarkDFSPath;
     context->functions.fpHandleBlockedEmbedIteration = _K4Search_HandleBlockedEmbedIteration;
     context->functions.fpHandleBlockedDescendantBicomp = _K4Search_HandleBlockedDescendantBicomp;
     context->functions.fpEmbedPostprocess = _K4Search_EmbedPostprocess;
     context->functions.fpCheckEmbeddingIntegrity = _K4Search_CheckEmbeddingIntegrity;
     context->functions.fpCheckObstructionIntegrity = _K4Search_CheckObstructionIntegrity;

     context->functions.fpInitGraphNode = _K4Search_InitGraphNode;
     context->functions.fpInitVertexRec = _K4Search_InitVertexRec;

     context->functions.fpInitGraph = _K4Search_InitGraph;
     context->functions.fpReinitializeGraph = _K4Search_ReinitializeGraph;
     context->functions.fpEnsureArcCapacity = _K4Search_EnsureArcCapacity;

     _K4Search_ClearStructures(context);

     // Store the K33 search context, including the data structure and the
     // function pointers, as an extension of the graph
     if (gp_AddExtension(theGraph, &K4SEARCH_ID, (void *) context,
                         _K4Search_DupContext, _K4Search_FreeContext,
                         &context->functions) != OK)
     {
         _K4Search_FreeContext(context);
         return NOTOK;
     }

     // Create the K33-specific structures if the size of the graph is known
     // Attach functions are always invoked after gp_New(), but if a graph
     // extension must be attached before gp_Read(), then the attachment
     // also happens before gp_InitGraph(), which means N==0.
     // However, sometimes a feature is attached after gp_InitGraph(), in
     // which case N > 0
     if (theGraph->N > 0)
     {
         if (_K4Search_CreateStructures(context) != OK ||
             _K4Search_InitStructures(context) != OK)
         {
             _K4Search_FreeContext(context);
             return NOTOK;
         }
     }

     return OK;
}

/********************************************************************
 gp_DetachK4Search()
 ********************************************************************/

int gp_DetachK4Search(graphP theGraph)
{
    return gp_RemoveExtension(theGraph, K4SEARCH_ID);
}

/********************************************************************
 _K4Search_ClearStructures()
 ********************************************************************/

void _K4Search_ClearStructures(K4SearchContext *context)
{
    if (!context->initialized)
    {
        // Before initialization, the pointers are stray, not NULL
        // Once NULL or allocated, free() or LCFree() can do the job
        context->sortedDFSChildLists = NULL;
        context->G = NULL;
        context->V = NULL;

        context->initialized = 1;
    }
    else
    {
        LCFree(&context->sortedDFSChildLists);
        if (context->G != NULL)
        {
            free(context->G);
            context->G = NULL;
        }
        if (context->V != NULL)
        {
            free(context->V);
            context->V = NULL;
        }
    }
}

/********************************************************************
 _K4Search_CreateStructures()
 Create uninitialized structures for the vertex and graph node
 levels, and initialized structures for the graph level
 ********************************************************************/
int  _K4Search_CreateStructures(K4SearchContext *context)
{
     int N = context->theGraph->N;
     int Gsize = context->theGraph->edgeOffset + context->theGraph->arcCapacity;

     if (N <= 0)
         return NOTOK;

     if ((context->sortedDFSChildLists = LCNew(context->theGraph->N)) == NULL ||
         (context->G = (K4Search_GraphNodeP) malloc(Gsize*sizeof(K4Search_GraphNode))) == NULL ||
         (context->V = (K4Search_VertexRecP) malloc(N*sizeof(K4Search_VertexRec))) == NULL
        )
     {
         return NOTOK;
     }

     return OK;
}

/********************************************************************
 _K4Search_InitStructures()
 ********************************************************************/
int  _K4Search_InitStructures(K4SearchContext *context)
{
     int I, N = context->theGraph->N;
     int Gsize = context->theGraph->edgeOffset + context->theGraph->arcCapacity;

     if (N <= 0)
         return OK;

     for (I = 0; I < Gsize; I++)
          _InitK4SearchGraphNode(context, I);

     for (I = 0; I < N; I++)
          _InitK4SearchVertexRec(context, I);

     return OK;
}

/********************************************************************
 ********************************************************************/

int  _K4Search_InitGraph(graphP theGraph, int N)
{
    K4SearchContext *context = NULL;
    gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

    if (context == NULL)
        return NOTOK;
    {
        theGraph->N = N;
        theGraph->edgeOffset = 2*N;
        if (theGraph->arcCapacity == 0)
        	theGraph->arcCapacity = 2*DEFAULT_EDGE_LIMIT*N;

        if (_K4Search_CreateStructures(context) != OK)
            return NOTOK;

        // This call initializes the base graph structures, but it also
        // initializes the custom graphnode and vertex level structures
        // due to the overloads of InitGraphNode and InitVertexRec
        context->functions.fpInitGraph(theGraph, N);
    }

    return OK;
}

/********************************************************************
 ********************************************************************/

void _K4Search_ReinitializeGraph(graphP theGraph)
{
    K4SearchContext *context = NULL;
    gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

    if (context != NULL)
    {
        // Reinitialization can go much faster if the underlying
        // init graph node and vertex rec functions are called,
        // rather than the overloads of this module, because it
        // avoids lots of unnecessary gp_FindExtension() calls.
        if (theGraph->functions.fpInitGraphNode == _K4Search_InitGraphNode &&
            theGraph->functions.fpInitVertexRec == _K4Search_InitVertexRec)
        {
            // Restore the graph function pointers
            theGraph->functions.fpInitGraphNode = context->functions.fpInitGraphNode;
            theGraph->functions.fpInitVertexRec = context->functions.fpInitVertexRec;

            // Reinitialize the graph
            context->functions.fpReinitializeGraph(theGraph);

            // Restore the function pointers that attach this feature
            theGraph->functions.fpInitGraphNode = _K4Search_InitGraphNode;
            theGraph->functions.fpInitVertexRec = _K4Search_InitVertexRec;

            // Do the reinitialization that is specific to this module
            _K4Search_InitStructures(context);
            LCReset(context->sortedDFSChildLists);
        }

        // If optimization is not possible, then just stick with what works.
        // Reinitialize the graph-level structure and then invoke the
        // reinitialize function.
        else
        {
            LCReset(context->sortedDFSChildLists);

            // The underlying function fpReinitializeGraph() implicitly initializes the K33
            // structures due to the overloads of fpInitGraphNode() and fpInitVertexRec().
            // It just does so less efficiently because each invocation of InitGraphNode
            // and InitVertexRec has to look up the extension again.
            //// _K4Search_InitStructures(context);
            context->functions.fpReinitializeGraph(theGraph);
        }
    }
}

/********************************************************************
 The current implementation does not support an increase of arc
 (edge record) capacity once the extension is attached to the graph
 data structure.  This is only due to not being necessary to support.
 For now, it is easy to ensure the correct capacity before attaching
 the extension, but support could be added later if there is some
 reason to do so.
 ********************************************************************/

int  _K4Search_EnsureArcCapacity(graphP theGraph, int requiredArcCapacity)
{
	return NOTOK;
}

/********************************************************************
 _K4Search_DupContext()
 ********************************************************************/

void *_K4Search_DupContext(void *pContext, void *theGraph)
{
     K4SearchContext *context = (K4SearchContext *) pContext;
     K4SearchContext *newContext = (K4SearchContext *) malloc(sizeof(K4SearchContext));

     if (newContext != NULL)
     {
         int N = ((graphP) theGraph)->N;
         int Gsize = ((graphP) theGraph)->edgeOffset + ((graphP) theGraph)->arcCapacity;

         *newContext = *context;

         newContext->theGraph = (graphP) theGraph;

         newContext->initialized = 0;
         _K4Search_ClearStructures(newContext);
         if (N > 0)
         {
             if (_K4Search_CreateStructures(newContext) != OK)
             {
                 _K4Search_FreeContext(newContext);
                 return NULL;
             }

             LCCopy(newContext->sortedDFSChildLists, context->sortedDFSChildLists);
             memcpy(newContext->G, context->G, Gsize*sizeof(K4Search_GraphNode));
             memcpy(newContext->V, context->V, N*sizeof(K4Search_VertexRec));
         }
     }

     return newContext;
}

/********************************************************************
 _K4Search_FreeContext()
 ********************************************************************/

void _K4Search_FreeContext(void *pContext)
{
     K4SearchContext *context = (K4SearchContext *) pContext;

     _K4Search_ClearStructures(context);
     free(pContext);
}

/********************************************************************
 _K4Search_CreateFwdArcLists()

 Puts the forward arcs (back edges from a vertex to its descendants)
 into a circular list indicated by the fwdArcList member, a task
 simplified by the fact that they have already been placed in
 succession at the end of the adjacency list.

 For K4 search, the forward edges must be sorted by DFS number of
 the descendant endpoint.  The sort is linear time, but it is a little
 slower, so we avoid this cost for the other planarity-related algorithms.

  Returns OK on success, NOTOK on internal code failure
 ********************************************************************/

int _K4Search_CreateFwdArcLists(graphP theGraph)
{
    K4SearchContext *context = NULL;

    gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);
    if (context == NULL)
        return NOTOK;

    // For isolating a K_4 homeomorph, we need the forward edges
    // of each vertex to be in sorted order by DFI of descendants.
    // Otherwise we just drop through to the normal processing...

    if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
    {
    int I, Jcur, Jnext, ancestor;

        // for each vertex v in order, we follow each of its back edges
        // to the twin forward edge in an ancestor u, then we move
        // the forward edge to the fwdArcList of u.  Since this loop
        // processes vertices in order, the fwdArcList of each vertex
        // will be in order by the neighbor indicated by the forward edges.

        for (I=0; I < theGraph->N; I++)
        {
        	// Skip this vertex if it has no edges
        	Jnext = gp_GetLastArc(theGraph, I);
        	if (!gp_IsArc(theGraph, Jnext))
        		continue;

            // Skip the forward edges, which are in succession at the
        	// end of the arc list (last and its predecessors)
            while (theGraph->G[Jnext].type == EDGE_FORWARD)
                Jnext = gp_GetPrevArc(theGraph, Jnext);

            // Now we want to put all the back arcs in a backArcList, too.
            // Since we've already skipped past the forward arcs, we continue
            // with the predecessor arcs until we either run out of arcs or
            // we find a DFS child arc (the DFS child arcs are in succession
            // at the beginning of the arc list, so when a child arc is
            // encountered in the predecessor direction, then there won't be
            // any more back arcs.
            while (gp_IsArc(theGraph, Jnext) &&
                   theGraph->G[Jnext].type != EDGE_DFSCHILD)
            {
                Jcur = Jnext;
                Jnext = gp_GetPrevArc(theGraph, Jnext);

                if (theGraph->G[Jcur].type == EDGE_BACK)
                {
                    // Remove the back arc from I's adjacency list
                	gp_DetachArc(theGraph, Jcur);
                    gp_SetPrevArc(theGraph, Jcur, NIL);
                    gp_SetNextArc(theGraph, Jcur, NIL);

                    // Determine the ancestor of vertex I to which Jcur connects
                    ancestor = theGraph->G[Jcur].v;

                    // Go to the forward arc in the ancestor
                    Jcur = gp_GetTwinArc(theGraph, Jcur);

                    // Remove the forward arc from the ancestor's adjacency list
                	gp_DetachArc(theGraph, Jcur);

                    // Add the forward arc to the end of the fwdArcList.
                    if (theGraph->V[ancestor].fwdArcList == NIL)
                    {
                        theGraph->V[ancestor].fwdArcList = Jcur;
                        gp_SetPrevArc(theGraph, Jcur, Jcur);
                        gp_SetNextArc(theGraph, Jcur, Jcur);
                    }
                    else
                    {
                    	gp_AttachArc(theGraph, NIL, theGraph->V[ancestor].fwdArcList, 1, Jcur);
                    }
                }
            }
        }

        // Since the fwdArcLists have been created, we do not fall through
        // to run the superclass implementation
        return OK;
    }

    // If we're not actually running a K4 search, then we just run the
    // superclass implementation
    return context->functions.fpCreateFwdArcLists(theGraph);
}

/********************************************************************
 _K4Search_CreateDFSTreeEmbedding()

 This function overloads the basic planarity version in the manner
 explained below.  Once the extra behavior is performed, the basic
 planarity version is invoked.

 First, the sortedDFSChildList of each vertex is computed. Each vertex
 receives the list of its DFS children sorted by their DFS numbers.
 This is linear time overall. The core planarity/outerplanarity
 algorithm computes a different list of the DFS children, the
 separatedDFSChildList, in which the DFS children are stored in order
 of their lowpoint values.

 Second, the p2dFwdArcCount of each vertex is computed. This is the
 number of forward arcs from the DFS parent of the vertex to DFS
 descendants of the vertex. This is computed based on a simultaneous
 traversal through the sortedDFSChildList and the sorted fwdArcList.

 Third, the subtree value of each forward arc (V, D) is determined. This
 value indicates the DFS child C of V whose DFS subtree contains the DFS
 descendant endpoint D of the forward arc. This can be computed during
 the setting of the p2dFwdArcCount values.

 Each DFS child is listed in DFI order in V[I].sortedDFSChildList.
 In V[I].fwdArcList, the forward arcs of all back edges are in order
 by DFI of the descendant endpoint of the edge.

 DFS descendants have a higher DFI than ancestors, so given two
 successive children C1 and C2, if a forward arc leads to a
 vertex D such that DFI(C1) < DFI(D) < DFI(C2), then the
 forward arc contributes to the count of C1 and has C1 as subtree.
 ********************************************************************/

void _K4Search_CreateDFSTreeEmbedding(graphP theGraph)
{
    K4SearchContext *context = NULL;
    gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

    if (context != NULL)
    {
        if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
        {
            int I, J, C1, C2, D, e;
            int N = theGraph->N;

            // First compute the sortedDFSChildList of each vertex
            for (I=0; I<N; I++)
            {
                J = gp_GetFirstArc(theGraph, I);

                // If a vertex has any DFS children, the edges
                // to them are stored in descending order of
                // the DFI's along the successor arc pointers, so
                // we traverse them and prepend each to the
                // ascending order sortedDFSChildList
                while (theGraph->G[J].type == EDGE_DFSCHILD)
                {
                    context->V[I].sortedDFSChildList =
                        LCPrepend(context->sortedDFSChildLists,
                                    context->V[I].sortedDFSChildList,
                                    theGraph->G[J].v);

                    J = gp_GetNextArc(theGraph, J);
                }
            }

            // Next compute the p2dFwdArcCount of each vertex and the
            // subtree of each forward arc.
            for (I=0; I<N; I++)
            {
            	// For each DFS child of the vertex I, ...
                C1 = context->V[I].sortedDFSChildList;
                e = theGraph->V[I].fwdArcList;
                while (C1 != NIL && gp_IsArc(theGraph, e))
                {
                	// Get the next higher numbered child C2
                    C2 = LCGetNext(context->sortedDFSChildLists,
                                   context->V[I].sortedDFSChildList, C1);

                    // If there is a next child C2, then we can restrict attention
                    // to the forward arcs with DFI less than C2
                    if (C2 != NIL)
                    {
    					D = theGraph->G[e].v;
    					while (D < C2)
    					{
                    		context->V[C1].p2dFwdArcCount++;
                    		context->G[e].subtree = C1;

                    		// Go to the next forward arc
							e = gp_GetNextArc(theGraph, e);
							if (e == theGraph->V[I].fwdArcList)
							{
								e = NIL;
								break;
							}
							D = theGraph->G[e].v;
    					}
                    }

                    // If C1 is the last DFS child (C2==NIL), then all remaining
                    // forward edges must connect to descendants of C1.
                    else
                    {
                    	while (gp_IsArc(theGraph, e))
                    	{
                    		context->V[C1].p2dFwdArcCount++;
                    		context->G[e].subtree = C1;

                    		// Go to the next forward arc
							e = gp_GetNextArc(theGraph, e);
							if (e == theGraph->V[I].fwdArcList)
								e = NIL;
                    	}
                    }

					// Move the DFS child context to C2
					C1 = C2;
                }
            }
       }

        // Invoke the superclass version of the function
        // Each DFS tree child arc is moved to the root copy of the vertex
        context->functions.fpCreateDFSTreeEmbedding(theGraph);
    }
}

/********************************************************************
 _K4Search_EmbedBackEdgeToDescendant()

 The forward and back arcs of the cycle edge are embedded by the planarity
 version of this function.
 However, for K_4 subgraph homeomorphism, we also maintain a forward
 arc counter in a DFS child C of each vertex V to indicate how many
 forward arcs there are from V to descendants of C.  Each forward arc
 has an indicator, 'subtree', of C.  When we embed the edge, we decrement
 the counter so that when the WalkDown resolves as much pertinence as
 possible along the external face of the bicomp rooted by R=C+N, then
 we can easily determine whether there is more unresolved pertinence
 by testing whether the forward arc count has dropped to zero.
 If not, then we either find a K4 or perform a reduction that enables
 the WalkDown to make more progress when reinvoked.
 ********************************************************************/

void _K4Search_EmbedBackEdgeToDescendant(graphP theGraph, int RootSide, int RootVertex, int W, int WPrevLink)
{
    K4SearchContext *context = NULL;
    gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

    if (context != NULL)
    {
        // K4 search may have been attached, but not enabled
        if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
        {
        	int fwdArc = theGraph->V[W].adjacentTo;
        	context->V[context->G[fwdArc].subtree].p2dFwdArcCount--;
        }

        // Invoke the superclass version of the function
        context->functions.fpEmbedBackEdgeToDescendant(theGraph, RootSide, RootVertex, W, WPrevLink);
    }
}

/********************************************************************
 ********************************************************************/

void _K4Search_InitGraphNode(graphP theGraph, int I)
{
    K4SearchContext *context = NULL;
    gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

    if (context != NULL)
    {
        context->functions.fpInitGraphNode(theGraph, I);
        _InitK4SearchGraphNode(context, I);
    }
}

void _InitK4SearchGraphNode(K4SearchContext *context, int I)
{
    context->G[I].pathConnector = NIL;
    context->G[I].subtree = NIL;
}

/********************************************************************
 ********************************************************************/

void _K4Search_InitVertexRec(graphP theGraph, int I)
{
    K4SearchContext *context = NULL;
    gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

    if (context != NULL)
    {
        context->functions.fpInitVertexRec(theGraph, I);
        _InitK4SearchVertexRec(context, I);
    }
}

void _InitK4SearchVertexRec(K4SearchContext *context, int I)
{
    context->V[I].p2dFwdArcCount = 0;
    context->V[I].sortedDFSChildList = NIL;
}

/********************************************************************
// This function used to be implemented by going to each
// vertex, asking for its DFS parent, then finding the
// edge that lead to that DFS parent and marking it.
// However, the K4 search performs certain bicomp reductions
// that are required to preserve the essential structure of
// the DFS tree.  As a result, sometimes a DFSParent has been
// reduced out of the graph, but the tree edge leads to the nearest
// ancestor still in the graph.  So, when we want to leave a vertex,
// we search for the DFS tree edge to the "parent" (nearest ancestor),
// then we mark the edge and use it to go up to the "parent".
 ********************************************************************/

int  _K4Search_MarkDFSPath(graphP theGraph, int ancestor, int descendant)
{
int  J, parent, N;

     N = theGraph->N;

     /* If we are marking from a root vertex upward, then go up to the parent
        copy before starting the loop */

     if (descendant >= N)
         descendant = theGraph->V[descendant-N].DFSParent;

     /* Mark the lowest vertex (the one with the highest number). */

     theGraph->G[descendant].visited = 1;

     /* Mark all ancestors of the lowest vertex, and the edges used to reach
        them, up to the given ancestor vertex. */

     while (descendant != ancestor)
     {
          if (descendant == NIL)
              return NOTOK;

          /* If we are at a bicomp root, then ascend to its parent copy and
                mark it as visited. */

          if (descendant >= N)
          {
              parent = theGraph->V[descendant-N].DFSParent;
          }

          /* If we are on a regular, non-virtual vertex then get the edge to
                the parent, mark the edge, then fall through to the code
                that marks the parent vertex. */
          else
          {

              /* Scan the edges for the one marked as the DFS parent */

              parent = NIL;
              J = gp_GetFirstArc(theGraph, descendant);
              while (gp_IsArc(theGraph, J))
              {
                  if (theGraph->G[J].type == EDGE_DFSPARENT)
                  {
                      parent = theGraph->G[J].v;
                      break;
                  }
                  J = gp_GetNextArc(theGraph, J);
              }

              /* If the desired edge was not found, then the data structure is
                    corrupt, so bail out. */

              if (parent == NIL)
                  return NOTOK;

              /* Mark the edge */

              theGraph->G[J].visited = 1;
              theGraph->G[gp_GetTwinArc(theGraph, J)].visited = 1;
          }

          /* Mark the parent, then hop to the parent and reiterate */

          theGraph->G[parent].visited = 1;
          descendant = parent;
     }

     return OK;
}

/********************************************************************
 * This function is called if the outerplanarity algorithm fails to
 * embed all back edges for a vertex I.  This means that an obstruction
 * to outerplanarity has occurred, so we determine if it is a subgraph
 * homeomorphic to K4.  If so, then NONEMBEDDABLE is returned.  If not,
 * then a reduction is performed that unobstructs outerplanarity and
 * OK is returned, which allows the outerplanarity algorithm to
 * proceed with iteration I-1 (or to stop if I==0).
 ********************************************************************/

int  _K4Search_HandleBlockedEmbedIteration(graphP theGraph, int I)
{
    if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
    {
    	// If the fwdArcList is empty, then the K4 was already isolated
    	// by _K4Search_HandleBlockedDescendantBicomp(), and we just
    	// return the NONEMBEDDABLE result in order to stop the embedding
    	// iteration loop.
		if (theGraph->V[I].fwdArcList == NIL)
			return NONEMBEDDABLE;

        return _SearchForK4InBicomps(theGraph, I);
    }
    else
    {
        K4SearchContext *context = NULL;
        gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpHandleBlockedEmbedIteration(theGraph, I);
        }
    }

    return NOTOK;
}

/********************************************************************
 This function is called when outerplanarity obstruction minor A is
 encountered by the WalkDown.  In the implementation for the core
 planarity/outerplanarity algorithm, this method simply pushes the
 blocked bicomp root onto the stack and returns NONEMBEDDABLE, which
 causes the WalkDown to terminate.  The embed postprocessing would
 then go on to isolate the obstruction.

 However, outerplanarity obstruction minor A corresponds to a K_{2,3}
 homeomorph.  This method invokes a search for a K_4 homeomorph that
 may be entangled with the K_{2,3} homeomorph.  If an entangled K_4
 homeomorph is found, then _SearchForK4() returns NONEMBEDDABLE, which
 causes the WalkDown to terminate as above.  This is correct since a
 K_4 homeomorph has been found and isolated, and the K4Search overload
 of EmbedPostprocess() does no additional work.

 On the other hand, if minor A is found but there is no entangled K_4
 homeomorph, then the blocked descendant was reduced to a single edge
 so that it no longer obstructs outerplanarity. Then, OK was returned
 to indicate that the WalkDown should proceed.  This function then
 sets the vertex W and directional information that must be returned
 so that WalkDown can proceed.

 Returns OK to proceed with WalkDown at W,
         NONEMBEDDABLE to terminate WalkDown of Root Vertex
         NOTOK for internal error
 ********************************************************************/

int  _K4Search_HandleBlockedDescendantBicomp(graphP theGraph, int I, int RootVertex, int R, int *pRout, int *pW, int *pWPrevLink)
{
    if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
    {
    	int RetVal = _SearchForK4InBicomp(theGraph, NULL, I, R);

    	// On internal error (NOTOK) or K4 found (NONEMBEDDABLE), we return.
    	if (RetVal != OK)
    		return RetVal;

    	// Since the bicomp rooted by R is now a singleton edge, either direction
    	// out of R and into W can be selected, as long as they are consistent
    	// We just choose the settings associated with selecting W as the next
    	// vertex from R on the external face.
    	*pRout = 0;
    	*pWPrevLink = 1;
    	*pW = _GetNextVertexOnExternalFace(theGraph, R, pWPrevLink);

        // Now return OK so the Walkdown can continue at W (i.e. *pW)
        return OK;
    }
    else
    {
        K4SearchContext *context = NULL;
        gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpHandleBlockedDescendantBicomp(theGraph, I, RootVertex, R, pRout, pW, pWPrevLink);
        }
    }

    return NOTOK;
}

/********************************************************************
 ********************************************************************/

int  _K4Search_EmbedPostprocess(graphP theGraph, int I, int edgeEmbeddingResult)
{
     // For K4 search, we just return the edge embedding result because the
     // search result has been obtained already.
     if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
     {
         return edgeEmbeddingResult;
     }

     // When not searching for K4, we let the superclass do the work
     else
     {
        K4SearchContext *context = NULL;
        gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpEmbedPostprocess(theGraph, I, edgeEmbeddingResult);
        }
     }

     return NOTOK;
}

/********************************************************************
 ********************************************************************/

int  _K4Search_CheckEmbeddingIntegrity(graphP theGraph, graphP origGraph)
{
     if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
     {
         return OK;
     }

     // When not searching for K4, we let the superclass do the work
     else
     {
        K4SearchContext *context = NULL;
        gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpCheckEmbeddingIntegrity(theGraph, origGraph);
        }
     }

     return NOTOK;
}

/********************************************************************
 ********************************************************************/

int  _K4Search_CheckObstructionIntegrity(graphP theGraph, graphP origGraph)
{
     // When searching for K4, we ensure that theGraph is a subgraph of
     // the original graph and that it contains a K4 homeomorph
     if (theGraph->embedFlags == EMBEDFLAGS_SEARCHFORK4)
     {
		int  degrees[4], imageVerts[4];

        if (_TestSubgraph(theGraph, origGraph) != TRUE)
            return NOTOK;

		if (_getImageVertices(theGraph, degrees, 3, imageVerts, 4) != OK)
			return NOTOK;

		if (_TestForCompleteGraphObstruction(theGraph, 4, degrees, imageVerts) == TRUE)
		{
			return OK;
		}

		return NOTOK;
     }

     // When not searching for K4, we let the superclass do the work
     else
     {
        K4SearchContext *context = NULL;
        gp_FindExtension(theGraph, K4SEARCH_ID, (void *)&context);

        if (context != NULL)
        {
            return context->functions.fpCheckObstructionIntegrity(theGraph, origGraph);
        }
     }

     return NOTOK;
}
