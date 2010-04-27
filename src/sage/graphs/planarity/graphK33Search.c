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

#include "graphK33Search.h"
#include "graphK33Search.private.h"

extern int K33SEARCH_ID;

#include "graph.h"

/* Imported functions */

extern void _FillVisitedFlags(graphP, int);
extern int  _FillVisitedFlagsInBicomp(graphP theGraph, int BicompRoot, int FillValue);
extern int  _FillVisitedFlagsInOtherBicomps(graphP theGraph, int BicompRoot, int FillValue);
extern void _FillVisitedFlagsInUnembeddedEdges(graphP theGraph, int FillValue);
extern int  _GetBicompSize(graphP theGraph, int BicompRoot);
extern int  _HideInternalEdges(graphP theGraph, int vertex);
extern int  _RestoreInternalEdges(graphP theGraph, int stackBottom);
extern int  _DeleteUnmarkedEdgesInBicomp(graphP theGraph, int BicompRoot);
extern int  _ClearInvertedFlagsInBicomp(graphP theGraph, int BicompRoot);
extern int  _ComputeArcType(graphP theGraph, int a, int b, int edgeType);
extern int  _SetEdgeType(graphP theGraph, int u, int v);

extern int  _GetNextVertexOnExternalFace(graphP theGraph, int curVertex, int *pPrevLink);
extern int  _JoinBicomps(graphP theGraph);
extern int  _OrientVerticesInBicomp(graphP theGraph, int BicompRoot, int PreserveSigns);
extern int  _OrientVerticesInEmbedding(graphP theGraph);
extern void _InvertVertex(graphP theGraph, int V);
extern int  _SetVisitedOnPath(graphP theGraph, int u, int v, int w, int x, int visited);
extern int  _OrientExternalFacePath(graphP theGraph, int u, int v, int w, int x);

extern int  _ChooseTypeOfNonplanarityMinor(graphP theGraph, int I, int R);
extern int  _MarkHighestXYPath(graphP theGraph);

extern int  _IsolateKuratowskiSubgraph(graphP theGraph, int I, int R);

extern int  _FindUnembeddedEdgeToCurVertex(graphP theGraph, int cutVertex, int *pDescendant);
extern int  _FindUnembeddedEdgeToSubtree(graphP theGraph, int ancestor, int SubtreeRoot, int *pDescendant);

extern int  _MarkPathAlongBicompExtFace(graphP theGraph, int startVert, int endVert);

extern int  _AddAndMarkEdge(graphP theGraph, int ancestor, int descendant);

extern int  _DeleteUnmarkedVerticesAndEdges(graphP theGraph);

extern int  _IsolateMinorE1(graphP theGraph);
extern int  _IsolateMinorE2(graphP theGraph);
extern int  _IsolateMinorE3(graphP theGraph);
extern int  _IsolateMinorE4(graphP theGraph);

extern int  _GetLeastAncestorConnection(graphP theGraph, int cutVertex);
extern int  _MarkDFSPathsToDescendants(graphP theGraph);
extern int  _AddAndMarkUnembeddedEdges(graphP theGraph);

/* Private functions for K_{3,3} searching. */

int  _SearchForK33(graphP theGraph, int I);

int  _SearchForK33InBicomp(graphP theGraph, K33SearchContext *context, int I, int R);

int  _RunExtraK33Tests(graphP theGraph, K33SearchContext *context);
int  _SearchForMinorE1(graphP theGraph);
int  _FinishIsolatorContextInitialization(graphP theGraph, K33SearchContext *context);
int  _SearchForDescendantExternalConnection(graphP theGraph, K33SearchContext *context, int cutVertex, int u_max);
int  _GetAdjacentAncestorInRange(graphP theGraph, K33SearchContext *context, int vertex,
                                int closerAncestor, int fartherAncestor);
int  _FindExternalConnectionDescendantEndpoint(graphP theGraph, int ancestor,
                                               int cutVertex, int *pDescendant);
int  _SearchForMergeBlocker(graphP theGraph, K33SearchContext *context, int I, int *pMergeBlocker);
int  _FindK33WithMergeBlocker(graphP theGraph, K33SearchContext *context, int I, int mergeBlocker);

int  _TestForLowXYPath(graphP theGraph);
int  _TestForZtoWPath(graphP theGraph);
int  _TestForStraddlingBridge(graphP theGraph, K33SearchContext *context, int u_max);
int  _ReduceBicomp(graphP theGraph, K33SearchContext *context, int R);
int  _ReduceExternalFacePathToEdge(graphP theGraph, K33SearchContext *context, int u, int x, int edgeType);
int  _ReduceXYPathToEdge(graphP theGraph, K33SearchContext *context, int u, int x, int edgeType);
int  _RestoreReducedPath(graphP theGraph, K33SearchContext *context, int J);
int  _RestoreAndOrientReducedPaths(graphP theGraph, K33SearchContext *context);

int  _IsolateMinorE5(graphP theGraph);
int  _IsolateMinorE6(graphP theGraph, K33SearchContext *context);
int  _IsolateMinorE7(graphP theGraph, K33SearchContext *context);

/****************************************************************************
 _SearchForK33()

  The strategy involves one special case in which, to achieve a linear time
  bound, we must delay the discovery of a K_{3,3} that caused a Walkdown
  failure prior to step I.  In such cases, vertex I was an ancestor with
  a connection to the bicomp on which the Walkdown failed, but it would
  have been too costly to find I at the time.  So, the bicomp was marked
  as non-mergeable prior to some ancestor of I.  If this function is
  invoked for step I, then we have found the connection from that bicomp
  prior to reaching the limiting ancestor of I. The bicomp and I are
  therefore part of a K_{3,3} that can now be isolated.

  Otherwise, a Walkdown failure in step I with a non-empty merge stack
  would have already resulted in an identified K_{3,3} by minor A, so
  we must have an empty merge stack now.

  We must first find all bicomp roots on which the Walkdown has failed
  in step I.  The fwdArcList of I contains the forward arcs of the
  back edges for I that we failed to embed.  Each forward arc leads to
  a descendant of I that is in a DFS subtree rooted by a child of I,
  where the child of I has the greatest DFI that is less than the DFI
  of the descendant indicated by the forward arc.  Each bicomp root
  that represents a vertex is uniquely associated with a DFS child
  of the vertex, so once we know the child of I whose subtree contains
  a descendant of I that the Walkdown couldn't reach, we can immediately
  deduce the root copy of I on which the Walkdown failed.

  For each such root copy of I, we test whether a K_{3,3} homemorph
  can be isolated based on that bicomp.  If so, then we return it.
  Otherwise, each bicomp can be reduced to a 4-cycle and the edges
  that the Walkdown failed to embed can be discarded.
 ****************************************************************************/

int  _SearchForK33(graphP theGraph, int I)
{
int  C1, C2, D, e, RetVal=OK, FoundOne;
K33SearchContext *context = NULL;

    gp_FindExtension(theGraph, K33SEARCH_ID, (void *)&context);
    if (context == NULL)
        return NOTOK;

/* Before we begin with the standard array of K_{3,3} tests, we handle
    one optimization case that may be left over from a prior step
    of the embedding algorithm.  If the embedding stack is non-empty,
    then the Walkdown either halted due to non-planarity minor A or
    because of the merge blocking optimization (see CASE 3 in the
    function RunExtraK33Tests()).  We test for the latter condition,
    and if it is found, then we isolate a K_{3,3} and return. */

     if (!sp_IsEmpty(theGraph->theStack))
     {
     int mergeBlocker;

         if (_SearchForMergeBlocker(theGraph, context, I, &mergeBlocker) != OK)
             return NOTOK;

         if (mergeBlocker != NIL)
         {
             if (_FindK33WithMergeBlocker(theGraph, context, I, mergeBlocker) != OK)
                 return NOTOK;

             return NONEMBEDDABLE;
         }
     }

     /* Each DFS child is listed in DFI order in V[I].sortedDFSChildList.
        In V[I].fwdArcList, the forward arcs of all unembedded back edges are
        in order by DFI of the descendant endpoint of the edge.

        DFS descendants have a higher DFI than ancestors, so given two
        successive children C1 and C2, if any forward arcs lead to a
        vertex D such that DFI(C1) < DFI(D) < DFI(C2), then the Walkdown
        failed to embed a back edge from I to a descendant D of C1. */

     e = theGraph->V[I].fwdArcList;
     D = theGraph->G[e].v;

     C1 = context->V[I].sortedDFSChildList;

     FoundOne = FALSE;

     while (C1 != NIL && e != NIL)
     {
        C2 = LCGetNext(context->sortedDFSChildLists,
                       context->V[I].sortedDFSChildList, C1);

        // If the edge e leads from I to a descendant D of C1,
        // then D will be less than C2 (as explained above),
        // so we search for a K_{3,3} in the bicomp rooted
        // by the root copy of I associated with C1.
        // (If C2 is NIL, then C1 is the last child)

        if (D < C2 || C2 == NIL)
        {
        	FoundOne = TRUE;
            RetVal = _SearchForK33InBicomp(theGraph, context, I, C1+theGraph->N);

            // If something went wrong, NOTOK was returned;
            // If a K_{3,3} was found, NONEMBEDDABLE was returned;
            // If OK was returned, then only a K5 was found, so
            // we continue searching any other bicomps on which
            // the Walkdown failed.

            if (RetVal != OK)
             break;
        }

        // Skip the edges that lead to descendants of C1 to get to those
        // edges that lead to descendants of C2.

        if (C2 != NIL)
        {
            while (D < C2 && gp_IsArc(theGraph, e))
            {
                e = gp_GetNextArc(theGraph, e);
                if (e == theGraph->V[I].fwdArcList)
                     e = NIL;
                else D = theGraph->G[e].v;
            }
        }

        // Move the DFS child context to C2
        C1 = C2;
     }

/* If we got through the loop with an OK value for each bicomp on
     which the Walkdown failed, then we return OK to indicate that only
     K5's were found (or there is a special case K_{3,3} that may be discovered
     later based on a setting we made in the data structure).
     The OK return allows the embedder to continue.

     If a K_{3,3} is ever found (or if an error occured), then RetVal
     will not be OK, and the loop terminates immediately so we can
     return the appropriate value.  If a K_{3,3} is found, then we must
     also handle the fact that some paths of the input graph may have
     been reduced to single edges by prior _ReduceBicomp() calls.

     NOTE: The variable FoundOne helps ensure we detect at least one
        bicomp on which the Walkdown failed (this should always be
        the case in an error-free implementation like this one!). */

     return FoundOne ? RetVal : NOTOK;
}

/****************************************************************************
 _SearchForK33InBicomp()
 ****************************************************************************/

int  _SearchForK33InBicomp(graphP theGraph, K33SearchContext *context, int I, int R)
{
isolatorContextP IC = &theGraph->IC;
int tempResult;

/* Begin by determining which non-planarity minor is detected */

     if (_ChooseTypeOfNonplanarityMinor(theGraph, I, R) != OK)
         return NOTOK;

/* If minor A is selected, then the root of the oriented bicomp has been changed */

     else R = IC->r;

/* Minors A to D result in the desired K_{3,3} homeomorph,
    so we isolate it and return NONEMBEDDABLE. */

     if (theGraph->IC.minorType & (MINORTYPE_A|MINORTYPE_B|MINORTYPE_C|MINORTYPE_D))
     {
        /* First we restore the orientations of the vertices in the
            one bicomp we have messed with so that there is no confusion. */

        if (_OrientVerticesInBicomp(theGraph, R, 1) != OK)
        	return NOTOK;

        /* Next we restore the orientation of the embedding so we
           can restore the reduced paths (because we avoid modifying
           the Kuratowski subgraph isolator to restore reduced paths,
           which are a construct of the K_{3,3} search). */

        if (_OrientVerticesInEmbedding(theGraph) != OK ||
        	_RestoreAndOrientReducedPaths(theGraph, context) != OK)
            return NOTOK;

        /* Next we simply call the Kuratowski subgraph isolation since
            we know now that it will isolate a K_{3,3}.
            For minor A, we need to set up the stack that would be
            available immediately after a Walkdown failure. */

        if (theGraph->IC.minorType & MINORTYPE_A)
        {
            sp_ClearStack(theGraph->theStack);
            sp_Push2(theGraph->theStack, R, NIL);
        }

        if (_IsolateKuratowskiSubgraph(theGraph, I, R) != OK)
            return NOTOK;

        return NONEMBEDDABLE;
     }

/* For minor E (a K5 minor), we run the additional tests to see if
    minors E1 to E4 apply since these minors isolate a K_{3,3} entangled
    with the K5. */

     IC->ux = _GetLeastAncestorConnection(theGraph, IC->x);
     IC->uy = _GetLeastAncestorConnection(theGraph, IC->y);
     IC->uz = _GetLeastAncestorConnection(theGraph, IC->z);

     if (IC->z != IC->w ||
         IC->uz > MAX(IC->ux, IC->uy) ||
         (IC->uz < MAX(IC->ux, IC->uy) && IC->ux != IC->uy) ||
         (IC->x != IC->px || IC->y != IC->py))
     {
        if (_OrientVerticesInBicomp(theGraph, R, 1) != OK)
        	return NOTOK;

        if (_OrientVerticesInEmbedding(theGraph) != OK ||
        	_RestoreAndOrientReducedPaths(theGraph, context) != OK)
            return NOTOK;

        if (_IsolateKuratowskiSubgraph(theGraph, I, R) != OK)
            return NOTOK;

        return NONEMBEDDABLE;
     }

/* If the Kuratowski subgraph isolator will not isolate a K_{3,3} based on minor E,
    then a K5 homeomorph could be isolated.  However, a K_{3,3} may still be tangled
    with the K5, so we now run the additional tests of the K_{3,3} search algorithm.

    If the search finds a K_{3,3} (tempResult of NONEMBEDDABLE), then we remove unwanted
    edges from the graph and return NONEMBEDDABLE.  If the search has a fault (NOTOK),
    then we return.  If the result is OK, then a K_{3,3} was not found at this time
    and we proceed with some clean-up work below. */

     if ((tempResult = _RunExtraK33Tests(theGraph, context)) != OK)
     {
         if (tempResult == NONEMBEDDABLE)
            if (_DeleteUnmarkedVerticesAndEdges(theGraph) != OK)
                return NOTOK;

         return tempResult;
     }

/* The extra cases for finding a K_{3,3} did not succeed, so the bicomp rooted by R
    is either a K5 homeomorph (with at most a superficially entangled K_{3,3}) or
    we have made the special setting that allows us to detect the one case that
    would be too costly to try now.  Either way, we can safely reduce the bicomp
    to the 4-cycle (R, X, W, Y, R) and proceed with the planarity algorithm.
    We also restore the mixed orientation of the bicomp (i.e. the proper
    orientation in the context of the edge signs) because this code can work
    when ReduceBicomp doesn't do any actual work. */

     if (_OrientVerticesInBicomp(theGraph, R, 1) != OK)
    	 return NOTOK;

     if (_ReduceBicomp(theGraph, context, R) != OK)
         return NOTOK;

/* Set visited flags to a high number so planarity algorithm
    can properly do Walkup procedure in future steps */

     if (_FillVisitedFlagsInBicomp(theGraph, IC->r, theGraph->N) != OK)
    	 return NOTOK;

/* We now intend to ignore the pertinence of W (conceptually eliminating
    the connection from W to the current vertex).  Note that none of the
    bicomps listed in the pertinentBicompList (nor their respective subtrees)
    will be visited again by the planarity algorithm because they must've
    been internally active.  If they were externally active and pertinent,
    then we would've found a K_{3,3} by non-planarity minor B. Thus, the original
    Walkup costs that identified the pertinent bicomps we intend to ignore are
    one-time costs, preserving linear time. */

     theGraph->V[IC->w].adjacentTo = NIL;
     theGraph->V[IC->w].pertinentBicompList = NIL;

     return OK;
}

/****************************************************************************
 _RunExtraK33Tests()
 ****************************************************************************/

int  _RunExtraK33Tests(graphP theGraph, K33SearchContext *context)
{
isolatorContextP IC = &theGraph->IC;
int u_max = MAX3(IC->ux, IC->uy, IC->uz), u;

/* Case 1: If there is a pertinent or externally active vertex other than W
            on the lower external face path between X and Y (the points of
            attachment of the x-y path), then we can isolate a K_{3,3} homeomorph
            by Minor E1. */

    if (_SearchForMinorE1(theGraph) != OK)
        return NOTOK;

    if (IC->w != IC->z)
    {
        if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
            _IsolateMinorE1(theGraph) != OK)
            return NOTOK;

        return NONEMBEDDABLE;
    }

/* Case 2: If W/Z can make an external connection to an ancestor of V
            that is descendant to u_{max}, then a K_{3,3} homeomorph can
            be isolated with Minor E2.

            NOTE: It may seem costly to check the entire subtree, but
            if it succeeds then we're done, and if the only connection
            is to u_{max} then we are sure to never make this check
            again on this subtree (if all the other K_{3,3} tests fail).

            OPTIMIZATION: We do not check for the connection if the
            least ancestor connection from W/Z leads to an ancestor
            of u_max.  The reason is that it could ultimately be too
            costly if the connection doesn't exist, and if the highest
            numbered ancestor H of the current vertex with an external
            connection from W is a descendant u_{max} (and if no other
            tests in this function succeed), then we will discover a
            K_{3,3} by Minor A or B in step H.

            OPTIMIZATION: When we do test for an external connection to
            an ancestor of V descendant to u_{max}, the result may be that
            only a connection to u_{max} exists.  The result is cached
            so that the subtrees of the vertex need not be traversed again
            should we need to make the test again.
            See _SearchForDescendantExternalConnection() */

    if (IC->uz == u_max)
    {
        u = _SearchForDescendantExternalConnection(theGraph, context, IC->w, u_max);
        if (u > u_max)
        {
            IC->uz = u;
            if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
                _IsolateMinorE2(theGraph) != OK)
                return NOTOK;

            return NONEMBEDDABLE;
        }
    }

/* Case 3: If X or Y can make an external connection to an ancestor of V
            that is descendant to u_{max}, then a K_{3,3} homeomorph
            can be isolated with Minor E3.

            NOTE: It may seem costly to check the entire subtree, but
            if it succeeds then we're done, and if the only connection
            is to u_{max} then we are sure to never make this check
            again on this subtree (if all the other K_{3,3} tests fail).

            OPTIMIZATION:  Due to the prior use of the Kuratowski subgraph
            isolator, we know that at most one of X, Y or W/Z could have an
            external connection to an ancestor of u_{max} = MAX(ux, uy, uz).
            If it is X or Y, then we do not check for the lower connection
            required to find Minor E3 because it might ultimately be too
            costly.  Instead, we mark the vertex with a 'merge barrier'
            of u_{max}.  If the planar embedder attempts to merge the vertex
            prior to step u_{max}, then the embedder has found the desired
            connection and a K_{3,3} is isolated at that time.

            OPTIMIZATION: When we can test for an external connection to
            an ancestor of V descendant to u_{max}, the result may be that
            only a connection to u_{max} exists.  The result is cached
            so that the subtrees of the vertex need not be traversed again
            should we need to make the test again.
            See _SearchForDescendantExternalConnection() */

    if (IC->ux < u_max)
        context->V[IC->x].mergeBlocker = u_max;
    else
    {
        u = _SearchForDescendantExternalConnection(theGraph, context, IC->x, u_max);
        if (u > u_max)
        {
            IC->ux = u;
            if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
                _IsolateMinorE3(theGraph) != OK)
                return NOTOK;

            return NONEMBEDDABLE;
        }
    }

    if (IC->uy < u_max)
        context->V[IC->y].mergeBlocker = u_max;
    else
    {
        u = _SearchForDescendantExternalConnection(theGraph, context, IC->y, u_max);
        if (u > u_max)
        {
            IC->uy = u;
            if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
                _IsolateMinorE3(theGraph) != OK)
                return NOTOK;

            return NONEMBEDDABLE;
        }
    }

/* Case 4: If there exists any x-y path with points of attachment px and py
            such that px!=x or py!=y, then a K_{3,3} homeomorph can be isolated
            with Minor E4. */

    if (_TestForLowXYPath(theGraph) != OK)
        return NOTOK;

    if (IC->px != IC->x || IC->py != IC->y)
    {
        if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
            _IsolateMinorE4(theGraph) != OK)
            return NOTOK;

        return NONEMBEDDABLE;
    }

/* Case 5: If the x-y path contains an internal vertex that starts a second
            internal path from the internal vertex to W/Z, then a K_{3,3} homeomorph
            can be isolated with Minor E5. */

    if (_TestForZtoWPath(theGraph) != OK)
        return NOTOK;

    if (theGraph->G[IC->w].visited)
    {
        if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
            _IsolateMinorE5(theGraph) != OK)
            return NOTOK;

        return NONEMBEDDABLE;
    }

/* Case 6: If uz < u_{max} and there is an external connection (other than external
            connections involving X, Y and W/Z) between an ancestor of u_{max} and a
            vertex in the range [V...u_{max}), then a K_{3,3} homeomorph can be
            isolated with Minor E6.

            OPTIMIZATION:  See _TestForStraddlingBridge() */

    if (IC->uz < u_max)
    {
        if (_TestForStraddlingBridge(theGraph, context, u_max) != NIL)
        {
            if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
                _IsolateMinorE6(theGraph, context) != OK)
                return NOTOK;

            return NONEMBEDDABLE;
        }
    }

/* Case 7: If ux < u_{max} or uy < u_{max} and there is an external connection
            between an ancestor of u_{max} and a vertex in the range [V...u_{max})
            (except for external connections involving X, Y and W/Z), then a K_{3,3}
            homeomorph can be isolated with Minor E7.

            OPTIMIZATION: Same as Case 6.*/

    if (IC->ux < u_max || IC->uy < u_max)
    {
        if (_TestForStraddlingBridge(theGraph, context, u_max) != NIL)
        {
            if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
                _IsolateMinorE7(theGraph, context) != OK)
                return NOTOK;

            return NONEMBEDDABLE;
        }
    }

/* If none of the tests found a K_{3,3}, then we return OK to indicate that nothing
    went wrong, but a K_{3,3} was not found. */

    return OK;
}

/****************************************************************************
 _SearchForMinorE1()
 Search along the external face below the x-y path for a vertex Z other
 than W that is externally active or pertinent.
 ****************************************************************************/

int _SearchForMinorE1(graphP theGraph)
{
int  Z=theGraph->IC.px, ZPrevLink=1;

     Z = _GetNextVertexOnExternalFace(theGraph, Z, &ZPrevLink);

     while (Z != theGraph->IC.py)
     {
         if (Z != theGraph->IC.w)
         {
            if (_VertexActiveStatus(theGraph, Z, theGraph->IC.v) == VAS_EXTERNAL)
            {
                theGraph->IC.z = Z;
                theGraph->IC.uz = _GetLeastAncestorConnection(theGraph, Z);
                return OK;
            }
            else if (theGraph->V[Z].pertinentBicompList != NIL ||
                     theGraph->V[Z].adjacentTo == theGraph->IC.v)
            {
                /* Swap the roles of W and Z */

                theGraph->IC.z = theGraph->IC.w;
                theGraph->IC.w = Z;

                /* If the new W (indicated by Z) was on the path (R, X, old W) then
                    the new Z (the old W, which has no type mark) is on the path
                    (X, new W, new Z, Y) so we change the type new Z to being on the
                    RYW path. Otherwise, the order is (X, new Z, new W, Y), so the
                    new Z (old W with no type) is type changed to be on the RXW path.*/

                if (theGraph->G[Z].type == VERTEX_LOW_RXW)
                     theGraph->G[theGraph->IC.z].type = VERTEX_LOW_RYW;
                else theGraph->G[theGraph->IC.z].type = VERTEX_LOW_RXW;

                /* For completeness, we change the new W to type unknown */

                theGraph->G[theGraph->IC.w].type = TYPE_UNKNOWN;

                /* The external activity ancestor connection of the new Z must be obtained */

                theGraph->IC.uz = _GetLeastAncestorConnection(theGraph, theGraph->IC.z);

                return OK;
            }
         }

         Z = _GetNextVertexOnExternalFace(theGraph, Z, &ZPrevLink);
     }

     return OK;
}

/****************************************************************************
 _FinishIsolatorContextInitialization()
 Once it has been decided that a desired subgraph can be isolated, it
 becomes safe to finish the isolator context initialization.
 ****************************************************************************/

int  _FinishIsolatorContextInitialization(graphP theGraph, K33SearchContext *context)
{
isolatorContextP IC = &theGraph->IC;

/* Restore the orientation of the bicomp on which we're working, then
    perform orientation of all vertices in graph. (An unnecessary but
    polite step that simplifies the description of key states of the
    data structures). */

     if (_OrientVerticesInBicomp(theGraph, IC->r, 1) != OK)
    	 return NOTOK;

     if (_OrientVerticesInEmbedding(theGraph) != OK)
    	 return NOTOK;

/* Restore any paths that were reduced to single edges */

     if (_RestoreAndOrientReducedPaths(theGraph, context) != OK)
         return NOTOK;

/* We assume that the current bicomp has been marked appropriately,
     but we must now clear the visitation flags of all other bicomps. */

     if (_FillVisitedFlagsInOtherBicomps(theGraph, IC->r, 0) != OK)
    	 return NOTOK;

/* To complete the normal behavior of _FillVisitedFlags() in the
    normal isolator context initialization, we also have to clear
    the visited flags on all edges that have not yet been embedded */

     _FillVisitedFlagsInUnembeddedEdges(theGraph, 0);

/* Now we can find the descendant ends of unembedded back edges based on
     the ancestor settings ux, uy and uz. */

     if (_FindExternalConnectionDescendantEndpoint(theGraph, IC->ux, IC->x, &IC->dx) != OK ||
         _FindExternalConnectionDescendantEndpoint(theGraph, IC->uy, IC->y, &IC->dy) != OK ||
         _FindExternalConnectionDescendantEndpoint(theGraph, IC->uz, IC->z, &IC->dz) != OK)
         return NOTOK;

/* Finally, we obtain the descendant end of an unembedded back edge to
     the current vertex. */

     if (_FindUnembeddedEdgeToCurVertex(theGraph, IC->w, &IC->dw) != TRUE)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _GetAdjacentAncestorInRange()
 Returns the ancestor of theVertex that is adjacent to theVertex by an
 unembedded back edge and has a DFI strictly between closerAncestor and
 fartherAncestor.
 Returns NIL if theVertex has no such neighboring ancestor.
 ****************************************************************************/

int _GetAdjacentAncestorInRange(graphP theGraph, K33SearchContext *context, int theVertex,
                                int closerAncestor, int fartherAncestor)
{
int J = context->V[theVertex].backArcList;

    while (gp_IsArc(theGraph, J))
    {
        if (theGraph->G[J].v < closerAncestor &&
            theGraph->G[J].v > fartherAncestor)
            return theGraph->G[J].v;

        J = gp_GetNextArc(theGraph, J);
        if (J == context->V[theVertex].backArcList)
            J = NIL;
    }
    return NIL;
}

/****************************************************************************
 _SearchForDescendantExternalConnection()
 Search the cutVertex and each subtree rooted by a vertex in the
 separatedDFSChildList of the cut vertex for an external connection
 to a vertex ancestor to the current vertex V and descendant to u_max.

 The function returns the descendant of u_max found to have an external
 connection to the given cut vertex.

 OPTIMIZATION: The subtrees are processed by preorder visitation.  If
 a vertex is visited and has a lowpoint indicating that it and its
 descendants make no external connections, then we prune the subtree,
 eliminating the vertex and its descendants from consideration.
 Otherwise, if the vertex has an externalConnectionAncestor setting,
 which must have been made by a prior invocation of this function,
 then we use that setting.  If both of these tests fail, then
 we examine the vertex and push its children for consideration.
 This ensures that this procedure never explores a subtree more than
 once during the life of the K_{3,3} search algorithm.

 Note that if a subtree is explored and the desired external connection
 descendant to u_{max} is found, then a K_{3,3} will be found, so we only
 have to concern ourselves with subtrees that connect only to u_{max}.
 Between steps v and u_{max}, the subtree search is optimized by setting
 'externalConnectionAncestor', and steps after u_{max} process ancestors
 of u_{max}.  Since this routine is only called if the cutVertex makes
 no external connections to ancestors of u_{max}, the lowpoint test
 prevents this subtree from being searched after step u_{max}.
 ****************************************************************************/

int  _SearchForDescendantExternalConnection(graphP theGraph, K33SearchContext *context, int cutVertex, int u_max)
{
isolatorContextP IC = &theGraph->IC;
int  u2 = _GetAdjacentAncestorInRange(theGraph, context, cutVertex, IC->v, u_max);
int  listHead, child, descendant;

/* Test cut vertex for external connection to descendant of u_max */

     if (u2 != NIL)
         return u2;

/* If there is no direct back edge connection from the cut vertex
        to a vertex on the path between V and u_max, then we will
        look for such a connection in the DFS subtrees rooted by
        separated DFS children of the vertex (ignoring those whose
        lowpoint indicates that they make no external connections) */

     /* Begin by pushing the separated DFS children of the cut vertex,
        except stop when the lowpoint is no longer less than V since
        external connections are not being made. */

     sp_ClearStack(theGraph->theStack);
     listHead = child = theGraph->V[cutVertex].separatedDFSChildList;
     while (child != NIL)
     {
         if (theGraph->V[child].Lowpoint >= IC->v)
             break;
         sp_Push(theGraph->theStack, child);
         child = LCGetNext(theGraph->DFSChildLists, listHead, child);
     }

     /* Now process the stack until it is empty or until we've found the desired connection. */

     while (!sp_IsEmpty(theGraph->theStack))
     {
         sp_Pop(theGraph->theStack, descendant);

         /* If the vertex has a lowpoint indicating that it makes no external
            connections, then skip the subtree rooted by the vertex */

         if (theGraph->V[descendant].Lowpoint < IC->v)
         {
             /* If a prior invocation has precalculated the result, use it. */

             if (context->V[descendant].externalConnectionAncestor != NIL)
             {
                 /* If the result is in the range we need, return it.  Otherwise,
                    skip the subtree rooted by the vertex. */

                 if (context->V[descendant].externalConnectionAncestor < IC->v &&
                     context->V[descendant].externalConnectionAncestor > u_max)
                     return context->V[descendant].externalConnectionAncestor;
             }

             /* If the subtree has not been explored, then explore it. */

             else
             {
                 /* Check the subtree root for the desired connection. */

                 u2 = _GetAdjacentAncestorInRange(theGraph, context, descendant, IC->v, u_max);
                 if (u2 != NIL)
                     return u2;

                 /* Push each child as a new subtree root to be considered,
                    except skip those whose lowpoint is too great. */

                 child = context->V[descendant].sortedDFSChildList;
                 while (child != NIL)
                 {
                     if (theGraph->V[child].Lowpoint < IC->v)
                         sp_Push(theGraph->theStack, child);

                     child = LCGetNext(context->sortedDFSChildLists,
                                       context->V[descendant].sortedDFSChildList, child);
                 }
             }
         }
     }

/* The only external connections from the cutVertex lead to u_max,
    so cache the result and return it. */

     context->V[cutVertex].externalConnectionAncestor = u_max;
     return u_max;
}

/****************************************************************************
 _FindExternalConnectionDescendantEndpoint()

 This operation is similar to _FindUnembeddedEdgeToAncestor() except that
 we need to be more precise in this case, finding an external connection
 from a given cut vertex to a *particular* given ancestor.

 NOTE: By external we don't mean externall active so much as not embedded in
       the bicomp containing the cut vertex.

 Returns OK if it finds that either the given cutVertex or one of its
    descendants in a separated bicomp has an unembedded back edge
    connection to the given ancestor vertex.
 Returns NOTOK otherwise (it is an error to not find the descendant because
    this function is only called if _SearchForDescendantExternalConnection()
    has already determined the existence of the descendant).
 ****************************************************************************/

int  _FindExternalConnectionDescendantEndpoint(graphP theGraph, int ancestor,
                                               int cutVertex, int *pDescendant)
{
int  listHead, child, J;

     // Check whether the cutVertex is directly adjacent to the ancestor
     // by an unembedded back edge.

     J = theGraph->V[ancestor].fwdArcList;
     while (gp_IsArc(theGraph, J))
     {
         if (theGraph->G[J].v == cutVertex)
         {
             *pDescendant = cutVertex;
             return OK;
         }

         J = gp_GetNextArc(theGraph, J);
         if (J == theGraph->V[ancestor].fwdArcList)
             J = NIL;
     }

     // Now check the descendants of the cut vertex to see if any make
     // a connection to the ancestor.
     listHead = child = theGraph->V[cutVertex].separatedDFSChildList;
     while (child != NIL)
     {
         if (theGraph->V[child].Lowpoint >= theGraph->IC.v)
             break;

         if (_FindUnembeddedEdgeToSubtree(theGraph, ancestor, child, pDescendant) == TRUE)
             return OK;

         child = LCGetNext(theGraph->DFSChildLists, listHead, child);
     }

     return NOTOK;
}

/****************************************************************************
 _SearchForMergeBlocker()

 This function helps to implement the merge blocking optimization of
 _SearchForDescendantExternalConnection().  The function RunExtraK33Tests()
 sets a mergeBlocker rather than run _SearchForDescendantExternalConnection()
 in certain cases.  This procedure is called by MergeBicomps to test the
 embedding stack for a merge blocker before merging any biconnected components.
 If a merge blocker is found, then the embedder's Walkdown function is
 terminated and SearchForK33() is subsequently called.  The blocked merge
 point is then used as the basis for isolating a K_{3,3}.

 Returns OK on success (whether or not the search found a merge blocker)
         NOTOK on internal function failure
         pMergeBlocker is set to NIL unless a merge blocker is found.
 ****************************************************************************/

int  _SearchForMergeBlocker(graphP theGraph, K33SearchContext *context, int I, int *pMergeBlocker)
{
stackP tempStack;
int  R, Rout, Z, ZPrevLink;

/* Set return result to 'not found' then return if there is no stack to inspect */

     *pMergeBlocker = NIL;

     if (sp_IsEmpty(theGraph->theStack))
         return OK;

/* Create a copy of the embedding stack */

     tempStack = sp_Duplicate(theGraph->theStack);
     if (tempStack == NULL)
         return NOTOK;

/* Search the copy of the embedding stack for a merge blocked vertex */

     while (!sp_IsEmpty(tempStack))
     {
         sp_Pop2(tempStack, R, Rout);
         sp_Pop2(tempStack, Z, ZPrevLink);

         if (context->V[Z].mergeBlocker != NIL &&
             context->V[Z].mergeBlocker < I)
         {
             *pMergeBlocker = Z;
             break;
         }
     }

     sp_Free(&tempStack);
     return OK;
}

/****************************************************************************
 _FindK33WithMergeBlocker()

 This function completes the merge blocking optimization by isolating a K_{3,3}
 based on minor E3 if a merge blocked vertex was previously found.

 Returns OK on success, NOTOK on internal function failure
 ****************************************************************************/

int  _FindK33WithMergeBlocker(graphP theGraph, K33SearchContext *context, int I, int mergeBlocker)
{
int  R, RPrevLink, u_max, u, J, W;
isolatorContextP IC = &theGraph->IC;

/* First, we orient the vertices so we can successfully restore all of the
    reduced paths.  This needs to be done before reconstructing the context
    for CASE 3 of RunExtraK33Tests() because the reconstruction involves
    using the Walkup to I from a descendant of I, which will not work if
    the descendant is in one of the reduced paths. */

     if (_OrientVerticesInEmbedding(theGraph) != OK ||
    	 _RestoreAndOrientReducedPaths(theGraph, context) != OK)
         return NOTOK;

/* Reconstruct the context that was present for CASE 3 of RunExtraK33Tests()
        when we decided to set a mergeBlocker rather than calling
        _SearchForDescendantExternalConnection() */

     /* Obtain the root of the bicomp containing the mergeBlocker. */

     RPrevLink = 1;
     R = mergeBlocker;
     while (R < theGraph->N)
        R = _GetNextVertexOnExternalFace(theGraph, R, &RPrevLink);

     /* Switch the 'current step' variable I to be equal to the
       non-virtual counterpart of the bicomp root. */

     I = theGraph->V[R - theGraph->N].DFSParent;

     /* Eliminate the visitation and pertinence settings for step u_max */

     _FillVisitedFlags(theGraph, I+1);

     for (J = 0; J < theGraph->N; J++)
     {
         theGraph->V[J].adjacentTo = NIL;
         theGraph->V[J].pertinentBicompList = NIL;
     }

     /* Restore the pertinence settings of step I by doing the Walkup for each
        back edge that was not embedded when step I was originally performed. */

     J = theGraph->V[I].fwdArcList;
     while (gp_IsArc(theGraph, J))
     {
        W = theGraph->G[J].v;
        theGraph->functions.fpWalkUp(theGraph, I, W);

        J = gp_GetNextArc(theGraph, J);
        if (J == theGraph->V[I].fwdArcList)
            J = NIL;
     }

/* Next, we make the standard initialization calls for when we have found
     a non-planarity condition. */

     sp_ClearStack(theGraph->theStack);

     if (_ChooseTypeOfNonplanarityMinor(theGraph, I, R) != OK)
         return NOTOK;

     IC->ux = _GetLeastAncestorConnection(theGraph, IC->x);
     IC->uy = _GetLeastAncestorConnection(theGraph, IC->y);
     IC->uz = _GetLeastAncestorConnection(theGraph, IC->z);

     u_max = MAX3(IC->ux,IC->uy,IC->uz);

/* Perform the remainder of CASE 3 of RunExtraK33Tests() */

     if (mergeBlocker == IC->x)
     {
         u = _SearchForDescendantExternalConnection(theGraph, context, IC->x, u_max);
         if (u > u_max)
         {
             IC->ux = u;
             if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
                 _IsolateMinorE3(theGraph) != OK)
                 return NOTOK;
         }
         else return NOTOK;
     }
     else if (mergeBlocker == IC->y)
     {
         u = _SearchForDescendantExternalConnection(theGraph, context, IC->y, u_max);
         if (u > u_max)
         {
             IC->uy = u;
             if (_FinishIsolatorContextInitialization(theGraph, context) != OK ||
                 _IsolateMinorE3(theGraph) != OK)
                 return NOTOK;
         }
         else return NOTOK;
     }
     else return NOTOK;

/* Do the final clean-up to obtain the K_{3,3} */

     if (_DeleteUnmarkedVerticesAndEdges(theGraph) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _TestForLowXYPath()
 Is there an x-y path that does not include X?
 If not, is there an x-y path that does not include Y?
 If not, then we restore the original x-y path.
 If such a low x-y path exists, then we adjust px or py accordingly,
    and we make sure that X or Y (whichever is excluded) and its edges are
    not marked visited.
 This method uses the stack, though it is called with an empty stack currently,
 it does happen to preserve any preceding stack content. This method pushes
 at most one integer per edge incident to the bicomp root plus two integers
 per vertex in the bicomp.
 ****************************************************************************/

int  _TestForLowXYPath(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;
int  result;
int  stackBottom;

/* Clear the previously marked X-Y path */

     if (_FillVisitedFlagsInBicomp(theGraph, IC->r, 0) != OK)
    	 return NOTOK;

/* Save the size of the stack before hiding any edges, so we will know
   how many edges to restore */

     stackBottom = sp_GetCurrentSize(theGraph->theStack);

/* Hide the internal edges of X */

     if (_HideInternalEdges(theGraph, IC->x) != OK)
    	 return NOTOK;

/* Try to find a low X-Y path that excludes X, then restore the
    internal edges of X. */

     result = _MarkHighestXYPath(theGraph);
     if (_RestoreInternalEdges(theGraph, stackBottom) != OK)
    	 return NOTOK;

/* If we found the low X-Y path, then return. */

     if (result == TRUE)
         return OK;

/* Hide the internal edges of Y */

     if (_HideInternalEdges(theGraph, IC->y) != OK)
    	 return NOTOK;

/* Try to find a low X-Y path that excludes Y, then restore the
    internal edges of Y. */

     result = _MarkHighestXYPath(theGraph);
     if (_RestoreInternalEdges(theGraph, stackBottom) != OK)
    	 return NOTOK;

/* If we found the low X-Y path, then return. */

     if (result == TRUE)
         return OK;

/* Restore the original X-Y path and return with no error
        (the search failure is reflected by no change to px and py */

     if (_MarkHighestXYPath(theGraph) != TRUE)
    	 return NOTOK;

     return OK;
}

/****************************************************************************
 _TestForZtoWPath()
 This function tests whether there is a path inside the bicomp leading from W
 to some internal node of the x-y path.  If there is, the path is marked.

 Upon function return, the marking of W distinguishes whether the path was found.
 The function returns NOTOK on internal error, OK otherwise.

 All internal vertices are marked as type unknown, as are W and the bicomp
 root.  There is an X-Y path marked visited.  So, we start a depth first
 search from W to find a visited vertex, except we prune the search to
 ignore vertices whose type is not unknown.

 The depth first search has to mark the vertices it has seen as visited,
 but we do not want to conflict with the visited/non-visited settings
 that have so far been used to isolate the X-Y path.  So, each vertex
 visited is marked with a NIL and pushed onto the resetList.  At the end,
 all vertices on the resetList have their visited flags reset to 0.

 For each vertex we visit, if it is an internal vertex on the X-Y path
 (i.e. visited=1 and type unknown), then we want to stop and unroll the
 stack to obtain the desired path (described below). If the vertex is type
 unknown, then we want to visit its unvisited neighbors.

 We want to manage the stack so that it when the desired vertex is found,
 the stack contains the desired path.  So, we do not simply push the
 neighbors of the vertex being visited.  First, we only push 'eligible'
 vertices (i.e. vertices with a type of unknown and visited not equal to
 NIL).  Second, when we decide a vertex v is eligible, we push (v, NIL).
 When we pop (v, NIL), we know that its type is unknown so we test
 whether it is the desired vertex by checking if its visited member is
 equal to 1.  If so, then we can stop the depth first search, process
 the resetList, then use the vertices and edges remaining on the
 stack to mark the desired path.

 If we pop (v, NIL) and find that the visited of v equals 0, then we
 set its visited to NIL.  Then we find the first edge record e leading
 to an eligible vertex w (i.e. a vertex with type unknown and visited
 not equal to NIL), and we push both (v, e) and (w, NIL).  Eventually all
 paths leading from w will be explored, and if none find the desired vertex,
 then (v, e) is popped.  Now we search the adjacency list of v starting
 after e to find the edge record that indicates the next eligible vertex
 to visit.  If none are found, then we simply go to the next iteration,
 which pops a 2-tuple containing the vertex u and an edge record e that
 indicates v as the neighbor of u.  Finally, if the stack empties without
 finding the desired vertex, then we simply process the resetStack and return.
 ****************************************************************************/

int  _TestForZtoWPath(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;
stackP resetList = sp_New(_GetBicompSize(theGraph, IC->r));
int  v, e, w;

     if (resetList == NULL) return NOTOK;

     sp_ClearStack(theGraph->theStack);
     sp_Push2(theGraph->theStack, IC->w, NIL);

     while (!sp_IsEmpty(theGraph->theStack))
     {
          sp_Pop2(theGraph->theStack, v, e);

          if (e == NIL)
          {
              if (theGraph->G[v].visited)
                  break;

              theGraph->G[v].visited = NIL;
              sp_Push(resetList, v);

              e = gp_GetFirstArc(theGraph, v);
          }
          else
              e = gp_GetNextArc(theGraph, e);

          while (gp_IsArc(theGraph, e))
          {
              w = theGraph->G[e].v;
              if (theGraph->G[w].visited != NIL &&
                  theGraph->G[w].type == TYPE_UNKNOWN)
              {
                  sp_Push2(theGraph->theStack, v, e);
                  sp_Push2(theGraph->theStack, w, NIL);
                  break;
              }
              e = gp_GetNextArc(theGraph, e);
          }
     }

     while (!sp_IsEmpty(resetList))
     {
         sp_Pop(resetList, v);
         theGraph->G[v].visited = 0;
     }
     sp_Free(&resetList);

     while (!sp_IsEmpty(theGraph->theStack))
     {
         sp_Pop2(theGraph->theStack, v, e);
         theGraph->G[v].visited = 1;
         theGraph->G[e].visited = 1;
         theGraph->G[gp_GetTwinArc(theGraph, e)].visited = 1;
     }

     return OK;
}

/****************************************************************************
 _TestForStraddlingBridge()
 We proceed on the path [V...u_{max}) from the current vertex V up to and
 excluding u_{max}.  For each vertex p, we test whether p has a least
 ancestor less than u_{max} and whether p has a DFS child c that is not an
 ancestor of X, Y and W and that has a connection to an ancestor of u_{max}
 (in other words, whether the child C has a lowpoint less than u_{max}).

 The separatedDFSChildList of each vertex already contains a list of
 the DFS children sorted by their lowpoint, and the list has not been
 reduced by bicomp merging because the vertices are not descendants of V.
 So, we can process a vertex by examining its leastAncestor and the
 lowpoint of one of the first two elements in its separatedDFSChildList.
 If the first child is an ancestor of X, Y and W, then we look at the
 second child.

 If no bridge straddling u_{max} is found, the function returns NIL.
 If a straddling bridge is found, the function returns a descendant d
 of p in the subtree rooted by c such that d has a leastAncestor less
 than u_{max}.  Given the vertex d, the path through the straddling
 bridge required in Minors E6 and E7 is easy to identify:  Mark the
 DFS tree path from d to p, and add and mark the edge from d to its
 least ancestor.

 OPTIMIZATION: If a straddling bridge is not found, then in each tree edge of
        the path [V...u_{max}) we set the member noStraddle equal to u_{max}.
        Then, we modify the above stated routine so that if it is testing
        for a straddling bridge of u_{max} along this path, it will stop
        if it encounters an edge with noStraddle equal to u_{max} then it
        will stop.  Also, the optimization will only set noStraddle equal to
        u_{max} on the portion of the path that is traversed.  Finally, if
        noStraddle is set to a value other than NIL, the setting will be
        ignored and it will not be changed.

        Due to this optimization, we do not traverse a path more than once
        to find out whether a vertex on the path has a bridge that straddles
        u_{max}.  This leaves two questions:
            1) What if a future step must determine whether there is a
                straddling bridge of an ancestor of u_{max}?
            2) What if a future step must determine whether there is a
                straddling bridge of a descendant of u_{max}?

        The condition described in the first question cannot occur because it
        would imply the ability to detect a straddling bridge now.
        The condition described by the second question may occur, but in the
        future step, the bicomp now being tested for a K_{3,3} will be part of
        a straddling bridge in that future step.  Thus, the straddling
        bridge query is asked at most twice along any DFS tree path.
 ****************************************************************************/

int  _TestForStraddlingBridge(graphP theGraph, K33SearchContext *context, int u_max)
{
isolatorContextP IC = &theGraph->IC;
int  p, c, d, excludedChild, e;

     p = IC->v;
     excludedChild = IC->r - theGraph->N;
     d = NIL;

/* Starting at V, traverse the ancestor path to u_max looking for a straddling bridge */

     while (p > u_max)
     {
         /* If we find a direct edge from p to an ancestor of u_max, the break. */

         if (theGraph->V[p].leastAncestor < u_max)
         {
             d = p;
             break;
         }

         /* Check for a path from p to an ancestor of u_max using the child
            of p with the least Lowpoint, except the child that is an
            ancestor of X, Y and W. */

         c = theGraph->V[p].separatedDFSChildList;
         if (c == excludedChild)
             c = LCGetNext(theGraph->DFSChildLists, c, c);

         if (c != NIL && theGraph->V[c].Lowpoint < u_max)
         {
             _FindUnembeddedEdgeToSubtree(theGraph, theGraph->V[c].Lowpoint, c, &d);
             break;
         }

         /* Check for noStraddle of u_max, break if found */

         e = gp_GetFirstArc(theGraph, p);
         if (context->G[e].noStraddle == u_max)
             break;

         /* Go to the next ancestor */

         excludedChild = p;
         p = theGraph->V[p].DFSParent;
     }

/* If d is NIL, then no straddling bridge was found, so we do the
        noStraddle optimization. */

     if (d == NIL)
     {
         c = IC->v;
         while (c != p)
         {
             e = gp_GetFirstArc(theGraph, c);
             if (context->G[e].noStraddle != NIL)
                 break;

             context->G[e].noStraddle = u_max;

             c = theGraph->V[c].DFSParent;
         }
     }

/* Return either NIL indicating no bridge straddling u_max or the descendant d
         used to help mark a straddling bridge that was found by this test. */

     return d;
}

/****************************************************************************
 _ReduceBicomp()

 We want to reduce the given biconnected component to a 4-cycle plus an
 internal edge connecting X and Y.  Each edge is to be associated with a
 path from the original graph, preserving the depth first search tree
 paths that help connect the vertices R, X, Y, and W.  If a K_{3,3} is later found,
 the paths are restored, but it is necessary to preserve the DFS tree so that
 functions like MarkDFSPath() will be able to pass through the restored bicomp.
 Also, if a K_{3,3} is later found due to the merge blocker optimization, then the
 internal X-Y path may be needed and, once the bicomp reduction is reversed,
 a full DFS subtree connecting all vertices in the bicomp will need to be
 restored or else functions that traverse the bicomp will not work.

 For example, _FindK33WithMergeBlocker() invokes ChooseTypeOfNonplanarityMinor()
 to help reconstruct the context under which the mergeBlocker was set.
 ChooseTypeOfNonplanarityMinor() calls _FillVisitedFlagsInBicomp(), which
 depends on the DFS tree.

 NOTE: The following are some general steps taken in this method:
       1) All edges in the bicomp are marked unvisited
       2) selected paths are marked visited
       3) unvisited edges are deleted
       4) the edges of the bicomp are marked unvisited again
       5) the remaining paths of the bicomp are reduced
       Some of the edges that get deleted in step 3 above may represent
       paths that were reduced in prior embedder iterations.  We delete
       the reduction edge but not the path it represents.
       If a K_{3,3} is ever found, then the edges of these reduced paths
       are still in the graph, though not connected to anything important.
       The desired K_{3,3} is marked visited, but step 4 above ensures that
       these reduction paths are not marked visited. Hence, they will be
       deleted when the K_{3,3} is isolated, and this routine does not
       need to restore any reduced paths on the edges it deletes.
       We also don't (and don't have the time to) restore any reduction
       edges along the paths we intend to keep.
 ****************************************************************************/

int  _ReduceBicomp(graphP theGraph, K33SearchContext *context, int R)
{
isolatorContextP IC = &theGraph->IC;
int  min, mid, max, A, A_edge, B, B_edge;
int  rxType, xwType, wyType, yrType, xyType;

/* The vertices in the bicomp need to be oriented so that functions
    like MarkPathAlongBicompExtFace() will work. */

     if (_OrientVerticesInBicomp(theGraph, R, 0) != OK)
    	 return NOTOK;

/* The reduced edges start with a default type of 'tree' edge. The
     tests below, which identify the additional non-tree paths
     needed to complete the reduced bicomp, also identify which
     reduced edges need to be cycle edges.*/

     rxType = xwType = wyType = yrType = xyType = EDGE_DFSPARENT;

/* Now we calculate some values that help figure out the shape of the
    DFS subtree whose structure will be retained in the bicomp. */

     min = MIN3(IC->x, IC->y, IC->w);
     max = MAX3(IC->x, IC->y, IC->w);
     mid = MAX3(MIN(IC->x, IC->y), MIN(IC->x, IC->w), MIN(IC->y, IC->w));

/* If the order of descendendancy from V goes first to X, then it can
    proceed either to W then Y or to Y then W */

     if (min == IC->x)
     {
         /* A is a descendant adjacent to the current vertex by a cycle edge
            whose DFS tree path to either mid or max is combined with the
            cycle edge to form the path that will be reduced to the
            external face cycle edge (V, max). */

         A_edge = gp_GetLastArc(theGraph, IC->r);
         A = theGraph->G[A_edge].v;
         yrType = EDGE_BACK;

         /* If Y is max, then a path parallel to the X-Y path will be a
            second path reduced to a cycle edge.  We find the neighbor B
            of min=X on the X-Y path.  The edge (B, min) is a cycle edge
            that, along with the DFS tree path (B, ..., max), will be
            retained and reduced to a cycle edge. */

         if (max == IC->y)
         {
             B_edge = gp_GetLastArc(theGraph, IC->x);
             while (B_edge != gp_GetFirstArc(theGraph, IC->x))
             {
                 if (theGraph->G[B_edge].visited) break;
                 B_edge = gp_GetPrevArc(theGraph, B_edge);
             }

             if (!theGraph->G[B_edge].visited)
                 return NOTOK;

             B = theGraph->G[B_edge].v;
             xyType = EDGE_BACK;
         }

         /* Otherwise, W is max so we find the neighbor B of min=X on the
            lower external face path (X, ..., W), which excludes V.  The
            cycle edge (B, min) and the DFS tree path (B, max) will be
            retained and reduced to a cycle edge.*/

         else if (max == IC->w)
         {
             B_edge = gp_GetFirstArc(theGraph, IC->x);
             B = theGraph->G[B_edge].v;
             xwType = EDGE_BACK;
         }

         else return NOTOK;
     }

/* Otherwise, the order of descendancy from V goes first to Y, then it
     proceeds to either W then X or to X then W. The */

     else
     {
         A_edge = gp_GetFirstArc(theGraph, IC->r);
         A = theGraph->G[A_edge].v;
         rxType = EDGE_BACK;

         if (max == IC->x)
         {
             B_edge = gp_GetFirstArc(theGraph, IC->y);
             while (B_edge != gp_GetLastArc(theGraph, IC->y))
             {
                 if (theGraph->G[B_edge].visited) break;
                 B_edge = gp_GetNextArc(theGraph, B_edge);
             }

             if (!theGraph->G[B_edge].visited)
                 return NOTOK;

             B = theGraph->G[B_edge].v;
             xyType = EDGE_BACK;
         }

         else if (max == IC->w)
         {
             B_edge = gp_GetLastArc(theGraph, IC->y);
             B = theGraph->G[B_edge].v;
             wyType = EDGE_BACK;
         }

         else return NOTOK;
     }

/* Now that we have collected the information on which cycle edge and
    which tree paths will actually be retained, we clear the visited
    flags so the current X-Y path will not be retained (an X-Y path
    formed mostly or entirely from DFS tree edges is retained). */

     if (_FillVisitedFlagsInBicomp(theGraph, R, 0) != OK)
    	 return NOTOK;

/* Now we mark the tree path from the maximum numbered vertex up
      to the bicomp root. This marks one of the following four paths:
      Case 1. (V, ..., X=min, ..., W=mid, ..., Y=max)
      Case 2. (V, ..., X=min, ..., Y=mid, ..., W=max)
      Case 3. (V, ..., Y=min, ..., W=mid, ..., X=max)
      Case 4. (V, ..., Y=min, ..., X=mid, ..., W=max) */

     if (theGraph->functions.fpMarkDFSPath(theGraph, R, max) != OK)
         return NOTOK;

/* Now we use A to mark a path on the external face corresponding to:
      Case 1. (V, ..., Y=max)
      Case 2. (V, ..., Y=mid)
      Case 3. (V, ..., X=max)
      Case 4. (V, ..., X=mid) */

     if (theGraph->functions.fpMarkDFSPath(theGraph, min==IC->x ? IC->y : IC->x, A) != OK)
         return NOTOK;

     theGraph->G[A_edge].visited = 1;
     theGraph->G[gp_GetTwinArc(theGraph, A_edge)].visited = 1;

/* Now we use B to mark either an X-Y path or a path of the external face
      corresponding to:
      Case 1. (X=min, ..., B, ..., Y=max)
      Case 2. (X=min, ..., B, ..., W=max)
      Case 3. (Y=min, ..., B, ..., X=max)
      Case 4. (Y=min, ..., B, ..., W=max) */

     if (theGraph->functions.fpMarkDFSPath(theGraph, max, B) != OK)
         return NOTOK;

     theGraph->G[B_edge].visited = 1;
     theGraph->G[gp_GetTwinArc(theGraph, B_edge)].visited = 1;

/* Delete the unmarked edges in the bicomp. Note that if an unmarked edge
 * represents a reduced path, then only the reduction edge is deleted here.
 * The path it represents is only deleted later (see NOTE above) */

     if (_DeleteUnmarkedEdgesInBicomp(theGraph, R) != OK)
    	 return NOTOK;

/* Clear all visited flags in the bicomp.
     This is the important "step 4" mentioned in the NOTE above */

     if (_FillVisitedFlagsInBicomp(theGraph, R, 0) != OK)
    	 return NOTOK;

/* Clear all orientation signs in the bicomp.
	Note that the whole bicomp may not be properly oriented at this point
	because we may have exchanged external face paths for internal
	DFS tree paths.  However, the reduced bicomp will be properly
	oriented, and the paths of degree 2 vertices will have their
	orientations fixed if/when reduction edges are restored. */

     if (_ClearInvertedFlagsInBicomp(theGraph, R) != OK)
    	 return NOTOK;

/* Reduce the paths to single edges.
	 Note that although the whole bicomp may not be properly oriented at this
	 point (as noted above), the four principal vertices R, X, W and Y still
	 are consistently oriented with one another, e.g. R's link[0] indicates
	 the external face path toward X that excludes W and Y, and X's link[1]
	 indicates that same path. */

     if (_ReduceExternalFacePathToEdge(theGraph, context, R, IC->x, rxType) != OK ||
         _ReduceExternalFacePathToEdge(theGraph, context, IC->x, IC->w, xwType) != OK ||
         _ReduceExternalFacePathToEdge(theGraph, context, IC->w, IC->y, wyType) != OK ||
         _ReduceExternalFacePathToEdge(theGraph, context, IC->y, R, yrType) != OK)
         return NOTOK;

     if (_ReduceXYPathToEdge(theGraph, context, IC->x, IC->y, xyType) != OK)
         return NOTOK;

/* The core planarity method used vertex visited flags in the Walkup, so we have to
   set the vertex visited flags so the remaining vertices will behave as though they
   are unvisited by Walkup when the embedder moves to the next vertex. */

     theGraph->G[R].visited =
     theGraph->G[IC->x].visited =
     theGraph->G[IC->y].visited =
     theGraph->G[IC->w].visited = IC->v;

     return OK;
}

/****************************************************************************
 _ReduceExternalFacePathToEdge()
 ****************************************************************************/

int  _ReduceExternalFacePathToEdge(graphP theGraph, K33SearchContext *context, int u, int x, int edgeType)
{
int  prevLink, v, w, e;

     /* If the path is a single edge, then no need for a reduction */

     prevLink = 1;
     v = _GetNextVertexOnExternalFace(theGraph, u, &prevLink);
     if (v == x)
     {
         theGraph->extFace[u].vertex[0] = x;
         theGraph->extFace[x].vertex[1] = u;
         return OK;
     }

     /* We have the endpoints u and x of the path, and we just computed the
        first vertex internal to the path and a neighbor of u.  Now we
        compute the vertex internal to the path and a neighbor of x. */

     prevLink = 0;
     w = _GetNextVertexOnExternalFace(theGraph, x, &prevLink);

     /* Delete the two edges that connect the path to the bicomp.
        If either edge is a reduction edge, then we have to restore
        the path it represents. We can only afford to visit the
        endpoints of the path.
        Note that in the restored path, the edge incident to each
        endpoint of the original path is a newly added edge,
        not a reduction edge. */

     e = gp_GetFirstArc(theGraph, u);
     if (context->G[e].pathConnector != NIL)
     {
         if (_RestoreReducedPath(theGraph, context, e) != OK)
             return NOTOK;
         e = gp_GetFirstArc(theGraph, u);
         v = theGraph->G[e].v;
     }
     gp_DeleteEdge(theGraph, e, 0);

     e = gp_GetLastArc(theGraph, x);
     if (context->G[e].pathConnector != NIL)
     {
         if (_RestoreReducedPath(theGraph, context, e) != OK)
             return NOTOK;
         e = gp_GetLastArc(theGraph, x);
         w = theGraph->G[e].v;
     }
     gp_DeleteEdge(theGraph, e, 0);

     /* Add the reduction edge, then set its path connectors so the original
        path can be recovered and set the edge type so the essential structure
        of the DFS tree can be maintained (The 'Do X to Bicomp' functions
        and functions like MarkDFSPath(0 depend on this). */

     gp_AddEdge(theGraph, u, 0, x, 1);

     e = gp_GetFirstArc(theGraph, u);
     context->G[e].pathConnector = v;
     theGraph->G[e].type = _ComputeArcType(theGraph, u, x, edgeType);

     e = gp_GetLastArc(theGraph, x);
     context->G[e].pathConnector = w;
     theGraph->G[e].type = _ComputeArcType(theGraph, x, u, edgeType);

     /* Set the external face info */

     theGraph->extFace[u].vertex[0] = x;
     theGraph->extFace[x].vertex[1] = u;

     return OK;
}

/****************************************************************************
 _ReduceXYPathToEdge()
 ****************************************************************************/

int  _ReduceXYPathToEdge(graphP theGraph, K33SearchContext *context, int u, int x, int edgeType)
{
int  e, v, w;

     e = gp_GetFirstArc(theGraph, u);
     e = gp_GetNextArc(theGraph, e);
     v = theGraph->G[e].v;

     /* If the XY-path is a single edge, then no reduction is needed */

     if (v == x)
         return OK;

     /* Otherwise, remove the two edges that join the XY-path to the bicomp */

     if (context->G[e].pathConnector != NIL)
     {
         if (_RestoreReducedPath(theGraph, context, e) != OK)
             return NOTOK;
         e = gp_GetFirstArc(theGraph, u);
         e = gp_GetNextArc(theGraph, e);
         v = theGraph->G[e].v;
     }
     gp_DeleteEdge(theGraph, e, 0);

     e = gp_GetFirstArc(theGraph, x);
     e = gp_GetNextArc(theGraph, e);
     w = theGraph->G[e].v;
     if (context->G[e].pathConnector != NIL)
     {
         if (_RestoreReducedPath(theGraph, context, e) != OK)
             return NOTOK;
         e = gp_GetFirstArc(theGraph, x);
         e = gp_GetNextArc(theGraph, e);
         w = theGraph->G[e].v;
     }
     gp_DeleteEdge(theGraph, e, 0);

     /* Now add a single edge to represent the XY-path */
     gp_InsertEdge(theGraph, u, gp_GetFirstArc(theGraph, u), 0,
    		                 x, gp_GetFirstArc(theGraph, x), 0);

     /* Now set up the path connectors so the original XY-path can be recovered if needed.
        Also, set the reduction edge's type to preserve the DFS tree structure */

     e = gp_GetFirstArc(theGraph, u);
     e = gp_GetNextArc(theGraph, e);
     context->G[e].pathConnector = v;
     theGraph->G[e].type = _ComputeArcType(theGraph, u, x, edgeType);

     e = gp_GetFirstArc(theGraph, x);
     e = gp_GetNextArc(theGraph, e);
     context->G[e].pathConnector = w;
     theGraph->G[e].type = _ComputeArcType(theGraph, x, u, edgeType);

     return OK;
}

/****************************************************************************
 _RestoreReducedPath()
 Given an edge record of an edge used to reduce a path, we want to restore
 the path in constant time.
 The path may contain more reduction edges internally, but we do not
 search for and process those since it would violate the constant time
 bound required of this function.
 return OK on success, NOTOK on failure
 ****************************************************************************/

int  _RestoreReducedPath(graphP theGraph, K33SearchContext *context, int J)
{
int  JTwin, u, v, w, x;
int  J0, J1, JTwin0, JTwin1;

     if (context->G[J].pathConnector == NIL)
         return OK;

     JTwin = gp_GetTwinArc(theGraph, J);

     u = theGraph->G[JTwin].v;
     v = context->G[J].pathConnector;
     w = context->G[JTwin].pathConnector;
     x = theGraph->G[J].v;

     /* Get the locations of the graph nodes between which the new
        graph nodes must be added in order to reconnect the path
        parallel to the edge. */

     J0 = gp_GetNextArc(theGraph, J);
     J1 = gp_GetPrevArc(theGraph, J);
     JTwin0 = gp_GetNextArc(theGraph, JTwin);
     JTwin1 = gp_GetPrevArc(theGraph, JTwin);

     /* We first delete the edge represented by J and JTwin. We do so before
        restoring the path to ensure we do not exceed the maximum arc capacity. */

     gp_DeleteEdge(theGraph, J, 0);

     /* Now we add the two edges to reconnect the reduced path represented
        by the edge [J, JTwin].  The edge record in u is added between J0 and J1.
        Likewise, the new edge record in x is added between JTwin0 and JTwin1. */

     if (gp_IsArc(theGraph, J0))
     {
    	 if (gp_InsertEdge(theGraph, u, J0, 1, v, gp_AdjacencyListEndMark(v), 0) != OK)
    		 return NOTOK;
     }
     else
     {
    	 if (gp_InsertEdge(theGraph, u, J1, 0, v, gp_AdjacencyListEndMark(v), 0) != OK)
    		 return NOTOK;
     }

     if (gp_IsArc(theGraph, JTwin0))
     {
    	 if (gp_InsertEdge(theGraph, x, JTwin0, 1, w, gp_AdjacencyListEndMark(w), 0) != OK)
    		 return NOTOK;
     }
     else
     {
    	 if (gp_InsertEdge(theGraph, x, JTwin1, 0, w, gp_AdjacencyListEndMark(w), 0) != OK)
    		 return NOTOK;
     }

     // Set the types of the newly added edges. In both cases, the first of the two
     // vertex parameters is known to be degree 2 because they are internal to the
     // path being restored, so this operation is constant time.
     if (_SetEdgeType(theGraph, v, u) != OK ||
         _SetEdgeType(theGraph, w, x) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _RestoreAndOrientReducedPaths()
 This function searches the embedding for any edges that are specially marked
 as being representative of a path that was previously reduced to a
 single edge by _ReduceBicomp().  The edge is replaced by the path.
 Note that the new path may contain more reduction edges, and these will be
 iteratively expanded by the outer for loop.

 If the edge records of an edge being expanded are the first or last arcs
 of the edge's vertex endpoints, then the edge may be along the external face.
 If so, then the vertices along the path being restored must be given a
 consistent orientation with the endpoints.  It is expected that the embedding
 will have been oriented prior to this operation.
 ****************************************************************************/

int  _RestoreAndOrientReducedPaths(graphP theGraph, K33SearchContext *context)
{
int  e, J, JTwin, u, v, w, x, visited;
int  J0, JTwin0, J1, JTwin1;

     for (e = 0; e < theGraph->M + sp_GetCurrentSize(theGraph->edgeHoles);)
     {
         J = theGraph->edgeOffset + 2*e;
         if (context->G[J].pathConnector != NIL)
         {
             visited = theGraph->G[J].visited;

             JTwin = gp_GetTwinArc(theGraph, J);
             u = theGraph->G[JTwin].v;
             v = context->G[J].pathConnector;
             w = context->G[JTwin].pathConnector;
             x = theGraph->G[J].v;

             /* Now we need the predecessor and successor edge records
                of J and JTwin.  The edge (u, v) will be inserted so
                that the record in u's adjacency list that indicates v
                will be between J0 and J1.  Likewise, the edge record
                (x -> w) will be placed between JTwin0 and JTwin1. */

             J0 = gp_GetNextArc(theGraph, J);
             J1 = gp_GetPrevArc(theGraph, J);
             JTwin0 = gp_GetNextArc(theGraph, JTwin);
             JTwin1 = gp_GetPrevArc(theGraph, JTwin);

             /* We first delete the edge represented by J and JTwin. We do so before
                restoring the path to ensure we do not exceed the maximum arc capacity. */

             gp_DeleteEdge(theGraph, J, 0);

             /* Now we add the two edges to reconnect the reduced path represented
                by the edge [J, JTwin].  The edge record in u is added between J0 and J1.
                Likewise, the new edge record in x is added between JTwin0 and JTwin1. */

             if (gp_IsArc(theGraph, J0))
             {
            	 if (gp_InsertEdge(theGraph, u, J0, 1, v, gp_AdjacencyListEndMark(v), 0) != OK)
            		 return NOTOK;
             }
             else
             {
            	 if (gp_InsertEdge(theGraph, u, J1, 0, v, gp_AdjacencyListEndMark(v), 0) != OK)
            		 return NOTOK;
             }

             if (gp_IsArc(theGraph, JTwin0))
             {
            	 if (gp_InsertEdge(theGraph, x, JTwin0, 1, w, gp_AdjacencyListEndMark(w), 0) != OK)
            		 return NOTOK;
             }
             else
             {
            	 if (gp_InsertEdge(theGraph, x, JTwin1, 0, w, gp_AdjacencyListEndMark(w), 0) != OK)
            		 return NOTOK;
             }

             /* Set the types of the newly added edges */

             if (_SetEdgeType(theGraph, u, v) != OK ||
                 _SetEdgeType(theGraph, w, x) != OK)
                 return NOTOK;

             /* We determine whether the reduction edge may be on the external face,
                in which case we will need to ensure that the vertices on the path
                being restored are consistently oriented.  This will accommodate
                future invocations of MarkPathAlongBicompExtFace().
                Note: If J0, J1, JTwin0 or JTwin1 is not an edge, then it is
                      because we've walked off the end of the edge record list,
                      which happens when J and JTwin are either the first or
                      last edge of the containing vertex.  In turn, the first
                      and last edges of a vertex are the ones that hold it onto
                      the external face, if it is on the external face. */

             if ((!gp_IsArc(theGraph, J0) && !gp_IsArc(theGraph, JTwin1)) ||
                 (!gp_IsArc(theGraph, J1) && !gp_IsArc(theGraph, JTwin0)))
             {
                 if (_OrientExternalFacePath(theGraph, u, v, w, x) != OK)
                     return NOTOK;
             }

             /* The internal XY path was already marked as part of the decision logic
                that made us decide we could find a K_{3,3} and hence that we should
                reverse all of the reductions.  Subsequent code counts on the fact
                that the X-Y path is already marked, so if we replace a marked edge
                with a path, then we need to mark the path. Similarly, for an unmarked
                edge, the replacement path should be unmarked. */

             if (_SetVisitedOnPath(theGraph, u, v, w, x, visited) != OK)
            	 return NOTOK;
         }
         else e++;
     }

     return OK;
}

/****************************************************************************
 _MarkStraddlingBridgePath()
 ****************************************************************************/

int  _MarkStraddlingBridgePath(graphP theGraph, int u_min, int u_max, int u_d, int d)
{
isolatorContextP IC = &theGraph->IC;
int p, J;

/* Find the point of intersection p between the path (v ... u_max)
       and the path (d ... u_max). */

     if (theGraph->functions.fpMarkDFSPath(theGraph, u_max, IC->r) != OK)
         return NOTOK;

     p = d;
     while (!theGraph->G[p].visited)
     {
         theGraph->G[p].visited = 1;

         J = gp_GetFirstArc(theGraph, p);
         while (gp_IsArc(theGraph, J))
         {
              if (theGraph->G[J].type == EDGE_DFSPARENT)
                  break;

              J = gp_GetNextArc(theGraph, J);
         }

         theGraph->G[J].visited = 1;
         theGraph->G[gp_GetTwinArc(theGraph, J)].visited = 1;

         p = theGraph->G[J].v;

         /* If p is a root copy, mark it visited and skip to the parent copy */
         if (p >= theGraph->N)
         {
             theGraph->G[p].visited = 1;
             p = theGraph->V[p-theGraph->N].DFSParent;
         }
     }

/* Unmark the path (p ... u_max), which was marked to help find p.
    The path from v to u_{max} is not needed to form a K_{3,3} except
    for the portion of the path up to p that, with the straddling
    bridge path, comprises part of the connection to u_d. In the
    minor, the path between v and p is edge contracted. */

     while (p != u_max)
     {
         J = gp_GetFirstArc(theGraph, p);
         while (gp_IsArc(theGraph, J))
         {
              if (theGraph->G[J].type == EDGE_DFSPARENT)
                  break;

              J = gp_GetNextArc(theGraph, J);
         }

         theGraph->G[J].visited = 0;
         theGraph->G[gp_GetTwinArc(theGraph, J)].visited = 0;

         p = theGraph->G[J].v;
         theGraph->G[p].visited = 0;

         /* If p is a root copy, clear its visited flag and skip to the
                parent copy */

         if (p >= theGraph->N)
         {
             p = theGraph->V[p-theGraph->N].DFSParent;
             theGraph->G[p].visited = 0;
         }
     }

/* The straddling bridge must join the path (u_max ... u_min).  If u_d is an
    ancestor of u_min, then mark the path that joins u_d to u_min. */

     if (u_d < u_min)
        if (theGraph->functions.fpMarkDFSPath(theGraph, u_d, u_min) != OK)
            return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorE5()
 The paths (x, w), (y, w) and (v, u_{max}) are not needed.
 The x-y path and the internal w-z path are already marked.
 ****************************************************************************/

int  _IsolateMinorE5(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;

     if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->x) != OK ||
         _MarkPathAlongBicompExtFace(theGraph, IC->y, IC->r) != OK ||
         theGraph->functions.fpMarkDFSPath(theGraph, MIN3(IC->ux,IC->uy,IC->uz),
                                                     MAX3(IC->ux,IC->uy,IC->uz)) != OK ||
         _MarkDFSPathsToDescendants(theGraph) != OK ||
         _JoinBicomps(theGraph) != OK ||
         _AddAndMarkUnembeddedEdges(theGraph) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorE6()
 The paths (x, y), (v, w) and (v, u_{max}) are not needed.
 The path through the straddling bridge that connects from an ancestor of
 u_{max} to v is required, but it may connnect to an ancestor p of v.
 In such a case, the path (v, p) is required, while (p, u_{max}) is not.
 ****************************************************************************/

int  _IsolateMinorE6(graphP theGraph, K33SearchContext *context)
{
isolatorContextP IC = &theGraph->IC;
int u_min, u_max, d, u_d;

/* Clear the previously marked x-y path */

     if (_FillVisitedFlagsInBicomp(theGraph, IC->r, 0) != OK)
    	 return NOTOK;

/* Clear dw to stop the marking of path (v, w) */

     IC->dw = NIL;

/* Mark (v, ..., x, ..., w, ..., y, ... v) */

     if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->r) != OK)
         return NOTOK;

/* Mark the path through the straddling bridge (except for the final
     edge (u_d, d) which is added last by convention). */

     u_min = MIN3(IC->ux,IC->uy,IC->uz);
     u_max = MAX3(IC->ux,IC->uy,IC->uz);
     d = _TestForStraddlingBridge(theGraph, context, u_max);
     u_d = theGraph->V[d].leastAncestor;

     if (_MarkStraddlingBridgePath(theGraph, u_min, u_max, u_d, d) != OK)
         return NOTOK;

/* Make the final markings and edge additions */

     if (theGraph->functions.fpMarkDFSPath(theGraph, u_min, u_max) != OK ||
         _MarkDFSPathsToDescendants(theGraph) != OK ||
         _JoinBicomps(theGraph) != OK ||
         _AddAndMarkUnembeddedEdges(theGraph) != OK ||
         _AddAndMarkEdge(theGraph, u_d, d) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorE7()
 ****************************************************************************/

int  _IsolateMinorE7(graphP theGraph, K33SearchContext *context)
{
isolatorContextP IC = &theGraph->IC;
int u_min, u_max, d, u_d;

/* Mark the appropriate two portions of the external face depending on
    symmetry condition */

     if (IC->uy < IC->ux)
     {
         if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->x) != OK ||
             _MarkPathAlongBicompExtFace(theGraph, IC->w, IC->y) != OK)
             return NOTOK;
     }
     else
     {
         if (_MarkPathAlongBicompExtFace(theGraph, IC->x, IC->w) != OK ||
             _MarkPathAlongBicompExtFace(theGraph, IC->y, IC->r) != OK)
             return NOTOK;
     }

/* Mark the path through the straddling bridge (except for the final
     edge (u_d, d) which is added last by convention). */

     u_min = MIN3(IC->ux,IC->uy,IC->uz);
     u_max = MAX3(IC->ux,IC->uy,IC->uz);
     d = _TestForStraddlingBridge(theGraph, context, u_max);
     u_d = theGraph->V[d].leastAncestor;

     if (_MarkStraddlingBridgePath(theGraph, u_min, u_max, u_d, d) != OK)
         return NOTOK;

/* Make the final markings and edge additions */

     if (theGraph->functions.fpMarkDFSPath(theGraph, u_min, u_max) != OK ||
         _MarkDFSPathsToDescendants(theGraph) != OK ||
         _JoinBicomps(theGraph) != OK ||
         _AddAndMarkUnembeddedEdges(theGraph) != OK ||
         _AddAndMarkEdge(theGraph, u_d, d) != OK)
         return NOTOK;

     return OK;
}
