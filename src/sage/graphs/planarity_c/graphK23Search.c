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

#include "graph.h"

/* Imported functions */

extern void _FillVisitedFlags(graphP, int);

extern int  _GetNextVertexOnExternalFace(graphP theGraph, int curVertex, int *pPrevLink);
extern int  _OrientVerticesInBicomp(graphP theGraph, int BicompRoot, int PreserveSigns);
extern int  _JoinBicomps(graphP theGraph);

extern int  _MarkHighestXYPath(graphP theGraph);

extern int  _FindUnembeddedEdgeToAncestor(graphP theGraph, int cutVertex, int *pAncestor, int *pDescendant);
extern int  _FindUnembeddedEdgeToCurVertex(graphP theGraph, int cutVertex, int *pDescendant);
extern int  _FindUnembeddedEdgeToSubtree(graphP theGraph, int ancestor, int SubtreeRoot, int *pDescendant);

extern int  _MarkPathAlongBicompExtFace(graphP theGraph, int startVert, int endVert);

extern int  _AddAndMarkEdge(graphP theGraph, int ancestor, int descendant);

extern int  _DeleteUnmarkedVerticesAndEdges(graphP theGraph);

extern int  _ChooseTypeOfNonOuterplanarityMinor(graphP theGraph, int I, int R);
extern int  _IsolateOuterplanarityObstructionA(graphP theGraph);
extern int  _IsolateOuterplanarityObstructionB(graphP theGraph);

/* Private function declarations for K_{2,3} searching */

int  _SearchForK23(graphP theGraph, int I);

int  _SearchForK23InBicomp(graphP theGraph, int I, int R);
int  _IsolateOuterplanarityObstructionE1orE2(graphP theGraph);
int  _IsolateOuterplanarityObstructionE3orE4(graphP theGraph);

/****************************************************************************
 _SearchForK23()
    We begin by identifying all root copies of I on which the
    Walkdown failed.  We can do this via straightforward traversal of
    descendant to ancestor DFS tree paths.  This is an O(n) total cost in the
    worst case.  If we find a K_{2,3}, then we can afford the time.  However,
    if we find only a K_4, then the outerplanarity algorithm must continue so
    we could not afford the worst case performance.  Fortunately, this
    operation is worst-case constant time if there are no K_{2,3} homeomorphs
    to be found in step I because a non-outerplanar biconnected graph either
    contains a K_{2,3} or is \emph{isomorphic} to K_4.
 ****************************************************************************/

int  _SearchForK23(graphP theGraph, int I)
{
int J, W, C, RetVal=OK;

/* Traverse the edges of I to find the unembedded forward edges to
    descendants.  For each such edge (I, W), traverse the DFS tree
    path up to mark the child of I that is the root of the subtree
    containing W.  Optimize with visitation flag. */

/* Traverse each unembedded back edge to the descendant endpoint... */

    J = theGraph->V[I].fwdArcList;

    // Ensure we have at least one bicomp on which Walkdown failed, which
    // should always be the case in an error free implementation
    if (!gp_IsArc(theGraph, J))
    	return NOTOK;

    while (J != NIL)
    {
        W = theGraph->G[J].v;

        /* Go from the descendant endpoint to find the ancestor that
            is a child of I, which in turn indicates the root of a
            bicomp on which the Walkdown failed to embed all back edges */

        C = W;
        while (theGraph->V[C].DFSParent != I)
            C = theGraph->V[C].DFSParent;

        RetVal = _SearchForK23InBicomp(theGraph, I, C+theGraph->N);

        /* If something went wrong, NOTOK was returned;
        If a K_{2,3} was found, NONEMBEDDABLE was returned;
        If OK was returned, then an isolated K_4 was found, so
        we continue searching any other bicomps on which the
        Walkdown failed. */

        if (RetVal != OK)
            break;

        /* Get the next unembedded back edge from I */

        J = gp_GetNextArc(theGraph, J);
        if (J == theGraph->V[I].fwdArcList)
            J = NIL;
    }

/* If we got through the loop with an OK value for each bicomp on
     which the Walkdown failed, then we return OK to indicate that only
     isolated K_4's were found.  This allows the embedder to continue.
     If a K_{2,3} is ever found (or if an error occurred), then RetVal
     will not be OK, and the loop terminates immediately so we can
     return the appropriate value. */

     return RetVal;
}

/****************************************************************************
 _SearchForK23InBicomp()
 ****************************************************************************/

int  _SearchForK23InBicomp(graphP theGraph, int I, int R)
{
isolatorContextP IC = &theGraph->IC;
int X, Y, XPrevLink, YPrevLink;

/* Begin by determining whether minor A, B or E is detected */

     if (_ChooseTypeOfNonOuterplanarityMinor(theGraph, I, R) != OK)
         return NOTOK;

/* Minors A and B result in the desired K_{2,3} homeomorph,
    so we isolate it and return NONEMBEDDABLE. */

     if (theGraph->IC.minorType & (MINORTYPE_A|MINORTYPE_B))
     {
         _FillVisitedFlags(theGraph, 0);

         if (theGraph->IC.minorType & MINORTYPE_A)
         {
             if (_FindUnembeddedEdgeToCurVertex(theGraph, IC->w, &IC->dw) != TRUE)
                 return NOTOK;

             if (_IsolateOuterplanarityObstructionA(theGraph) != OK)
                 return NOTOK;
         }
         else if (theGraph->IC.minorType & MINORTYPE_B)
         {
         int SubtreeRoot = LCGetPrev(theGraph->BicompLists,
                                     theGraph->V[IC->w].pertinentBicompList, NIL);

             if (_FindUnembeddedEdgeToSubtree(theGraph, IC->v, SubtreeRoot, &IC->dw) != TRUE)
                 return NOTOK;

             if (_IsolateOuterplanarityObstructionB(theGraph) != OK)
                 return NOTOK;
         }

         if (_DeleteUnmarkedVerticesAndEdges(theGraph) != OK)
             return NOTOK;

         return NONEMBEDDABLE;
     }

/* For minor E (a K_4) , we run the additional tests to see if a K_{2,3} is
    entangled with the K_4.  If not, then we return OK to indicate that
    the outerplanarity embedder should proceed as if the K_4 had not
    been found. */

    /* If any vertices other than R, X, Y and W exist along the
        external face, then we can obtain a K_{2,3} by minor E1 or E2 */

     X = IC->x;
     Y = IC->y;
     XPrevLink = 1;
     YPrevLink = 0;
     if (IC->w != _GetNextVertexOnExternalFace(theGraph, X, &XPrevLink) ||
         IC->w != _GetNextVertexOnExternalFace(theGraph, Y, &YPrevLink))
     {
         _FillVisitedFlags(theGraph, 0);

         if (_IsolateOuterplanarityObstructionE1orE2(theGraph) != OK)
             return NOTOK;

         if (_DeleteUnmarkedVerticesAndEdges(theGraph) != OK)
             return NOTOK;

         return NONEMBEDDABLE;
     }

 /* If X, Y or W make either a direct back edge connection or a
        connection through a separated child bicomp to an ancestor of
        the current vertex I, then we can obtain a K_{2,3} by minor
        E3 or E4. Note that this question is query on X, Y and W is
        equivalent to the planarity version of external activity. */

     if (FUTUREPERTINENT(theGraph, X, I) ||
         FUTUREPERTINENT(theGraph, Y, I) ||
         FUTUREPERTINENT(theGraph, IC->w, I))
     {
         _FillVisitedFlags(theGraph, 0);

         if (_IsolateOuterplanarityObstructionE3orE4(theGraph) != OK)
             return NOTOK;

         if (_DeleteUnmarkedVerticesAndEdges(theGraph) != OK)
             return NOTOK;

         return NONEMBEDDABLE;
     }

/* The extra cases for finding a K_{2,3} failed, so the bicomp rooted
    by R is a separable subgraph of the input that is isomorphic
    to K_4.  So, we restore the original vertex orientation of
    the bicomp (because it's polite, not because we really have to).
    Then, we return OK to tell the outerplanarity embedder that it
    can ignore this K_4 and keep processing. */

    if (_OrientVerticesInBicomp(theGraph, R, 1) != OK)
    	return NOTOK;

    return OK;
}

/****************************************************************************
 _IsolateOuterplanarityObstructionE1orE2()
 ****************************************************************************/

int  _IsolateOuterplanarityObstructionE1orE2(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;
int XPrevLink = 1;

     if (_MarkHighestXYPath(theGraph) != TRUE)
         return NOTOK;

/* Isolate E1 */

     if (theGraph->IC.px != theGraph->IC.x)
     {
         if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->w) != OK ||
             _MarkPathAlongBicompExtFace(theGraph, IC->py, IC->r) != OK)
             return NOTOK;
     }
     else if (theGraph->IC.py != theGraph->IC.y)
     {
         if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->x) != OK ||
             _MarkPathAlongBicompExtFace(theGraph, IC->w, IC->r) != OK)
             return NOTOK;
     }

/* Isolate E2 */

     else if (IC->w != _GetNextVertexOnExternalFace(theGraph, IC->x, &XPrevLink))
     {
         if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->y) != OK)
             return NOTOK;
     }

     else
     {
         if (_MarkPathAlongBicompExtFace(theGraph, IC->x, IC->r) != OK)
             return NOTOK;
     }

/* Final bits are in common */

     if (_FindUnembeddedEdgeToCurVertex(theGraph, IC->w, &IC->dw) != TRUE ||
         theGraph->functions.fpMarkDFSPath(theGraph, IC->w, IC->dw) != OK ||
         _JoinBicomps(theGraph) != OK ||
         _AddAndMarkEdge(theGraph, IC->v, IC->dw) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateOuterplanarityObstructionE3orE4()
 ****************************************************************************/

int  _IsolateOuterplanarityObstructionE3orE4(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;
int u, d, XorY;

/* Minor E3 */

     if (FUTUREPERTINENT(theGraph, theGraph->IC.x, theGraph->IC.v) ||
         FUTUREPERTINENT(theGraph, theGraph->IC.y, theGraph->IC.v))
     {
         if (_MarkHighestXYPath(theGraph) != TRUE)
             return NOTOK;

         if (FUTUREPERTINENT(theGraph, theGraph->IC.x, theGraph->IC.v))
              XorY = theGraph->IC.x;
         else XorY = theGraph->IC.y;

         /* The cases of X externally active and Y externally active
                are the same except for the bicomp external face marking
                (because parameter order is important) */

         if (XorY == theGraph->IC.x)
         {
             if (_MarkPathAlongBicompExtFace(theGraph, IC->x, IC->w) != OK ||
                 _MarkPathAlongBicompExtFace(theGraph, IC->y, IC->r) != OK)
                 return NOTOK;
         }
         else
         {
             if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->x) != OK ||
                 _MarkPathAlongBicompExtFace(theGraph, IC->w, IC->y) != OK)
                 return NOTOK;
         }

         if (_FindUnembeddedEdgeToCurVertex(theGraph, IC->w, &IC->dw) != TRUE)
             return NOTOK;

         if (_FindUnembeddedEdgeToAncestor(theGraph, XorY, &u, &d) != TRUE)
             return NOTOK;

         if (theGraph->functions.fpMarkDFSPath(theGraph, u, IC->v) != OK ||
             theGraph->functions.fpMarkDFSPath(theGraph, XorY, d) != OK ||
             theGraph->functions.fpMarkDFSPath(theGraph, IC->w, IC->dw) != OK ||
             _JoinBicomps(theGraph) != OK ||
             _AddAndMarkEdge(theGraph, u, d) != OK ||
             _AddAndMarkEdge(theGraph, IC->v, IC->dw) != OK)
             return NOTOK;

         return OK;
     }

/* Otherwise, isolate Minor E4 (reduce to minor A) */

     if (_FindUnembeddedEdgeToAncestor(theGraph, IC->w, &u, &d) != TRUE)
         return NOTOK;

     IC->v = u;
     IC->dw = d;
     return _IsolateOuterplanarityObstructionA(theGraph);
}
