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
extern void _ClearIsolatorContext(graphP theGraph);

extern int  _GetNextVertexOnExternalFace(graphP theGraph, int curVertex, int *pPrevLink);
extern int  _JoinBicomps(graphP theGraph);

extern int  _InitializeNonplanarityContext(graphP theGraph, int I, int R);
extern int  _MarkHighestXYPath(graphP theGraph);

extern int  _FindUnembeddedEdgeToAncestor(graphP theGraph, int cutVertex, int *pAncestor, int *pDescendant);
extern int  _FindUnembeddedEdgeToCurVertex(graphP theGraph, int cutVertex, int *pDescendant);
extern int  _FindUnembeddedEdgeToSubtree(graphP theGraph, int ancestor, int SubtreeRoot, int *pDescendant);

extern int  _MarkPathAlongBicompExtFace(graphP theGraph, int startVert, int endVert);

extern int  _AddAndMarkEdge(graphP theGraph, int ancestor, int descendant);

extern int  _DeleteUnmarkedVerticesAndEdges(graphP theGraph);

/* Private function declarations (exported to system) */

int  _IsolateOuterplanarObstruction(graphP theGraph, int I, int R);

int  _ChooseTypeOfNonOuterplanarityMinor(graphP theGraph, int I, int R);

int  _IsolateOuterplanarityObstructionA(graphP theGraph);
int  _IsolateOuterplanarityObstructionB(graphP theGraph);
int  _IsolateOuterplanarityObstructionE(graphP theGraph);

/****************************************************************************
 _ChooseTypeOfNonOuterplanarityMinor()
 A constant time implementation is easily feasible but only constant amortized
 time is needed for the outerplanarity obstruction isolation, which also
 benefits from having the bicomp rooted by R oriented.
 If an extension algorithm requires constant actual time, then this function
 should not be used and instead the minor should be decided without orienting
 the bicomp.
 ****************************************************************************/

int  _ChooseTypeOfNonOuterplanarityMinor(graphP theGraph, int I, int R)
{
int  N, X, Y, W;

	 // Create the initial non-outerplanarity obstruction isolator state.
     if (_InitializeNonplanarityContext(theGraph, I, R) != OK)
         return NOTOK;

     N = theGraph->N;
     R = theGraph->IC.r;
     X = theGraph->IC.x;
     Y = theGraph->IC.y;
     W = theGraph->IC.w;

     // If the root copy is not a root copy of the current vertex I,
     // then the Walkdown terminated on a descendant bicomp, which is Minor A.
     if (theGraph->V[R - N].DFSParent != I)
     {
         theGraph->IC.minorType |= MINORTYPE_A;
         return OK;
     }

     // If W has a pertinent child bicomp, then we've found Minor B.
     // Notice this is different from planarity, in which minor B is indicated
     // only if the pertinent child bicomp is also externally active under the
     // planarity processing model (i.e. future pertinent).
     if (theGraph->V[W].pertinentBicompList != NIL)
     {
         theGraph->IC.minorType |= MINORTYPE_B;
         return OK;
     }

     // The only other result is minor E (we will search for the X-Y path later)
     theGraph->IC.minorType |= MINORTYPE_E;
     return OK;
}

/****************************************************************************
 _IsolateOuterplanarObstruction()
 ****************************************************************************/

int  _IsolateOuterplanarObstruction(graphP theGraph, int I, int R)
{
int  RetVal;

/* A subgraph homeomorphic to K_{2,3} or K_4 will be isolated by using the visited
   flags, 1=keep edge/vertex and 0=omit. Here we initialize to omit all, then we
   subsequently set visited to 1 on all edges and vertices in the homeomorph. */

	 _FillVisitedFlags(theGraph, 0);

/* Next we determineg which of the non-outerplanarity Minors was encountered
        and the principal bicomp on which the isolator will focus attention. */

     if (_ChooseTypeOfNonOuterplanarityMinor(theGraph, I, R) != OK)
         return NOTOK;

/* Find the path connecting the pertinent vertex w with the current vertex v */

     if (theGraph->IC.minorType & MINORTYPE_B)
     {
     isolatorContextP IC = &theGraph->IC;
     int SubtreeRoot = LCGetPrev(theGraph->BicompLists,
                                 theGraph->V[IC->w].pertinentBicompList, NIL);

         if (_FindUnembeddedEdgeToSubtree(theGraph, IC->v, SubtreeRoot, &IC->dw) != TRUE)
             return NOTOK;
     }
     else
     {
     isolatorContextP IC = &theGraph->IC;

         if (_FindUnembeddedEdgeToCurVertex(theGraph, IC->w, &IC->dw) != TRUE)
             return NOTOK;
     }

/* For minor E, we need to find and mark an X-Y path */

     if (theGraph->IC.minorType & MINORTYPE_E)
     {
        if (_MarkHighestXYPath(theGraph) != TRUE)
             return NOTOK;
     }

/* Call the appropriate isolator */

     if (theGraph->IC.minorType & MINORTYPE_A)
         RetVal = _IsolateOuterplanarityObstructionA(theGraph);
     else if (theGraph->IC.minorType & MINORTYPE_B)
         RetVal = _IsolateOuterplanarityObstructionB(theGraph);
     else if (theGraph->IC.minorType & MINORTYPE_E)
         RetVal = _IsolateOuterplanarityObstructionE(theGraph);
     else
    	 RetVal = NOTOK;

/* Delete the unmarked edges and vertices, and return */

     if (RetVal == OK)
         RetVal = _DeleteUnmarkedVerticesAndEdges(theGraph);

     return RetVal;
}

/****************************************************************************
 _IsolateOuterplanarityObstructionA(): Isolate a K2,3 homeomorph
 ****************************************************************************/

int  _IsolateOuterplanarityObstructionA(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;

     if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->r) != OK ||
         theGraph->functions.fpMarkDFSPath(theGraph, IC->v, IC->r) != OK ||
         theGraph->functions.fpMarkDFSPath(theGraph, IC->w, IC->dw) != OK ||
         _JoinBicomps(theGraph) != OK ||
         _AddAndMarkEdge(theGraph, IC->v, IC->dw) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateOuterplanarityObstructionB(): Isolate a K2,3 homeomorph
 ****************************************************************************/

int  _IsolateOuterplanarityObstructionB(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;

     if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->r) != OK ||
         theGraph->functions.fpMarkDFSPath(theGraph, IC->w, IC->dw) != OK ||
         _JoinBicomps(theGraph) != OK ||
         _AddAndMarkEdge(theGraph, IC->v, IC->dw) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateOuterplanarityObstructionE(): Isolate a K4 homeomorph
 ****************************************************************************/

int  _IsolateOuterplanarityObstructionE(graphP theGraph)
{
isolatorContextP IC = &theGraph->IC;

     if (_MarkPathAlongBicompExtFace(theGraph, IC->r, IC->r) != OK ||
         theGraph->functions.fpMarkDFSPath(theGraph, IC->w, IC->dw) != OK ||
         _JoinBicomps(theGraph) != OK ||
         _AddAndMarkEdge(theGraph, IC->v, IC->dw) != OK)
         return NOTOK;

     return OK;
}
