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

#define GRAPHISOLATOR_C

#include "graph.h"

/* Imported functions */

extern void _FillVisitedFlags(graphP, int);

extern int  _GetNextVertexOnExternalFace(graphP theEmbedding, int curVertex, int *pPrevLink);
extern int  _JoinBicomps(graphP theEmbedding);
extern void _RecordPertinentChildBicomp(graphP theEmbedding, int I, int RootVertex);
extern int  _GetPertinentChildBicomp(graphP theEmbedding, int W);

extern int _ChooseTypeOfNonplanarityMinor(graphP theEmbedding, int I, int R);

/* Private function declarations (exported within system) */

int _IsolateKuratowskiSubgraph(graphP theEmbedding, int I);

int  _FindUnembeddedEdgeToAncestor(graphP theEmbedding, int cutVertex,
                                   int *pAncestor, int *pDescendant);
int  _FindUnembeddedEdgeToCurVertex(graphP theEmbedding, int cutVertex,
                                    int *pDescendant);
int  _FindUnembeddedEdgeToSubtree(graphP theEmbedding, int ancestor,
                                  int SubtreeRoot, int *pDescendant);

int  _MarkDFSPath(graphP theEmbedding, int ancestor, int descendant);
int  _MarkPathAlongBicompExtFace(graphP theEmbedding, int startVert, int endVert);

int  _AddAndMarkEdge(graphP theEmbedding, int ancestor, int descendant);
void _AddBackEdge(graphP theEmbedding, int ancestor, int descendant);
int  _DeleteUnmarkedVerticesAndEdges(graphP theEmbedding);

int  _InitializeIsolatorContext(graphP theEmbedding);

int  _IsolateMinorA(graphP theEmbedding);
int  _IsolateMinorB(graphP theEmbedding);
int  _IsolateMinorC(graphP theEmbedding);
int  _IsolateMinorD(graphP theEmbedding);
int  _IsolateMinorE(graphP theEmbedding);

int  _IsolateMinorE1(graphP theEmbedding);
int  _IsolateMinorE2(graphP theEmbedding);
int  _IsolateMinorE3(graphP theEmbedding);
int  _IsolateMinorE4(graphP theEmbedding);

int  _GetLeastAncestorConnection(graphP theEmbedding, int cutVertex);
int  _MarkDFSPathsToDescendants(graphP theEmbedding);
int  _AddAndMarkUnembeddedEdges(graphP theEmbedding);

/****************************************************************************
 gp_IsolateKuratowskiSubgraph()
 ****************************************************************************/

int  _IsolateKuratowskiSubgraph(graphP theEmbedding, int I)
{
int  RetVal;

/* Begin by determining which of the non-planarity Minors was encountered
        and the principal bicomp on which the isolator will focus attention. */

     if (_ChooseTypeOfNonplanarityMinor(theEmbedding, I, NIL) != OK)
         return NOTOK;

     if (_InitializeIsolatorContext(theEmbedding) != OK)
         return NOTOK;

/* Call the appropriate isolator */

     if (theEmbedding->IC.minorType & FLAGS_MINOR_A)
         RetVal = _IsolateMinorA(theEmbedding);
     else if (theEmbedding->IC.minorType & FLAGS_MINOR_B)
         RetVal = _IsolateMinorB(theEmbedding);
     else if (theEmbedding->IC.minorType & FLAGS_MINOR_C)
         RetVal = _IsolateMinorC(theEmbedding);
     else if (theEmbedding->IC.minorType & FLAGS_MINOR_D)
         RetVal = _IsolateMinorD(theEmbedding);
     else if (theEmbedding->IC.minorType & FLAGS_MINOR_E)
         RetVal = _IsolateMinorE(theEmbedding);

/* Delete the unmarked edges and vertices, and return */

     if (RetVal == OK)
         RetVal = _DeleteUnmarkedVerticesAndEdges(theEmbedding);

     return RetVal;
}

/****************************************************************************
 _InitializeIsolatorContext()
 ****************************************************************************/

int  _InitializeIsolatorContext(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

/* Obtains the edges connecting X and Y to ancestors of the current vertex */

     if (_FindUnembeddedEdgeToAncestor(theEmbedding, IC->x, &IC->ux, &IC->dx) != OK ||
         _FindUnembeddedEdgeToAncestor(theEmbedding, IC->y, &IC->uy, &IC->dy) != OK)
         return NOTOK;

/* For Minor B, we seek the last pertinent child biconnected component, which
     is externally active, and obtain the DFS child in its root edge.
     This child is the subtree root containing vertices with connections to
     both the current vertex and an ancestor of the current vertex. */

     if (theEmbedding->IC.minorType & FLAGS_MINOR_B)
     {
     int SubtreeRoot = LCGetPrev(theEmbedding->BicompLists,
                                 theEmbedding->V[IC->w].pertinentBicompList, NIL);

         IC->uz = theEmbedding->V[SubtreeRoot].Lowpoint;

         if (_FindUnembeddedEdgeToSubtree(theEmbedding, IC->v, SubtreeRoot, &IC->dw) != OK ||
             _FindUnembeddedEdgeToSubtree(theEmbedding, IC->uz, SubtreeRoot, &IC->dz) != OK)
             return NOTOK;
     }

/* For all other minors, we obtain  */

     else
     {
         if (_FindUnembeddedEdgeToCurVertex(theEmbedding, IC->w, &IC->dw) != OK)
             return NOTOK;

         if (theEmbedding->IC.minorType & FLAGS_MINOR_E)
             if (_FindUnembeddedEdgeToAncestor(theEmbedding, IC->z, &IC->uz, &IC->dz) != OK)
                 return NOTOK;
     }

     return OK;
}

/****************************************************************************
 _IsolateMinorA(): Isolate a K3,3 homeomorph
 ****************************************************************************/

int  _IsolateMinorA(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (_MarkPathAlongBicompExtFace(theEmbedding, IC->r, IC->r) != OK ||
         _MarkDFSPath(theEmbedding, MIN(IC->ux, IC->uy), IC->r) != OK ||
         _MarkDFSPathsToDescendants(theEmbedding) != OK ||
         _JoinBicomps(theEmbedding) != OK ||
         _AddAndMarkUnembeddedEdges(theEmbedding) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorB(): Isolate a K3,3 homeomorph
 ****************************************************************************/

int  _IsolateMinorB(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (_MarkPathAlongBicompExtFace(theEmbedding, IC->r, IC->r) != OK ||
         _MarkDFSPath(theEmbedding, MIN3(IC->ux,IC->uy,IC->uz),
                                    MAX3(IC->ux,IC->uy,IC->uz)) != OK ||
         _MarkDFSPathsToDescendants(theEmbedding) != OK ||
         _JoinBicomps(theEmbedding) != OK ||
         _AddAndMarkUnembeddedEdges(theEmbedding) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorC(): Isolate a K3,3 homeomorph
 ****************************************************************************/

int  _IsolateMinorC(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (theEmbedding->G[IC->px].type == VERTEX_HIGH_RXW)
     {
     int highY = theEmbedding->G[IC->py].type == VERTEX_HIGH_RYW
                 ? IC->py : IC->y;
         if (_MarkPathAlongBicompExtFace(theEmbedding, IC->r, highY) != OK)
             return NOTOK;
     }
     else
     {
         if (_MarkPathAlongBicompExtFace(theEmbedding, IC->x, IC->r) != OK)
             return NOTOK;
     }

     if (_MarkDFSPathsToDescendants(theEmbedding) != OK ||
         _MarkDFSPath(theEmbedding, MIN(IC->ux, IC->uy), IC->r) != OK ||
         _JoinBicomps(theEmbedding) != OK ||
         _AddAndMarkUnembeddedEdges(theEmbedding) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorD(): Isolate a K3,3 homeomorph
 ****************************************************************************/

int  _IsolateMinorD(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (_MarkPathAlongBicompExtFace(theEmbedding, IC->x, IC->y) != OK ||
         _MarkDFSPath(theEmbedding, MIN(IC->ux, IC->uy), IC->r) != OK ||
         _MarkDFSPathsToDescendants(theEmbedding) != OK ||
         _JoinBicomps(theEmbedding) != OK ||
         _AddAndMarkUnembeddedEdges(theEmbedding) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorE()
 ****************************************************************************/

int  _IsolateMinorE(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

/* Minor E1: Isolate a K3,3 homeomorph */

     if (IC->z != IC->w)
         return _IsolateMinorE1(theEmbedding);

/* Minor E2: Isolate a K3,3 homeomorph */

     if (IC->uz > MAX(IC->ux, IC->uy))
         return _IsolateMinorE2(theEmbedding);

/* Minor E3: Isolate a K3,3 homeomorph */

     if (IC->uz < MAX(IC->ux, IC->uy) && IC->ux != IC->uy)
         return _IsolateMinorE3(theEmbedding);

/* Minor E4: Isolate a K3,3 homeomorph */

     else if (IC->x != IC->px || IC->y != IC->py)
         return _IsolateMinorE4(theEmbedding);

/* Minor E: Isolate a K5 homeomorph */

     if (_MarkPathAlongBicompExtFace(theEmbedding, IC->r, IC->r) != OK ||
         _MarkDFSPath(theEmbedding, MIN3(IC->ux, IC->uy, IC->uz), IC->r) != OK ||
         _MarkDFSPathsToDescendants(theEmbedding) != OK ||
         _JoinBicomps(theEmbedding) != OK ||
         _AddAndMarkUnembeddedEdges(theEmbedding) != OK)
         return NOTOK;

     return OK;
}

/****************************************************************************
 _IsolateMinorE1()

 Reduce to Minor C if the vertex Z responsible for external activity
 below the X-Y path does not equal the pertinent vertex W.
 ****************************************************************************/

int  _IsolateMinorE1(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (theEmbedding->G[IC->z].type == VERTEX_LOW_RXW)
     {
         theEmbedding->G[IC->px].type = VERTEX_HIGH_RXW;
         IC->x=IC->z; IC->ux=IC->uz; IC->dx=IC->dz;
     }
     else if (theEmbedding->G[IC->z].type == VERTEX_LOW_RYW)
     {
         theEmbedding->G[IC->py].type = VERTEX_HIGH_RYW;
         IC->y=IC->z; IC->uy=IC->uz; IC->dy=IC->dz;
     }
     else return NOTOK;

     IC->z = IC->uz = IC->dz = NIL;
     theEmbedding->IC.minorType ^= FLAGS_MINOR_E;
     theEmbedding->IC.minorType |= (FLAGS_MINOR_C|FLAGS_MINOR_E1);
     return _IsolateMinorC(theEmbedding);
}

/****************************************************************************
 _IsolateMinorE2()

 If uZ (which is the ancestor of I that is adjacent to Z) is a
 descendant of both uY and uX, then we reduce to Minor A
 ****************************************************************************/

int  _IsolateMinorE2(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     _FillVisitedFlags(theEmbedding, 0);

     IC->v = IC->uz;
     IC->dw = IC->dz;
     IC->z = IC->uz = IC->dz = NIL;

     theEmbedding->IC.minorType ^= FLAGS_MINOR_E;
     theEmbedding->IC.minorType |= (FLAGS_MINOR_A|FLAGS_MINOR_E2);
     return _IsolateMinorA(theEmbedding);
}

/****************************************************************************
 _IsolateMinorE3()
 ****************************************************************************/

int  _IsolateMinorE3(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (IC->ux < IC->uy)
     {
         if (_MarkPathAlongBicompExtFace(theEmbedding, IC->r, IC->px) != OK ||
             _MarkPathAlongBicompExtFace(theEmbedding, IC->w, IC->y) != OK)
             return NOTOK;
     }
     else
     {
         if (_MarkPathAlongBicompExtFace(theEmbedding, IC->x, IC->w) != OK ||
             _MarkPathAlongBicompExtFace(theEmbedding, IC->py, IC->r) != OK)
             return NOTOK;
     }

     if (_MarkDFSPath(theEmbedding, MIN3(IC->ux, IC->uy, IC->uz), IC->r) != OK ||
         _MarkDFSPathsToDescendants(theEmbedding) != OK ||
         _JoinBicomps(theEmbedding) != OK ||
         _AddAndMarkUnembeddedEdges(theEmbedding) != OK)
         return NOTOK;

     theEmbedding->IC.minorType |= FLAGS_MINOR_E3;
     return OK;
}

/****************************************************************************
 _IsolateMinorE4()
 ****************************************************************************/

int  _IsolateMinorE4(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (IC->px != IC->x)
     {
         if (_MarkPathAlongBicompExtFace(theEmbedding, IC->r, IC->w) != OK ||
             _MarkPathAlongBicompExtFace(theEmbedding, IC->py, IC->r) != OK)
             return NOTOK;
     }
     else
     {
         if (_MarkPathAlongBicompExtFace(theEmbedding, IC->r, IC->px) != OK ||
             _MarkPathAlongBicompExtFace(theEmbedding, IC->w, IC->r) != OK)
             return NOTOK;
     }

     if (_MarkDFSPath(theEmbedding, MIN3(IC->ux, IC->uy, IC->uz),
                                    MAX3(IC->ux, IC->uy, IC->uz)) != OK ||
         _MarkDFSPathsToDescendants(theEmbedding) != OK ||
         _JoinBicomps(theEmbedding) != OK ||
         _AddAndMarkUnembeddedEdges(theEmbedding) != OK)
         return NOTOK;

     theEmbedding->IC.minorType |= FLAGS_MINOR_E4;
     return OK;
}

/****************************************************************************
 _GetLeastAncestorConnection()

 This function searches for an ancestor of the current vertex I adjacent by a
 cycle edge to the given cutVertex or one of its DFS descendants appearing in
 a separated bicomp. The given cutVertex is assumed to be externally active
 such that either the leastAncestor or the lowpoint of a separated DFS child
 is less than I.  We obtain the minimum possible connection from the cutVertex
 to an ancestor of I.
 ****************************************************************************/

int  _GetLeastAncestorConnection(graphP theEmbedding, int cutVertex)
{
int  subtreeRoot = theEmbedding->V[cutVertex].separatedDFSChildList;
int  ancestor = theEmbedding->V[cutVertex].leastAncestor;

     if (subtreeRoot != NIL &&
         ancestor > theEmbedding->V[subtreeRoot].Lowpoint)
         ancestor = theEmbedding->V[subtreeRoot].Lowpoint;

     return ancestor;
}

/****************************************************************************
 _FindUnembeddedEdgeToAncestor()

 This function searches for an ancestor of the current vertex I adjacent by a
 cycle edge to the given cutVertex or one of its DFS descendants appearing in
 a separated bicomp.

 The given cutVertex is assumed to be externally active such that either the
 leastAncestor or the lowpoint of a separated DFS child is less than I.
 We obtain the minimum possible connection from the cutVertex to an ancestor
 of I, then compute the descendant accordingly.
 ****************************************************************************/

int  _FindUnembeddedEdgeToAncestor(graphP theEmbedding, int cutVertex,
                                   int *pAncestor, int *pDescendant)
{
     *pAncestor = _GetLeastAncestorConnection(theEmbedding, cutVertex);

     if (*pAncestor == theEmbedding->V[cutVertex].leastAncestor)
     {
         *pDescendant = cutVertex;
         return OK;
     }
     else
     {
     int subtreeRoot = theEmbedding->V[cutVertex].separatedDFSChildList;

         return _FindUnembeddedEdgeToSubtree(theEmbedding, *pAncestor,
                                             subtreeRoot, pDescendant);
     }
}

/****************************************************************************
 _FindUnembeddedEdgeToCurVertex()

 Given the current vertex I, we search for an edge connecting I to either
 a given pertinent vertex W or one of its DFS descendants in the subtree
 indicated by the the last pertinent child biconnected component.
 ****************************************************************************/

int  _FindUnembeddedEdgeToCurVertex(graphP theEmbedding, int cutVertex, int *pDescendant)
{
int  RetVal = OK, I = theEmbedding->IC.v;

     if (theEmbedding->V[cutVertex].adjacentTo != NIL)
         *pDescendant = cutVertex;
     else
     {
     int subtreeRoot = theEmbedding->V[cutVertex].pertinentBicompList;

         RetVal = _FindUnembeddedEdgeToSubtree(theEmbedding, I,
                                               subtreeRoot, pDescendant);
     }

     return RetVal;
}

/****************************************************************************
 _FindUnembeddedEdgeToSubtree()

 Given the root vertex of a DFS subtree and an ancestor of that subtree,
 find a vertex in the subtree that is adjacent to the ancestor by a
 cycle edge.
 ****************************************************************************/

int  _FindUnembeddedEdgeToSubtree(graphP theEmbedding, int ancestor,
                                  int SubtreeRoot, int *pDescendant)
{
int  J, Z, ZNew;

     *pDescendant = NIL;

/* If SubtreeRoot is a root copy, then we change to the DFS child in the
        DFS tree root edge of the bicomp rooted by SubtreeRoot. */

     if (SubtreeRoot >= theEmbedding->N)
         SubtreeRoot -= theEmbedding->N;

/* Find the least descendant of the cut vertex incident to the ancestor. */

     J = theEmbedding->V[ancestor].fwdArcList;
     while (J != NIL)
     {
          if (theEmbedding->G[J].v >= SubtreeRoot)
          {
              if (*pDescendant == NIL || *pDescendant > theEmbedding->G[J].v)
                  *pDescendant = theEmbedding->G[J].v;
          }

          J = theEmbedding->G[J].link[0];
          if (J == theEmbedding->V[ancestor].fwdArcList)
              J = NIL;
     }

     if (*pDescendant == NIL) return NOTOK;

/* Make sure the identified descendant actually descends from the cut vertex */

     Z = *pDescendant;
     while (Z != SubtreeRoot)
     {
         ZNew = theEmbedding->V[Z].DFSParent;
         if (ZNew == NIL || ZNew == Z)
             return NOTOK;
         Z = ZNew;
     }

/* Return successfully */

     return OK;
}


/****************************************************************************
 _MarkPathAlongBicompExtFace()

 Sets the visited flags of vertices and edges on the external face of a
 bicomp from startVert to endVert, inclusive, by following the link[0] out of
 each visited vertex.
 ****************************************************************************/

int  _MarkPathAlongBicompExtFace(graphP theEmbedding, int startVert, int endVert)
{
int  Z, ZPrevLink, ZPrevArc;

/* Mark the start vertex (and if it is a root copy, mark the parent copy too. */

     theEmbedding->G[startVert].visited = 1;

/* For each vertex visited after the start vertex, mark the vertex and the
        edge used to get there.  Stop after marking the ending vertex. */

     Z = startVert;
     do {
        ZPrevLink = 1;
        Z = _GetNextVertexOnExternalFace(theEmbedding, Z, &ZPrevLink);

        ZPrevArc = theEmbedding->G[Z].link[ZPrevLink];

        theEmbedding->G[ZPrevArc].visited = 1;
        theEmbedding->G[gp_GetTwinArc(theEmbedding, ZPrevArc)].visited = 1;
        theEmbedding->G[Z].visited = 1;

     } while (Z != endVert);

     return OK;
}

/****************************************************************************
 _MarkDFSPath()

 Sets visited flags of vertices and edges from descendant to ancestor,
 including root copy vertices, and including the step of hopping from
 a root copy to its parent copy.

 The DFSParent of a vertex indicates the next vertex to visit, so we
 search the adjacency list looking for the edge leading either to the
 DFSParent or to a root copy of the DFSParent.  We start by marking the
 descendant, then traverse up the DFS tree, stopping after we mark
 the ancestor.
 ****************************************************************************/

int  _MarkDFSPath(graphP theEmbedding, int ancestor, int descendant)
{
int  J, parent, Z, N;

     N = theEmbedding->N;

     /* If we are marking from a root vertex upward, then go up to the parent
        copy before starting the loop */

     if (descendant >= N)
         descendant = theEmbedding->V[descendant-N].DFSParent;

     /* Mark the lowest vertex */
     theEmbedding->G[descendant].visited = 1;

     /* Mark all ancestors of the lowest vertex, and the edges used to reach
        them, up to the given ancestor vertex. */

     while (descendant != ancestor)
     {
          /* Get the parent vertex */

          parent = theEmbedding->V[descendant].DFSParent;

          /* If the descendant was a DFS tree root, then obviously
                we aren't going to find the ancestor, so something is wrong.*/

          if (parent == NIL || parent == descendant)
              return NOTOK;

          /* Find the edge from descendant that leads either to
                parent or to a root copy of the parent.
                When the edge is found, mark it and break the loop */

          J = theEmbedding->G[descendant].link[0];
          while (J >= 2*theEmbedding->N)
          {
              Z = theEmbedding->G[J].v;
              if (Z < N && Z == parent ||
                  Z >= N && theEmbedding->V[Z-N].DFSParent == parent)
              {
                  theEmbedding->G[J].visited = 1;
                  theEmbedding->G[gp_GetTwinArc(theEmbedding, J)].visited = 1;
                  break;
              }
              J = theEmbedding->G[J].link[0];
          }

          /* Mark the parent copy of the DFS parent */
          theEmbedding->G[parent].visited = 1;

          /* Hop to the parent */
          descendant = parent;
     }

     return OK;
}

/****************************************************************************
 _MarkDFSPathsToDescendants()
 ****************************************************************************/

int  _MarkDFSPathsToDescendants(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (_MarkDFSPath(theEmbedding, IC->x, IC->dx) != OK ||
         _MarkDFSPath(theEmbedding, IC->y, IC->dy) != OK)
         return NOTOK;

     if (IC->dw != NIL)
         if (_MarkDFSPath(theEmbedding, IC->w, IC->dw) != OK)
             return NOTOK;

     if (IC->dz != NIL)
         if (_MarkDFSPath(theEmbedding, IC->w, IC->dz) != OK)
             return NOTOK;

     return OK;
}

/****************************************************************************
 _AddAndMarkUnembeddedEdges()
 ****************************************************************************/

int  _AddAndMarkUnembeddedEdges(graphP theEmbedding)
{
isolatorContextP IC = &theEmbedding->IC;

     if (_AddAndMarkEdge(theEmbedding, IC->ux, IC->dx) != OK ||
         _AddAndMarkEdge(theEmbedding, IC->uy, IC->dy) != OK)
         return NOTOK;

     if (IC->dw != NIL)
         if (_AddAndMarkEdge(theEmbedding, IC->v, IC->dw) != OK)
             return NOTOK;

     if (IC->dz != NIL)
         if (_AddAndMarkEdge(theEmbedding, IC->uz, IC->dz) != OK)
             return NOTOK;

     return OK;
}

/****************************************************************************
 _AddAndMarkEdge()

 Adds edge records for the edge (ancestor, descendant) and marks the edge
 records and vertex structures that represent the edge.
 ****************************************************************************/

int _AddAndMarkEdge(graphP theEmbedding, int ancestor, int descendant)
{
    _AddBackEdge(theEmbedding, ancestor, descendant);

    /* Mark the edge so it is not deleted */

    theEmbedding->G[ancestor].visited = 1;
    theEmbedding->G[theEmbedding->G[ancestor].link[0]].visited = 1;
    theEmbedding->G[theEmbedding->G[descendant].link[0]].visited = 1;
    theEmbedding->G[descendant].visited = 1;

    return OK;
}

/****************************************************************************
 _AddBackEdge()

 This function transfers the edge records for the edge between the ancestor
 and descendant from the forward edge list of the ancestor to the adjacency
 lists of the ancestor and descendant.
 ****************************************************************************/

void _AddBackEdge(graphP theEmbedding, int ancestor, int descendant)
{
int fwdArc, backArc;

    /* We get the two edge records of the back edge to embed. */

     fwdArc = theEmbedding->V[ancestor].fwdArcList;
     while (fwdArc != NIL)
     {
          if (theEmbedding->G[fwdArc].v == descendant)
              break;

          fwdArc = theEmbedding->G[fwdArc].link[0];
          if (fwdArc == theEmbedding->V[ancestor].fwdArcList)
              fwdArc = NIL;
     }

     if (fwdArc == NIL)
         return;

    backArc = gp_GetTwinArc(theEmbedding, fwdArc);

    /* The forward arc is removed from the fwdArcList of the ancestor. */

    if (theEmbedding->V[ancestor].fwdArcList == fwdArc)
    {
        if (theEmbedding->G[fwdArc].link[0] == fwdArc)
             theEmbedding->V[ancestor].fwdArcList = NIL;
        else theEmbedding->V[ancestor].fwdArcList = theEmbedding->G[fwdArc].link[0];
    }

    theEmbedding->G[theEmbedding->G[fwdArc].link[0]].link[1] = theEmbedding->G[fwdArc].link[1];
    theEmbedding->G[theEmbedding->G[fwdArc].link[1]].link[0] = theEmbedding->G[fwdArc].link[0];

    /* The forward arc is added to the adjacency list of the ancestor. */

    theEmbedding->G[fwdArc].link[1] = ancestor;
    theEmbedding->G[fwdArc].link[0] = theEmbedding->G[ancestor].link[0];
    theEmbedding->G[theEmbedding->G[ancestor].link[0]].link[1] = fwdArc;
    theEmbedding->G[ancestor].link[0] = fwdArc;

    /* The back arc is added to the adjacency list of the descendant. */

    theEmbedding->G[backArc].v = ancestor;

    theEmbedding->G[backArc].link[1] = descendant;
    theEmbedding->G[backArc].link[0] = theEmbedding->G[descendant].link[0];
    theEmbedding->G[theEmbedding->G[descendant].link[0]].link[1] = backArc;
    theEmbedding->G[descendant].link[0] = backArc;
}

/****************************************************************************
 _DeleteUnmarkedVerticesAndEdges()

 For each vertex, traverse its adjacency list and delete all unvisited edges.
 ****************************************************************************/

int  _DeleteUnmarkedVerticesAndEdges(graphP theEmbedding)
{
int  I, J, fwdArc, descendant;

     /* All of the forward and back arcs of all of the edge records
        were removed from the adjacency lists in the planarity algorithm
        preprocessing.  We now put them back into the adjacency lists
        (and we do not mark them), so they can be properly deleted below. */

     for (I = 0; I < theEmbedding->N; I++)
     {
         while (theEmbedding->V[I].fwdArcList != NIL)
         {
             fwdArc = theEmbedding->V[I].fwdArcList;
             descendant = theEmbedding->G[fwdArc].v;
             _AddBackEdge(theEmbedding, I, descendant);
         }
     }

     /* Now we delete all unmarked edges.  We don't delete vertices
        from the embedding, but the ones we should delete will become
        degree zero. */

     for (I = 0; I < theEmbedding->N; I++)
     {
          J = theEmbedding->G[I].link[0];
          while (J >= 2*theEmbedding->N)
          {
                if (theEmbedding->G[J].visited)
                     J = theEmbedding->G[J].link[0];
                else J = gp_DeleteEdge(theEmbedding, J, 0);
          }
     }

     return OK;
}
