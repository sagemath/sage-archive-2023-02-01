#ifndef GRAPH_DRAWPLANAR_PRIVATE_H
#define GRAPH_DRAWPLANAR_PRIVATE_H

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

#ifdef __cplusplus
extern "C" {
#endif

// Additional equipment for each graph node (edge arc or vertex)
/*
        pos, start, end: used to store a visibility representation, or
                horvert diagram of a planar graph.
                For vertices, vertical position, horizontal range
                For edges, horizontal position, vertical range
*/
typedef struct
{
     int  pos, start, end;
} DrawPlanar_GraphNode;

typedef DrawPlanar_GraphNode * DrawPlanar_GraphNodeP;

// Additional equipment for each vertex
/*
        drawingFlag, ancestor, ancestorChild: used to collect information needed
                to help 'draw' a visibility representation.  During planar
                embedding, a vertex is determined to be between its DFS parent and
                a given ancestor (the vertex being processed) or beyond the parent
                relative to the ancestor. In post processing, the relative
                orientation of the parent and ancestor are determined,
                then the notion of between/beyond resolves to above/below or
                below/above depending on whether the ancestor is above or below,
                respectively, the parent.  The ancestorChild are used to help r
                esolve this latter question.
        tie[2]  stores information along the external face during embedding
                that is pertinent to helping break ties in the decisions about
                vertical vertex positioning.  When vertices are first merged
                together into a bicomp, we cannot always decide right away which
                vertices will be above or below others.  But as we traverse the
                external face removing inactive vertices, these positional ties
                can be resolved.
*/
typedef struct
{
        int drawingFlag, ancestor, ancestorChild;
        int tie[2];
} DrawPlanar_VertexRec;

typedef DrawPlanar_VertexRec * DrawPlanar_VertexRecP;

#define DRAWINGFLAG_BEYOND     0
#define DRAWINGFLAG_TIE        1
#define DRAWINGFLAG_BETWEEN    2
#define DRAWINGFLAG_BELOW      3
#define DRAWINGFLAG_ABOVE      4

typedef struct
{
    // Helps distinguish initialize from re-initialize
    int initialized;

    // The graph that this context augments
    graphP theGraph;

    // Parallel array for additional graph node level equipment
    DrawPlanar_GraphNodeP G;

    // Parallel array for additional vertex level equipment
    DrawPlanar_VertexRecP V;

    // Overloaded function pointers
    graphFunctionTable functions;

} DrawPlanarContext;

#ifdef __cplusplus
}
#endif

#endif

