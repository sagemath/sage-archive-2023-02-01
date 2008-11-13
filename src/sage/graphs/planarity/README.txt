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

Overview
========

The main header file, graph.h, contains declarations of all the
functions available.  Here are the main ones to consider:

gp_New() - allocates an empty graph structure and returns a pointer to it.
                Example: graphP theGraph = gp_New();

gp_Free() - frees a graph data structure and nulls out the pointer.
                Take care to pass the address of the pointer returned by gp_New()
                Example: gp_Free(&theGraph);

gp_Read() - initializes a given graph, then adds edges to it according
                to content of a given file.
                Example: if (gp_Read(theGraph, "file.txt") != OK) error

gp_Write() - writes the graph to a file in an adjacency list and adjacency
                matrix formats
                Example: gp_Write(theGraph, "out.txt", WRITE_ADJLIST);

gp_Embed() - the main function that receives a graph and rearranges its
                adjacency lists to produce either a combinatorial planar
                embedding or a minimal non-planar subgraph.
                Example: Result = gp_Embed(theGraph, EMBEDFLAGS_PLANAR);

gp_SortVertices() - can be used after gp_Embed() to recover the original
                        numbering of the graph that appeared, for example, in
                        the input file.  By default, gp_Embed() assumes that the
                        graph should remain with its depth first search numbering,
                        not the original numbering.
                        Example: gp_SortVertices(theGraph);

To create graphs without reading them from a file, the following additional
functions are useful:

gp_InitGraph() - given N, this function allocates within a graph structure
                        enough memory for N vertices and 3N edges.
                        Example:

gp_AddEdge() - allows the addition of a single edge to a previously created
                and initialized graph.
                Example: if (gp_AddEdge(theGraph, u, 0, v, 0) != OK) error


We can find the original distribution at:

http://www.cs.brown.edu/sites/jgaa/volume08.html (http://www.cs.brown.edu/sites/jgaa/accepted/2004/BoyerMyrvold2004.8.3/planarity.zip)

A Dr. Dobbs article that gives a brief overview of the software is at http://www.ddj.com/architect/184406070

The following was in the NOTICE file of the distributed version with
the paper.

Planar Graph Embedder and Non-Planar Subgraph Isolator
Copyright 1999-2005 by John M. Boyer

This Work includes a reference implementation for the following journal paper:

John M. Boyer and Wendy J. Myrvold, On the Cutting Edge: Simplified O(n)
Planarity by Edge Addition.  Journal of Graph Algorithms and Applications,
Vol. 8, No. 3, pp. 241-273, 2004.
