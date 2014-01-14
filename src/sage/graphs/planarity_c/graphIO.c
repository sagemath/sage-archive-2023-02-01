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
#include <string.h>

#include "graph.h"

/* Private functions (exported to system) */

int  _ReadAdjMatrix(graphP theGraph, FILE *Infile);
int  _ReadAdjList(graphP theGraph, FILE *Infile);
int  _WriteAdjList(graphP theGraph, FILE *Outfile);
int  _WriteAdjMatrix(graphP theGraph, FILE *Outfile);
int  _WriteDebugInfo(graphP theGraph, FILE *Outfile);

/********************************************************************
 _ReadAdjMatrix()
 This function reads the undirected graph in upper triangular matrix format.
 Though O(N^2) time is required, this routine is useful during
 reliability testing due to the wealth of graph generating software
 that uses this format for output.
 Returns: OK, NOTOK on internal error, NONEMBEDDABLE if too many edges
 ********************************************************************/

int _ReadAdjMatrix(graphP theGraph, FILE *Infile)
{
int N, I, W, Flag, ErrorCode;

     if (Infile == NULL) return NOTOK;
     fscanf(Infile, " %d ", &N);
     if (gp_InitGraph(theGraph, N) != OK)
          return NOTOK;

     for (I = 0, ErrorCode = OK; I < N-1 && ErrorCode==OK; I++)
     {
          theGraph->G[I].v = I;
          for (W = I+1; W < N; W++)
          {
               fscanf(Infile, " %1d", &Flag);
               if (Flag)
               {
                   ErrorCode = gp_AddEdge(theGraph, I, 0, W, 0);
                   if (ErrorCode != OK) break;
               }
          }
     }

     return ErrorCode;
}

/********************************************************************
 _ReadAdjList()
 This function reads the graph in adjacency list format.

 The file format is
 On the first line    : N= number of vertices
 On N subsequent lines: #: a b c ... -1
 where # is a vertex number and a, b, c, ... are its neighbors.

 NOTE:  The vertex number is for file documentation only.  It is an
        error if the vertices are not in sorted order in the file.

 NOTE:  If a loop edge is found, it is ignored without error.

 NOTE:  This routine supports digraphs.  For a directed arc (I -> W),
        an edge record is created in both vertices, I and W, and the
        edge record in I's adjacency list is marked OUTONLY while the
        edge record in W's list is marked INONLY.
        This makes it easy to used edge directedness when appropriate
        but also seamlessly process the corresponding undirected graph.

 Returns: OK, NOTOK on internal error, NONEMBEDDABLE if too many edges
 ********************************************************************/

int  _ReadAdjList(graphP theGraph, FILE *Infile)
{
int N, I, W, ErrorCode, adjList, J;

     if (Infile == NULL) return NOTOK;
     fgetc(Infile);                             /* Skip the N= */
     fgetc(Infile);
     fscanf(Infile, " %d ", &N);                /* Read N */
     if (gp_InitGraph(theGraph, N) != OK)
     {
    	  printf("Failed to init graph");
          return NOTOK;
     }

     // Clear the visited members of the vertices so they can be used
     // during the adjacency list read operation
     for (I=0; I < N; I++)
          theGraph->G[I].visited = 0;

     // Do the adjacency list read operation for each vertex in order
     for (I = 0, ErrorCode = OK; I < N && ErrorCode==OK; I++)
     {
          // Read the vertex number
          fscanf(Infile, "%d", &theGraph->G[I].v);

          // The vertices are expected to be in numeric ascending order
          if (theGraph->G[I].v != I)
        	  return NOTOK;

          // Skip the colon after the vertex number
          fgetc(Infile);

          // If the vertex already has a non-empty adjacency list, then it is
          // the result of adding edges during processing of preceding vertices.
          // The list is removed from the current vertex I and saved for use
          // during the read operation for I.  Adjacencies to preceding vertices
          // are pulled from this list, if present, or added as directed edges
          // if not.  Adjacencies to succeeding vertices are added as undirected
          // edges, and will be corrected later if the succeeding vertex does not
          // have the matching adjacency using the following mechanism.  After the
          // read operation for a vertex I, any adjacency nodes left in the saved
          // list are converted to directed edges from the preceding vertex to I.
          adjList = gp_GetFirstArc(theGraph, I);
          if (gp_IsArc(theGraph, adjList))
          {
        	  // Store the adjacency node location in the visited member of each
        	  // of the preceding vertices to which I is adjacent so that we can
        	  // efficiently detect the adjacency during the read operation and
        	  // efficiently find the adjacency node.
        	  J = gp_GetFirstArc(theGraph, I);
			  while (gp_IsArc(theGraph, J))
			  {
				  theGraph->G[theGraph->G[J].v].visited = J;
				  J = gp_GetNextArc(theGraph, J);
			  }

        	  // Make the adjacency list circular, for later ease of processing
			  gp_SetPrevArc(theGraph, adjList, gp_GetLastArc(theGraph, I));
			  gp_SetNextArc(theGraph, gp_GetLastArc(theGraph, I), adjList);

        	  // Remove the list from the vertex
			  gp_SetFirstArc(theGraph, I, gp_AdjacencyListEndMark(I));
			  gp_SetLastArc(theGraph, I, gp_AdjacencyListEndMark(I));
          }

          // Read the adjacency list.
          while (1)
          {
        	 // Read the next adjacent vertex, with NIL indicating the list end
             fscanf(Infile, " %d ", &W);
             if (W < 0) break;

             // Vertex numbers must be less than N
             if (W >= N)
                  ErrorCode = NOTOK;

             // Loop edges are not supported, but no reason to throw an error if they occur
             // If a loop occurs, we just do like the ostrich and ignore it
             else if (W == I)
            	 ErrorCode = OK;

             // If the adjacency is to a succeeding, higher numbered vertex,
             // then we'll add an undirected edge for now
             else if (I < W)
             {
             	 ErrorCode = gp_AddEdge(theGraph, I, 0, W, 0);
             }

             // If the adjacency is to a preceding, lower numbered vertex, then
             // we have to pull the adjacency node from the preexisting adjList,
             // if it is there, and if not then we have to add a directed edge.
             else
             {
            	 // If the adjacency node (arc) already exists, then we add it
            	 // as the new first arc of the vertex and delete it from adjList
            	 if (theGraph->G[W].visited)
            	 {
            		 J = theGraph->G[W].visited;

            		 // Remove the arc J from the adjList construct
            		 theGraph->G[W].visited = 0;
            		 if (adjList == J)
            		 {
            			 if ((adjList = gp_GetNextArc(theGraph, J)) == J)
            				 adjList = NIL;
            		 }
            		 gp_SetPrevArc(theGraph, gp_GetNextArc(theGraph, J), gp_GetPrevArc(theGraph, J));
            		 gp_SetNextArc(theGraph, gp_GetPrevArc(theGraph, J), gp_GetNextArc(theGraph, J));

            		 gp_AttachFirstArc(theGraph, I, J);
            	 }

            	 // If an adjacency node to the lower numbered vertex W does not
            	 // already exist, then we make a new directed arc from the current
            	 // vertex I to W.
            	 else
            	 {
            		 // It is added as the new first arc in both vertices
                	 ErrorCode = gp_AddEdge(theGraph, I, 0, W, 0);
                	 if (ErrorCode == OK)
                		 // Note that this call also sets OUTONLY on the twin arc
                		 gp_SetDirection(theGraph, gp_GetFirstArc(theGraph, W), EDGEFLAG_DIRECTION_INONLY);
            	 }
             }

             if (ErrorCode != OK) break;
          }

          // If there are still adjList entries after the read operation
          // then those entries are not representative of full undirected edges.
          // Rather, they represent incoming directed arcs from other vertices
          // into vertex I. They need to be added back into I's adjacency list but
          // marked as "INONLY", while the twin is marked "OUTONLY" (by the same function).
          while (gp_IsArc(theGraph, adjList))
          {
        	  J = adjList;

			  theGraph->G[theGraph->G[J].v].visited = 0;

 			  if ((adjList = gp_GetNextArc(theGraph, J)) == J)
 				  adjList = NIL;

     		  gp_SetPrevArc(theGraph, gp_GetNextArc(theGraph, J), gp_GetPrevArc(theGraph, J));
     		  gp_SetNextArc(theGraph, gp_GetPrevArc(theGraph, J), gp_GetNextArc(theGraph, J));

     		  gp_AttachFirstArc(theGraph, I, J);
     		  gp_SetDirection(theGraph, J, EDGEFLAG_DIRECTION_INONLY);
          }
     }

     return ErrorCode;
}

/********************************************************************
 _ReadLEDAGraph()
 Reads the edge list from a LEDA file containing a simple undirected graph.
 ********************************************************************/

int  _ReadLEDAGraph(graphP theGraph, FILE *Infile)
{
char Line[256];
int N, I, M, J, u, v;

    /* Skip the lines that say LEDA.GRAPH and give the node and edge types */
    fgets(Line, 255, Infile);
    fgets(Line, 255, Infile);
    fgets(Line, 255, Infile);

    /* Read the number of vertices, then skip that many more lines. */
    fgets(Line, 255, Infile);
    sscanf(Line, " %d", &N);
    for (I = 0; I < N; I++)
        fgets(Line, 255, Infile);

    /* Initialize the graph */
     if (gp_InitGraph(theGraph, N) != OK)
          return NOTOK;

    /* Read the number of edges */
    fgets(Line, 255, Infile);
    sscanf(Line, " %d", &M);

    /* Read and add each edge, omitting duplicates */
    for (J = 0; J < M; J++)
    {
        fgets(Line, 255, Infile);
        sscanf(Line, " %d %d", &u, &v);
        if (u != v && !gp_IsNeighbor(theGraph, u-1, v-1))
        {
             if (gp_AddEdge(theGraph, u-1, 0, v-1, 0) != OK)
                 return NOTOK;
        }
    }

    return OK;
}

/********************************************************************
 gp_Read()
 Opens the given file, determines whether it is in adjacency list or
 matrix format based on whether the file start with N or just a number,
 calls the appropriate read function, then closes the file and returns
 the graph.

 Digraphs and loop edges are not supported in the adjacency matrix format,
 which is upper triangular.

 In the adjacency list format, digraphs are supported.  Loop edges are
 ignored without producing an error.

 Pass "stdin" for the FileName to read from the stdin stream

 Returns: OK, NOTOK on internal error, NONEMBEDDABLE if too many edges
 ********************************************************************/

int gp_Read(graphP theGraph, char *FileName)
{
FILE *Infile;
char Ch;
int RetVal;

     if (strcmp(FileName, "stdin") == 0)
          Infile = stdin;
     else if ((Infile = fopen(FileName, READTEXT)) == NULL)
          return NOTOK;

     Ch = (char) fgetc(Infile);
     ungetc(Ch, Infile);
     if (Ch == 'N')
          RetVal = _ReadAdjList(theGraph, Infile);
     else if (Ch == 'L')
          RetVal = _ReadLEDAGraph(theGraph, Infile);
     else RetVal = _ReadAdjMatrix(theGraph, Infile);

     if (RetVal == OK)
     {
         void *extraData = NULL;
         long filePos = ftell(Infile);
         long fileSize;

         fseek(Infile, 0, SEEK_END);
         fileSize = ftell(Infile);
         fseek(Infile, filePos, SEEK_SET);

         if (filePos < fileSize)
         {
            extraData = malloc(fileSize - filePos + 1);
            fread(extraData, fileSize - filePos, 1, Infile);
         }
/*// Useful for quick debugging of IO extensibility
         if (extraData == NULL)
             printf("extraData == NULL\n");
         else
         {
             char *extraDataString = (char *) extraData;
             extraDataString[fileSize - filePos] = '\0';
             printf("extraData = '%s'\n", extraDataString);
         }
*/

         if (extraData != NULL)
         {
             RetVal = theGraph->functions.fpReadPostprocess(theGraph, extraData, fileSize - filePos);
             free((void *) extraData);
         }
     }

     if (strcmp(FileName, "stdin") != 0)
         fclose(Infile);

     return RetVal;
}

int  _ReadPostprocess(graphP theGraph, void *extraData, long extraDataSize)
{
     return OK;
}

/********************************************************************
 _WriteAdjList()
 For each vertex, we write its number, a colon, the list of adjacent vertices,
 then a NIL.  The vertices occupy the first N positions of theGraph.  Each
 vertex is also has indicators of the first and last adjacency nodes (arcs)
 in its adjacency list.

 Returns: NOTOK if either param is NULL; OK otherwise (after printing
                adjacency list representation to Outfile).
 ********************************************************************/

int  _WriteAdjList(graphP theGraph, FILE *Outfile)
{
int I, J;

     if (theGraph==NULL || Outfile==NULL) return NOTOK;

     fprintf(Outfile, "N=%d\n", theGraph->N);
     for (I=0; I < theGraph->N; I++)
     {
          fprintf(Outfile, "%d:", I);

          J = gp_GetLastArc(theGraph, I);
          while (gp_IsArc(theGraph, J))
          {
        	  if (!gp_GetDirection(theGraph, J, EDGEFLAG_DIRECTION_INONLY))
                  fprintf(Outfile, " %d", theGraph->G[J].v);

              J = gp_GetPrevArc(theGraph, J);
          }
          fprintf(Outfile, " %d\n", NIL);
     }
     return OK;
}

/********************************************************************
 _WriteAdjMatrix()
 Outputs upper triangular matrix representation capable of being
 read by _ReadAdjMatrix()

 Note: This routine does not support digraphs and will return an
       error if a directed edge is found.

 returns OK for success, NOTOK for failure
 ********************************************************************/

int  _WriteAdjMatrix(graphP theGraph, FILE *Outfile)
{
int  I, J, K;
char *Row = NULL;

     if (theGraph != NULL)
         Row = (char *) malloc((theGraph->N+1)*sizeof(char));

     if (Row==NULL || theGraph==NULL || Outfile==NULL)
     {
         if (Row != NULL) free(Row);
         return NOTOK;
     }

     fprintf(Outfile, "%d\n", theGraph->N);
     for (I = 0; I < theGraph->N; I++)
     {
          for (K = 0; K <= I; K++)
               Row[K] = ' ';
          for (K = I+1; K < theGraph->N; K++)
               Row[K] = '0';

          J = gp_GetFirstArc(theGraph, I);
          while (gp_IsArc(theGraph, J))
          {
        	  if (gp_GetDirection(theGraph, J, EDGEFLAG_DIRECTION_INONLY))
        		  return NOTOK;

              if (theGraph->G[J].v > I)
                  Row[theGraph->G[J].v] = '1';

              J = gp_GetNextArc(theGraph, J);
          }

          Row[theGraph->N] = '\0';
          fprintf(Outfile, "%s\n", Row);
     }

     free(Row);
     return OK;
}

/********************************************************************
 _WriteDebugInfo()
 Writes adjacency list, but also includes the type value of each
 edge (e.g. is it DFS child  arc, forward arc or back arc?), and
 the L, A and DFSParent of each vertex.
 ********************************************************************/

int  _WriteDebugInfo(graphP theGraph, FILE *Outfile)
{
int I, J, Gsize;

     if (theGraph==NULL || Outfile==NULL) return NOTOK;

     /* Print parent copy vertices and their adjacency lists */

     fprintf(Outfile, "DEBUG N=%d M=%d\n", theGraph->N, theGraph->M);
     for (I=0; I < theGraph->N; I++)
     {
          fprintf(Outfile, "%d(P=%d,lA=%d,LowPt=%d,v=%d):",
                             I, theGraph->V[I].DFSParent,
                                theGraph->V[I].leastAncestor,
                                theGraph->V[I].Lowpoint,
                                theGraph->G[I].v);

          J = gp_GetFirstArc(theGraph, I);
          while (gp_IsArc(theGraph, J))
          {
              fprintf(Outfile, " %d(J=%d)", theGraph->G[J].v, J);
              J = gp_GetNextArc(theGraph, J);
          }

          fprintf(Outfile, " %d\n", NIL);
     }

     /* Print any root copy vertices and their adjacency lists */

     for (I = theGraph->N; I < 2*theGraph->N; I++)
     {
          if (theGraph->G[I].v == NIL)
              continue;

          fprintf(Outfile, "%d(copy of=%d, DFS child=%d):",
                           I, theGraph->G[I].v, I-theGraph->N);

          J = gp_GetFirstArc(theGraph, I);
          while (gp_IsArc(theGraph, J))
          {
              fprintf(Outfile, " %d(J=%d)", theGraph->G[J].v, J);
              J = gp_GetNextArc(theGraph, J);
          }

          fprintf(Outfile, " %d\n", NIL);
     }

     /* Print information about vertices and root copy (virtual) vertices */
     fprintf(Outfile, "\nVERTEX INFORMATION\n");
     for (I=0; I < 2*theGraph->N; I++)
     {
         if (theGraph->G[I].v == NIL)
             continue;

         fprintf(Outfile, "V[%3d] v=%3d, type=%c, first arc=%3d, last arc=%3d\n",
                          I,
                          theGraph->G[I].v,
                          theGraph->G[I].type,
                          gp_GetFirstArc(theGraph, I),
                          gp_GetLastArc(theGraph, I));
     }

     /* Print information about edges */

     fprintf(Outfile, "\nEDGE INFORMATION\n");
     Gsize = theGraph->edgeOffset + theGraph->arcCapacity;
     for (J=theGraph->edgeOffset; J < Gsize; J++)
     {
          if (theGraph->G[J].v == NIL)
              continue;

          fprintf(Outfile, "E[%3d] v=%3d, type=%c, next arc=%3d, prev arc=%3d\n",
                           J,
                           theGraph->G[J].v,
                           theGraph->G[J].type,
                           gp_GetNextArc(theGraph, J),
                           gp_GetPrevArc(theGraph, J));
     }

     return OK;
}

/********************************************************************
 gp_Write()
 Writes theGraph into the file.
 Pass "stdout" or "stderr" to FileName to write to the corresponding stream
 Pass WRITE_ADJLIST, WRITE_ADJMATRIX or WRITE_DEBUGINFO for the Mode

 NOTE: For digraphs, it is an error to use a mode other than WRITE_ADJLIST

 Returns NOTOK on error, OK on success.
 ********************************************************************/

int  gp_Write(graphP theGraph, char *FileName, int Mode)
{
FILE *Outfile;
int RetVal;

     if (theGraph == NULL || FileName == NULL)
    	 return NOTOK;

     if (strcmp(FileName, "nullwrite") == 0)
    	  return OK;

     if (strcmp(FileName, "stdout") == 0)
          Outfile = stdout;
     else if (strcmp(FileName, "stderr") == 0)
          Outfile = stderr;
     else if ((Outfile = fopen(FileName, WRITETEXT)) == NULL)
          return NOTOK;

     switch (Mode)
     {
         case WRITE_ADJLIST   :
        	 RetVal = _WriteAdjList(theGraph, Outfile);
             break;
         case WRITE_ADJMATRIX :
        	 RetVal = _WriteAdjMatrix(theGraph, Outfile);
             break;
         case WRITE_DEBUGINFO :
        	 RetVal = _WriteDebugInfo(theGraph, Outfile);
             break;
         default :
        	 RetVal = NOTOK;
        	 break;
     }

     if (RetVal == OK)
     {
         void *extraData = NULL;
         long extraDataSize;

         RetVal = theGraph->functions.fpWritePostprocess(theGraph, &extraData, &extraDataSize);

         if (extraData != NULL)
         {
             if (!fwrite(extraData, extraDataSize, 1, Outfile))
                 RetVal = NOTOK;
             free(extraData);
         }
     }

     if (strcmp(FileName, "stdout") == 0 || strcmp(FileName, "stderr") == 0)
         fflush(Outfile);

     else if (fclose(Outfile) != 0)
         RetVal = NOTOK;

     return RetVal;
}

/********************************************************************
 _WritePostprocess()

 By default, no additional information is written.
 ********************************************************************/

int  _WritePostprocess(graphP theGraph, void **pExtraData, long *pExtraDataSize)
{
     return OK;
}

/********************************************************************
 _Log()

 When the project is compiled with LOGGING enabled, this method writes
 a string to the file PLANARITY.LOG in the current working directory.
 On first write, the file is created or cleared.
 Call this method with NULL to close the log file.
 ********************************************************************/

void _Log(char *Str)
{
static FILE *logfile = NULL;

    if (logfile == NULL)
    {
        if ((logfile = fopen("PLANARITY.LOG", WRITETEXT)) == NULL)
        	return;
    }

    if (Str != NULL)
    {
        fprintf(logfile, "%s", Str);
        fflush(logfile);
    }
    else
        fclose(logfile);
}

void _LogLine(char *Str)
{
	_Log(Str);
	_Log("\n");
}

static char LogStr[512];

char *_MakeLogStr1(char *format, int one)
{
	sprintf(LogStr, format, one);
	return LogStr;
}

char *_MakeLogStr2(char *format, int one, int two)
{
	sprintf(LogStr, format, one, two);
	return LogStr;
}

char *_MakeLogStr3(char *format, int one, int two, int three)
{
	sprintf(LogStr, format, one, two, three);
	return LogStr;
}

char *_MakeLogStr4(char *format, int one, int two, int three, int four)
{
	sprintf(LogStr, format, one, two, three, four);
	return LogStr;
}

char *_MakeLogStr5(char *format, int one, int two, int three, int four, int five)
{
	sprintf(LogStr, format, one, two, three, four, five);
	return LogStr;
}
