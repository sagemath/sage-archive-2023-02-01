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

#define GRAPHIO_C

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
 Returns: OK, NOTOK on internal error, NONPLANAR if too many edges
 ********************************************************************/

int _ReadAdjMatrix(graphP theGraph, FILE *Infile)
{
int N, I, J, Flag, ErrorCode;

     if (Infile == NULL) return NOTOK;
     fscanf(Infile, " %d ", &N);
     if (gp_InitGraph(theGraph, N) != OK)
          return NOTOK;

     for (I = 0, ErrorCode = OK; I < N-1 && ErrorCode==OK; I++)
     {
          theGraph->G[I].v = I;
          for (J = I+1; J < N; J++)
          {
               fscanf(Infile, " %1d", &Flag);
               if (Flag)
               {
                   ErrorCode = gp_AddEdge(theGraph, I, 0, J, 0);
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

 NOTE:  The vertex number is skipped; the vertices are expected to be
        in sorted order in the file.  The vertex number is for file
        documentation only.

 NOTE:  This reader will not read all edges in the same order as
        they exist in your input file.  When a low numbered vertex u
        is being read, an edge to a higher numbered vertex v causes
        this reader to add both (u,v) and (v,u) to the graph.  When
        v is processed, the edge to u is simply ignored.  This ensures
        that the graph from this reader is undirected, but it is an
        error to try feeding a digraph to this reader because arcs
        from high to low numbered vertices are ignored even if the
        arc from low to high numbered vertex was absent.

 Returns: OK, NOTOK on internal error, NONPLANAR if too many edges
 ********************************************************************/

int  _ReadAdjList(graphP theGraph, FILE *Infile)
{
int N, I, J, ErrorCode;

     if (Infile == NULL) return NOTOK;
     fgetc(Infile);                             /* Skip the N= */
     fgetc(Infile);
     fscanf(Infile, " %d ", &N);                /* Read N */
     if (gp_InitGraph(theGraph, N) != OK)
          return NOTOK;

     for (I = 0, ErrorCode = OK; I < N && ErrorCode==OK; I++)
     {
          theGraph->G[I].v = I;
          fscanf(Infile, "%*d");                /* Skip vertex # and colon */
          fgetc(Infile);
          while (1)                             /* Read Adj List */
          {
             fscanf(Infile, " %d ", &J);

             if (J < 0) break;                  /* If NIL, then we are done */
             if (J >= N)
                  ErrorCode = NOTOK;
             else if (I > J)
                  ErrorCode = OK;
             else ErrorCode = gp_AddEdge(theGraph, I, 0, J, 0);

             if (ErrorCode != OK) break;
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
 Pass "stdin" for the FileName to read from the stdin stream
 Returns: OK, NOTOK on internal error, NONPLANAR if too many edges
 ********************************************************************/

int gp_Read(graphP theGraph, char *FileName)
{
FILE *Infile;
char Ch;
int RetVal = NOTOK;

     if (strcmp(FileName, "stdin")==OK)
          Infile = stdin;
     else if ((Infile = fopen(FileName, READTEXT)) == NULL)
          return RetVal;

     Ch = (char) fgetc(Infile);
     ungetc(Ch, Infile);
     if (Ch == 'N')
          RetVal = _ReadAdjList(theGraph, Infile);
     else if (Ch == 'L')
          RetVal = _ReadLEDAGraph(theGraph, Infile);
     else RetVal = _ReadAdjMatrix(theGraph, Infile);

     if (strcmp(FileName, "stdin") != OK)
         fclose(Infile);

     return RetVal;
}

/********************************************************************
 _WriteAdjList()
 For each vertex, we write its number, a colon, the list of adjacent vertices,
 then a NIL.  The vertices occupy the first N positions of theGraph.  Each
 vertex is also the head of a circular list kept by link[0] and link[1].
 Aside from the vertex itself, the other elements in the list represent
 the edges between the vertex and its neighbors, and these edge records
 reside at or above position 2N in theGraph.

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

          J = theGraph->G[I].link[1];
          while (J >= theGraph->N)
          {
              fprintf(Outfile, " %d", theGraph->G[J].v);
              J = theGraph->G[J].link[1];
          }
          fprintf(Outfile, " %d\n", NIL);
     }
     return OK;
}

/********************************************************************
 _WriteAdjMatrix()
 Outputs upper triangular matrix representation capable of being
 read by _ReadAdjMatrix()
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

          J = theGraph->G[I].link[0];
          while (J >= theGraph->N)
          {
              if (theGraph->G[J].v > I)
                  Row[theGraph->G[J].v] = '1';

              J = theGraph->G[J].link[0];
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
int I, J;

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

          J = theGraph->G[I].link[0];
          while (J >= theGraph->N)
          {
              fprintf(Outfile, " %d(J=%d)", theGraph->G[J].v, J);
              J = theGraph->G[J].link[0];
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

          J = theGraph->G[I].link[0];
          while (J >= 2*theGraph->N)
          {
              fprintf(Outfile, " %d(J=%d)", theGraph->G[J].v, J);
              J = theGraph->G[J].link[0];
          }

          fprintf(Outfile, " %d\n", NIL);
     }

     /* Print all graph node information for vertices (0 to N-1),
        root copy vertices (N to 2N-1), and edges (2N to 8N-1) */

     fprintf(Outfile, "\nGRAPH NODES\n", NIL);
     for (I=0; I < 8*theGraph->N; I++)
     {
          if (theGraph->G[I].v == NIL)
              continue;

          fprintf(Outfile, "G[%3d] v=%3d, type=%c, link[0]=%3d, link[1]=%3d\n",
                           I,
                           theGraph->G[I].v,
                           theGraph->G[I].type,
                           theGraph->G[I].link[0],
                           theGraph->G[I].link[1]);
     }

     return OK;
}

/********************************************************************
 gp_Write()
 Writes theGraph into the file.
 Pass "stdout" or "stderr" to FileName to write to the corresponding stream
 Pass WRITE_ADJLIST, WRITE_ADJMATRIX or WRITE_DEBUGINFO for the Mode
 Returns NOTOK on parameter error, OK on success.
 ********************************************************************/

int  gp_Write(graphP theGraph, char *FileName, int Mode)
{
FILE *Outfile;
int RetVal;

     if (theGraph == NULL || FileName == NULL) return NOTOK;

     if (strcmp(FileName, "stdout") == OK)
          Outfile = stdout;
     else if (strcmp(FileName, "stderr") == OK)
          Outfile = stderr;
     else if ((Outfile = fopen(FileName, WRITETEXT)) == NULL)
          return NOTOK;

     switch (Mode)
     {
         case WRITE_ADJLIST   : RetVal = _WriteAdjList(theGraph, Outfile);
                                break;
         case WRITE_ADJMATRIX : RetVal = _WriteAdjMatrix(theGraph, Outfile);
                                break;
         case WRITE_DEBUGINFO : RetVal = _WriteDebugInfo(theGraph, Outfile);
                                break;
     }

     if (strcmp(FileName, "stdout")==OK || strcmp(FileName, "stderr")==OK)
         fflush(Outfile);

     else if (fclose(Outfile) != OK)
         RetVal = NOTOK;

     return RetVal;
}
