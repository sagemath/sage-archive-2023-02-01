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

#include "planarity.h"

void GetNumberIfZero(int *pNum, char *prompt, int min, int max);
void ReinitializeGraph(graphP *pGraph, int ReuseGraphs, char command);
graphP MakeGraph(int Size, char command);

/****************************************************************************
 RandomGraphs()
 Top-level method to randomly generate graphs to test the algorithm given by
 the command parameter.
 The number of graphs to generate, and the number of vertices for each graph,
 can be sent as the second and third params.  For each that is sent as zero,
 this method will prompt the user for a value.
 ****************************************************************************/

#define NUM_MINORS  9

int  RandomGraphs(char command, int NumGraphs, int SizeOfGraphs)
{
char theFileName[256];
int  I, countUpdateFreq;
int Result=OK, MainStatistic=0;
int  ObstructionMinorFreqs[NUM_MINORS];
graphP theGraph=NULL, origGraph=NULL;
platform_time start, end;
int embedFlags = GetEmbedFlags(command);
int ReuseGraphs = TRUE;

     GetNumberIfZero(&NumGraphs, "Enter number of graphs to generate:", 1, 1000000000);
     GetNumberIfZero(&SizeOfGraphs, "Enter size of graphs:", 1, 10000);

   	 theGraph = MakeGraph(SizeOfGraphs, command);
   	 origGraph = MakeGraph(SizeOfGraphs, command);
   	 if (theGraph == NULL || origGraph == NULL)
   	 {
   		 gp_Free(&theGraph);
   		 return NOTOK;
   	 }

     // Initialize a secondary statistics array
     for (I=0; I<NUM_MINORS; I++)
          ObstructionMinorFreqs[I] = 0;

   	 // Seed the random number generator with "now". Do it after any prompting
   	 // to tie randomness to human process of answering the prompt.
   	 srand(time(NULL));

   	 // Select a counter update frequency that updates more frequently with larger graphs
   	 // and which is relatively prime with 10 so that all digits of the count will change
   	 // even though we aren't showing the count value on every iteration
   	 countUpdateFreq = 3579 / SizeOfGraphs;
   	 countUpdateFreq = countUpdateFreq < 1 ? 1 : countUpdateFreq;
   	 countUpdateFreq = countUpdateFreq % 2 == 0 ? countUpdateFreq+1 : countUpdateFreq;
   	 countUpdateFreq = countUpdateFreq % 5 == 0 ? countUpdateFreq+2 : countUpdateFreq;

   	 // Start the count
     fprintf(stdout, "0\r");
     fflush(stdout);

     // Start the timer
     platform_GetTime(start);

     // Generate and process the number of graphs requested
     for (I=0; I < NumGraphs; I++)
     {
          if ((Result = gp_CreateRandomGraph(theGraph)) == OK)
          {
              if (tolower(OrigOut)=='y')
              {
                  sprintf(theFileName, "random\\%d.txt", I%10);
                  gp_Write(theGraph, theFileName, WRITE_ADJLIST);
              }

              gp_CopyGraph(origGraph, theGraph);

              if (strchr("pdo234", command))
              {
                  Result = gp_Embed(theGraph, embedFlags);

                  if (gp_TestEmbedResultIntegrity(theGraph, origGraph, Result) != Result)
                      Result = NOTOK;

                  if (Result == OK)
                  {
                       MainStatistic++;

                       if (tolower(EmbeddableOut) == 'y')
                       {
                           sprintf(theFileName, "embedded\\%d.txt", I%10);
                           gp_Write(theGraph, theFileName, WRITE_ADJMATRIX);
                       }

                       if (tolower(AdjListsForEmbeddingsOut) == 'y')
                       {
                           sprintf(theFileName, "adjlist\\%d.txt", I%10);
                           gp_Write(theGraph, theFileName, WRITE_ADJLIST);
                       }
                  }
                  else if (Result == NONEMBEDDABLE)
                  {
                       if (embedFlags == EMBEDFLAGS_PLANAR || embedFlags == EMBEDFLAGS_OUTERPLANAR)
                       {
                           if (theGraph->IC.minorType & MINORTYPE_A)
                                ObstructionMinorFreqs[0] ++;
                           else if (theGraph->IC.minorType & MINORTYPE_B)
                                ObstructionMinorFreqs[1] ++;
                           else if (theGraph->IC.minorType & MINORTYPE_C)
                                ObstructionMinorFreqs[2] ++;
                           else if (theGraph->IC.minorType & MINORTYPE_D)
                                ObstructionMinorFreqs[3] ++;
                           else if (theGraph->IC.minorType & MINORTYPE_E)
                                ObstructionMinorFreqs[4] ++;

                           if (theGraph->IC.minorType & MINORTYPE_E1)
                                ObstructionMinorFreqs[5] ++;
                           else if (theGraph->IC.minorType & MINORTYPE_E2)
                                ObstructionMinorFreqs[6] ++;
                           else if (theGraph->IC.minorType & MINORTYPE_E3)
                                ObstructionMinorFreqs[7] ++;
                           else if (theGraph->IC.minorType & MINORTYPE_E4)
                                ObstructionMinorFreqs[8] ++;

                           if (tolower(ObstructedOut) == 'y')
                           {
                               sprintf(theFileName, "obstructed\\%d.txt", I%10);
                               gp_Write(theGraph, theFileName, WRITE_ADJMATRIX);
                           }
                       }
                  }
              }
              else if (command == 'c')
              {
      			if ((Result = gp_ColorVertices(theGraph)) == OK)
      				 Result = gp_ColorVerticesIntegrityCheck(theGraph, origGraph);
				if (Result == OK && gp_GetNumColorsUsed(theGraph) <= 5)
					MainStatistic++;
              }

              // If there is an error in processing, then write the file for debugging
              if (Result != OK && Result != NONEMBEDDABLE)
              {
                   sprintf(theFileName, "error\\%d.txt", I%10);
                   gp_Write(origGraph, theFileName, WRITE_ADJLIST);
              }
          }

          // Reinitialize or recreate graphs for next iteration
          ReinitializeGraph(&theGraph, ReuseGraphs, command);
          ReinitializeGraph(&origGraph, ReuseGraphs, command);

          // Show progress, but not so often that it bogs down progress
          if (quietMode == 'n' && (I+1) % countUpdateFreq == 0)
          {
              fprintf(stdout, "%d\r", I+1);
              fflush(stdout);
          }

          // Terminate loop on error
          if (Result != OK && Result != NONEMBEDDABLE)
          {
        	  ErrorMessage("\nError found\n");
              Result = NOTOK;
              break;
          }
     }

     // Stop the timer
     platform_GetTime(end);

     // Finish the count
     fprintf(stdout, "%d\n", NumGraphs);
     fflush(stdout);

     // Free the graph structures created before the loop
     gp_Free(&theGraph);
     gp_Free(&origGraph);

     // Print some demographic results
     if (Result == OK || Result == NONEMBEDDABLE)
         Message("\nNo Errors Found.");
     sprintf(Line, "\nDone (%.3lf seconds).\n", platform_GetDuration(start,end));
     Message(Line);

     // Report statistics for planar or outerplanar embedding
     if (embedFlags == EMBEDFLAGS_PLANAR || embedFlags == EMBEDFLAGS_OUTERPLANAR)
     {
         sprintf(Line, "Num Embedded=%d.\n", MainStatistic);
         Message(Line);

         for (I=0; I<5; I++)
         {
        	  // Outerplanarity does not produces minors C and D
        	  if (embedFlags == EMBEDFLAGS_OUTERPLANAR && (I==2 || I==3))
        		  continue;

              sprintf(Line, "Minor %c = %d\n", I+'A', ObstructionMinorFreqs[I]);
              Message(Line);
         }

         if (!(embedFlags & ~EMBEDFLAGS_PLANAR))
         {
             sprintf(Line, "\nNote: E1 are added to C, E2 are added to A, and E=E3+E4+K5 homeomorphs.\n");
             Message(Line);

             for (I=5; I<NUM_MINORS; I++)
             {
                  sprintf(Line, "Minor E%d = %d\n", I-4, ObstructionMinorFreqs[I]);
                  Message(Line);
             }
         }
     }

     // Report statistics for graph drawing
     else if (embedFlags == EMBEDFLAGS_DRAWPLANAR)
     {
         sprintf(Line, "Num Graphs Embedded and Drawn=%d.\n", MainStatistic);
         Message(Line);
     }

     // Report statistics for subgraph homeomorphism algorithms
     else if (embedFlags == EMBEDFLAGS_SEARCHFORK23)
     {
         sprintf(Line, "Of the generated graphs, %d did not contain a K_{2,3} homeomorph as a subgraph.\n", MainStatistic);
         Message(Line);
     }
     else if (embedFlags == EMBEDFLAGS_SEARCHFORK33)
     {
         sprintf(Line, "Of the generated graphs, %d did not contain a K_{3,3} homeomorph as a subgraph.\n", MainStatistic);
         Message(Line);
     }
     else if (embedFlags == EMBEDFLAGS_SEARCHFORK4)
     {
         sprintf(Line, "Of the generated graphs, %d did not contain a K_4 homeomorph as a subgraph.\n", MainStatistic);
         Message(Line);
     }

     // Report statistics for vertex coloring
     else if (command == 'c')
     {
         sprintf(Line, "Num Graphs colored with 5 or fewer colors=%d.\n", MainStatistic);
         Message(Line);
     }

     FlushConsole(stdout);

     return Result==OK || Result==NONEMBEDDABLE ? OK : NOTOK;
}

/****************************************************************************
 GetNumberIfZero()
 Internal function that gets a number if the given *pNum is zero.
 The prompt is displayed if the number must be obtained from the user.
 Whether the given number is used or obtained from the user, the function
 ensures it is in the range [min, max] and assigns the midpoint value if
 it is not.
 ****************************************************************************/

void GetNumberIfZero(int *pNum, char *prompt, int min, int max)
{
	if (*pNum == 0)
	{
	    Prompt(prompt);
	    scanf(" %d", pNum);
	}

	if (min < 1) min = 1;
	if (max < min) max = min;

	if (*pNum < min || *pNum > max)
	{
		*pNum = (max + min) / 2;
        sprintf(Line, "Number out of range [%d, %d]; changed to %d\n", min, max, *pNum);
        ErrorMessage(Line);
	}
}

/****************************************************************************
 MakeGraph()
 Internal function that makes a new graph, initializes it, and attaches an
 algorithm to it based on the command.
 ****************************************************************************/

graphP MakeGraph(int Size, char command)
{
	graphP theGraph;
    if ((theGraph = gp_New()) == NULL || gp_InitGraph(theGraph, Size) != OK)
    {
    	ErrorMessage("Error creating space for a graph of the given size.\n");
    	gp_Free(&theGraph);
    	return NULL;
    }

// Enable the appropriate feature. Although the same code appears in SpecificGraph,
// it is deliberately not separated to a common utility because SpecificGraph is
// used as a self-contained tutorial.  It is not that hard to update both locations
// when new algorithms are added.

	switch (command)
	{
		case 'd' : gp_AttachDrawPlanar(theGraph); break;
		case '2' : gp_AttachK23Search(theGraph); break;
		case '3' : gp_AttachK33Search(theGraph); break;
		case '4' : gp_AttachK4Search(theGraph); break;
		case 'c' : gp_AttachColorVertices(theGraph); break;
	}

	return theGraph;
}

/****************************************************************************
 ReinitializeGraph()
 Internal function that will either reinitialize the given graph or free it
 and make a new one just like it.
 ****************************************************************************/

void ReinitializeGraph(graphP *pGraph, int ReuseGraphs, char command)
{
	if (ReuseGraphs)
		gp_ReinitializeGraph(*pGraph);
	else
	{
		graphP newGraph = MakeGraph((*pGraph)->N, command);
		gp_Free(pGraph);
		*pGraph = newGraph;
	}
}

/****************************************************************************
 Creates a random maximal planar graph, then adds 'extraEdges' edges to it.
 ****************************************************************************/

int RandomGraph(char command, int extraEdges, int numVertices, char *outfileName, char *outfile2Name)
{
int  Result;
platform_time start, end;
graphP theGraph=NULL, origGraph;
int embedFlags = GetEmbedFlags(command);
char saveEdgeListFormat;

     GetNumberIfZero(&numVertices, "Enter number of vertices:", 1, 1000000);
     if ((theGraph = MakeGraph(numVertices, command)) == NULL)
    	 return NOTOK;

     srand(time(NULL));

     Message("Creating the random graph...\n");
     platform_GetTime(start);
     if (gp_CreateRandomGraphEx(theGraph, 3*numVertices-6+extraEdges) != OK)
     {
         ErrorMessage("gp_CreateRandomGraphEx() failed\n");
         return NOTOK;
     }
     platform_GetTime(end);

     sprintf(Line, "Created random graph with %d edges in %.3lf seconds. ", theGraph->M, platform_GetDuration(start,end));
     Message(Line);
     FlushConsole(stdout);

     // The user may have requested a copy of the random graph before processing
     if (outfile2Name != NULL)
     {
         gp_Write(theGraph, outfile2Name, WRITE_ADJLIST);
     }

     origGraph = gp_DupGraph(theGraph);

     // Do the requested algorithm on the randomly generated graph
     Message("Now processing\n");
     FlushConsole(stdout);

     if (strchr("pdo234", command))
     {
         platform_GetTime(start);
         Result = gp_Embed(theGraph, embedFlags);
         platform_GetTime(end);

    	 gp_SortVertices(theGraph);

         if (gp_TestEmbedResultIntegrity(theGraph, origGraph, Result) != Result)
             Result = NOTOK;
     }
     else if (command == 'c')
     {
         platform_GetTime(start);
    	 Result = gp_ColorVertices(theGraph);
         platform_GetTime(end);
     }
     else
    	 Result = NOTOK;

     // Write what the algorithm determined and how long it took
     WriteAlgorithmResults(theGraph, Result, command, start, end, NULL);

     // On successful algorithm result, write the output file and see if the
     // user wants the edge list formatted file.
     if (Result == OK || Result == NONEMBEDDABLE)
     {
    	 if (outfileName != NULL)
    		 gp_Write(theGraph, outfileName, WRITE_ADJLIST);

         Prompt("Do you want to save the generated graph in edge list format (y/n)? ");
         fflush(stdin);
         scanf(" %c", &saveEdgeListFormat);
         if (tolower(saveEdgeListFormat) == 'y')
         {
        	 char *fileName = "maxPlanarEdgeList.txt";
             if (extraEdges > 0)
            	 fileName = "nonPlanarEdgeList.txt";

             SaveAsciiGraph(theGraph, fileName);
             sprintf(Line, "Edge list format saved to '%s'\n", fileName);
        	 Message(Line);
         }
     }
     else ErrorMessage("Failure occurred");

     gp_Free(&theGraph);
     gp_Free(&origGraph);

     FlushConsole(stdout);
     return Result;
}
