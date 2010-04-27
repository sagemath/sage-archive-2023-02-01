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

/****************************************************************************
 SpecificGraph()
 ****************************************************************************/

int SpecificGraph(char command, char *infileName, char *outfileName, char *outfile2Name)
{
graphP theGraph, origGraph;
platform_time start, end;
int Result;

    // Get the filename of the graph to test
    if ((infileName = ConstructInputFilename(infileName)) == NULL)
	    return NOTOK;

    // Create the graph and, if needed, attach the correct algorithm to it
    theGraph = gp_New();

	switch (command)
	{
		case 'd' : gp_AttachDrawPlanar(theGraph); break;
		case '2' : gp_AttachK23Search(theGraph); break;
		case '3' : gp_AttachK33Search(theGraph); break;
		case '4' : gp_AttachK4Search(theGraph); break;
		case 'c' : gp_AttachColorVertices(theGraph); break;
	}

    // Read the graph into memory
	Result = gp_Read(theGraph, infileName);
	if (Result == NONEMBEDDABLE)
	{
		Message("The graph contains too many edges.\n");
		// Some of the algorithms will still run correctly with some edges removed.
		if (strchr("pdo234", command))
		{
			Message("Some edges were removed, but the algorithm will still run correctly.\n");
			Result = OK;
		}
	}

	// If there was an unrecoverable error, report it
	if (Result != OK)
		ErrorMessage("Failed to read graph\n");

	// Otherwise, call the correct algorithm on it
	else
	{
		// Copy the graph for integrity checking
        origGraph = gp_DupGraph(theGraph);

        // Run the algorithm
        if (strchr("pdo234", command))
        {
    		int embedFlags = GetEmbedFlags(command);
	        platform_GetTime(start);
			Result = gp_Embed(theGraph, embedFlags);
	        platform_GetTime(end);
	        Result = gp_TestEmbedResultIntegrity(theGraph, origGraph, Result);
        }
        else
        {
	        platform_GetTime(start);
        	if (command == 'c')
        	{
    			if ((Result = gp_ColorVertices(theGraph)) == OK)
    				 Result = gp_ColorVerticesIntegrityCheck(theGraph, origGraph);
        	}
        	else
    			Result = NOTOK;
   	        platform_GetTime(end);
        }

        // Write what the algorithm determined and how long it took
        WriteAlgorithmResults(theGraph, Result, command, start, end, infileName);

        // Free the graph obtained for integrity checking.
        gp_Free(&origGraph);
	}

	// Report an error, if there was one, free the graph, and return
	if (Result != OK && Result != NONEMBEDDABLE)
	{
		ErrorMessage("AN ERROR HAS BEEN DETECTED\n");
		Result = NOTOK;
	}

	// Provide the output file(s)
	else
	{
        // Restore the vertex ordering of the original graph (undo DFS numbering)
        if (strchr("pdo234", command))
            gp_SortVertices(theGraph);

        // Determine the name of the primary output file
        outfileName = ConstructPrimaryOutputFilename(infileName, outfileName, command);

        // For some algorithms, the primary output file is not always written
        if ((strchr("pdo", command) && Result == NONEMBEDDABLE) ||
        	(strchr("234", command) && Result == OK))
        	;

        // Write the primary output file, if appropriate to do so
        else
			gp_Write(theGraph, outfileName, WRITE_ADJLIST);

        // NOW WE WANT TO WRITE THE SECONDARY OUTPUT FILE

		// When called from the menu system, we want to write the planar or outerplanar
		// obstruction, if one exists. For planar graph drawing, we want the character
        // art rendition.  An empty but non-NULL string is passed to indicate the necessity
        // of selecting a default name for the second output file.
		if (outfile2Name != NULL)
		{
			if ((command == 'p' || command == 'o') && Result == NONEMBEDDABLE)
			{
				// By default, use the same name as the primary output filename
				if (strlen(outfile2Name) == 0)
				    outfile2Name = outfileName;
				gp_Write(theGraph, outfile2Name, WRITE_ADJLIST);
			}
			else if (command == 'd' && Result == OK)
			{
				// By default, add ".render.txt" to the primary output filename
				if (strlen(outfile2Name) == 0)
   				    strcat((outfile2Name = outfileName), ".render.txt");
				gp_DrawPlanar_RenderToFile(theGraph, outfile2Name);
			}
		}
	}

	// Free the graph
	gp_Free(&theGraph);

	// Flush any remaining message content to the user, and return the result
    FlushConsole(stdout);
	return Result;
}

/****************************************************************************
 WriteAlgorithmResults()
 ****************************************************************************/

void WriteAlgorithmResults(graphP theGraph, int Result, char command, platform_time start, platform_time end, char *infileName)
{
	if (infileName)
		 sprintf(Line, "The graph '%s' ", infileName);
	else sprintf(Line, "The graph ");
	Message(Line);

	switch (command)
	{
		case 'p' : sprintf(Line, "is%s planar.\n", Result==OK ? "" : " not"); break;
		case 'd' : sprintf(Line, "is%s planar.\n", Result==OK ? "" : " not"); break;
		case 'o' : sprintf(Line, "is%s outerplanar.\n", Result==OK ? "" : " not"); break;
		case '2' : sprintf(Line, "has %s subgraph homeomorphic to K_{2,3}.\n", Result==OK ? "no" : "a"); break;
		case '3' : sprintf(Line, "has %s subgraph homeomorphic to K_{3,3}.\n", Result==OK ? "no" : "a"); break;
		case '4' : sprintf(Line, "has %s subgraph homeomorphic to K_4.\n", Result==OK ? "no" : "a"); break;
		case 'c' : sprintf(Line, "has been %d-colored.\n", gp_GetNumColorsUsed(theGraph)); break;
		default  : sprintf(Line, "nas not been processed due to unrecognized command.\n"); break;
	}
	Message(Line);

	sprintf(Line, "Algorithm '%s' executed in %.3lf seconds.\n",
			GetAlgorithmName(command), platform_GetDuration(start,end));
	Message(Line);
}
