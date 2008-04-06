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

#include <string.h>
#include <stdlib.h>
#include <time.h>
#include <ctype.h>
#include "graph.h"

/****************************************************************************
 MAIN
 ****************************************************************************/

int main(int argc, char *argv[])
{
    if (argc >= 3)
    {
    graphP theGraph = gp_New();
    int Result;

        Result = gp_Read(theGraph, argv[1]);
        if (Result == NONPLANAR)
            return 0;
        if (Result != OK)
        {
            fprintf(stderr, "Failed to read graph %s\n", argv[1]);
            return NOTOK;
        }

        Result = gp_Embed(theGraph, EMBEDFLAGS_PLANAR);

        if (Result == OK)
        {
            gp_SortVertices(theGraph);
            gp_Write(theGraph, argv[2], WRITE_ADJLIST);
        }

        else if (Result == NONPLANAR)
        {
            if (argc >= 5 && strcmp(argv[3], "-n")==0)
            {
                gp_SortVertices(theGraph);
                gp_Write(theGraph, argv[4], WRITE_ADJLIST);
            }

            Result = OK;
        }
        else
            Result = NOTOK;

        gp_Free(&theGraph);

        return Result;
    }
    else
    {
        printf("Planarity 1.0\n");
        printf("Copyright (c) 2005 by John M. Boyer\n\n");

        printf("This program is provided to you as-is with no warranty.\n");
        printf("You are licensed to use this program for any purpose\n");
        printf("provided the copyright message in License.txt appears\n");
        printf("in the acknowledgements of derivative works that include\n");
        printf("this program or its parts.\n\n");

        printf("Send feedback to jboyer@acm.org\n");

        printf("Usage: planarity input.txt embedding.txt [-n kuratowskiSubgraph.txt]\n\n");

        printf("The input graph can be in an adjacency list format,\n");
        printf("adjacency matrix format, or a simple LEDA graph.\n");
        printf("The resulting graph, a combinatorial planar embedding or\n");
        printf("a Kuratowski subgraph, is in the adjacency list format.\n\n");

        printf("Adjacency list format:\n");
        printf("N=5\n");
        printf("0: 2 1 4 3 -1\n");
        printf("1: 2 4 0 3 -1\n");
        printf("2: 0 1 4 -1\n");
        printf("3: 4 0 1 -1\n");
        printf("4: 1 0 3 2 -1\n\n");

        printf("Adjacency matrix format:\n");
        printf("5\n");
        printf("  1 1 1 1\n");
        printf("    1 1 1\n");
        printf("      0 1\n");
        printf("        1\n\n");

        printf("Loops and duplicate edges are not supported and not removed.\n");

        return OK;
    }
}

