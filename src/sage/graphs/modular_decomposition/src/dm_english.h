/******************************************************

Copyright 2004, 2010 Fabien de Montgolfier
fm@liafa.jussieu.fr

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

**********************************************************/

/********************************************************

MODULAR DECOMPOSITION OF UNDIRECTED GRAPHS
by Fabien de Montgolfier

The program dm.c offer a single function,
decomposition_modulaire().

The input is a graph, its output is its modular decomposition tree.

Input graph is stored using adjacency lists, and trees using a pointers representation (see below)
********************************************************/

#include <stdio.h>
#include <stdlib.h>

/**********************************
Definition of graphs (for the input)
The graph has n vertices. Each vertex is numbered from 0 to n-1.
A graph is a structure (not an object-oriented program !).
If you declare "graphe Gr", Gr.n is the numer of vertices of Gr.
Gr.G is an array of size n. Gr.G[i] points the first element
of adjaceny list of vertex i (NULL if i is isolated)
An adjency list is an usual linked list (vertex;pointer to next).
Adjacency lists may be unsorted.

WARNING : the input graph *MUST* be symmetric : if j belongs to adjacency list of i then i belongs to adjacency list of j. Graphs that do not respect this produce unpredictible and false results.
**********************************/


/* structure for creating an adjacency list */
typedef struct Adj {
  int s; // number of the vertex
  struct Adj *suiv; // adress of next pair in the list, NULL if last
} adj;

typedef struct {
  int n; //number of vertices of the graph
  adj **G; // array of size n. G[i] points the first pair of the adjaceny list of vertex i

} graphe;

/********************************
Output : definition of modular decomposition tree.
Each internal node is labelled SERIE (for series), PARALLELE (for parallel) or PREMIER (for prime) depending of the quotient's type.
Each leaf is labelled FEUILLE and also contains the vertex number of the leaf.
As the tree is an inclusion tree, the vertex-set corresponding to an internal node correspond to the vertices numbers of the leaves that descend from that tree. The function decomposition_modulaire() return a pointer to the root of the tree.



/* define the type of nodes. UNKN,MODULE,ARTEFACT are for internal use*/

#define FEUILLE 0  // the node is a leaf
#define UNKN 1
#define MODULE 2
#define ARTEFACT 3
#define SERIE 4    // series composition
#define PARALLELE 5  // parallel composition
#define PREMIER 6  // prime composition

/* defines a node of the tree */

typedef struct Noeud {
  int type;	// is FEUILLE, SERIE, PARALLELE or PREMIER
  struct Noeud *pere;	// adress of parent node, NULL if root
  struct Fils *fpere;	// points the head of the linked list of sons (if type is not FEUILLE, else is NULL)
  int ps;	// internal use
  int bg;	// internal use
  int ds;       // internal use
  int bd;	// internal use
  int sommet;	// internal use
  int nom;	// if type=FEUILLE, number of the corresponding vertex of the graph
  struct Fils *fils; // points the head of the linked list of sons
  struct Fils *lastfils;  // internal use (points the last item in the listed list of sons)
  int id;	// internal use (node unique ID)
} noeud;

/* linked list that strore the sons of an internal node (in any order) */

typedef struct Fils {
  struct Noeud *pointe; // adress of the node in the tree
  struct Fils *suiv; // adress of the next pair in the list, NULL if last
} fils;

/* prototype of the function.
   Input is a graph, output the root of the modular decomposition tree */

noeud *decomposition_modulaire(graphe G);








