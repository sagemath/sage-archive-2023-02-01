// A simple program for calculating rank-width and rank-decompositions.
// Compile using gcc -O2 --std=c99 -pedantic rw.c simplerw.c -ligraph

// Philipp Klaus Krause, philipp@informatik.uni-frankfurt.de, pkk@spth.de, 2009 - 2011
// Copyright (c) 2009-2011 Philipp Klaus Krause

// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 2 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc.,
// 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <igraph/igraph.h>

#include "rw.h"

static int num_vertices = 18;

static uint_fast32_t bitmask(int i)
{
	return(1ul << i);
}

static void set_am(int i, int j, int val)
{
	adjacency_matrix[i] &= ~bitmask(j);
	adjacency_matrix[j] &= ~bitmask(i);
	if(val)
	{
		adjacency_matrix[i] |= bitmask(j);
		adjacency_matrix[j] |= bitmask(i);
	}
}

static void tab(int l)
{
	while(l--)
		printf("\t");
}

void print_rank_dec(subset_t s, int l)
{
	tab(l);
	printf("cslot: %x\n", (unsigned)s);
	if(cslots[s] == 0)
		return;
	print_rank_dec(cslots[s], l + 1);
	print_rank_dec(s & ~cslots[s], l + 1);
}

/*int random_graph(void)
{
	int i, j;

	if(init_rw_dec(num_vertices))
		return(-1);

	srand(23);
	for(i = 0; i < num_vertices; i++)
		for(j = 0; j < num_vertices; j++)
			set_am(i, j, rand() % 2);

	for(i = 0; i < num_vertices; i++)
		printf("Adj. mat.: %x\n", (unsigned)(adjacency_matrix[i]));

	return(0);
}*/

int read_graph(const char *format, const char * filename)
{
	int ret, i, j;
	FILE *file;
	igraph_t igraph;
	igraph_matrix_t imatrix;
	igraph_bool_t isimple;

	// Open file
	if(!(file = fopen(filename, "r")))
	{
		printf("Failed to open file %s.\n", filename);
		return(-1);
	}

	// Read graph from file
	if(!strcmp("--edgelist", format))
		ret = igraph_read_graph_edgelist(&igraph, file, 0, 0);
	else if(!strcmp("--graphdb", format))
		ret = igraph_read_graph_graphdb(&igraph, file, 0);
	else if(!strcmp("--graphml", format))
		ret = igraph_read_graph_graphml(&igraph, file, 0);
	else if(!strcmp("--gml", format))
		ret = igraph_read_graph_gml(&igraph, file);
	else if(!strcmp("--pajek", format))
		ret = igraph_read_graph_pajek(&igraph, file);
	else
	{
		printf("Unknown file format %s.\n", format);
		fclose(file);
		return(-1);
	}

	fclose(file);

	if(ret)
	{
		printf("Failed to read a graph from file %s.\n", filename);
		return(-1);
	}

	// Check for undirectedness
	if(igraph_is_directed(&igraph) || igraph_is_simple(&igraph, &isimple) || !isimple)
	{
		printf("Input is not a simple, undirected graph from file %s.\n", filename);
		igraph_destroy(&igraph);
		return(-1);
	}

	// Get adjacency matrix
	if(igraph_matrix_init(&imatrix, 1, 1))
	{
		igraph_destroy(&igraph);
		return(-1);
	}
	igraph_get_adjacency(&igraph, &imatrix, IGRAPH_GET_ADJACENCY_BOTH);
	igraph_destroy(&igraph);
	if(igraph_matrix_nrow(&imatrix) > MAX_VERTICES)
	{
		igraph_matrix_destroy(&imatrix);
		printf("Graph too large.\n");
		return(-1);
	}

	num_vertices = igraph_matrix_nrow(&imatrix);

	if(init_rw_dec(num_vertices))
	{
		destroy_rw();
		printf("Not enough memory.\n");
		return(-1);
	}

	for(i = 0; i < num_vertices; i++)
		for(j = 0; j < num_vertices; j++)
			set_am(i, j, MATRIX(imatrix, i, j));

	igraph_matrix_destroy(&imatrix);

	return(ret ? -1 : 0);
}

int main(int argc, char *argv[])
{
	int i;

	if(argc /*==*/ <= 2 || argc > 4)
	{
		printf("Usage: rw [--format filename]\n");
		printf("Supported formats: edgelist, graphdb, graphml, gml, pajek.\n");
		return(-1);
	}

	if(/*argc <= 1 ? random_graph() :*/ read_graph(argv[1], argv[2]))
	{
		printf("Failed to create input graph.\n");
		return(-1);
	}

	printf("%d vertices.\n", (int)num_vertices);

	for(i = 0; i <= num_vertices; i++)
	{
		printf("Calculating for subsets of size %d.\n", i);
		calculate_level(i);
	}

	printf("rank-width: %d\n", (int)get_rw());

	print_rank_dec(0x7ffffffful >> (31 - num_vertices), 0);

	destroy_rw();

	return(0);
}

