// Calculate rank-width and rank-decompositions of graphs.

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

// This is a program that calculates rank-width and rank-decompositions. It is based on ideas from "Computing rank-width exactly" by Sang-il Oum, "Sopra una formula numerica" by Ernesto Pascal and "Generation of a Vector from the Lexicographical Index" by B.P. Buckles and M. Lybanon.
// On 2009's computers it works quite well up to graph sizes of about 28 nodes.
// It is an implementation of the trivial algorithm from "Computing rank-width exactly". For larger graphs (more than about 40 nodes) the more algorithm based on fast subset convolution from the same paper would probably be faster.

#include <stdint.h>

// Use data type uint_leastN_t. N is an upper limit on the size of the graphs that can be handled. N=32 seems to be a good compromise for now (the code works well with other values of N).
// uint_leastN_t is faster than uint_fastN_t here, since the bottleneck is cache misses.
#ifndef __RANKWIDTH_H_SUBSET_T__
#define __RANKWIDTH_H_SUBSET_T__
typedef uint_least32_t subset_t;
#endif

#define MAX_VERTICES 32

// Input graph.
//extern subset_t adjacency_matrix[NUM_VERTICES];
subset_t *adjacency_matrix;

// Output rank-decomposition
// This array is meant to contain max{f(Y), w(f, Y)} at the position corresponding to Y.
static uint_fast8_t *slots;
subset_t *cslots;

// Initialization (for getting rank-width only). Returns 0 on success.
//int init_rw(uint_fast8_t n);

// Initialization (for getting both rank-width and rank-decomposition). Returns 0 on success.
int init_rw_dec(uint_fast8_t n);

// Free memory allocated during initialization.
void destroy_rw(void);

// Calculate everything. May take some time.
void calculate_all(void);

// Calculate a single level only
void calculate_level(uint_fast8_t subset_size);

// Get the rank-width.
uint_fast8_t get_rw(void);

static uint_fast8_t subset_size;

static uint_fast8_t num_vertices;
