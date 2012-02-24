// Calculate rank-width and rank-decompositions of graphs.

// Philipp Klaus Krause, philipp@informatik.uni-frankfurt.de, pkk@spth.de, 2009 - 2012
// Copyright (c) 2009-2012 Philipp Klaus Krause

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

#include <stdlib.h>

#include "rw.h"

static uint_fast8_t subset_size;

static uint_fast8_t num_vertices;

subset_t *adjacency_matrix;

// This array is meant to contain max{f(Y), w(f, Y)} at the position corresponding to Y.
static uint_fast8_t *slots;
subset_t *cslots;

uint_fast8_t cut_rank(const subset_t s)
{
	subset_t am[num_vertices];
	subset_t x, y;
	uint_fast8_t i, j;
	uint_fast8_t rank = 0;

	for(i = 0; i < num_vertices; i++)
		am[i] =  (s & (1ul << i)) ? 0 : (adjacency_matrix[i] & s);

	for(i = 0; i < num_vertices; i++)
	{
		y = 0;
		for(j = rank; j < num_vertices; j++)
		{
			x = am[j];
			if(x & 1)
			{
				if(!y)
				{
					y = x;
					am[j] = am[rank++];
					continue;
				}
				x ^= y;
			}
			am[j] = x >> 1;
		}
	}

	return(rank);
}

// Compute the binomial coefficient C(n, k).
// Time complexity: O(n) * complexity of integer division.
// Could be replaced by a lookup table of size O(n^2), but since this is not the bopttleneck, it doesn't matter.
subset_t binomial_coefficient(uint_fast8_t n, uint_fast8_t k)
{
	uint_fast8_t i, delta, max;
	subset_t c;

	if(n < k)
		return(0);

	if(n == k)
		return(1);

	if(k < n - k)
	{
		delta = n - k;
		max = k;
	}
	else
	{
		delta = k;
		max = n - k;
	}

	c = delta + 1;

	for(i = 2; i <= max; i++)
		c = (c * (delta + i)) / i;

	return(c);
}

// Convert unsigned integer from combinadic to machine representation.
subset_t comb_to_int(subset_t c)
{
	uint_fast8_t k, n;
	subset_t i;

	i = 0;
	for(k = 0, n = 1; k < num_vertices; k++, c >>= 1)
		if(c & 1)
			i += binomial_coefficient(k, n++);
	return(i);
}

// Return largest value v where v < a and  Choose(v,b) <= x.
static uint_fast8_t largest_v(uint_fast8_t a, uint_fast8_t b, subset_t x)
{
	uint_fast8_t v = a - 1;

	while(binomial_coefficient(v,b) > x)
		v--;

	return(v);
}

// Convert unsigned integer from machine representation to combinadic.
static subset_t int_to_comb(subset_t i)
{
	uint_fast8_t j, v;
	uint_fast8_t a = num_vertices;
	uint_fast8_t b = subset_size;
	subset_t c = 0;

	for(j = 0; j < subset_size; j++, b--)
	{
		v = largest_v(a, b, i);
		i = i - binomial_coefficient(v, b);
		a = v;
		c |= (1ul << v);
	}

	return(c);
}

// Masked increment.
static subset_t subset_inc(subset_t v, subset_t mask)
{
	return((v - mask) & mask);
}

// Returns rank-width for a subset of size at least 2 given that slots already contains correct values for all nonempty subsets of sizes less than the size of s.
// This is where most of the time is spent.
uint_fast8_t width(subset_t s)
{
	uint_fast8_t w = UINT_FAST8_MAX, v, v1, v2;
	subset_t ss;
	subset_t cs;
	for(ss = subset_inc(0, s); ss != s; ss = subset_inc(ss, s))
	{
		v1 = slots[ss];
		v2 = slots[s & ~ss];
		v = v1 > v2 ? v1 : v2;
		if(v < w)
		{
			w = v;
			cs = ss;
		}
	}
	cslots[s] = cs;
	return(w);
}

void fill_slot(subset_t i)
{
	uint_fast8_t v, w;
	subset_t s = int_to_comb(i);
	v = cut_rank(s);
	w = width(s);
	slots[s] = v > w ? v : w;
}

void calculate_level(uint_fast8_t ss)
{
	uint_fast8_t i;

	subset_size = ss;

	if(subset_size == 0)
		slots[0] = 0;
	else if(subset_size == 1)
		for(i = 0; i < num_vertices; i++)
		{
			slots[1ul << i] = cut_rank(1ul << i);
			cslots[1ul << i] = 0x0;
		}
	else
	{
		subset_t i;
		const subset_t end = binomial_coefficient(num_vertices, subset_size);
		for(i = 0; i < end; i++)
			fill_slot(i);
	}
}

void calculate_all(void)
{
	uint_fast8_t i;

	for(i = 0; i < num_vertices; i++)
	{
		slots[1ul << i] = cut_rank(1ul << i);
		cslots[1ul << i] = 0x0;
	}

	for(subset_size = 2; subset_size <= num_vertices; subset_size++)
	{
		subset_t i;
		const subset_t end = binomial_coefficient(num_vertices, subset_size);
		for(i = 0; i < end; i++)
			fill_slot(i);
	}
}

int init_rw(uint_fast8_t n)
{
	// If sizeof(uint_fast8_t) * (1ul << n) overflows, it wraps around to 0, since size_t and unsigned long are unsigned integer types.
        if(  (n > MAX_VERTICES) || ( (n>0) && !(sizeof(uint_fast8_t) * (1ul << n)) ) )
		return(-1);

	num_vertices = n;
	adjacency_matrix = malloc(sizeof(subset_t) * n);
	slots = malloc(sizeof(uint_fast8_t) * (1ul << n));
	cslots = 0;
	return((adjacency_matrix && slots) ? 0 : -1);
}

int init_rw_dec(uint_fast8_t n)
{
	// If sizeof(subset_t) * (1ul << n) overflows, it wraps around to 0, since size_t and unsigned long are unsigned integer types.
	if(n && !(sizeof(subset_t) * (1ul << n)))
		return(-1);

	if(init_rw(n))
		return(-1);
	cslots = malloc(sizeof(subset_t) * (1ul << n));
	return(cslots ? 0 : -1);
}

void destroy_rw(void)
{
	free(cslots);
	free(slots);
	free(adjacency_matrix);
}

uint_fast8_t get_rw(void)
{
	return(slots[0xfffffffful >> (32 - num_vertices)]);
}

