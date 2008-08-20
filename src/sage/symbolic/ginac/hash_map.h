/** @file hash_map.h
 *
 *  Replacement for map<> using hash tables. */

/*
 *  GiNaC Copyright (C) 1999-2008 Johannes Gutenberg University Mainz, Germany
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef __GINAC_HASH_MAP_H__
#define __GINAC_HASH_MAP_H__

#include <list>
#include <iterator>
#include <algorithm>
#include <functional>
#include <utility>


namespace GiNaC {

/*
 *  "Hashmap Light" - buckets only contain one value, quadratic probing,
 *  grows automatically
 */

namespace internal {

// List of prime numbers shamelessly stolen from GCC STL
enum { num_primes = 29 };

static const unsigned long prime_list[num_primes] =
{
	31ul,        53ul,         97ul,         193ul,       389ul,
	769ul,       1543ul,       3079ul,       6151ul,      12289ul,
	24593ul,     49157ul,      98317ul,      196613ul,    393241ul,
	786433ul,    1572869ul,    3145739ul,    6291469ul,   12582917ul,
	25165843ul,  50331653ul,   100663319ul,  201326611ul, 402653189ul,
	805306457ul, 1610612741ul, 3221225473ul, 4294967291ul
};

inline unsigned long next_prime(unsigned long n)
{
	const unsigned long *first = prime_list;
	const unsigned long *last = prime_list + num_primes;
	const unsigned long *pos = std::lower_bound(first, last, n);
	return pos == last ? *(last - 1) : *pos;
}

} // namespace internal


// Define default arguments
template <typename T, template <class> class A = std::allocator>
class exhashmap;


/** Pair Associative Container with 'ex' objects as keys, that is implemented
 *  with a hash table and can be used as a replacement for map<> in many cases.
 *
 *  Differences to map<>:
 *   - no lower_bound()/upper_bound()
 *   - no reverse iterators, no rbegin()/rend()
 *   - no operator<()
 *   - comparison functor is hardcoded to ex_is_less
 *   - bucket_count() returns the number of buckets allocated in the hash table
 *   - insert() and erase() invalidate all iterators
 *   - average complexity of find() is constant time, worst case is O(n) */
template <typename T, template <class> class A>
class exhashmap {
public:
	static const unsigned min_num_buckets = 31; // must be prime

	// Standard types
	typedef ex key_type;
	typedef T mapped_type;
	typedef std::pair<key_type, T> value_type;
	typedef ex_is_less key_compare;
	typedef ex_is_equal key_equal;
	typedef value_type & reference;
	typedef const value_type & const_reference;
	typedef value_type * pointer;
	typedef const value_type * const_pointer;

protected:
	// Private types
	enum bucket_state {
		EMPTY,  ///< bucket empty (never used)
		USED,   ///< bucket in use
		ERASED  ///< bucket empty (element deleted), but may be part of a search chain
	};
	typedef std::pair<bucket_state, value_type> Bucket;

public:
	// More standard types
	typedef A<Bucket> allocator_type;

protected:
	// More private types
	typedef std::vector<Bucket, allocator_type> Table;

	typedef typename Table::iterator table_iterator;
	typedef typename Table::const_iterator table_const_iterator;

public:
	// Iterators
	template <typename Pointer, typename Reference, class TableIterator>
	class exhashmap_iterator : public std::iterator<std::forward_iterator_tag, value_type, typename Table::difference_type, Pointer, Reference> {
	protected:
		friend class exhashmap;

	public:
		exhashmap_iterator() {}
		exhashmap_iterator(TableIterator t, TableIterator te)
		 : it(t), table_end(te) {}

		// Allow iterator to const_iterator conversion
		template <typename P, typename R, class TI>
		exhashmap_iterator(const exhashmap_iterator<P, R, TI> &other)
		 : it(other.get_it_()), table_end(other.get_table_end_()) {}

		typename exhashmap_iterator::reference operator*() const
		{
			return it->second;
		}

		typename exhashmap_iterator::pointer operator->() const
		{
			return &(it->second);
		}

		exhashmap_iterator &operator++()
		{
			increment();
			return *this;
		}

		exhashmap_iterator operator++(int)
		{
			exhashmap_iterator tmp = *this;
			increment();
			return tmp;
		}

		template <typename P, typename R, class TI>
		bool operator==(const exhashmap_iterator<P, R, TI> &other) const
		{
			return it == other.get_it_();
		}

		template <typename P, typename R, class TI>
		bool operator!=(const exhashmap_iterator<P, R, TI> &other) const
		{
			return it != other.get_it_();
		}

		// Private access function
		TableIterator get_it_() const { return it; }
		TableIterator get_table_end_() const { return table_end; }

	protected:
		TableIterator it;        ///< Pointer to current bucket
		TableIterator table_end; ///< Pointer to one-past-last bucket

		void increment()
		{
			if (it != table_end)
				++it;

			// Skip empty and erased buckets
			while (it != table_end && it->first != USED)
				++it;
		}
	};

	typedef exhashmap_iterator<value_type*, value_type&, table_iterator> iterator;
	typedef exhashmap_iterator<const value_type*, const value_type&, table_const_iterator> const_iterator;

	// More standard types
	typedef typename Table::size_type size_type;
	typedef typename Table::difference_type difference_type;

	class value_compare : public std::binary_function<value_type, value_type, bool>, private key_compare {
		friend class exhashmap;
	public:
		bool operator()(const value_type &lhs, const value_type &rhs) const
		{
			return key_compare::operator()(lhs.first, rhs.first);
		}

		bool operator()(const key_type &lhs, const value_type &rhs) const
		{
			return key_compare::operator()(lhs, rhs.first);
		}

		bool operator()(const value_type &lhs, const key_type &rhs) const
		{
			return key_compare::operator()(lhs.first, rhs);
		}
	};

protected:
	// Private data
	size_type num_entries; ///< Number of values stored in container (cached for faster operation of size())
	size_type num_buckets; ///< Number of buckets (= hashtab.size())
	Table hashtab;         ///< Vector of buckets, each bucket is kept sorted

	/** Return index of key in hash table. */
	static size_type hash_index(const key_type &x, size_type nbuckets)
	{
		return x.gethash() % nbuckets;
	}

	static table_iterator find_bucket(const key_type &x, table_iterator tab, size_type nbuckets);
	static table_const_iterator find_bucket(const key_type &x, table_const_iterator tab, size_type nbuckets);
	static table_iterator find_bucket_for_insertion(const key_type &x, table_iterator tab, size_type nbuckets);

	/** Return pointer to bucket corresponding to key (or first empty bucket). */
	table_iterator find_bucket(const key_type &x)
	{
		return find_bucket(x, hashtab.begin(), num_buckets);
	}

	/** Return pointer to bucket corresponding to key (or first empty bucket). */
	table_const_iterator find_bucket(const key_type &x) const
	{
		return find_bucket(x, hashtab.begin(), num_buckets);
	}

	/** Return pointer to bucket corresponding to key (or first empty or erased bucket). */
	table_iterator find_bucket_for_insertion(const key_type &x)
	{
		return find_bucket_for_insertion(x, hashtab.begin(), num_buckets);
	}

	/** Return number of entries above which the table will grow. */
	size_type hwm() const
	{
		// Try to keep at least 25% of the buckets free
		return num_buckets - (num_buckets >> 2);
	}

	void grow();

public:
	// 23.3.1.1 Construct/copy/destroy
	exhashmap()
	 : num_entries(0), num_buckets(min_num_buckets), hashtab(num_buckets, std::make_pair(EMPTY, std::make_pair(0, mapped_type()))) {}

	explicit exhashmap(size_type nbuckets)
	 : num_entries(0), num_buckets(internal::next_prime(nbuckets)), hashtab(num_buckets, std::make_pair(EMPTY, std::make_pair(0, mapped_type()))) {}

	template <class InputIterator>
	exhashmap(InputIterator first, InputIterator last)
	 : num_entries(0), num_buckets(min_num_buckets), hashtab(num_buckets, std::make_pair(EMPTY, std::make_pair(0, mapped_type())))
	{
		insert(first, last);
	}

	exhashmap &operator=(const exhashmap &other)
	{
		exhashmap(other).swap(*this);
		return *this;
	}

	// Iterators
	iterator begin()
	{
		// Find first used bucket
		table_iterator bucket = hashtab.begin();
		while (bucket != hashtab.end() && bucket->first != USED)
			++bucket;
		return iterator(bucket, hashtab.end());
	}

	const_iterator begin() const
	{
		// Find first used bucket
		table_const_iterator bucket = hashtab.begin();
		while (bucket != hashtab.end() && bucket->first != USED)
			++bucket;
		return const_iterator(bucket, hashtab.end());
	}

	iterator end()
	{
		return iterator(hashtab.end(), hashtab.end());
	}

	const_iterator end() const
	{
		return const_iterator(hashtab.end(), hashtab.end());
	}

	// Capacity
	bool empty() const
	{
		return num_entries == 0;
	}

	size_type size() const
	{
		return num_entries;
	}

	size_type max_size() const
	{
		return hashtab.max_size();
	}

	size_type bucket_count() const
	{
		return num_buckets;
	}

	// 23.3.1.2 Element access
	T &operator[](const key_type &x)
	{
		return insert(value_type(x, mapped_type())).first->second;
	}

	// Modifiers
	std::pair<iterator, bool> insert(const value_type &x);

	iterator insert(iterator pos, const value_type &x)
	{
		return insert(x).first;
	}

	template <class InputIterator>
	void insert(InputIterator first, InputIterator last)
	{
		for (; first != last; ++first)
			insert(*first);
	}

	void erase(iterator position)
	{
		table_iterator bucket = position.get_it_();
		bucket->first = ERASED;
		bucket->second.first = 0;
		--num_entries;
	}

	size_type erase(const key_type &x);

	void swap(exhashmap &other)
	{
		hashtab.swap(other.hashtab);
		std::swap(num_buckets, other.num_buckets);
		std::swap(num_entries, other.num_entries);
	}

	void clear();

	// Observers
	key_compare key_comp() const
	{
		return key_compare();
	}

	value_compare value_comp() const
	{
		return value_compare();
	}

	// 23.3.1.3 Map operations
	iterator find(const key_type &x);
	const_iterator find(const key_type &x) const;

	size_type count(const key_type &x) const
	{
		return find(x) == end() ? 0 : 1;
	}

	std::pair<iterator, iterator> equal_range(const key_type &x)
	{
		iterator i = find(x);
		if (i == end())
			return std::make_pair(i, i);
		else {
			iterator j = ++i;
			return std::make_pair(i, j);
		}
	}

	std::pair<const_iterator, const_iterator> equal_range(const key_type &x) const
	{
		const_iterator i = find(x);
		if (i == end())
			return std::make_pair(i, i);
		else {
			const_iterator j = ++i;
			return std::make_pair(i, j);
		}
	}

	friend bool operator==(const exhashmap &lhs, const exhashmap &rhs)
	{
		if (lhs.num_entries != rhs.num_entries || lhs.num_buckets != rhs.num_buckets)
			return false;

		// We can't compare the tables directly as the elements may be
		// in different order due to the collision handling. We therefore
		// look up each value from the lhs map in the rhs map separately.
		for (const_iterator itl = lhs.begin(); itl != lhs.end(); ++itl) {
			const_iterator itr = rhs.find(itl->first);
			if (itr == rhs.end())
				return false;
			if (itl->second != itr->second)
				return false;
		}
		return true;
	}

	friend bool operator!=(const exhashmap &lhs, const exhashmap &rhs)
	{
		return !(lhs == rhs);
	}

#if 0
	void dump() const
	{
		std::clog << "num_entries = " << num_entries << std::endl;
		std::clog << "num_buckets = " << num_buckets << std::endl;
		size_type t = 0;
		for (table_const_iterator it = hashtab.begin(); it != hashtab.end(); ++it, ++t) {
			std::clog << " bucket " << t << ": ";
			std::clog << (it->first == EMPTY ? "free" : (it->first == USED ? "used" : "erased")) << ", " << it->second.first << " -> " << it->second.second << std::endl;
		}
	}
#endif
};

/** Return pointer to bucket corresponding to key (or first empty bucket). */
template <typename T, template <class> class A>
inline typename exhashmap<T, A>::table_iterator exhashmap<T, A>::find_bucket(const key_type &x, table_iterator tab, size_type nbuckets)
{
	// Quadratic probing
	size_type h = hash_index(x, nbuckets);
	size_type d = 1;
	table_iterator it = tab + h;
	while (it->first != EMPTY && !(it->first == USED && key_equal()(it->second.first, x))) {
		h = (h + d) % nbuckets;
		d += 2;
		it = tab + h;
	}
	return it;
}

/** Return pointer to bucket corresponding to key (or first empty bucket). */
template <typename T, template <class> class A>
inline typename exhashmap<T, A>::table_const_iterator exhashmap<T, A>::find_bucket(const key_type &x, table_const_iterator tab, size_type nbuckets)
{
	// Quadratic probing
	size_type h = hash_index(x, nbuckets);
	size_type d = 1;
	table_const_iterator it = tab + h;
	while (it->first != EMPTY && !(it->first == USED && key_equal()(it->second.first, x))) {
		h = (h + d) % nbuckets;
		d += 2;
		it = tab + h;
	}
	return it;
}

/** Return pointer to bucket corresponding to key (or first empty or erased bucket). */
template <typename T, template <class> class A>
inline typename exhashmap<T, A>::table_iterator exhashmap<T, A>::find_bucket_for_insertion(const key_type &x, table_iterator tab, size_type nbuckets)
{
	// Quadratic probing
	size_type h = hash_index(x, nbuckets);
	size_type d = 1;
	table_iterator it = tab + h;
	while (it->first != EMPTY && !key_equal()(it->second.first, x)) {
		h = (h + d) % nbuckets;
		d += 2;
		it = tab + h;
	}
	return it;
}

/** Grow hash table */
template <typename T, template <class> class A>
void exhashmap<T, A>::grow()
{
	// Allocate new empty hash table
	size_type new_num_buckets = internal::next_prime(num_buckets + 1);
	Table new_hashtab;
	new_hashtab.resize(new_num_buckets);
	for (table_iterator it = new_hashtab.begin(); it != new_hashtab.end(); ++it)
		it->first = EMPTY;

	// Re-insert all elements into new table
	for (table_const_iterator it = hashtab.begin(); it != hashtab.end(); ++it) {
		if (it->first == USED) {
			table_iterator bucket = find_bucket(it->second.first, new_hashtab.begin(), new_num_buckets);
			*bucket = *it;
		}
	}

	// Swap with the old table
	hashtab.swap(new_hashtab);
	num_buckets = new_num_buckets;
}

template <typename T, template <class> class A>
std::pair<typename exhashmap<T, A>::iterator, bool> exhashmap<T, A>::insert(const value_type &x)
{
	table_iterator bucket = find_bucket_for_insertion(x.first);
	if (bucket->first == USED) {
		// Value already in map
		return std::make_pair(iterator(bucket, hashtab.end()), false);
	} else {
		// Insert new value
		bucket->first = USED;
		bucket->second = x;
		++num_entries;
		if (num_entries >= hwm()) {
			grow();
			bucket = find_bucket(x.first);
		}
		return std::make_pair(iterator(bucket, hashtab.end()), true);
	}
}

template <typename T, template <class> class A>
typename exhashmap<T, A>::size_type exhashmap<T, A>::erase(const key_type &x)
{
	iterator i = find(x);
	if (i != end()) {
		erase(i);
		return 1;
	} else
		return 0;
}

template <typename T, template <class> class A>
typename exhashmap<T, A>::iterator exhashmap<T, A>::find(const key_type &x)
{
	table_iterator bucket = find_bucket(x);
	if (bucket->first == USED)
		return iterator(bucket, hashtab.end());
	else
		return end();
}

template <typename T, template <class> class A>
typename exhashmap<T, A>::const_iterator exhashmap<T, A>::find(const key_type &x) const
{
	table_const_iterator bucket = find_bucket(x);
	if (bucket->first == USED)
		return const_iterator(bucket, hashtab.end());
	else
		return end();
}

template <typename T, template <class> class A>
void exhashmap<T, A>::clear()
{
	for (table_iterator i = hashtab.begin(); i != hashtab.end(); ++i) {
		i->first = EMPTY;
		i->second.first = 0;
		i->second.second = mapped_type();
	}
	num_entries = 0;
}

} // namespace GiNaC


// Specializations of Standard Library algorithms
namespace std {

/** Specialization of std::swap() for exhashmap. */
template <typename T, template <class> class A>
inline void swap(GiNaC::exhashmap<T, A> &lhs, GiNaC::exhashmap<T, A> &rhs)
{
	lhs.swap(rhs);
}

} // namespace std

#endif // ndef __GINAC_HASH_MAP_H__
