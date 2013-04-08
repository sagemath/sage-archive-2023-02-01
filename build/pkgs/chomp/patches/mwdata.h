/// @addtogroup multiwork
/// @{

/////////////////////////////////////////////////////////////////////////////
///
/// @file multiwork/mwdata.h
///
/// This file contains the definition of the MultiWork data class.
///
/// @author Pawel Pilarczyk
///
/////////////////////////////////////////////////////////////////////////////

// Copyright (C) 1997-2010 by Pawel Pilarczyk.
//
// This file is part of the Homology Library.  This library is free software;
// you can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation;
// either version 2 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this software; see the file "license.txt".  If not, write to the
// Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
// MA 02111-1307, USA.

// Started on August 10, 2004. Last revision: January 25, 2010.


#ifndef _CHOMP_MULTIWORK_MWDATA_H_
#define _CHOMP_MULTIWORK_MWDATA_H_

#include "chomp/multiwork/mwconfig.h"

// include an appropriate header file for 'string' if necessary
#if mwSTRING
	#include <string>
#endif // mwSTRING

// include an appropriate header for 'istream' and 'ostream' if necessary
#if mwSTREAMS
	#include <iostream>
#endif // use mwSTREAMS
#include <stdio.h>


namespace chomp {
namespace multiwork {

// --------------------------------------------------
// --------------------- mwData ---------------------
// --------------------------------------------------

/// This class is used to convert data structures into
/// a single sequence of bytes and to retrieve this data
/// for the purpose of communication between a coordinator and workers.
/// It is assumed that 'int' has at least 32 bits and 'double' has 64 bits.
/// For 'bool', 'long long' and 'string', see the corresponding constants.
/// In the buffer: 'char' = 8-bit, 'short' = 16-bit, 'int' = 32-bit,
/// 'long' and 'long long' = 64-bit, 'float' = 32-bit, 'double' = 64-bit
/// (all big endian).
/// NOTE: Until January 21, 2010 'long' was supposed to be 32 bits long.
class mwData
{
public:
	/// The default constructor of empty data.
	mwData ();

	/// The destructor.
	~mwData ();

	/// The copy constructor.
	mwData (const mwData &x);

	/// The assignment operator.
	mwData &operator = (const mwData &x);

	/// Swaps the data with another data structure.
	mwData &Swap (mwData &x);

	/// Takes the data from another data structure
	/// and makes the other data empty.
	mwData &Take (mwData &x);

	/// Takes raw data that was allocated with the "new" operator.
	mwData &Take (unsigned char *buffer, int length);

	/// Takes raw data that was allocated with the "new" operator.
	mwData &Take (char *buffer, int length);

	/// Returns a pointer to the data buffer.
	const char *Buffer () const;

	/// Returns the length of all the data in the buffer.
	int Length () const;

	/// Rewinds the buffer to the beginning for reading.
	mwData &Rewind ();

	/// Forgets the buffer and makes the data empty.
	mwData &Reset ();

	/// Sets the current reading position in the buffer.
	int Position (int newpos);

	/// Returns the current reading position in the buffer.
	int Position () const;

	/// Appends raw data to the buffer.
	mwData &Append (const char *buffer, int length);

	/// Appends raw data to the buffer.
	mwData &Append (const unsigned char *buffer, int length);

	/// Appends the entire data buffer to the data.
	mwData &Append (const mwData &x);

	/// Returns the currently pointed data in the buffer.
	const char *Current () const;

	// append and retrieve an object to the buffer
	mwData &Append (const int &x);
	mwData &Retrieve (int &x);
	mwData &Append (const unsigned int &x);
	mwData &Retrieve (unsigned int &x);
	mwData &Append (const short &x);
	mwData &Retrieve (short &x);
	mwData &Append (const unsigned short &x);
	mwData &Retrieve (unsigned short &x);
	mwData &Append (const long &x);
	mwData &Retrieve (long &x);
	mwData &Append (const unsigned long &x);
	mwData &Retrieve (unsigned long &x);
	#if mwLONGLONG
	mwData &Append (const long long &x);
	mwData &Retrieve (long long &x);
	mwData &Append (const unsigned long long &x);
	mwData &Retrieve (unsigned long long &x);
	#endif // mwLONGLONG
	mwData &Append (const char &x);
	mwData &Retrieve (char &x);
	mwData &Append (const unsigned char &x);
	mwData &Retrieve (unsigned char &x);
	#if mwBOOL
	mwData &Append (const bool &x);
	mwData &Retrieve (bool &x);
	#endif // mwBOOL
	mwData &Append (const float &x);
	mwData &Retrieve (float &x);
	mwData &Append (const double &x);
	mwData &Retrieve (double &x);
	#if mwSTRING
	mwData &Append (const std::string &x);
	mwData &Retrieve (std::string &x);
	#endif // mwSTRING
	mwData &Append (const char *x);
	mwData &Retrieve (char *x);
	mwData &Append (const unsigned char *x);
	mwData &Retrieve (unsigned char *x);

	/// Skips a zero-terminated string in the buffer.
	mwData &SkipString ();

private:
	/// The data buffer.
	unsigned char *buf;

	/// The length of the data in the buffer and the write position.
	int len;

	/// The length of the allocated buffer memory space.
	int allocated;

	/// The current read position in the buffer.
	int pos;

	/// Increases the buffer by the given number of bytes
	/// (or more) beyond the current position.
	void IncreaseBuffer (int n);

	/// Appends a data piece to the buffer and swaps bytes if necessary.
	mwData &AppendBytes (const unsigned char *x, int n);

	/// Retrieve a data piece from the buffer and swaps bytes
	/// if necessary.
	mwData &RetrieveBytes (unsigned char *x, int n);

}; /* class mwData */

// --------------------------------------------------

inline mwData::mwData ()
{
	buf = (unsigned char *) 0;
	len = 0;
	allocated = 0;
	pos = 0;
	return;
} /* mwData::mwData */

inline mwData::mwData (const mwData &x)
{
	len = x. len;
	allocated = x. allocated;
	pos = x. pos;
	if (allocated)
	{
		buf = new unsigned char [allocated];
		if (!buf)
			throw "No memory for mwData copying constructor.";
		for (int i = 0; i < len; ++ i)
			buf [i] = x. buf [i];
	}
	else
		buf = (unsigned char *) 0;
	return;
} /* mwData::mwData */

inline mwData &mwData::operator = (const mwData &x)
{
	if (this == &x)
		return *this;
	if (buf)
		delete [] buf;
	len = x. len;
	allocated = x. allocated;
	pos = x. pos;
	if (allocated)
	{
		buf = new unsigned char [allocated];
		if (!buf)
			throw "No memory for mwData assignment operator.";
		for (int i = 0; i < len; ++ i)
			buf [i] = x. buf [i];
	}
	else
		buf = (unsigned char *) 0;
	return *this;
} /* mwData::operator = */

inline mwData &mwData::Swap (mwData &x)
{
	// swap the buffers
	unsigned char *buf0;
	buf0 = buf;
	buf = x. buf;
	x. buf = buf0;

	// swap the numbers
	int number;
	number = pos;
	pos = x. pos;
	x. pos = number;
	number = len;
	len = x. len;
	x. len = number;
	number = allocated;
	allocated = x. allocated;
	x. allocated = number;

	return *this;
} /* mwData::Swap */

inline void swap (mwData &x, mwData &y)
{
	x. Swap (y);
	return;
} /* swap */

inline mwData &mwData::Take (mwData &x)
{
	if (buf)
		delete [] buf;
	buf = x. buf;
	x. buf = 0;

	pos = x. pos;
	x. pos = 0;

	len = x. len;
	x. len = 0;

	allocated = x. allocated;
	x. allocated = 0;

	return *this;
} /* mwData::Take */

inline mwData &mwData::Take (unsigned char *buffer, int length)
{
	if (buf)
		delete [] buf;
	buf = buffer;
	pos = 0;
	len = length;
	allocated = length;
	return *this;
} /* mwData::Take */

inline mwData &mwData::Take (char *buffer, int length)
{
	return Take ((unsigned char *) buffer, length);
} /* mwData::Take */

inline mwData::~mwData ()
{
	if (buf)
		delete [] buf;
	return;
} /* mwData::~mwData */

inline const char *mwData::Buffer () const
{
	return reinterpret_cast<const char *> (buf);
} /* mwData::Buffer */

inline int mwData::Length () const
{
	return len;
} /* mwData::Length */

inline void mwData::IncreaseBuffer (int n)
{
	// if it is not necessary to increase the buffer, do nothing
	if (len + n <= allocated)
		return;

	// allocate a new buffer
	int allocated1 = allocated + allocated + n;
	unsigned char *buf1 = new unsigned char [allocated1];
	if (!buf1)
		throw "Not enough memory to increase mwData buffer.";

	// copy the previous buffer to the new one
	for (int i = 0; i < len; ++ i)
		buf1 [i] = buf [i];

	// replace the old buffer with the new one
	delete [] buf;
	buf = buf1;
	allocated = allocated1;
	return;
} /* mwData::IncreaseBuffer */

inline mwData &mwData::AppendBytes (const unsigned char *x, int n)
{
	// increase the buffer if necessary
	IncreaseBuffer (n);

	// if the system is big-endian, the bytes must be copied directly
	const int testnumber = 1;
	if (!*((char *) &testnumber))
	{
		for (int i = 0; i < n; ++ i)
			buf [len ++] = *(x ++);
	}
	// otherwise the bytes must be copied in the reverse order
	else
	{
		x += n;
		for (int i = 0; i < n; ++ i)
			buf [len ++] = *(-- x);
	}
	return *this;
} /* mwData::AppendBytes */

inline mwData &mwData::RetrieveBytes (unsigned char *x, int n)
{
	// if there is not enough data in the buffer, do nothing
	if (len - pos < n)
		return *this;

	// if the system is big-endian, the bytes must be just copied
	const int testnumber = 1;
	if (!*((char *) &testnumber))
	{
		for (int i = 0; i < n; ++ i)
			*(x ++) = buf [pos ++];
	}
	// otherwise the bytes must be swapped
	else
	{
		x += n;
		for (int i = 0; i < n; ++ i)
			*(-- x) = buf [pos ++];
	}
	return *this;
} /* mwData::RetrieveBytes */

inline mwData &mwData::Append (const unsigned char *x, int n)
{
	IncreaseBuffer (n);
	for (int i = 0; i < n; ++ i)
		buf [len ++] = *(x ++);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Append (const char *x, int n)
{
	return Append ((const unsigned char *) x, n);
} /* mwData::Append */

/*
inline mwData::operator const unsigned char * () const
{
	return buf;
}*/ /* mwData::operator const unsigned char * */

/*
inline mwData::operator const char * () const
{
	return (const char *) buf;
}*/ /* mwData::operator const char * */

/*
inline mwData::operator int () const
{
	return len;
}*/ /* mwData::operator int */

inline const char *mwData::Current () const
{
	return (const char *) (buf + pos);
} /* mwData::Current */

inline mwData &mwData::Append (const mwData &x)
{
	return Append (x. Buffer (), x. Length ());
} /* mwData::Append */

inline mwData &mwData::Rewind ()
{
	pos = 0;
	return *this;
} /* mwData::Rewind */

inline mwData &mwData::Reset ()
{
	len = 0;
	pos = 0;
	allocated = 0;
	if (buf)
		delete [] buf;
	buf = 0;
	return *this;
} /* mwData::Reset */

inline int mwData::Position (int newpos)
{
	if ((newpos >= 0) && (newpos <= len))
		pos = newpos;
	return pos;
} /* mwData::Position */

inline int mwData::Position () const
{
	return pos;
} /* mwData::Position */

// --------------------------------------------------

inline mwData &mwData::Append (const int &x)
{
	IncreaseBuffer (4);
	buf [len ++] = (unsigned char) ((x >> 24) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 16) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (int &x)
{
	if (len - pos < 4)
		return *this;
	x = (int) (buf [pos ++]) << 24;
	x |= (int) (buf [pos ++]) << 16;
	x |= (int) (buf [pos ++]) << 8;
	x |= buf [pos ++];
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const unsigned int &x)
{
	IncreaseBuffer (4);
	buf [len ++] = (unsigned char) ((x >> 24) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 16) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (unsigned int &x)
{
	if (len - pos < 4)
		return *this;
	x = (int) (buf [pos ++]) << 24;
	x |= (int) (buf [pos ++]) << 16;
	x |= (int) (buf [pos ++]) << 8;
	x |= buf [pos ++];
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const short &x)
{
	IncreaseBuffer (2);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (short &x)
{
	if (len - pos < 2)
		return *this;
	x = (short) ((short) (buf [pos ++]) << 8);
	x |= buf [pos ++];
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const unsigned short &x)
{
	IncreaseBuffer (2);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (unsigned short &x)
{
	if (len - pos < 2)
		return *this;
	x = (unsigned short) ((unsigned short) (buf [pos ++]) << 8);
	x |= buf [pos ++];
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const long &x)
{
	IncreaseBuffer (8);
#if (__LONG_MAX__ > 2147483647)
	buf [len ++] = (unsigned char) (((unsigned long) x >> 56) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 48) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 40) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 32) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 24) & 0xFF);
#else
	buf [len ++] = 0;
	buf [len ++] = 0;
	buf [len ++] = 0;
	buf [len ++] = 0;
	buf [len ++] = (unsigned char) (((unsigned long) x >> 24) & 0xFF);
#endif
	buf [len ++] = (unsigned char) ((x >> 16) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (long &x)
{
	if (len - pos < 8)
		return *this;
#if (__LONG_MAX__ > 2147483647)
	x = (long) (buf [pos ++]) << 56;
	x |= (long) (buf [pos ++]) << 48;
	x |= (long) (buf [pos ++]) << 40;
	x |= (long) (buf [pos ++]) << 32;
	x |= (long) (buf [pos ++]) << 24;
#else
	pos += 4;
	x = (long) (buf [pos ++]) << 24;
#endif
	x |= (long) (buf [pos ++]) << 16;
	x |= (long) (buf [pos ++]) << 8;
	x |= (long) (buf [pos ++]);
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const unsigned long &x)
{
	IncreaseBuffer (8);
#if (__LONG_MAX__ > 2147483647)
	buf [len ++] = (unsigned char) ((x >> 56) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 48) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 40) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 32) & 0xFF);
#else
	buf [len ++] = 0;
	buf [len ++] = 0;
	buf [len ++] = 0;
	buf [len ++] = 0;
#endif
	buf [len ++] = (unsigned char) ((x >> 24) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 16) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (unsigned long &x)
{
	if (len - pos < 8)
		return *this;
#if (__LONG_MAX__ > 2147483647)
	x = (long) (buf [pos ++]) << 56;
	x |= (long) (buf [pos ++]) << 48;
	x |= (long) (buf [pos ++]) << 40;
	x |= (long) (buf [pos ++]) << 32;
	x |= (long) (buf [pos ++]) << 24;
#else
	pos += 4;
	x = (long) (buf [pos ++]) << 24;
#endif
	x |= (long) (buf [pos ++]) << 16;
	x |= (long) (buf [pos ++]) << 8;
	x |= (long) (buf [pos ++]);
	return *this;
} /* mwData::Retrieve */

#if mwLONGLONG

inline mwData &mwData::Append (const long long &x)
{
	IncreaseBuffer (8);
	buf [len ++] = (unsigned char) ((x >> 56) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 48) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 40) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 32) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 24) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 16) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (long long &x)
{
	if (len - pos < 8)
		return *this;
	x = (long long) (buf [pos ++]) << 56;
	x |= (long long) (buf [pos ++]) << 48;
	x |= (long long) (buf [pos ++]) << 40;
	x |= (long long) (buf [pos ++]) << 32;
	x |= (long long) (buf [pos ++]) << 24;
	x |= (long long) (buf [pos ++]) << 16;
	x |= (long long) (buf [pos ++]) << 8;
	x |= buf [pos ++];
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const unsigned long long &x)
{
	IncreaseBuffer (8);
	buf [len ++] = (unsigned char) ((x >> 56) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 48) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 40) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 32) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 24) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 16) & 0xFF);
	buf [len ++] = (unsigned char) ((x >> 8) & 0xFF);
	buf [len ++] = (unsigned char) (x & 0xFF);
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (unsigned long long &x)
{
	if (len - pos < 8)
		return *this;
	x = (long long) (buf [pos ++]) << 56;
	x |= (long long) (buf [pos ++]) << 48;
	x |= (long long) (buf [pos ++]) << 40;
	x |= (long long) (buf [pos ++]) << 32;
	x |= (long long) (buf [pos ++]) << 24;
	x |= (long long) (buf [pos ++]) << 16;
	x |= (long long) (buf [pos ++]) << 8;
	x |= buf [pos ++];
	return *this;
} /* mwData::Retrieve */

#endif // mwLONGLONG

inline mwData &mwData::Append (const char &x)
{
	IncreaseBuffer (1);
	buf [len ++] = (unsigned char) x;
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (char &x)
{
	if (len - pos < 1)
		return *this;
	x = (char) (buf [pos ++]);
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const unsigned char &x)
{
	IncreaseBuffer (1);
	buf [len ++] = x;
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (unsigned char &x)
{
	if (len - pos < 1)
		return *this;
	x = buf [pos ++];
	return *this;
} /* mwData::Retrieve */

#if mwBOOL

inline mwData &mwData::Append (const bool &x)
{
	IncreaseBuffer (1);
	buf [len ++] = x ? '\1' : '\0';
	return *this;
} /* mwData::Append */

inline mwData &mwData::Retrieve (bool &x)
{
	if (len - pos < 1)
		return *this;
	x = buf [pos ++] ? true : false;
	return *this;
} /* mwData::Retrieve */

#endif // mwBOOL

inline mwData &mwData::Append (const float &x)
{
	IncreaseBuffer (4);
	return AppendBytes ((const unsigned char *) &x, 4);
} /* mwData::Append */

inline mwData &mwData::Retrieve (float &x)
{
	if (len - pos < 4)
		return *this;
	return RetrieveBytes ((unsigned char *) &x, 4);
} /* mwData::Retrieve */

inline mwData &mwData::Append (const double &x)
{
	IncreaseBuffer (8);
	return AppendBytes ((const unsigned char *) &x, 8);
} /* mwData::Append */

inline mwData &mwData::Retrieve (double &x)
{
	if (len - pos < 8)
		return *this;
	return RetrieveBytes ((unsigned char *) &x, 8);
} /* mwData::Retrieve */

#if mwSTRING

inline mwData &mwData::Append (const std::string &x)
{
	const char *str = x. c_str ();
	int length = 0;
	while (str [length ++]);
	IncreaseBuffer (length);
	return Append (reinterpret_cast<const unsigned char *> (str),
		length);
} /* mwData::Append */

inline mwData &mwData::Retrieve (std::string &x)
{
	int pos0 = pos;
	while ((pos0 < len) && buf [pos0])
		++ pos0;
	if (pos0 >= len)
		return *this;
	x = std::string (reinterpret_cast<const char *> (buf + pos));
	pos = pos0 + 1;
	return *this;
} /* mwData::Retrieve */

#endif // mwSTRING

inline mwData &mwData::Append (const unsigned char *x)
{
	int length = 0;
	while (x [length ++]);
	IncreaseBuffer (length);
	return Append (x, length);
} /* mwData::Append */

inline mwData &mwData::Retrieve (unsigned char *x)
{
	int pos0 = pos;
	while ((pos0 < len) && buf [pos0])
		++ pos0;
	if (pos0 >= len)
		return *this;
	while (pos <= pos0)
		*(x ++) = buf [pos ++];
	return *this;
} /* mwData::Retrieve */

inline mwData &mwData::Append (const char *x)
{
	return Append (reinterpret_cast<const unsigned char *> (x));
} /* mwData::Append */

inline mwData &mwData::Retrieve (char *x)
{
	return Retrieve ((unsigned char *) x);
} /* mwData::Retrieve */

// --------------------------------------------------

inline mwData &mwData::SkipString ()
{
	while ((pos < len) && buf [pos])
		++ pos;
	if (pos < len)
		++ pos;
	return *this;
} /* mwData::SkipString */

// --------------------------------------------------

template <class type>
inline mwData &operator << (mwData &m, const type &x)
{
	return m. Append (x);
} /* operator << */

template <class type>
inline mwData &operator >> (mwData &m, type &x)
{
	return m. Retrieve (x);
} /* operator >> */

inline mwData &operator << (mwData &m, const char *x)
{
	return m. Append (x);
} /* operator << */

inline mwData &operator << (mwData &m, const unsigned char *x)
{
	return m. Append (x);
} /* operator << */

inline mwData &operator >> (mwData &m, char *x)
{
	return m. Retrieve (x);
} /* operator >> */

inline mwData &operator >> (mwData &m, unsigned char *x)
{
	return m. Retrieve (x);
} /* operator >> */

#if mwSTREAMS
inline std::ostream &operator << (std::ostream &s, const mwData &m)
{
	if (m. Length ())
		s. write (m. Buffer (), m. Length ());
	return s;
} /* operator << */

inline std::istream &operator >> (std::istream &s, mwData &m)
{
	char *buf = NULL;
	int pos = 0, len = 0;

	// read the entire stream to the given buffer
	int ch = s. get ();
	while (ch != EOF)
	{
		if (len <= pos)
		{
			len = pos + pos + 3;
			char *newbuf = new char [len];
			if (!newbuf)
				break;
			for (int i = 0; i < pos; ++ i)
				newbuf [i] = buf [i];
			delete [] buf;
			buf = newbuf;
		}
		buf [pos ++] = (unsigned char) (ch);
		ch = s. get ();
	}

	// allocate a new buffer of the exact size and take the data
	if (pos)
	{
		char *newbuf = new char [pos];
		for (int i = 0; i < pos; ++ i)
			newbuf [i] = buf [i];
		delete [] buf;
		m. Take (newbuf, pos);
	}

	return s;
} /* operator >> */
#endif


} // namespace multiwork
} // namespace chomp

#endif // _CHOMP_MULTIWORK_MWDATA_H_

/// @}

