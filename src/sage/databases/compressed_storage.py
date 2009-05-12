"""
Compression for ZODB.
"""

#*****************************************************************************
#
#   Sage: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
###############################################################################
#
# This code is a modified version of the code I got here:
#
#      http://www.zope.org/Members/tsarna/CompressedStorage
#
# I modified it (very slightly) to work with the new ZODB...
#*****************************************************************************


# This is the original copyright notice:
#
#       $Id$
#
# Copyright (c) 1999 Tyler C. Sarna (tsarna@endicor.com)
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions
# are met:
# 1. Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
# 3. The name of the author may not be used to endorse or promote products
#    derived from this software without specific prior written permission
#
# THIS SOFTWARE IS PROVIDED BY THE AUTHOR ``AS IS'' AND ANY EXPRESS OR
# IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
# OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
# IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT
# NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF
# THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#

"""
Compressed Layered Storage

This is a  storage that layers over another storage, compressing
objects larger than a certain threshold (default 2K).

The format of compressed pickles is an optional prefix byte followed by
data.  If the prefix byte is 1, the data is compressed.  If the
prefix byte is 0, the data is uncompressed.  If the first byte is not 0
or 1, the data is assumed uncompressed.  Since python pickles never
start with 0 or 1, they are stored directly without a prefix byte.  This
makes compatibility with existing storages easier -- you can simply
stick a CompressedStorage atop an existing storage and it will work.

However, if some other layered storage sits atop this one, say an
encrypting storage, the "pickle" this layer receives may start with a 0
or 1.  In that event, if this layer won't be compressing it, a 0 byte
will be prefixed to assure that CompressedStorage doesn't try to
decompress data that isn't compressed.

"""

import bz2


# The following would only be used for pack, and pack doesn't work...
class ZWrap:
    """
    Wrap a function to produce a new one that can take a compressed
    string as the first argument.
    """

    def __init__(self, f):
        self.func = f

    def __call__(self, *args, **kw):
        if type(args) != type(()):
            args = tuple(args)
        if args[0][0] == '\1':
            args = (bz2.decompress(args[0][1:]),) + args[1:]
        elif args[0][1] == '\0':
            args = (args[0][1:],) + args[1:]
        return apply(self.func, args, kw)


class CompressedStorage:
    def __init__(self, base, thresh=2048):
        self._base = base
        self._zthresh = thresh

    def __getattr__(self, attr):
        """
        Pseudo-acquisition for base's stuff that we don't
        otherwise override.
        """
        return getattr(self._base, attr)

    def load(self, oid, version):
        p, serial = self._base.load(oid, version)
        if p[0] == '\1':
            p = bz2.decompress(p[1:])
        elif p[0] == '\0':
            p = p[1:]
        return p, serial

    def store(self, oid, serial, data, version, transaction):
        if self._is_read_only:
            raise POSException.ReadOnlyError()
        datalen = len(data)
        compressed = 0
        if datalen >= self._zthresh:
            ndata = '\1' + bz2.compress(data)
            # print datalen - len(ndata), datalen, len(ndata)
            if len(ndata) < datalen:
                data = ndata
                compressed = 1

        if not compressed and data[0] in '\0\1':
            data = '\0' + data

        return self._base.store(oid, serial, data, version, transaction)

    def pack(self, t, referencesf):
        return self._base.pack(t, ZWrap(referencesf))
