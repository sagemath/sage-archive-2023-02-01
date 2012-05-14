#!/usr/bin/env python
"""
This file provides functions to convert Mercurial 1.0 bundles to and
from JSON dumps.
"""

import sys, struct, json, cStringIO
from mercurial.changegroup import readexactly, readbundle, writebundle

#
# bundle --> JSON
#

def unpack_groups(fh):
    """
    A generator of parsed groups from a bundle file. Decompress the bundle
    and discard the magic header before calling this function (i.e. call this
    function with a "changegroup" as input). Expects bundle format HG10
    (the bundle format introduced in Mercurial 1.0).
    """
    yield [chunk for chunk in unpack_chunks(fh)] # the changeloge group
    yield [chunk for chunk in unpack_chunks(fh)] # the manifest group
    while True:
        length, = struct.unpack('>l', readexactly(fh, 4))
        if length <= 4:
            # found a "null meta chunk", which ends the changegroup
            break
        filename = readexactly(fh, length - 4).encode('string_escape')
        yield (filename, [chunk for chunk in unpack_chunks(fh)]) # a file group

def unpack_chunks(fh):
    """
    A generator of parsed chunks of a "group" in a bundle file. Place the
    input head at the beginning of a group, and this function will yield
    parsed chunks until the end of the group (a "null chunk").
    """
    while True:
        length, = struct.unpack('>l', readexactly(fh, 4))
        if length <= 4:
            # found a "null chunk", which ends the group
            break
        if length < 84:
            raise Exception("negative data length")
        node, p1, p2, cs = struct.unpack( '20s20s20s20s',
                readexactly(fh, 80) )
        yield { 'node': node.encode('hex')
              , 'p1': p1.encode('hex')
              , 'p2': p2.encode('hex')
              , 'cs': cs.encode('hex')
              , 'data': [patch for patch in unpack_patches(fh, length - 84)] }

def unpack_patches(fh, remaining):
    """
    A generator of patches from the data field in a chunk. As there is
    no delimiter for this data field, we require a length argument.
    """
    while remaining >= 12:
        start, end, blocklen = struct.unpack('>lll', readexactly(fh, 12))
        remaining -= 12
        if blocklen > remaining:
            raise Exception("unexpected end of patch stream")
        block = readexactly(fh, blocklen)
        remaining -= blocklen
        yield { 'start': start
              , 'end': end
              , 'blocklen': blocklen
              , 'block': block.encode('string_escape')
              }

    if remaining > 0:
        print remaining
        raise Exception("unexpected end of patch stream")

def to_json(ifilename, ofilename, compression='none'):
    """
    Given an HG10xx file (Mercurial 1.0 bundle) ``ifilename``, convert
    it to JSON and dump it to ``ofilename``
    """
    if ifilename:
        ifile = open(ifilename, 'rb')
    else:
        ifile = sys.stdin
    fh = readbundle(ifile, ifilename)

    oobj = [group for group in unpack_groups(fh)]
    # this is unimplemented in Mercurial 1.8.4, unfortunately (?)
    #fh.close()

    if ofilename:
        ofile = open(ofilename, 'w')
    else:
        ofile = sys.stdout
    json.dump(oobj, ofile, indent=4)
    ofile.close()

#
# JSON --> bundle
#

def pack_groups(cg_obj):
    """
    Given a JSON object ``cg_obj`` representing a changegroup, return the
    changegroup.
    """
    return ''.join([ pack_chunks(cg_obj[0])   # the changelog group
                   , pack_chunks(cg_obj[1]) ] # the manifest group
                 + [ ''.join([ struct.pack('>l', len(pair[0].decode('string_escape')) + 4)
                             , pair[0].decode('string_escape') # a filename
                             , pack_chunks(pair[1]) ]) # a file group
                   for pair in cg_obj[2:] ] # over all files
                 + [ struct.pack('>l', 0) ]) # a null meta chunk to end files

def pack_chunks(group_obj):
    """
    Given a JSON object ``group_obj`` representing a group, create the group
    and write it to ``fh``
    """
    def chunk_gen():
        for chunk in group_obj:
            chunk['data'] = pack_patches(chunk['data'])
            yield ''.join([ struct.pack('>l', len(chunk['data']) + 84)
                          , chunk['node'].decode('string_escape').decode('hex')
                          , chunk['p1'].decode('string_escape').decode('hex')
                          , chunk['p2'].decode('string_escape').decode('hex')
                          , chunk['cs'].decode('string_escape').decode('hex')
                          , chunk['data'] ])
    return ''.join([ packed_chunk for packed_chunk in chunk_gen() ]
                 + [ struct.pack('>l', 0) ]) # a null chunk to end the group

def pack_patches(data_obj):
    """
    Pack a given list of patch objects into a binary patch stream.
    """
    return ''.join([ ''.join([ struct.pack( '>lll'
                                          , patch['start']
                                          , patch['end']
                                          , patch['blocklen'] )
                             , patch['block'].decode('string_escape') ])
                   for patch in data_obj ])

def to_hg(ifilename, ofilename, compression='bzip2'):
    """
    Given a JSON file produced by this script, at ``ifilename``, convert
    it into an HG10UN file (Mercurial 1.0 uncompressed bundle) at
    ``ofilename``.
    """
    magics = { 'none': 'HG10UN', 'bzip2': 'HG10BZ', 'gzip': 'HG10GZ' }
    try:
        magic = magics[compression]
    except KeyError:
        raise ValueError("unsupported compression type '{0}'".format(
                compression ))

    if ifilename:
        ifile = open(ifilename, 'r')
    else:
        ifile = sys.stdin
    iobj = json.load(ifile, encoding='ascii')

    writebundle(cStringIO.StringIO(pack_groups(iobj)), ofilename, magic)

#
# Command line usage
#

def alternative_cmdline_handler():
    import os
    if len(sys.argv) != 2:
        sys.exit(os.EX_USAGE)
    filename = sys.argv[1]
    if filename.endswith('.hg'):
        to_json(filename, filename[:-3] + '.json')
    elif filename.endswith('.json'):
        to_hg(filename, filename[:-5] + '.hg')
    else:
        sys.exit(os.EX_USAGE)
    sys.exit(os.EX_OK)

if __name__ == '__main__':
    try:
        import argparse
    except:
        alternative_cmdline_handler()

    parser = argparse.ArgumentParser( description='Convert Mercurial bundle '
            'format into JSON text files and back', epilog='''
There is a more concise alternative syntax as well; pass a .hg file to
get a .json file in the same directory and vice versa. Defaults to bz2
compression for bundles.
                    ''')

    parser.add_argument( 'filename', type=str, help='the filename of the '
            'bundle in HG10UN or JSON format' )
    parser.add_argument( '--to-json', dest='selected_converter',
            action='store_const', const=to_json,
            help='convert from HG10UN to JSON' )
    parser.add_argument( '--to-bundle', dest='selected_converter',
            action='store_const', const=to_hg,
            help='convert from JSON to HG10UN' )
    parser.add_argument( '-t', '--type', dest='compression',
            choices=['none', 'bzip2', 'gzip'], help='compression format for '
            'hg bundle', default='bzip2')
    parser.add_argument( '-o', '--output-file', help='file to write (default '
            'is stdout for json files, temp file for bundles')

    args = parser.parse_args()
    if not args.selected_converter:
        alternative_cmdline_handler()

    args.selected_converter( args.filename, args.output_file,
            compression=args.compression )
