#!/usr/bin/env python

# Run this script on a vanilla sage build on a clean OS X with -w to figure out
# what external libraries vanilla sage links to.  This forms a reference
# for figuring out build problems with non-vanilla configurations.
#
# dphilp 15/9/8

import platform, string, os, re, getopt, sys

def usage() :
  print '''Script for finding all external libraries used by Sage. (OS X only.)

Usage options, use one only:
  -o list_name    create a list of linked external libraries
  -w              create whitelist and save it in SAGE_ROOT/data/extcode
  -h              this message
  -c              check against whitelist (default)
'''

sage_root = os.environ.get('SAGE_ROOT')
mode = ''
output_filename = ''

# What is the major OS X version number?
def os_version_str() :
  result = platform.mac_ver()[0]
  result = result.split('.')
  result = result[:-1]
  result = string.join(result, '-')
  return result

def whitelist_filename() :
  return sage_root+'/data/extcode/sage-library-whitelist-osx-'+os_version_str()

# At startup, look for sage_root, check it's a Mac, check the
# options make sense, and define the output_filename
# On completion, mode = 'output' if whitelisting or spitting out a file,
# and output_filename = <something correct> if mode = 'output'.
# Otherwise, mode = '' and we are in checking mode.
def startup() :
  global sage_root, output_filename, mode
  # Check that sage_root is defined
  if (sage_root == None) :
    print 'Set $SAGE_ROOT before calling this script.'
    sys.exit(1)

  # Check this is a Mac
  if not(platform.system() == 'Darwin') :
    print "This script is for Mac OS X only"
    sys.exit(1)

  # Check the options are good; define mode and output_filename
  try:
    opts, args = getopt.getopt(sys.argv[1:], 'o:whc')
  except getopt.GetoptError, err:
    print str(err)
    usage()
    sys.exit(1)
  if len(opts) > 1 :
    usage()
    sys.exit(1)
  for o, a in opts :
    if o == '-w' :
      mode = 'output'
      output_filename = whitelist_filename()
    elif o == '-o' :
      mode = 'output'
      if os.path.exists(a) :
        print 'Output file (' + a + ') already exists, aborting.'
        sys.exit(1)
      else:
        output_filename = a
    elif o == '-c' :
      pass
    else :   # -h
      usage()
      sys.exit(1)


# Find all binaries that could be pulling in external libraries
def get_sage_library_list() :
  # At the moment, only look for .so and .dylib files
  found_sage_libraries = os.popen('find '+sage_root+' -name *.so -or -name *.dylib').readlines()
  return [i.strip() for i in found_sage_libraries]


# For removing the trailing (compatiblity version .....) stuff from otool -L
brace_pattern = re.compile('\(.*\)')
def trim_linked_file(x) :
  return re.sub(brace_pattern, '', x).strip()


# Is a file external? Yes, unless it starts with SAGE_ROOT, or is a relative path
nonlocal_file_pattern = re.compile('^(?!'+sage_root+')/')
def is_external(filename) :
  return re.match(nonlocal_file_pattern, filename)

def external_libraries_for_file(filename) :
  # Find what the binary links to
  linked_file_list = os.popen('otool -L ' + filename).readlines()
  # Trim, and ignore first result (it's the binary name)
  linked_file_list = [trim_linked_file(f) for f in linked_file_list[1:]]
  # Return non-local files only
  return [f for f in linked_file_list if is_external(f)]


# -w or -o file end up here.
def run_output_mode() :
  sage_libraries = get_sage_library_list()
  ext_linked_libs = [ external_libraries_for_file(i) for i in sage_libraries ]
  # Remove duplicates by adding them to a set
  ext_linked_libs = reduce(lambda p,q: p.union(q), ext_linked_libs, set())
  # Add \n to each
  ext_linked_libs = [ i + '\n' for i in ext_linked_libs ]
  ext_linked_libs.sort()
  # Write it out. output_filename is already set.
  output_file = open(output_filename, 'w')
  output_file.writelines(ext_linked_libs)
  output_file.close()

# default mode is here
def run_check_mode() :
  found_non_whitelisted = False
  if not(os.path.exists(whitelist_filename())) :
    print 'Whitelist file (' + whitelist_filename() + ') not found.  Use -w to generate one.'
    sys.exit(1)
  whitelist = open(whitelist_filename()).readlines()
  whitelist = [i.strip() for i in whitelist]
  for lib in get_sage_library_list() :
    ext_linked_libs = external_libraries_for_file(lib)
    for extlib in ext_linked_libs :
      if not(extlib in whitelist) :
        print lib + ' links to non-whitelisted file ' + extlib
        found_non_whitelisted = True
  if found_non_whitelisted : sys.exit(1)

def main() :
  startup()
  if mode == 'output':
    run_output_mode()
  else:
    run_check_mode()
  sys.exit(0)

if __name__ == "__main__":
    main()
