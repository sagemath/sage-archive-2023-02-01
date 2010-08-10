#/bin/sh

##############Find the location of f95##########
if [ `./fortran_type.pl` = "g95" ]; then
f95_dir=$SAGE_LOCAL/lib/gcc-lib
f95_dir=$f95_dir/`ls $f95_dir`
f95_dir=$f95_dir/`ls $f95_dir`
echo "Using f95, f95 lib in"
echo $f95_dir
fi
###############################################




if [ `uname` = "Darwin" ]; then
    cd "$SAGE_LOCAL"/lib
    atlas_command="ld  -L"$SAGE_LOCAL"/lib  -dylib  -o libatlas.dylib -lm -lc -all_load libatlas.a"
    cblas_command="ld  -L"$SAGE_LOCAL"/lib -dylib -o libcblas.dylib -lm -lc -latlas -all_load libcblas.a"
    echo $atlas_command
    $atlas_command
    echo $cblas_command
    $cblas_command

fi



if [ `uname` = "Linux" ] || [ `uname` = "FreeBSD" ]; then
    if [ `./fortran_type.pl` =  "g95" ]; then
	cd "$SAGE_LOCAL"/lib
	lapack_command="ld -L"$f95_dir" -L"$SAGE_LOCAL"/lib  -shared -soname liblapack.so -o liblapack.so  --whole-archive liblapack.a --no-whole-archive -lc -lm -lf95"
	f77blas_command="ld -L"$f95_dir" -L"$SAGE_LOCAL"/lib  -shared -soname libf77blas.so -o libf77blas.so  --whole-archive libf77blas.a --no-whole-archive -lc -lm -lf95"

	echo $lapack_command
	$lapack_command
	echo $f77blas_command
	$f77blas_command

    else
	cd "$SAGE_LOCAL"/lib
	lapack_command="ld -L"$SAGE_LOCAL"/lib  -shared -soname liblapack.so -o liblapack.so  --whole-archive liblapack.a --no-whole-archive -lc -lm -lgfortran"
	f77blas_command="ld -L"$SAGE_LOCAL"/lib -shared -soname libf77blas.so -o libf77blas.so  --whole-archive libf77blas.a --no-whole-archive -lc -lm -lgfortran"
	echo $lapack_command
	$lapack_command
	echo $f77blas_command
	$f77blas_command
    fi

fi

# Build dynamic libraries on Solaris, using the Sun equivalent of the GNU
# options above. Rather than call 'ld', the exact path to 'ld' is given
# since we know for 100% sure that 'ld' will reside in /usr/ccs/bin
# on Solaris - this is true of both Solaris 10 and OpenSolaris
# I'm not 100% sure if these options are optimal (there might be some
# unnecessary commands, but its the exact equivalent of the GNU
# commands above.

if [ "x`uname`" = xSunOS ]; then
   cd "$SAGE_LOCAL/lib"

   if [ "x$SAGE64" = xyes ] ; then
      # To create a 64-bit shared library, the linker flag
      # '-64' must be added. Note this is not the same as
      # the compiler flag '-m64'
      LINKER_FLAG=-64
   fi

   # Build libatlas.so libf77blas.so and libcblas.so
   # Actually, libatlas.so and libcblas.so gets built from spkg-install,
   # but this is not working correctly for 64-bit libraries. Hence I simply
   # remake those two libraries here. This is why 4 libraries are built
   # in this Solaris section, but only two are shown above for FreeBSD,
   # OS X and Linux.

   for ATLAS_LIBRARY in libatlas libf77blas libcblas ; do
      echo "Building shared library $ATLAS_LIBRARY.so from the static library $ATLAS_LIBRARY.a"
     /usr/ccs/bin/ld $LINKER_FLAG -L"$SAGE_LOCAL/lib"  -G -h $ATLAS_LIBRARY.so -o $ATLAS_LIBRARY.so  -zallextract  $ATLAS_LIBRARY.a -zdefaultextract -lc -lm -lgfortran
      if [ $? -ne 0 ]; then
         echo "Failed to build ATLAS library $ATLAS_LIBRARY.so"
         exit 1
      else
         echo "$ATLAS_LIBRARY.so has been built on Solaris."
      fi
   done
   rm liblapack.so # liblapack.so causes problems with R on Solaris.
fi

