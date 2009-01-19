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



if [ `uname` = "Linux" ]; then
    if [ `./fortran_type.pl` =  "g95" ]; then
	cd "$SAGE_LOCAL"/lib
	lapack_command="ld -L"$f95_dir" -L"$SAGE_LOCAL"/lib  -shared -soname liblapack.so -o liblapack.so  --whole-archive liblapack.a --no-whole-archive -lc -lm -lf95"
	f77blas_command="ld -L"$f95_dir" -L"$SAGE_LOCAL"/lib  -shared -soname libf77blas.so -o libf77blas.so  --whole-archive libf77blas.a --no-whole-archive -lc -lm -lf95"

	echo $lapack_command
	$lapack_command
	echo $f77blas_command
	$f77blas_command

    else
        echo "Not on Linux, so nuking liblapack.so to work around problems with Numpy and Scipy"
	cd "$SAGE_LOCAL"/lib
        rm -rf liblapack.so
    fi

fi

