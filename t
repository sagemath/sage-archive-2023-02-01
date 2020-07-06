if test x`locale -a | grep C\.UTF-8` != x; then
 export LC_ALL=C.UTF-8;
else
 export LC_ALL=C;
fi

echo $LC_ALL
