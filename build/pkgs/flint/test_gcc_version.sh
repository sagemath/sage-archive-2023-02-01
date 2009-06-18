#!/bin/sh

GCC_VERSION=`gcc -dumpversion`

case $GCC_VERSION in
    3.4*)
        echo "Found gcc 3.4.x"
	exit 0
	;;
    3.*)
        echo "WARNING: gcc version less than 3.4.0"
	exit 1
	;;
    2.*)
        echo "WARNING: gcc version less than 3.4.0"
	exit 1
	;;
    1.*)
        echo "WARNING: gcc version less than 3.4.0"
	exit 1
	;;
    *)
        echo "Found gcc 4 or later"
	exit 0
	;;
esac

