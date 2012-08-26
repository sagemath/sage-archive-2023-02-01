#/bin/bash

if [ $# -ge 1 ]; then
    VERSION=$1
else
    VERSION="UNCOMMITTED"
fi

if [ ! -e .git ]; then
    echo "Must be called from the git repository."
    echo
    echo "git clone git://github.com/robertwb/sage-patchbot.git"
    exit 1
fi

rm -rf workspace-*
TMP=$(mktemp -d workspace-XXXXXX)
ORIGINAL=$(pwd)

git diff --stat

if [ "$VERSION" == "UNCOMMITTED" ]; then
    git diff > $TMP/uncommitted.patch
    HEAD=$(git rev-parse HEAD)
else
    status=$(sage -hg status)
    if [ -n "$status" ]; then
	echo "Uncommitted hg changes."
	echo "$status"
	exit 1
    fi
    HEAD=$(git rev-parse $VERSION)
    if [ "$?" -ne "0" ]; then
	echo "Unknown tag or commit: $VERSION"
	exit 1
    fi
fi

# Clone the repo.
cd $TMP
git clone $ORIGINAL patchbot-$VERSION
cd patchbot-$VERSION
if [ -e "../uncommitted.patch" ]; then
    patch -p1 < ../uncommitted.patch
fi

# Format as an spkg.
rm -rf .git*
cp -r $ORIGINAL/.hg .
mv README.txt src
cd ..

sage -spkg patchbot-$VERSION
cp patchbot-$VERSION.spkg $ORIGINAL

cd $ORIGINAL
rm -rf $TMP
