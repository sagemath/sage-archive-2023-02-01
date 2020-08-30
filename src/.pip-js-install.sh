#! /bin/sh -x
set -e
python -m pip install "$@"
if [ -n "$JS_DEPS" ]; then
    nodeenv --python-virtualenv --verbose
    npm install -g --no-package-lock --no-save --verbose $JS_DEPS
fi
