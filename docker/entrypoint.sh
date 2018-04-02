#!/bin/bash
if [ x"$1" = jupyter ]; then
    # If "jupyter" is given as a first argument, we start a jupyter notebook
    # with reasonable default parameters for running it inside a container.
    shift
    exec sage -n jupyter --no-browser --ip='*' --port=8888 "$@"
else
    exec sage -sh -c "$*"
fi
