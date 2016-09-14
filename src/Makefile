.PHONY: all sage clean

all: sage

# SAGE_BUILD_DIR set by src/bin/sage-env
ifndef SAGE_BUILD_DIR
$(error SAGE_BUILD_DIR not set. Use top-level Makefile!)
endif
# setup.py also depends on SAGE_CYTHONIZED, defined in src/sage/env.py
# TODO: It would be better if it inferred the location from the build-base instead of using SAGE_CYTHONIZED.
sage:
	python -u setup.py build --build-base=$(SAGE_BUILD_DIR)/sagelib install

clean:
	@echo "Deleting Sage library build artifacts..."
	rm -rf c_lib build .cython_version # from old sage versions
	rm -rf $(SAGE_BUILD_DIR)/sagelib
	find . -name '*.pyc' | xargs rm -f
	rm -f sage/libs/pari/auto_*
	rm -rf sage/ext/interpreters


