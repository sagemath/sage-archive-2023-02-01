.PHONY: all sage clean

all: sage

# SAGE_BUILD_DIR set by src/bin/sage-env
ifndef SAGE_BUILD_DIR
$(error SAGE_BUILD_DIR not set. Use top-level Makefile!)
endif
# setup.py also depends on SAGE_CYTHONIZED, defined in src/sage/env.py
# TODO: It would be better if it inferred the location from the build-base instead of using SAGE_CYTHONIZED.
sage: sage/libs/pari/auto_gen.pxi sage/ext/interpreters/__init__.py
	python -u setup.py build --build-base=$(SAGE_BUILD_DIR)/sagelib install

clean:
	@echo "Deleting Sage library build artifacts..."
	rm -rf c_lib build .cython_version # from old sage versions
	rm -rf $(SAGE_BUILD_DIR)/sagelib
	find . -name '*.pyc' | xargs rm -f
	rm -f sage/libs/pari/auto_*
	rm -rf sage/ext/interpreters


# Auto-generated files
# TODO: These belong in SAGE_BUILD_DIR as well.
sage/libs/pari/auto_gen.pxi: $(SAGE_LOCAL)/share/pari/pari.desc \
        sage/libs/pari/decl.pxi sage_setup/autogen/pari/*.py
	python -c "from sage_setup.autogen.pari import rebuild; rebuild()"

sage/ext/interpreters/__init__.py: sage_setup/autogen/interpreters.py
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"
