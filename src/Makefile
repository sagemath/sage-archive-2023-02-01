.PHONY: all sage clean

all: sage

# SAGE_VERSION set by build/make/install.
# SAGE_BUILD_DIR set by build/bin/sage-spkg
sage: sage/libs/pari/auto_gen.pxi sage/ext/interpreters/__init__.py
	python -u setup.py --build-base=$SAGE_BUILD_DIR/sagelib-$SAGE_VERSION install

clean:
	@echo "Deleting Sage library build artifacts..."
	rm -rf c_lib build
	find . -name '*.pyc' | xargs rm -f
	rm -f sage/libs/pari/auto_*
	rm -rf sage/ext/interpreters


# Auto-generated files
sage/libs/pari/auto_gen.pxi: $(SAGE_LOCAL)/share/pari/pari.desc \
        sage/libs/pari/decl.pxi sage_setup/autogen/pari/*.py
	python -c "from sage_setup.autogen.pari import rebuild; rebuild()"

sage/ext/interpreters/__init__.py: sage_setup/autogen/interpreters.py
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"
