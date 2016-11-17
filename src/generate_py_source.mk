# This Makefile is called by setup.py
# to build some generated Python/Cython sources.

all: sage/libs/pari/auto_gen.pxi sage/ext/interpreters/__init__.py

# Auto-generated files
# TODO: Adjustments for VPATH builds will be necessary.
sage/libs/pari/auto_gen.pxi: $(SAGE_LOCAL)/share/pari/pari.desc \
        sage/libs/pari/decl.pxi sage_setup/autogen/pari/*.py
	python -c "from sage_setup.autogen.pari import rebuild; rebuild()"

sage/ext/interpreters/__init__.py: sage_setup/autogen/interpreters.py
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"
