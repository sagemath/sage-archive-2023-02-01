# This Makefile is called by setup.py
# to build some generated Python/Cython sources.

all: sage/libs/cypari2/auto_gen.pxi sage/ext/interpreters/__init__.py

# Auto-generated files
# TODO: Adjustments for VPATH builds will be necessary.
sage/libs/cypari2/auto_gen.pxi: $(SAGE_LOCAL)/share/pari/pari.desc \
        sage_setup/autogen/pari/*.py
	rm -f sage/libs/pari/auto_*.pxi
	python -c "from sage_setup.autogen.pari import rebuild; rebuild()"

sage/ext/interpreters/__init__.py: sage_setup/autogen/interpreters.py
	python -c "from sage_setup.autogen.interpreters import rebuild; rebuild('sage/ext/interpreters')"
