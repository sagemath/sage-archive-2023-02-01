At present Gambit is not an optional package (which is the first step to becoming a full part of Sage).

To install Gambit so that it is part of your Sage package follow the following:

    $ cd SAGE_ROOT/upstream
    $ wget http://sourceforge.net/projects/gambit/files/gambit13/13.1.2/gambit-13.1.2.tar.gz/download
    $ cd SAGE_ROOT
    $ ./sage -i gambit

The above steps should download the gambit 13.1.2 tarbell and then use Sage to unpacke as necessary.

To check that this is worked, in Sage simply type `import gambit`.
