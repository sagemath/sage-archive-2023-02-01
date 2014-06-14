At present Gambit is not an optional package (which is the first step to becoming a full part of Sage).

To install Gambit so that you can use NormalFormGame functionality, follow the steps below.
Make sure that you are using the most recent recent branch form this repository.

    git clone --single-branch --branch sage_integration https://github.com/tturocy/gambit
    cd gambit
    make dist
    mv gambit-14.0.2.tar.gz SAGE_ROOT/upstream
    cd SAGE_ROOT
    ./sage -i gambit

The above steps should download the relative Gambit branch, create a tarbell and then add it to sage.

To check that this is worked, in Sage simply type `import gambit`.
