At present Gambit is not an optional package (which is the first step to becoming a full part of Sage).
The good guys over at [Gambit]() have made some tweaks to their code to enable easier integration in to Sage: that all lives in this git branch: [https://github.com/tturocy/gambit/tree/sage_integration](https://github.com/tturocy/gambit/tree/sage_integration).

To install Gambit so that it is part of your Sage package follow the following instructions.

- Getting Gambit:

    - If you already have a clone of the gambit repository, make sure you have the `sage-integration` branch checked out:

        git remote add turocy https://github.com/tturocy/gambit.git
        git checkout -b sage_integration
        git pull turocy sage_integration

    - If you don't have a clone of the gambit repository:

        git clone https://github.com/tturocy/gambit.git
        git checkout origin/sage_integration

- Getting Gambit in to Sage:

    cd GAMBIT_ROOT
    make dist
    cd SAGE_ROOT/upstream
    mv GAMBIT_ROOT/gambit-14.0.2.tar.gz ./gambit-13.1.2.tar.gz  # Note here that we are changing the name of the file (this is ugly and should be fixed)
    cd SAGE_ROOT
    ./sage -i gambit

The above steps should download the gambit 13.1.2 tarbell and then use Sage to unpacke as necessary.

To check that this is worked, in Sage simply type `import gambit`.
