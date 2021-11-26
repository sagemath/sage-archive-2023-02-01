FROM gitpod/workspace-full

# Install system packages
RUN sudo apt-get update
## Needed for bootstrap
RUN sudo apt-get install -y gettext
## From dockerfile
RUN sudo apt-get install -y --no-install-recommends gfortran gcc g++ libstdc++-9-dev sudo openssl
RUN sudo apt-get install -y wget build-essential automake m4 dpkg-dev python libssl-dev rdfind
## Recommended by ./configure
### We do not install tox, since it pulls in javascript-common which does not install for some reason
RUN sudo apt-get install -y libflint-arb-dev libbrial-dev libbrial-groebner-dev libcdd-dev libcdd-tools ecl libec-dev eclib-tools fflas-ffpack flintqs libgc-dev gfan libgiac-dev xcas libgsl-dev libiml-dev lcalc liblfunction-dev libhomfly-dev libopenblas-dev palp pari-gp2c libpari-dev pari-doc pari-elldata pari-galdata pari-galpol pari-seadata libppl-dev ppl-dev python3 libpython3-dev python3-distutils r-base-dev r-cran-lattice librw-dev libsuitesparse-dev libsymmetrica2-dev sympow tachyon libzmq3-dev libzn-poly-dev python3-venv libgf2x-dev cliquer libcliquer-dev gmp-ecm libecm-dev glpk-utils libglpk-dev libbraiding-dev liblrcalc-dev libm4rie-dev libmpc-dev libmpfi-dev libmpfr-dev nauty libntl-dev libplanarity-dev planarity 
## Homebrew has more up-to-date package versions
### We do not install ecl from brew, since this breaks the build of maxima
RUN brew update && brew upgrade
RUN brew install arb flint fplll pari pari-elldata pari-galdata pari-galpol pari-seadata tox || true
### Give prio to brew over system packages
#Currently doesn't work as gitpod gp executable is then hidden by pari/gp
#ENV PATH="/home/linuxbrew/.linuxbrew/bin:$PATH"

# Configure 
## Gitpod sets PIP_USER: yes by default, which leads to problems during build (e.g pip not being installed in the venv)
RUN unset PIP_USER
## Gitpod installs pyenv by default, and sage's pip install targets the pyenv python for some reason
RUN pyenv global system
