.. _sec-installation-conda:

Install from conda-forge
========================

SageMath can be installed from `conda-forge <https://conda-forge.org>`_ on Linux
and macOS running x86-64 that most current desktops and laptops use.

To install SageMath, install `Miniconda <https://conda.io/miniconda.html>`_ and
then type in the following commands in a terminal:

* Add the conda-forge channel: ``conda config --add channels conda-forge``
* Update all packages: ``conda update --all``
* Create a new environment containing SageMath: ``conda create -n sage sage python=X``, where
  ``X`` is version of Python, e.g. ``2.7``
* Enter the new environment: ``conda activate sage``
* Start SageMath: ``sage``

.. note:: The dependency resolution process of conda is slow; a faster experimental resolver
   we tried with success is `mamba <https://github.com/QuantStack/mamba>`_. It can be used as follows.

   .. code-block:: shell

      conda install mamba -c conda-forge # installs mamba
      mamba create -n sage sage -c conda-forge # replaces "conda create..."

.. note:: At the time of this edit (2019/06/12) the newest SageMath available on conda-forge
   is 8.3, and it needs to be installed with the correct dependencies, as follows:

   .. code-block:: shell

      conda create -n sage sage python=2.7 -c conda-forge/label/cf201901 # or
      mamba create -n sage sage python=2.7 -c conda-forge/label/cf201901 # same, using mamba

