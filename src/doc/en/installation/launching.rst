.. _sec-launching:

Launching SageMath
==================

Now we assume that you installed SageMath properly on your system. This
section quickly explains how to start the Sage console and the Jupyter
Notebook from the command line.

If you did install the Windows version or the macOS application you
should have icons available on your desktops or launching menus. Otherwise
you are strongly advised to create shortcuts for Sage as indicated at the end
of the "Linux" Section in :ref:`sec-installation-from-binaries`. Assuming
that you have this shortcut, running

.. CODE-BLOCK:: bash

    sage

in a console starts a Sage session.  To quit the session enter ``quit`` and
then press ``<Enter>``.

To start a Jupyter Notebook instead of a Sage console, run the command

.. CODE-BLOCK:: bash

    sage -n jupyter

instead of just ``sage``. To quit the Jupyter Notebook press ``<Ctrl> + <c>``
twice in the console where you launched the command.

Using a Jupyter Notebook remotely
---------------------------------

If Sage is installed on a remote machine to which you have ``ssh`` access, you
can launch a Jupyter Notebook using a command such as

.. CODE-BLOCK:: bash

    ssh -L localhost:8888:localhost:8888 -t USER@REMOTE sage -n jupyter --no-browser --port=8888

where ``USER@REMOTE`` needs to be replaced by the login details to the remote
machine. This uses local port forwarding to connect your local machine to the
remote one. The command will print a URL to the console which you can copy and
paste in a web browser.

Note that this assumes that a firewall which might be present between server
and client allows connections on port 8888. See details on port forwarding on
the internet, e.g. https://www.ssh.com/ssh/tunneling/example.

------------------------------------------------------------------------

For further reading you can have a look at the other documents in the
SageMath documentation at http://doc.sagemath.org/.


Setting up SageMath as a Jupyter kernel in an existing Jupyter notebook or JupyterLab installation
--------------------------------------------------------------------------------------------------

You may already have a global installation of Jupyter.  For added convenience, it is possible to link your installation of SageMath into your Jupyter installation, adding it to the list of available kernels that can be selected in the notebook or JupyterLab interface.

If ``$SAGE_LOCAL`` is the installation prefix of your Sage installation (the default is ``$SAGE_ROOT/local``)
and you can start the Jupyter notebook by typing ``jupyter notebook``, then the following command will install SageMath as a new kernel.

.. CODE-BLOCK:: bash

    jupyter kernelspec install --user $SAGE_LOCAL/share/jupyter/kernels/sagemath

This installs the kernel under the name ``sagemath``.  If you wish to rename it to something more specific in order to distinguish between different installations of SageMath, you can use the additional option ``--name``, for example

.. CODE-BLOCK:: bash

    jupyter kernelspec install --user $SAGE_LOCAL/share/jupyter/kernels/sagemath --name sagemath-dev-worktree

For the full functionality of the SageMath kernel in your global Jupyter installation, additionally some Notebook Extension packages need to be installed (or linked) into the environment from which the Jupyter installation runs.  You can check the presence of some of these packages using the command

.. CODE-BLOCK:: bash

    jupyter nbextension list

 - For the Sage interacts, you will need the package ``widgetsnbextension`` installed in the Python environment of the Jupyter installation.  If your Jupyter installation is coming from the system package manager, it is best to install ``widgetsnbextension`` in the same way.  Otherwise, install it using ``pip``.

..  - For 3D graphics using Three.js, ...............

..  - For 3D graphics using jsmol, ..............


