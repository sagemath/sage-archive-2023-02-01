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

.. _sec-launching-wsl-post-installation:

WSL Post-installation steps
---------------------------

If you've installed Sage Math from source on WSL, there are a couple of extra steps you can do to make your life easier:

Create a notebook launch script
"""""""""""""""""""""""""""""""

If you plan to use JupyterLab, install that first.

Now create a script called ``~/sage_nb.sh`` containing the following lines, and fill in the correct paths for your desired starting directory and ``SAGE_ROOT``


.. CODE-BLOCK:: bash

    #!/bin/bash
    # Switch to desired windows directory
    cd /mnt/c/path/to/desired/starting/directory
    # Start the Jupyter notebook
    SAGE_ROOT/sage --notebook
    # Alternatively you can run JupyterLab - delete the line above, and uncomment the line below
    #SAGE_ROOT/sage --notebook jupyterlab

Make it executable:

.. CODE-BLOCK:: bash

    chmod ug+x ~/sage_nb.sh

Run it to test:

.. CODE-BLOCK:: bash

    cd ~
    ./sage_nb.sh

The Jupyter(Lab) server should start in the terminal window, and you windows browser should open a page showing the Jupyter or JupyterLab starting page, at the directory you specified.

Create a shortcut
"""""""""""""""""

This is a final nicety that lets you start the Jupyter or JupyterLab server in one click:

* Open Windows explorer, and type ``%APPDATA%\Microsoft\Windows\Start Menu\Programs`` in the address bar and press enter. This is the folder that contains you start menu shortcuts. If you want the sage shortcut somewhere else (like your desktop), open that folder instead.
* Open a separate window and go to ``%LOCALAPPDATA%\Microsoft\WindowsApps\``
* Right-click-drag the ``ubuntu.exe`` icon from the second window into the first, then choose ``Create shortcuts here`` from the context menu when you drop it. 
* To customize this shortcut, right-click on it and choose properties.

  * On the General tab:

    * Change the name to whatever you want, e.g. "Sage 9.2 JupyterLab"
  
  * On the Shortcut tab:

    * Change Target to: ``ubuntu.exe run ~/sage_nb.sh``
    * Change Start in to: ``%USERPROFILE%``
    * Change Run to: Minimised
    * Change the icon if you want

Now hit the start button or key and type the name you gave it. it should appear in the list, and should load the server and fire up your browser when you click on it.

------------------------------------------------------------------------

For further reading you can have a look at the other documents in the
SageMath documentation at http://doc.sagemath.org/.

