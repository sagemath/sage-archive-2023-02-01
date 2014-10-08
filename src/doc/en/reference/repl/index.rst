The Sage Command Line
=====================

The Sage Read-Eval-Print-Loop (REPL) is based on IPython. In this
document, you'll find how the IPython integration works. You should
also be familiar with the documentation for IPython. 

For more details about using the Sage command line, see the Sage
tutorial.


.. toctree::
   :maxdepth: 2

   options
   startup
   environ
   sage/misc/trace

   sage/repl/readline_extra_commands
   sage/repl/interpreter
   sage/repl/ipython_extension


Pretty Printing
---------------

In addition to making input nicer, we also modify how results are
printed. This again builds on how IPython formats output. Technically,
this works using a modified displayhook in Python.

.. toctree::
   :maxdepth: 2
   
   sage/repl/display/formatter
   sage/repl/display/pretty_print
   sage/repl/display/fancy_repr
   sage/repl/display/util
   sage/repl/display/python_hook
   



.. include:: ../footer.txt
