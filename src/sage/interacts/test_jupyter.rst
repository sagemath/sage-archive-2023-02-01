Test the Sage interact library in Jupyter

We need to setup a proper test environment for widgets::

    sage: from ipywidgets.widgets.tests import setup_test_comm
    sage: setup_test_comm()

Check that the interact library is available within Jupyter::

    sage: from sage.repl.ipython_kernel.all_jupyter import *
    sage: interacts.geometry.special_points()
    Interactive function <function special_points at ...> with 10 widgets
      title: HTML(value=u'<h2>Special points in triangle</h2>', description=u'title')
      a0: IntSlider(value=30, min=0, max=360, step=1, description=u'A')
      a1: IntSlider(value=180, min=0, max=360, step=1, description=u'B')
      a2: IntSlider(value=300, min=0, max=360, step=1, description=u'C')
      show_median: Checkbox(value=False, description=u'Medians')
      show_pb: Checkbox(value=False, description=u'Perpendicular Bisectors')
      show_alt: Checkbox(value=False, description=u'Altitudes')
      show_ab: Checkbox(value=False, description=u'Angle Bisectors')
      show_incircle: Checkbox(value=False, description=u'Incircle')
      show_euler: Checkbox(value=False, description=u"Euler's Line")
