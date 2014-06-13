from sage.misc.package import is_package_installed
from cooperative_game import CooperativeGame
if is_package_installed('gambit'):
    from normal_form_game import *
