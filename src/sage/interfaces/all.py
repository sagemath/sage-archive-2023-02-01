# interfaces to other interpreters

from frobby import frobby
from four_ti_2 import four_ti_2
from axiom import Axiom, axiom, axiom_console
from fricas import FriCAS, fricas, fricas_console

from expect import is_ExpectElement
from gap import gap, gap_reset_workspace, gap_console, gap_version, is_GapElement, Gap
from gap3 import gap3, gap3_console, gap3_version, Gap3
from genus2reduction import genus2reduction, Genus2reduction
from gfan import gfan,  Gfan
from giac import giac, giac_console, Giac
from gp import gp, gp_console, gp_version, is_GpElement, Gp
from gnuplot import gnuplot, gnuplot_console
from kash import  kash, kash_console, kash_version, is_KashElement, Kash
from lisp import lisp, lisp_console, Lisp
from magma import magma, magma_console, magma_version, Magma, is_MagmaElement
from magma_free import magma_free
from macaulay2 import macaulay2, macaulay2_console, Macaulay2
from maple import maple, maple_console, Maple
from maxima_abstract import maxima_console
from maxima import maxima, is_MaximaElement, Maxima
# import problems
#from maxima_lib import maxima_lib
from mathematica import mathematica, mathematica_console, Mathematica
from matlab import matlab, matlab_console, matlab_version, Matlab
from mupad import mupad, mupad_console, Mupad  # NOT functional yet
from mwrank import mwrank, Mwrank, mwrank_console
from octave import octave, octave_console, octave_version, Octave
from qepcad import qepcad, qepcad_console, qepcad_version, qepcad_formula
from qsieve import qsieve
from singular import singular, singular_console, singular_version, is_SingularElement, Singular
from sage0 import sage0 as sage0, sage0_console, sage0_version, Sage
from scilab import scilab
from tachyon import tachyon_rt
from psage import PSage
from ecm import ECM, ecm
from povray import povray
from lie import lie, lie_console, LiE
from r import r, r_console, R, r_version, is_RElement
from read_data import read_data

interfaces = ['gap', 'gap3', 'giac', 'gp', 'mathematica', 'gnuplot', \
              'kash', 'magma', 'macaulay2', 'maple', 'maxima', \
              'mathematica', 'mwrank', 'octave', 'r', \
              'singular', 'sage0', 'sage']
