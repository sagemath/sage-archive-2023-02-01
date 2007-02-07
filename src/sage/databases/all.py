"""nodoctest
This file gathers together all the tables in SAGE.

    * ConwayPolynomials -- database of selected Conway polynomials.

    * CremonaDatabase() - Cremona's tables of elliptic curves and related data.

    * Gamma0Wt2Database() -- table of arithmetic informationa about
                    newforms of weight 2 on Gamma_0(N).

    * JonesDatabase() -- returns the John Jones table of number fields
                with bounded ramification and degree <= 6.

    * SteinWatkinsDatabase() - The Stein-Watkins tables of elliptic curves
                and related data.

    * Sloane's tables -- sloane_sequence, sloane_find

    * Linear codes -- linear_code_bound

    * Symbolic Data -- benchmark and test ideals

EXAMPLES:
    sage: ConwayPolynomials()
    Frank Luebeck's database of Conway polynomials

    sage: CremonaDatabase()
    Cremona's database of elliptic curves

    sage: Gamma0Wt2Database()
    Table of arithmetic informationa about newforms of weight 2 on Gamma_0(N)

    sage: JonesDatabase()
    John Jones's table of number fields with bounded ramification and degree <= 6

    sage: SteinWatkinsDatabase()
    The Stein-Watkins database of elliptic curves

    sage: SymbolicData()
    SymbolicData with 372 ideals

"""

#*****************************************************************************
#       SAGE: System for Algebra and Geometry Experimentation
#
#       Copyright (C) 2005 William Stein <wstein@gmail.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#    This code is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
#    General Public License for more details.
#
#  The full text of the GPL is available at:
#
#                  http://www.gnu.org/licenses/
#*****************************************************************************



from conway import ConwayPolynomials

from cremona import CremonaDatabase, \
     cremona_letter_code, parse_cremona_label, \
     old_cremona_letter_code, is_optimal_id

from gamma0wt2 import Gamma0Wt2Database

from jones import JonesDatabase

from stein_watkins import SteinWatkinsAllData, SteinWatkinsPrimeData

from install import database_install

from sloane import sloane_sequence, sloane_find, SloaneEncyclopedia

from lincodes import linear_code_bound

from odlyzko import zeta_zeros

from db_modular_polynomials import \
     ClassicalModularPolynomialDatabase, \
     DedekindEtaModularPolynomialDatabase, \
     DedekindEtaModularCorrespondenceDatabase, \
     AtkinModularPolynomialDatabase, \
     AtkinModularCorrespondenceDatabase

from db_class_polynomials import \
     HilbertClassPolynomialDatabase

from symbolic_data import SymbolicData
