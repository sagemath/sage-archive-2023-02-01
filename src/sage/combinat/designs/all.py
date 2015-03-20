"""
Combinatorial design features that are imported by default in the interpreter namespace
"""
from block_design import (BlockDesign)

from ext_rep import (designs_from_XML, designs_from_XML_url)

from incidence_structures import (IncidenceStructure,
                          IncidenceStructureFromMatrix)

from incidence_structures import IncidenceStructure as Hypergraph

from covering_design import (CoveringDesign,
                             schonheim,
                             trivial_covering_design)

import design_catalog as designs

import sage.combinat.designs.steiner_quadruple_systems
import sage.combinat.designs.orthogonal_arrays
import sage.combinat.designs.latin_squares
import sage.combinat.designs.bibd
import sage.combinat.designs.group_divisible_designs
