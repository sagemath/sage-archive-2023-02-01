from block_design import (BlockDesign)

from ext_rep import (designs_from_XML, designs_from_XML_url)

from incidence_structures import (IncidenceStructure,
                          IncidenceStructureFromMatrix)

from covering_design import (CoveringDesign,
                             schonheim,
                             trivial_covering_design)

import design_catalog as designs

from sage.misc.superseded import deprecated_callable_import
deprecated_callable_import(14499,
                           'sage.combinat.designs.block_design',
                           globals(),
                           locals(),
                           ["ProjectiveGeometryDesign",
                            "AffineGeometryDesign",
                            "WittDesign",
                            "HadamardDesign",
                            "BlockDesign_generic"],
                           ("This method soon will not be available in that "
                            "way anymore. To use it, you can now call it by "
                            "typing designs.%(name)s"))
deprecated_callable_import(14499,
                           'sage.combinat.designs.covering_design',
                           globals(),
                           locals(),
                           ["best_known_covering_design_www"],
                           ("This method soon will not be available in that "
                            "way anymore. To use it, you can now call it by "
                            "typing designs.best_known_covering_design_from_LJCR"))

del deprecated_callable_import

import sage.combinat.designs.steiner_quadruple_systems
import sage.combinat.designs.latin_squares
