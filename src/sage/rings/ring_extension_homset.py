#############################################################################
#    Copyright (C) 2019 Xavier Caruso <xavier.caruso@normalesup.org>
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 2 of the License, or
#    (at your option) any later version.
#                  http://www.gnu.org/licenses/
#****************************************************************************

from sage.rings.homset import RingHomset_generic
from sage.rings.ring_extension import RingExtension_generic
from sage.rings.ring_extension_morphism import RingExtensionHomomorphism

class RingExtensionHomset(RingHomset_generic):
    def __call__(self, *args, **kwargs):
        return RingExtensionHomomorphism(self, *args, **kwargs)

    def _coerce_impl(self, x):
        if isinstance(x, RingExtensionHomomorphism):
            x = x._backend()
        domain = self.domain()
        if isinstance(domain, RingExtension_generic):
            domain = domain._backend()
        codomain = self.codomain()
        if isinstance(codomain, RingExtension_generic):
            codomain = codomain._backend()
        if domain is x.domain() and codomain is x.codomain():
            return RingExtensionHomomorphism(self, x)
        raise TypeError
