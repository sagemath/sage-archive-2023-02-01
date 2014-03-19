"""
Catalog Of Crystal Models For (Tensor Products Of) Kirillov-Reshetikhin
Crystals
"""
from kirillov_reshetikhin import KirillovReshetikhinCrystal as KashiwaraNakashimaTableaux
from kirillov_reshetikhin import KirillovReshetikhinCrystalFromLSPaths as LSPaths
from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux

def RiggedConfigurations(cartan_type, r, s):
    """
    Return the KR crystal `B^{r,s}` using rigged configurations.

    EXAMPLES::

        sage: K1 = crystals.kirillov_reshetikhin.RiggedConfigurations(['A',6,2], 2, 1)
        sage: K2 = crystals.kirillov_reshetikhin.LSPaths(['A',6,2], 2, 1)
        sage: K1.digraph().is_isomorphic(K2.digraph(), edge_labels=True)
        True
    """
    from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
    return RiggedConfigurations(cartan_type, [[r,s]])

