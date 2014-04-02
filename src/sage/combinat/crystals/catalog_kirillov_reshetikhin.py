"""
Catalog Of Crystal Models For Kirillov-Reshetikhin Crystals

We currently have the following models:

* :func:`~sage.combinat.crystals.kirillov_reshetikhin.KashiwaraNakashimaTableaux`
* :func:`LSPaths <sage.combinat.crystals.kirillov_reshetikhin.KirillovReshetikhinCrystalFromLSPaths>`
* :class:`~sage.combinat.rigged_configurations.kr_tableaux.KirillovReshetikhinTableaux`
* :func:`~sage.combinat.crystals.catalog_kirillov_reshetikhin.RiggedConfigurations`
"""
from kirillov_reshetikhin import KirillovReshetikhinCrystal as KashiwaraNakashimaTableaux
from kirillov_reshetikhin import KirillovReshetikhinCrystalFromLSPaths as LSPaths
from sage.combinat.rigged_configurations.kr_tableaux import KirillovReshetikhinTableaux

def RiggedConfigurations(cartan_type, r, s):
    """
    Return the KR crystal `B^{r,s}` using
    :class:`rigged configurations <RiggedConfigurations>`.

    This is the set `RC(B^{r,s})` or with the
    `L_i^{(a)} = \delta_{ar} \delta_{is}`.

    EXAMPLES::

        sage: K1 = crystals.kirillov_reshetikhin.RiggedConfigurations(['A',6,2], 2, 1)
        sage: K2 = crystals.kirillov_reshetikhin.LSPaths(['A',6,2], 2, 1)
        sage: K1.digraph().is_isomorphic(K2.digraph(), edge_labels=True)
        True
    """
    from sage.combinat.rigged_configurations.rigged_configurations import RiggedConfigurations
    return RiggedConfigurations(cartan_type, [[r,s]])

