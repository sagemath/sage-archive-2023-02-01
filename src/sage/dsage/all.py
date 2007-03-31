"""nodoctest
"""
from sage.dsage.dsage import dsage
from sage.dsage.dist_functions.all import *
# from sage.dsage.dist_functions.dist_factor import DistributedFactor

def DSage(server=None, port=8081, username=None,
          pubkey_file=None, privkey_file=None):
    from sage.dsage.interface.dsage_interface import BlockingDSage
    return BlockingDSage(server=server, port=port, username=username,
                         pubkey_file=pubkey_file, privkey_file=privkey_file)

# def DistributedFactor(dsage, n):
#     from sage.dsage.dist_functions.dist_factor import DistributedFactor
#     return DistributedFactor(dsage, n)

