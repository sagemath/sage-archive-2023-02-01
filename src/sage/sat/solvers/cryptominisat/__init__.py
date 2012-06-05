try:
    from cryptominisat import CryptoMiniSat
except ImportError:
    raise ImportError("Failed to import 'sage.sat.solvers.cryptominisat.CryptoMiniSat'. Run \"install_package('cryptominisat')\" to install it.")
