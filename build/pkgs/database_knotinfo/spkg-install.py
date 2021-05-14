import os
from sage.env import SAGE_SHARE
from sage.misc.misc import sage_makedirs

install_root = os.path.join(SAGE_SHARE, 'knotinfo')

if __name__ == '__main__':
    sage_makedirs(install_root)
    print("Creating the KnotInfo database.")
    from sage.databases.knotinfo_db import KnotInfoDataBase
    KnotInfoDataBase(install=True)
    os.system('cp package-version.txt %s' %install_root)
