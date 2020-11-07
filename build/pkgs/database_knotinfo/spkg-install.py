import os
from sage.all import save
from sage.env import SAGE_SHARE
from sage.misc.misc import sage_makedirs
from sage.databases.knotinfo_db import KnotInfoDataBase

install_root = os.path.join(SAGE_SHARE, 'knotinfo')

if __name__ == '__main__':
    sage_makedirs(install_root)
    print("Creating the KnotInfo database.")
    ki_db = KnotInfoDataBase()
    ki_db.create_col_dict_sobj()
    ki_db.create_data_sobj()
    os.system('cp package-version.txt %s' %install_root)
