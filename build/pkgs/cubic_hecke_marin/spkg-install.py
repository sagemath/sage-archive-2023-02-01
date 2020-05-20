import os
from sage.all import save
from sage.env import SAGE_SHARE
from sage.misc.misc import sage_makedirs
from sage.databases.cubic_hecke_db import CubicHeckeDataBase

install_root = os.path.join(SAGE_SHARE, 'cubic_hecke_marin')

if __name__ == '__main__':
    sage_makedirs(install_root)
    print("Creating Iwan Marin's Cubic Hecke database.")
    cha_db = CubicHeckeDataBase()
    cha_db.create_static_db_marin_basis()
    cha_db.create_static_db_marin_regular()
    cha_db.create_static_db_marin_regular(right=True)
    cha_db.create_static_db_marin_split()
    os.system('cp package-version.txt %s' %install_root)
