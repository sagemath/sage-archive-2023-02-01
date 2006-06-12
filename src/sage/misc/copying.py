import os
import pager


class License:
    def __call__(self):
        pager.pager()(open(os.environ['SAGE_ROOT'] + '/COPYING.txt').read())

    def __repr__(self):
        return "Type license() to see the full license text."


license = License()

