from sage.categories.map import Map
from sage.categories.homset import Hom
from sage.categories.sets_cat import Sets

class Valuation(Map):
    def __init__(self, domain, codomain):
        Map.__init__(self, Hom(domain, codomain, Sets()))
