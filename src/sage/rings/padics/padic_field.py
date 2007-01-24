class pAdicFieldGeneric(field.Field):
    def characteristic(self): pass

class pAdicField(pAdicFieldGeneric):
    def __call__(self, x, prec=infinity): pass
    def __cmp__(self, other): pass
    def __init__(self, p, default_prec = 20, series_print=True): pass
    def _coerce_(self, x): pass
    def _repr_(self): pass
    def absolute_discriminant(self): pass
    def automorphisms(self): pass
    def base_ring(self): pass
    def class_field(self, group=None, map=None, generators=None): pass
    def composite(self, subfield1, subfield2): pass
    def default_prec(self, n=None): pass
    def defining_polynomial(self): pass
    def degree(self): pass
    def different(self): pass
    def discriminant(self): pass
    def e(self): pass
    def f(self): pass
    def galois_group(self): pass
    def get_default_prec(self): pass
    def get_extension(self): pass
    def get_print_mode(self): pass
    def ground_ring(self): pass
    def has_GNB(self): pass
    def has_pth_root(self): pass
    def has_root_of_unity(self, n): pass
    def hom(self, field): pass
    def inertia_degree(self): pass
    def inertia_subfield(self): pass
    def integer_ring(self): pass
    def is_abelian(self): pass
    def is_isomorphic(self, K): pass
    def is_normal(self): pass
    def norm_equation(self): pass
    def norm_group(self): pass
    def norm_group_discriminant(self, group=None, map=None, generators=None): pass
    def number_of_extensions(self, degree, discriminant=None, e=None, f=None): pass
    def list_of_extensions(self, degree, discriminant=None, e=None, f=None): pass
    def prime(self): pass
    def principal_unit_group(self): pass
    def ramification_index(self): pass
    def ramification_filtration(self): pass
    def random(self): pass
    def residue_characteristic(self): pass
    def residue_class_field(self): pass
    def subfield(self, list): pass
    def subfield_lattice(self): pass
    def subfields_of_degree(self, n): pass
    def teichmuller(self, x): pass
    def teichmuller_system(self): pass
    def uniformizer(self): pass
    def unit_group(self): pass
    def unit_group_gens(self): pass

class pAdicFieldLazy(pAdicFieldGeneric):
    def __call__(self, x, prec=infinity): pass
    def __cmp__(self, other): pass
    def __init__(self, p, default_prec = 20, series_print=True): pass
    def _coerce_(self, x): pass
    def _repr_(self): pass
    def absolute_discriminant(self): pass
    def automorphisms(self): pass
    def base_ring(self): pass
    def class_field(self, group=None, map=None, generators=None): pass
    def composite(self, subfield1, subfield2): pass
    def default_prec(self, n=None): pass
    def defining_polynomial(self): pass
    def degree(self): pass
    def different(self): pass
    def discriminant(self): pass
    def e(self): pass
    def f(self): pass
    def galois_group(self): pass
    def get_default_prec(self): pass
    def get_extension(self): pass
    def get_print_mode(self): pass
    def ground_ring(self): pass
    def has_GNB(self): pass
    def has_pth_root(self): pass
    def has_root_of_unity(self, n): pass
    def hom(self, field): pass
    def inertia_degree(self): pass
    def inertia_subfield(self): pass
    def integer_ring(self): pass
    def is_abelian(self): pass
    def is_isomorphic(self, K): pass
    def is_normal(self): pass
    def norm_equation(self): pass
    def norm_group(self): pass
    def norm_group_discriminant(self, group=None, map=None, generators=None): pass
    def number_of_extensions(self, degree, discriminant=None, e=None, f=None): pass
    def list_of_extensions(self, degree, discriminant=None, e=None, f=None): pass
    def prime(self): pass
    def principal_unit_group(self): pass
    def ramification_index(self): pass
    def ramification_filtration(self): pass
    def random(self): pass
    def residue_characteristic(self): pass
    def residue_class_field(self): pass
    def subfield(self, list): pass
    def subfield_lattice(self): pass
    def subfields_of_degree(self, n): pass
    def teichmuller(self, x): pass
    def teichmuller_system(self): pass
    def uniformizer(self): pass
    def unit_group(self): pass
    def unit_group_gens(self): pass

Qp = pAdicField

def is_pAdicField(x):
    return isinstance(x, pAdicField_generic)

