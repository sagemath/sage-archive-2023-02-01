from sage.misc.cachefunc import cached_method
from sage.structure.parent import Parent
from sage.categories.graded_algebras import GradedAlgebras
from sage.structure.unique_representation import UniqueRepresentation
from sage.symbolic.ring import SR
from sage.manifolds.differentiable.mixed_form import MixedForm

class MixedFormAlgebra(Parent, UniqueRepresentation):
    # Declare elements:
    Element = MixedForm

    #####
    # __init__
    #
    # Constructor
    #####

    def __init__(self, vector_field_module):
        if vector_field_module is None:
            raise ValueError("underlying vector field module must be provided")
        domain = vector_field_module._domain
        dest_map = vector_field_module._dest_map
        # Set name and latex_name:
        name = "Omega^*(" + domain._name
        latex_name = r"\Omega^*\left(" + domain._latex_name
        if dest_map is not domain.identity_map():
            dm_name = dest_map._name
            dm_latex_name = dest_map._latex_name
            if dm_name is None:
                dm_name = "unnamed map"
            if dm_latex_name is None:
                dm_latex_name = r"\mathrm{unnamed\; map}"
            name += "," + dm_name
            latex_name += "," + dm_latex_name
        self._name = name + ")"
        self._latex_name = latex_name + r"\right)"
        # Add this algebra to the category of graded algebras:
        base_field = domain.base_field()
        if domain.base_field_type() in ['real', 'complex']:
            base_field = SR
        Parent.__init__(self, base=base_field,
                        category=GradedAlgebras(base_field))
        # Define attributes:
        self._domain = domain
        self._ambient_domain = vector_field_module._ambient_domain
        self._dest_map = dest_map
        self._vmodule = vector_field_module
        self._max_deg = vector_field_module._ambient_domain.dim()

    #####
    # _repr_
    #
    # Return description
    #####

    def _repr_(self):
        description = "Graded algebra " + self._name + " of mixed differential forms "
        if self._dest_map is self._domain.identity_map():
            description += "on the {}".format(self._domain)
        else:
            description += "along the {} mapped into the {} ".format(
                self._domain, self._ambient_domain)
            if self._dest_map._name is None:
                dm_name = "unnamed map"
            else:
                dm_name = self._dest_map._name
            description += "via " + dm_name
        return description

    #####
    # _element_constructor_
    #
    # Construct the element MixedForm with coercion
    #####

    def _element_constructor_(self, comp=None, name=None, latex_name=None):
        if comp is None:
            return self.element_class(self, comp, name, latex_name)
        # Prepare list:
        if isinstance(comp, list):
            if len(comp) != self._max_deg + 1:
                raise IndexError(
                    "input list must have length {}".format(self._max_deg + 1))
            comp_list = comp
        elif isinstance(comp, self.Element):
            comp_list = comp[:]
        elif comp in self._domain.scalar_field_algebra():
            comp_list = [0] * (self._max_deg + 1)
            comp_list[0] = comp
        else:
            comp_list = comp_list = [0] * (self._max_deg + 1)
            # Try...except?
            deg = comp.degree()
            if deg <= self._max_deg:
                comp_list[deg] = comp
        # Use already existing coercions:
        comp_list = [self._domain.diff_form_module(j, self._dest_map).coerce(comp_list[j])
                     for j in range(0, self._max_deg + 1)]
        try:
            if name is None and comp._name is not None:
                name = comp._name
            if latex_name is None and comp._latex_name is not None:
                latex_name = comp._latex_name
        except AttributeError:
            pass
        return self.element_class(self, comp_list, name, latex_name)

    #####
    # _coerce_map_from_
    #
    #
    #####

    def _coerce_map_from_(self, S):
        if isinstance(S, self.__class__):
            # coercion by domain restriction
            return self._domain.is_subset(
                S._domain) and self._ambient_domain.is_subset(S._ambient_domain)
        # Test scalar_field_algebra separately to ensure coercion from reals etc.
        if self._domain.scalar_field_algebra().has_coerce_map_from(S):
            return True
        # This is tricky: To check whether S is coerciable to a DiffFormModule,
        # we need to determine its degree first.
        # But if there is no degree-method, it will raise an error. Let's try to catch that!
        try:
            deg = S.degree()
            if self._domain.diff_form_module(deg,
                                             self._dest_map).has_coerce_map_from(S):
                return True
        except (NotImplementedError, AttributeError, TypeError):
            pass
        return False

    #####
    # _repr_
    #
    # Return latex name
    #####

    def _latex_(self):
        return self._latex_name

    #####
    # zero
    #
    # Return zero element
    #####

    @cached_method
    def zero(self):
        resu_comp = [0] * (self._max_deg + 1)
        return self._element_constructor_(comp=resu_comp, name='zero',
                                          latex_name='0')

    #####
    # one
    #
    # Return one element
    #####

    @cached_method
    def one(self):
        resu_comp = [0] * (self._max_deg + 1)
        resu_comp[0] = 1
        return self._element_constructor_(comp=resu_comp, name='one',
                                          latex_name='1')
