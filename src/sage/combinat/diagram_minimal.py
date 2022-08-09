from sage.categories.sets_cat import Sets
from sage.structure.unique_representation import UniqueRepresentation
from sage.structure.list_clone import ClonableArray
from sage.structure.parent import Parent
from sage.misc.inherit_comparison import InheritComparisonClasscallMetaclass

class Diagram(ClonableArray, metaclass=InheritComparisonClasscallMetaclass):
    @staticmethod
    def __classcall_private__(cls, cells):
        return Diagrams()(cells)

    def __init__(self, parent, cells):
        self._cells = {c: True for c in cells}
        ClonableArray.__init__(self, parent, cells)

    def check(self):
        pass

class Diagrams(UniqueRepresentation, Parent):
    def __init__(self):
        Parent.__init__(self, category=Sets())

    def _element_constructor_(self, cells):
        return self.element_class(self, cells)

    Element = Diagram