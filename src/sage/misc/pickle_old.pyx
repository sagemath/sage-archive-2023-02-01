from sage.misc.misc import deprecation

def make_element_old(_class, _dict, parent):
    """
    Used for unpickling old pickles of Element objects (and subclasses).

    This function is deprecated and is kept only to support old pickles.
    """
    deprecation("Your data is stored in an old format. Please use the save() function to store your data in a more recent format.")
    new_object = _class.__new__(_class)
    new_object._set_parent(parent)
    new_object.__dict__ = _dict
    return new_object
