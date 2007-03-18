def ModularAbelianVariety(X):
    try:
        return X.modular_abelian_variety()
    except AttributeError:
        raise ValueError, "No known way to associate a modular abelian variety to %s"%X
