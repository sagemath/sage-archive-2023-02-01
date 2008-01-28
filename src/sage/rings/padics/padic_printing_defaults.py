def mode(mode = None):
    import sage.rings.padics.padic_printing as printing
    if printing._printer_defaults is None:
        printing._printer_defaults = printing.pAdicPrinterDefaults()
    if mode is None:
        return printing._printer_defaults.mode
    else:
        if mode in ['val-unit','series','terse','digits','bars']:
            printing._printer_defaults.mode = mode
        else:
            raise ValueError, "invalid printing mode"

def allow_negatives(neg = None):
    import sage.rings.padics.padic_printing as printing
    if printing._printer_defaults is None:
        printing._printer_defaults = printing.pAdicPrinterDefaults()
    if neg is None:
        return not printing._printer_defaults.pos
    else:
        printing._printer_defaults.pos = not neg

def max_series_terms(max = None):
    import sage.rings.padics.padic_printing as printing
    if printing._printer_defaults is None:
        printing._printer_defaults = printing.pAdicPrinterDefaults()
    if max is None:
        return printing._printer_defaults.max_ram_terms
    else:
        from sage.rings.integer import Integer
        printing._printer_defaults.max_ram_terms = Integer(max)

def max_unram_terms(max = None):
    import sage.rings.padics.padic_printing as printing
    if printing._printer_defaults is None:
        printing._printer_defaults = printing.pAdicPrinterDefaults()
    if max is None:
        return printing._printer_defaults.max_unram_terms
    else:
        from sage.rings.integer import Integer
        printing._printer_defaults.max_unram_terms = Integer(max)

def max_poly_terms(max = None):
    import sage.rings.padics.padic_printing as printing
    if printing._printer_defaults is None:
        printing._printer_defaults = printing.pAdicPrinterDefaults()
    if max is None:
        return printing._printer_defaults.max_terse_terms
    else:
        from sage.rings.integer import Integer
        printing._printer_defaults.max_terse_terms = Integer(max)

def sep(sep = None):
    import sage.rings.padics.padic_printing as printing
    if printing._printer_defaults is None:
        printing._printer_defaults = printing.pAdicPrinterDefaults()
    if sep is None:
        return printing._printer_defaults.sep
    else:
        printing._printer_defaults.sep = sep

def alphabet(alphabet = None):
    import sage.rings.padics.padic_printing as printing
    if printing._printer_defaults is None:
        printing._printer_defaults = printing.pAdicPrinterDefaults()
    if alphabet is None:
        return printing._printer_defaults.alphabet
    else:
        printing._printer_defaults.alphabet = alphabet


