from sagedoc import search_src

def class_graph(classname=None, filename=None, base_classname=None):
    # Get
    z = search_src('^class', interact=False)
    classes = []
    for x in z.splitlines():
        if '|' in x or '-----------------------------' in x:
            continue   # banner
        i = x.find('#')
        if i != -1:
            x = x[:i]
        x = x.strip()
        if x[-1] != ':':
            continue
        i = x.find(':')
        if filename:
            if not (filename in x[:i]):
                continue
        x = x[i+1:].strip().lstrip('class').strip().strip(':')
        i = x.find('(')
        cls = x[:i].strip()
        if classname:
            if not classname in cls:
                continue
        b = x[i+1:].strip(')')

        bases = []
        for a in b.split():
            i = a.rfind('.')
            if i == -1:
                nm = a
            else:
                nm = a[i+1:]
            if base_classname:
                if not (base_classname in nm):
                    continue
            bases.append(nm)

        classes.append((cls, bases))

    classes = dict(classes)

    return classes
