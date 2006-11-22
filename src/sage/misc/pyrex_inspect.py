import inspect, os

SAGE_ROOT = os.environ["SAGE_ROOT"]

def getsource(obj, is_binary):
    if not is_binary:
        try:
            return inspect.getsource(obj)
        except Exception:
            return None
    try:
        # Perhaps it's binary and is defined in a Pyrex file.
        try:
            d = inspect.getdoc(obj)
        except Exception, msg:
            return None
        i = d.find('\n')
        line0 = d[:i]
        if line0[:5] == "File:" and '.pyx' in line0:
            # Pyrex file -- we know the file position.
            j = line0.find('(')
            filename = '%s/local/lib/python/site-packages/%s'%(SAGE_ROOT, line0[5:j].strip())
            if not os.path.exists(filename):
                return ''
            l = line0.rfind('line ')
            source_lines = open(filename).readlines()
            s = line0[l+5:].rstrip(')').strip()
            lineno = eval(s)
            while lineno > 0 and source_lines[lineno].lstrip()[:4] != 'def ':
                lineno -= 1
            F = source_lines[lineno:]
            # Go down the file until we find another line indented the same as the first line.
            line0 = F[0]
            indent = len(line0) - len(line0.lstrip())     # indention level

            src = 'File: ' + filename + '\n\n' + F[0]
            i = 1
            while i < len(F):
                line = F[i]
                if len(line.strip()) == 0:
                    src += line
                    i += 1
                    continue
                line_indent = len(line) - len(line.lstrip())

                if line_indent <= indent:
                    break
                src += line
                i += 1
            return src
        return None

    except Exception, msg:
        return None
