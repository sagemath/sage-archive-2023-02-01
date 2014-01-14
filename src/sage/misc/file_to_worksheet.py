import os

def file_to_worksheet(filename):
    """
    Convert a Python file to a worksheet suitable for editing and debugging.
    """
    r = open(filename).readlines()
    for i in range(len(r)):
        s = r[i].lstrip()
        if s.startswith('sage:'):
            r[i] = '{{{%s}}}\n'%(s[5:].rstrip())

    s = '<pre>' + ''.join(r) + '</pre>'
    open(os.path.abspath(filename) + '.ws', 'w').write(s)

##     try:
##         while True:
##             i = i + find_first_prompt(s[i:])
##             j = i + find_end_of_input(s[i:])
##             k = j + find_end_of_block(s[j:])
##     except ValueError:
##         pass

## def find_first_prompt(s):
##     j = s.find('sage:')
##     if j == -1:
##         raise ValueError
##     return j

## def find_end_of_block(s):
##     """
##     Next blank line, triple quotes, or sage:
##     """



