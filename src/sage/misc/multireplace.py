"multi_replace"

##########################################################################
#
# multi_replace function
#
# By Xavier Defrang.
#
# From the Python cookbook:
#
#  http://aspn.activestate.com/ASPN/Cookbook/Python/Recipe/81330
#
##########################################################################


import re


#
# The simplest, lambda-based implementation
#

def multiple_replace(dic, text):
    """
    Replace in 'text' all occurrences of any key in the given
    dictionary by its corresponding value.  Returns the new string.
    """
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dic)))

    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start():mo.end()]], text)
