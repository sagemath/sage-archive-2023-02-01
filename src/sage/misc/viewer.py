r"""
Determination of Programs for Viewing web pages, etc.
"""

import os

def cmd_exists(cmd):
    return os.system('which %s 2>/dev/null >/dev/null'%cmd) == 0

BROWSER    = None
DVI_VIEWER = None
PDF_VIEWER = None
PNG_VIEWER = None

def viewer():
    global BROWSER, DVI_VIEWER, PDF_VIEWER, PNG_VIEWER

    if not (BROWSER is None):
        return BROWSER

    if os.environ.has_key('SAGE_BROWSER'):
        BROWSER = os.environ['SAGE_BROWSER']
        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER
        return BROWSER

    if os.uname()[0] == 'Darwin':
        # Simple on OS X, since there is an open command that opens
        # anything, using the user's preferences.
        # sage-open -- a wrapper around OS X open that
        # turns off any of SAGE's special library stuff.

        BROWSER = 'sage-open'
        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER

    elif os.uname()[0][:6] == 'CYGWIN':
        # Windows is also easy, since it has a system for
        # determining what opens things.
        # Bobby Moreti provided the following.

        if not os.environ.has_key('BROWSER'):
            systemroot = os.environ['SYSTEMROOT'].replace(':','/').replace('\\','')
            systemroot = '/cygdrive/' + systemroot
            BROWSER = '%s/system32/rundll32.exe url.dll,FileProtocolHandler'%\
                      systemroot

        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER

    elif cmd_exists('xdg-open'):

        # On other OS'es try xdg-open if present.
        # See http://portland.freedesktop.org/xdg-utils-1.0.

            BROWSER = 'xdg-open'
            DVI_VIEWER = BROWSER
            PDF_VIEWER = BROWSER
            PNG_VIEWER = BROWSER

    else:

        # If all fails try to get something from the environment.

        try:
            BROWSER = os.environ['BROWSER']
        except KeyError:
            BROWSER = 'less'  # silly default; lets hope it doesn't come to this!
            for cmd in ['firefox', 'konqueror', 'mozilla', 'mozilla-firefox']:
                if cmd_exists(cmd):
                    BROWSER = cmd
                    break
        DVI_VIEWER = BROWSER
        PDF_VIEWER = BROWSER
        PNG_VIEWER = BROWSER

        # Alternatives, if they are set in the environment or available.
        try:
            DVI_VIEWER = os.environ['DVI_VIEWER']
        except KeyError:
            for cmd in ['xdvi', 'kdvi']:
                if cmd_exists(cmd):
                    DVI_VIEWER = cmd
                    break
        try:
            PDF_VIEWER = os.environ['PDF_VIEWER']
        except KeyError:
            for cmd in ['acroread', 'xpdf']:
                if cmd_exists(cmd):
                    PDF_VIEWER = cmd
                    break

    return BROWSER


def browser():
    viewer()
    return "sage-native-execute " + BROWSER

def dvi_viewer():
    viewer()
    return "sage-native-execute " + DVI_VIEWER

def pdf_viewer():
    viewer()
    return "sage-native-execute " + PDF_VIEWER

def png_viewer():
    viewer()
    return "sage-native-execute " + PNG_VIEWER
