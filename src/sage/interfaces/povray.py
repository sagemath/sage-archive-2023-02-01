r"""
POV-Ray, The Persistence of Vision Ray Tracer
"""


from sage.misc.pager import pager
import os

class POVRay:
    """
    POV-Ray The Persistence of Vision Ray Tracer

    INPUT:
        pov_file -- complete path to the .pov file you want to be rendered
        outfile -- the filename you want to save your result to
        **kwarg -- additionally keyword arguments you want to pass to POVRay

    OUTPUT:
        Image is written to the file you specified in outfile

    EXAMPLES:

    AUTHOR:
        SAGE Interface written by Yi Qiang (yqiang _atNOSPAM_ gmail.com)
        POVRay: http://www.povray.org
    """
    def __repr__(self):
        return 'POV-Ray The Persistence of Vision Ray Tracer'

    def __call__(self, pov_file, outfile='sage.ppm', block=True, **kwarg):
        if not os.path.isfile(pov_file):
            return "%s not found" % (pov_file)

        outfile = os.path.abspath(os.path.expanduser(outfile))
        ext = outfile[-4:].lower()
        try:
            width = kwarg['width']
            height = kwarg['height']
        except:
            return "You must specify a width and height."

        # test for start row and end row
        try:
            sr = kwarg['sr']
            er = kwarg['er']
        except:
            sr = 1
            er = height

        cmd = "povray +I%s +FP +SR%s +ER%s +W%s +H%s +O%s -D" % (pov_file,
                                                                 sr, er,
                                                                 width, height,
                                                                 outfile)
        if not block:
            cmd += ' &'
        os.system(cmd)

    def usage(self):
        r = os.popen('povray').read()
        pager()(r)

povray = POVRay()
