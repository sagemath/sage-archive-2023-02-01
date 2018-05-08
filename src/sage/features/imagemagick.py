# -*- coding: utf-8 -*-
r"""
Check for imagemagick
"""

from . import Executable

class ImageMagick(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    ``ImageMagick``

    EXAMPLES::

        sage: from sage.features.imagemagick import ImageMagick
        sage: ImageMagick().is_present()  # optional: imagemagick
        FeatureTestResult('convert', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.imagemagick import ImageMagick
            sage: isinstance(ImageMagick(), ImageMagick)
            True
        """
        Executable.__init__(self, "convert", executable="convert",
                url="https://www.imagemagick.org/")

