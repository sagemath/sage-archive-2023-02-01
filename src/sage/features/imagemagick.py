# -*- coding: utf-8 -*-
r"""
Check for imagemagick

Currently we only check for the presence of ``convert``. When needed other
commands like ``magick``, ``magick-script``, ``convert``, ``mogrify``,
``identify``, ``composite``, ``montage``, ``compare``, etc. could be also
checked in this module.
"""

from . import Executable

class ImageMagick(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of
    ``ImageMagick``

    Currently, only the availability of ``convert`` is checked.

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

