# -*- coding: utf-8 -*-
r"""
Check for FFmpeg
"""

from . import Executable

class FFmpeg(Executable):
    r"""
    A :class:`sage.features.Feature` describing the presence of ``FFmpeg``

    EXAMPLES::

        sage: from sage.features.ffmpeg import FFmpeg
        sage: FFmpeg().is_present()  # optional: ffmpeg
        FeatureTestResult('FFmpeg', True)
    """
    def __init__(self):
        r"""
        TESTS::

            sage: from sage.features.ffmpeg import FFmpeg
            sage: isinstance(FFmpeg(), FFmpeg)
            True
        """
        Executable.__init__(self, "FFmpeg", executable="ffmpeg",
                url="https://www.ffmpeg.org/")

