# -*- encoding: utf-8 -*-
r"""
Catalog of all available output container types.

If you define another output type then you must add it to the imports here.
"""

from .output_basic import (
    OutputPlainText,
    OutputAsciiArt,
    OutputMathJax,
)

from .output_graphics import (
    OutputImagePng,
    OutputImageGif,
    OutputImageJpg,
    OutputImageSvg,
    OutputImagePdf,
    OutputImageDvi,
)

from .output_graphics3d import (
    OutputSceneJmol,
    OutputSceneLightwave,
    OutputSceneCanvas3d,
)
