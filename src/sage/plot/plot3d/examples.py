r"""
Introduction

EXAMPLES::

    sage: x, y = var('x y')
    sage: W = plot3d(sin(pi*((x)^2+(y)^2))/2,(x,-1,1),(y,-1,1), frame=False, color='purple', opacity=0.8)
    sage: S = sphere((0,0,0),size=0.3, color='red', aspect_ratio=[1,1,1])
    sage: show(W + S, figsize=8)
"""



