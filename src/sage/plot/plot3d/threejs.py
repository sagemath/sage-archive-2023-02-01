r"""
Three.js Support for 3D Plots

"""

#*****************************************************************************
#       Copyright (C) 2016 Paul Masson <paulmasson@analyticphysics.com>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#  as published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#                  http://www.gnu.org/licenses/
#*****************************************************************************

def threejs(plot, **kwds):

    from sage.plot.plot3d.base import Graphics3d
    if not isinstance(plot, Graphics3d):
        raise TypeError('input plot must be an instance of Graphics3d')

    options = {}
    options['frame'] = kwds.get('frame', True)
    options['aspect_ratio'] = kwds.get('aspect_ratio', [1,1,1])
    options['axes_labels'] = kwds.get('axes_labels', ['x','y','z'])
    options['decimals'] = int(kwds.get('decimals', 0))
    options['opacity'] = float(kwds.get('opacity', 1))
    options['thickness'] = float(kwds.get('thickness', 1))

    lights = "[{x:0,y:0,z:10},{x:0,y:0,z:-10}]"

    b = plot.bounding_box()
    bounds = "[{{x:{},y:{},z:{}}},{{x:{},y:{},z:{}}}]".format(
             b[0][0],b[0][1],b[0][2],b[1][0],b[1][1],b[1][2])

    import json
    points, lines = [], []
    if not hasattr(plot, 'all'):
        plot += Graphics3d()
    for p in plot.all:
        if hasattr(p, 'loc'):
            color = p._extra_kwds.get('color', 'blue')
            points.append("{{point:{}, size:{}, color:'{}'}}".format(json.dumps(p.loc), p.size, color))
        if hasattr(p, 'points'):
            color = p._extra_kwds.get('color', 'blue')
            lines.append("{{points:{}, color:'{}'}}".format(json.dumps(p.points), color))

    from sage.plot.plot3d.base import flatten_list
    surfaces = plot.json_repr(plot.default_render_params())
    surfaces = flatten_list(surfaces)

    if len(points) == 0 and len(lines) == 0 and len(surfaces) == 0:
        raise ValueError('no data for this plot')

    import os
    from sage.env import SAGE_SRC
    filename = os.path.join(SAGE_SRC, 'sage',
                            'plot', 'plot3d', 'threejs_template.html')
    f = open(filename, 'r')
    html = f.read()
    f.close()

    html = html.replace('SAGE_OPTIONS', json.dumps(options))
    html = html.replace('SAGE_LIGHTS', lights)
    html = html.replace('SAGE_BOUNDS', bounds)
    html = html.replace('SAGE_POINTS', str(points))
    html = html.replace('SAGE_LINES', str(lines))
    html = html.replace('SAGE_SURFACES', str(surfaces))

    from sage.misc.temporary_file import tmp_filename
    temp_filename = tmp_filename(ext='.html')

    from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
    backend = BackendIPythonCommandline()

    f = open(temp_filename, 'w')
    f.write(html)
    f.close()

    backend.launch_viewer(temp_filename, 'Template')

