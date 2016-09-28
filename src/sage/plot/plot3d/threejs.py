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

def threejs(plot):

    from sage.plot.plot3d.base import Graphics3d
    if not isinstance(plot, Graphics3d):
        raise TypeError('input plot must be an instance of Graphics3d')

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

    html = threejs_template().format(lights, bounds, points, lines, surfaces)

    from sage.misc.temporary_file import tmp_filename
    temp_filename = tmp_filename(ext='.html')

    from sage.repl.rich_output.backend_ipython import BackendIPythonCommandline
    backend = BackendIPythonCommandline()

    f = open(temp_filename, 'w')
    f.write(html)
    f.close()

    backend.launch_viewer(temp_filename, 'Template')


def threejs_template():

    return r"""
<!DOCTYPE html>
<html>
<head>
<title></title>
<meta charset=utf-8>
<meta name=viewport content='width=device-width, user-scalable=no, minimum-scale=1.0, maximum-scale=1.0'>
<style>

    body {{ margin: 0px; overflow: hidden; }}
  
</style>
</head>

<body>

<script src=http://rawgit.com/mrdoob/three.js/r80/build/three.js></script>
<script src=http://rawgit.com/mrdoob/three.js/r80/examples/js/controls/OrbitControls.js></script>

<script>

    var scene = new THREE.Scene();

    var renderer = new THREE.WebGLRenderer( {{ antialias: true }} );
    renderer.setSize( window.innerWidth, window.innerHeight );
    renderer.setClearColor( 0xffffff, 1 );
    document.body.appendChild( renderer.domElement ); 

    var lights = {};
    for ( var i=0 ; i < lights.length ; i++ ) {{
        var light = new THREE.DirectionalLight( 0xdddddd, 1 );
        light.position.set( lights[i].x, lights[i].y, lights[i].z );
        scene.add( light );
    }}

    scene.add( new THREE.AmbientLight( 0x404040, 1 ) );

    var bounds = {};

    var points = {};
    for ( var i=0 ; i < points.length ; i++ ) {{
        eval( 'var json = ' + points[i] );
        moveBoundsOffPoint( json );
    }}

    function moveBoundsOffPoint() {{
        var v = json.point;
        if ( v[0] == bounds[0].x ) bounds[0].x -= 1;
        if ( v[1] == bounds[0].y ) bounds[0].y -= 1;
        if ( v[2] == bounds[0].z ) bounds[0].z -= 1;
        if ( v[0] == bounds[1].x ) bounds[1].x += 1;
        if ( v[1] == bounds[1].y ) bounds[1].y += 1;
        if ( v[2] == bounds[1].z ) bounds[1].z += 1;
    }}

    var rRange = Math.sqrt( Math.pow( bounds[1].x - bounds[0].x, 2 )
                            + Math.pow( bounds[1].x - bounds[0].x, 2 ) );
    var zRange = bounds[1].z - bounds[0].z;

    var a = [ 1, 1, 1 ]; // aspect multipliers
    var autoAspect = 2.5;
    if ( zRange > autoAspect * rRange ) a[2] = autoAspect * rRange / zRange; 

    var box = new THREE.Geometry();
    box.vertices.push( new THREE.Vector3( bounds[0].x, bounds[0].y, a[2]*bounds[0].z ) );
    box.vertices.push( new THREE.Vector3( bounds[1].x, bounds[1].y, a[2]*bounds[1].z ) );
    var boxMesh = new THREE.LineSegments( box );
    scene.add( new THREE.BoxHelper( boxMesh, 'black' ) );

    scene.position.set( -( bounds[0].x + bounds[1].x ) / 2,
                        -( bounds[0].y + bounds[1].y ) / 2,
                        -a[2]*( bounds[0].z + bounds[1].z ) / 2 )

    // optional axis helper 
    //scene.add( new THREE.AxisHelper( Math.min( [ bounds[1].x, bounds[1].y, a[2]*bounds[1].z ] ) ) );

    var camera = new THREE.PerspectiveCamera( 45, window.innerWidth / window.innerHeight, 0.1, 1000 ); 
    camera.up.set( 0, 0, 1 );
    var cameraOut = a[2]*zRange;
    camera.position.set( cameraOut, 1.3*cameraOut, .7*cameraOut );
    camera.lookAt( scene.position );

    var controls = new THREE.OrbitControls( camera, renderer.domElement );

    window.addEventListener( 'resize', function() {{
        
        renderer.setSize( window.innerWidth, window.innerHeight );
        camera.aspect = window.innerWidth / window.innerHeight;
        camera.updateProjectionMatrix();
        
    }} );
    
    for ( var i=0 ; i < points.length ; i++ ) {{
        eval( 'var json = ' + points[i] );
        addPoint( json );
    }}

    function addPoint( json ) {{
        var geometry = new THREE.Geometry();
        var v = json.point;
        geometry.vertices.push( new THREE.Vector3( v[0], v[1], a[2]*v[2] ) );
        var canvas = document.createElement( 'canvas' );
        canvas.width = 128;
        canvas.height = 128;
        var context = canvas.getContext( '2d' );
        context.arc( 64, 64, 64, 0, 2 * Math.PI );
        context.fillStyle = json.color;
        context.fill();
        var texture = new THREE.Texture( canvas );
        texture.needsUpdate = true;
        var material = new THREE.PointsMaterial( {{ size: json.size/100, map: texture,
                                                    transparent: true, alphaTest: .1 }} );
        scene.add( new THREE.Points( geometry, material ) );
    }}

    var lines = {};
    for ( var i=0 ; i < lines.length ; i++ ) {{
        eval( 'var json = ' + lines[i] );
        addLine( json );
    }}

    function addLine( json ) {{
        var geometry = new THREE.Geometry();
        for ( var i=0 ; i < json.points.length - 1 ; i++ ) {{
            var v = json.points[i];
            geometry.vertices.push( new THREE.Vector3( v[0], v[1], a[2]*v[2] ) );
            var v = json.points[i+1];
            geometry.vertices.push( new THREE.Vector3( v[0], v[1], a[2]*v[2] ) );
        }}
        var material = new THREE.LineBasicMaterial( {{ color: json.color }} );
        scene.add( new THREE.LineSegments( geometry, material ) );
    }}

    var surfaces = {};
    for ( var i=0 ; i < surfaces.length ; i++ ) {{
        eval( 'var json = ' + surfaces[i] );
        addSurface( json );
    }}

    function addSurface( json ) {{
        var geometry = new THREE.Geometry();
        for ( var i=0 ; i < json.vertices.length ; i++ ) {{
            var v = json.vertices[i];
            geometry.vertices.push( new THREE.Vector3( v.x, v.y, a[2]*v.z ) );
        }}
        for ( var i=0 ; i < json.faces.length ; i++ ) {{
            var f = json.faces[i]
            for ( var j=0 ; j < f.length - 2 ; j++ ) {{
                geometry.faces.push( new THREE.Face3( f[0], f[j+1], f[j+2] ) );
            }}
        }}
        geometry.computeVertexNormals();
        var material = new THREE.MeshPhongMaterial( {{ color: json.color , side: THREE.DoubleSide,
                                                       shininess: 20 }} );
        scene.add( new THREE.Mesh( geometry, material ) );
    }}

    function render() {{

        requestAnimationFrame( render ); 
        renderer.render( scene, camera );
        controls.update();

    }}
    
    render();

</script>

</body>
</html>
"""
