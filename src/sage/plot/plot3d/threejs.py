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

    data = plot.json_repr(plot.default_render_params())
    if len(data) == 0:
        raise ValueError('no json_repr for this plot')
    from sage.plot.plot3d.base import flatten_list
    data = flatten_list(data)

    b = plot.bounding_box()
    bounds = "[{{x:{},y:{},z:{}}},{{x:{},y:{},z:{}}}]".format(
             b[0][0],b[0][1],b[0][2],b[1][0],b[1][1],b[1][2])

    lights = "[{x:0,y:0,z:10},{x:0,y:0,z:-10}]"

    html = threejs_template().format(bounds, lights, data)

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

    // initialize common variables

    var scene = new THREE.Scene();
    var renderer = new THREE.WebGLRenderer( {{ antialias: true }} );
    renderer.setSize( window.innerWidth, window.innerHeight );
    renderer.setClearColor( 0xffffff, 1 );
    document.body.appendChild( renderer.domElement ); 

    var bounds = {};
    var box = new THREE.Geometry();
    box.vertices.push( new THREE.Vector3( bounds[0].x, bounds[0].y, bounds[0].z ) );
    box.vertices.push( new THREE.Vector3( bounds[1].x, bounds[1].y, bounds[1].z ) );
    var boxMesh = new THREE.LineSegments( box );
    scene.add( new THREE.BoxHelper( boxMesh, 'black' ) );

    scene.add( new THREE.AxisHelper( bounds[1].z ) );

    var camera = new THREE.PerspectiveCamera( 45, window.innerWidth / window.innerHeight, 0.1, 1000 ); 
    camera.up.set( 0, 0, 1 );
    camera.position.set( 1.5*bounds[1].x, 1.5*bounds[1].y, 1.5*bounds[1].z ); 
    var controls = new THREE.OrbitControls( camera, renderer.domElement );

    var lights = {};
    for ( var i=0 ; i < lights.length ; i++ ) {{
        var light = new THREE.DirectionalLight( 0xdddddd, 1 );
        light.position.set( lights[i].x, lights[i].y, lights[i].z );
        scene.add( light );
    }}

    scene.add( new THREE.AmbientLight( 0x404040, 1 ) );

    window.addEventListener( 'resize', function() {{
        
        renderer.setSize( window.innerWidth, window.innerHeight );
        camera.aspect = window.innerWidth / window.innerHeight;
        camera.updateProjectionMatrix();
        
    }} );
    
    // initialize geometry and material

    var data = {};
    for ( var i=0 ; i < data.length ; i++ ) {{
        eval( 'var json = ' + data[i] );
        addGeometry( json );
    }}

    function addGeometry( json ) {{
        var geometry = new THREE.Geometry();
        for ( var i=0 ; i < json.vertices.length ; i++ ) {{
            var v = json.vertices[i];
            geometry.vertices.push( new THREE.Vector3( v.x, v.y, v.z ) );
        }}
        for ( var i=0 ; i < json.faces.length ; i++ ) {{
            var f = json.faces[i]
            for ( var j=0 ; j < f.length - 2 ; j++ ) {{
                geometry.faces.push( new THREE.Face3( f[0], f[j+1], f[j+2] ) );
            }}
        }}
        geometry.computeVertexNormals();
        var material = new THREE.MeshPhongMaterial( {{ color: json.color , side: THREE.DoubleSide,
                                                       shininess: 20}} );
        scene.add( new THREE.Mesh( geometry, material ) );
    }}

    function render() {{

        requestAnimationFrame( render ); 
        renderer.render( scene, camera );
        controls.update();

    }}
    
    // initialize render
    render();

</script>

</body>
</html>
"""
