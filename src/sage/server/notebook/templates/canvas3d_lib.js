/* This JavaScript library implements a software 3D renderer drawing to a
 * canvas context, which makes up the 'canvas3d' backend used in 3D plotting. */

canvas3d = (function() {
    var CANVAS_SIZE = 400;
    // The focal length is used to compute the projection transformation.
    var FOCAL_LENGTH = 300;
    // This scaling factor is applied to the figure before it is displayed.
    var FIGURE_SCALE = 90;
    // The figure will be translated by this amount on the Z axis.
    var FIGURE_ZOFFSET = -250;
    // This coefficient affects how rapidly the figure is scaled when the user
    // holds SHIFT and drags.
    var SCALING_SENSITIVITY = 0.04;

    ///////////////////////////////////////////////////////////////////////
    //
    // 3D Math and Utilities
    //
    ///////////////////////////////////////////////////////////////////////

    /* 4x4 matrices are stored in arrays of size 12 in row-major order, where
     * the last row is implied to be [0, 0, 0, 1]. This constitutes an affine
     * transformation matrix. */

    /* Three-dimensional vectors are stored as objects with three fields, one of
     * each 'x', 'y', and 'z'. There is no constructor function for these objects,
     * so the preferred way is using object literal notation ({ ... }). */


    // For constructing special transformation matrices:

    function make_identity_affine() {
        return [1, 0, 0, 0,
                0, 1, 0, 0,
                0, 0, 1, 0];
    }
    function make_translation_affine(x, y, z) {
        return [1, 0, 0, x,
                0, 1, 0, y,
                0, 0, 1, z];
    }
    function make_dilation_affine(a) {
        return [a, 0, 0, 0,
                0, a, 0, 0,
                0, 0, a, 0];
    }

    function mult_matrix(a, b) {
        /* Multiply two matrices and return the product. This computation is
         * unrolled, and includes the implied last rows. */

        return [
            (a[0] * b[0]) + (a[1] * b[4]) + (a[2] * b[8]),
            (a[0] * b[1]) + (a[1] * b[5]) + (a[2] * b[9]),
            (a[0] * b[2]) + (a[1] * b[6]) + (a[2] * b[10]),
            (a[0] * b[3]) + (a[1] * b[7]) + (a[2] * b[11]) + a[3],
            (a[4] * b[0]) + (a[5] * b[4]) + (a[6] * b[8]),
            (a[4] * b[1]) + (a[5] * b[5]) + (a[6] * b[9]),
            (a[4] * b[2]) + (a[5] * b[6]) + (a[6] * b[10]),
            (a[4] * b[3]) + (a[5] * b[7]) + (a[6] * b[11]) + a[7],
            (a[8] * b[0]) + (a[9] * b[4]) + (a[10] * b[8]),
            (a[8] * b[1]) + (a[9] * b[5]) + (a[10] * b[9]),
            (a[8] * b[2]) + (a[9] * b[6]) + (a[10] * b[10]),
            (a[8] * b[3]) + (a[9] * b[7]) + (a[10] * b[11]) + a[11]
        ];
    }

    // Vector operations:

    function transform_point(t, p) {
        /* Transform the point by the given transformation matrix. */

        return {
            x: p.x * t[0] + p.y * t[1] + p.z * t[2] + t[3],
            y: p.x * t[4] + p.y * t[5] + p.z * t[6] + t[7],
            z: p.x * t[8] + p.y * t[9] + p.z * t[10] + t[11]
        };
    }

    function vec3_dot(a, b) {
        /* Take the dot product of two three-dimensional vectors. */

        return a.x * b.x + a.y * b.y + a.z * b.z;
    }

    function vec3_norm(v) {
        /* Turn the vector into a unit vector by scaling every component. */

        var mag = Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
        v.x /= mag;
        v.y /= mag;
        v.z /= mag;
    }

    function project_point(p) {
        /* Project a 3D point onto a 2D image plane using the pinhole camera model
         * (http://en.wikipedia.org/wiki/Pinhole_camera_model). */

        var coeff = FOCAL_LENGTH / p.z;
        return { x: coeff * p.x, y: coeff * p.y };
    }

    ///////////////////////////////////////////////////////////////////////
    //
    // Trackball Rotation
    //
    ///////////////////////////////////////////////////////////////////////

    /* This trackball code is based off of the paper "Virtual Trackballs Revisited"
     * (http://image.diku.dk/research/trackballs/index.html) and the accompanying
     * code in the OpenTissue framework (http://www.opentissue.org/). */

    function Trackball(radius) {
        /* Construct a trackball with the given radius. */

        this.radius = radius;
        this.transform = make_identity_affine();
        // This is a unit vector corresponding to the position on the sphere
        // where the user started dragging (as initialized in begin_drag()).
        this.anchor_vec = null;
        // This unit vector corresponds to the position on the sphere where
        // the user's mouse cursor is currently located (as updated during drag()).
        this.cursor_vec = null;
    }

    Trackball.prototype.begin_drag = function(p) {
        /* Indicate that the user has started dragging. */

        this.project_onto_surface(p);
        vec3_norm(p);
        this.anchor_vec = this.cursor_vec = p;
    }

    Trackball.prototype.end_drag = function() {
        /* Indicate that the user has stopped dragging. This function returns
         * the overall transformation matrix for the completed drag. */

        var old_transform = this.transform;
        this.transform = make_identity_affine();
        return old_transform;
    }

    Trackball.prototype.drag = function(p) {
        /* Indicate that the user is dragging. This will update the transformation
         * matrix stored in Trackball.transform. Call this function repeatedly
         * in-between calls to begin_drag() and end_drag(). */

        this.project_onto_surface(p);
        vec3_norm(p);
        this.cursor_vec = p;
        this.compute_transform();
    }

    Trackball.prototype.compute_transform = function() {
        /* Recompute the transformation matrix for the trackball, based on the
         * anchor position and the current position. */

        // We calculate the rotation quaternion as the negative of the product
        // of Q_a and Q_c, where Q_a is the unit quaternion in the direction of
        // this.anchor_vec and Q_c is the unit quaternion in the direction of
        // this.current_vec.
        var quat_s = vec3_dot(this.anchor_vec, this.cursor_vec);
        var quat_x = -(this.anchor_vec.y * this.cursor_vec.z - this.anchor_vec.z * this.cursor_vec.y);
        var quat_y = -(this.anchor_vec.z * this.cursor_vec.x - this.anchor_vec.x * this.cursor_vec.z);
        var quat_z = -(this.anchor_vec.x * this.cursor_vec.y - this.anchor_vec.y * this.cursor_vec.x);

        // Convert the quaternion into a matrix.
        this.transform[0]  = 1 - 2 * ( (quat_y * quat_y) + (quat_z * quat_z));
        this.transform[5]  = 1 - 2 * ( (quat_x * quat_x) + (quat_z * quat_z));
        this.transform[10] = 1 - 2 * ( (quat_y * quat_y) + (quat_x * quat_x));
        this.transform[4]  =     2 * ( (quat_x * quat_y) + (quat_s * quat_z));
        this.transform[1]  =     2 * ( (quat_x * quat_y) - (quat_s * quat_z));
        this.transform[8]  =     2 * (-(quat_s * quat_y) + (quat_x * quat_z));
        this.transform[2]  =     2 * ( (quat_s * quat_y) + (quat_x * quat_z));
        this.transform[9]  =     2 * ( (quat_z * quat_y) + (quat_s * quat_x));
        this.transform[6]  =     2 * ( (quat_z * quat_y) - (quat_s * quat_x));
    }

    Trackball.prototype.project_onto_surface = function(p) {
        /* This function implements Bell's function (described in the paper) for
         * projecting a 2D point onto the trackball sphere. Note that this modifies
         * the provided 2D point by added a 'z' field for the third-dimensional
         * coordinate. */

        var radius_squared = this.radius * this.radius;
        var length_squared = p.x * p.x + p.y * p.y;
        if(length_squared <= radius_squared / 2)
            p.z = Math.sqrt(radius_squared - length_squared);
        else
            p.z = radius_squared / (2 * Math.sqrt(length_squared));
    }

    ///////////////////////////////////////////////////////////////////////
    //
    // Rendering and Notebook Integration
    //
    ///////////////////////////////////////////////////////////////////////

    var browser_opera = navigator.userAgent.indexOf('Opera') > -1;
    function viewport_offset_for_element(e) {
        /* Returns the X and Y coordinates of the element relative to the viewport. */

        // This code was adapted from Element.viewportOffset in the Prototype
        // JavaScript framework (http://www.prototypejs.org/).

        var value_top = 0, value_left = 0;
        var element = e;
        // Traverse up the element heirarchy, keeping track of each element's
        // relative offset, until we have an absolute offset for the original
        // element.
        do {
            value_top += element.offsetTop || 0;
            value_left += element.offsetLeft || 0;
            // Safari fix?
            if(element.offsetParent == document.body &&
               element.style.position == 'absolute') break;
        } while(element = element.offsetParent);
        // Now we must compute the scroll offset, and adjust our result accordingly.
        element = e;
        do {
            if(!browser_opera ||
               (element.tagName && (element.tagName.toUpperCase() == 'BODY'))) {

                value_top -= element.scrollTop || 0;
                value_left -= element.scrollLeft || 0;
            }
        } while(element = element.parentNode);
        return { x: value_left, y: value_top };
    }

    function setup(canvas) {
        /* Setup the provided canvas with event listeners and state (stored in
         * closures) to be an interactive 3D model viewer. */

        var dragging = false;
        var last_mouse_x = 0, last_mouse_y = 0, initialized = false;
        var camera_scale = 1;
        var camera_transform = make_identity_affine();
        var trackball = new Trackball(CANVAS_SIZE / 2);
        var pending_update = null; // A timeout ID which may correspond to some
                                   // callback to redraw the canvas.
        function adapt_mouse_pos_for_trackball(evt) {
            var canvas_pos = viewport_offset_for_element(canvas);
            return {
                x: evt.clientX - canvas_pos.x - CANVAS_SIZE / 2,
                y: evt.clientY - canvas_pos.y - CANVAS_SIZE / 2
            };
        }
        function update() {
            var t = make_identity_affine();
            t = mult_matrix(camera_transform, t);
            t = mult_matrix(trackball.transform, t);
            t = mult_matrix(make_dilation_affine(FIGURE_SCALE * camera_scale), t);
            t = mult_matrix(make_translation_affine(0, 0, -FIGURE_ZOFFSET), t);
            draw(canvas, t);
            pending_update = null;
        }
        function schedule_update() {
            if(pending_update != null)
                clearTimeout(pending_update);
            pending_update = setTimeout(update, 0);
        }

        canvas.addEventListener("mousedown", function(evt) {
            trackball.begin_drag(adapt_mouse_pos_for_trackball(evt));
            dragging = true;
        }, false);
        canvas.addEventListener("mouseup", function(evt) {
            camera_transform = mult_matrix(trackball.end_drag(), camera_transform);
            dragging = false;
        }, false);
        canvas.addEventListener("mousemove", function(evt) {
            if(initialized) {
                if(dragging) {
                    if(evt.shiftKey) {
                        camera_scale += SCALING_SENSITIVITY * (evt.clientY - last_mouse_y);
                        if(camera_scale < 0.2)
                            camera_scale = 0.2;
                    } else
                        trackball.drag(adapt_mouse_pos_for_trackball(evt));
                    schedule_update();
                }
            } else
                initialized = true;
            last_mouse_x = evt.clientX;
            last_mouse_y = evt.clientY;
        }, false);

        update();
    }

    function render_model(ctx, transform, model) {
        if("color" in model)
            ctx.strokeStyle = ("#" + model.color);
        else
            ctx.strokeStyle = "black";

        for(var i = 0; i < model.faces.length; i++) {
            ctx.beginPath();
            var points = new Array(model.faces[i].length);
            var culled_points_count = 0;
            for(var j = 0; j < model.faces[i].length; j++) {
                var transformed_vertex =
                      transform_point(transform, model.vertices[model.faces[i][j]]);
                points[j] = project_point(transformed_vertex);
                if(points[j].x < CANVAS_SIZE / -2 || points[j].x > CANVAS_SIZE / 2 ||
                   points[j].y < CANVAS_SIZE / -2 || points[j].y > CANVAS_SIZE / 2 ||
                   transformed_vertex.z < 0) {

                    culled_points_count++;
                }
            }
            if(culled_points_count < model.faces[i].length) {
                for(var j = 0; j < model.faces[i].length; j++) {
                    if(j == 0)
                        ctx.moveTo(points[j].x, points[j].y);
                    else
                        ctx.lineTo(points[j].x, points[j].y);
                }
            }
            ctx.closePath();
            ctx.stroke();
        }
    }

    function draw(canvas, transform) {
        /* Redraw the specified canvas. Vertex and face data are stored as a
         * property of the canvas object. */

        var ctx = canvas.getContext('2d');
        ctx.save();
        ctx.clearRect(0, 0, CANVAS_SIZE, CANVAS_SIZE);
        ctx.translate(CANVAS_SIZE / 2, CANVAS_SIZE / 2);
        for(var i = 0; i < canvas.data.length; i++)
            render_model(ctx, transform, canvas.data[i]);
        ctx.restore();
    }

    var viewer_count = 0;
    function viewer(url) {
        var canvas_id = "canvas3d-viewer" + (viewer_count++);
        cell_writer.write(('<canvas style="border: 2px solid black" id="' + canvas_id +
                           '" width="' + CANVAS_SIZE + '" height="' + CANVAS_SIZE +'">') +
                          'Sorry, but you need a browser that supports the &lt;canvas&gt; tag.' +
                          '</canvas>');
        // Send an XHR to get the JSON model data stored at the URL
        var xhr = new XMLHttpRequest();
        xhr.onreadystatechange = function() {
            if(xhr.readyState == 4) {
                if(xhr.status == 200 && xhr.responseText != null) {
                    var canvas = document.getElementById(canvas_id);
                    canvas.data = eval('(' + xhr.responseText + ')');
                    setup(canvas);
                }
            }
        }
        xhr.open('GET', url, true);
        xhr.send(null);
    }

    return {
        viewer: viewer
    };
})();
