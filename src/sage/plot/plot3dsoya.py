

import sys, os, os.path

from sage.misc.functional import log




def plot3d_soya(A, fullscreen=False):
    """
    Use soya3d to draw a surface such that the point over the
    (i,j) position is A[i,j].

    EXAMPLES:
        sage: A = MatrixSpace(RR,32)([sin(i/2)*cos(j/2) for i in range(32) for j in range(32)])
        sage: from sage.plot.plot3d import plot3d_soya      # optional
        sage.: plot3d_soya(A)           # optional
    """

    try:
        import soya, soya.sdlconst
    except ImportError, msg:
        raise RuntimeError, "%s\nYou must install soya3d first."%msg

    if fullscreen:
        soya.init("Soya fullscreen tutorial", 800, 600, 1)
    else:
        soya.init()
    scene = soya.World()

    # Creates a new landscape in the scene. land_size is the
    # dimension of the landscape ;
    # it must be of the form (2 ** n) + 1.

    #land_size = 33
    nr = A.nrows()
    nc = A.ncols()
    n = max(nr, nc)
    land_size = 2**(log(n,2).ceil())  + 1
    land = soya.Land(scene, land_size, land_size)

    # Sets value for each height.
    # Other vertex-setting methods include:
    #  - Land.set_material     (i, j, material)
    #  - Land.set_vertex_color (i, j, color) where color is a (red, green, blue, alpha) tuple
    #  - Land.set_vertex_option(i, j, hidden, invisible, non_solid, force_presence)

    m = abs(A[0,0])
    for i in range(A.nrows()):
      for j in range(A.ncols()):
        #land.set_height(i, j, random.random())
        land.set_height(i, j, float(A[i,j]))
        m = max(m, abs(A[i,j]))

    # Multiplies all the heights by 4
    # land.multiply_height(4.0)

    # Adds a light.
    light = soya.Light(scene)
    light.set_xyz(0.0, 2*m, 0.0)

    # Add a camera and a loop to render
    class MovableCamera(soya.Camera):
      def __init__(self, parent):
        soya.Camera.__init__(self, parent)

        self.speed = soya.Vector(self)
        self.rotation_lateral_speed  = 0.0
        self.rotation_vertical_speed = 0.0

      def begin_round(self):
        soya.Camera.begin_round(self)

        for event in soya.process_event():
          if event[0] == soya.sdlconst.KEYDOWN:
            if   event[1] == soya.sdlconst.K_UP:     self.speed.z = -1.0
            elif event[1] == soya.sdlconst.K_DOWN:   self.speed.z =  1.0
            elif event[1] == soya.sdlconst.K_j:      self.speed.y = -1.0
            elif event[1] == soya.sdlconst.K_k:      self.speed.y =  1.0
            elif event[1] == soya.sdlconst.K_LEFT:   self.rotation_lateral_speed =  10.0
            elif event[1] == soya.sdlconst.K_RIGHT:  self.rotation_lateral_speed = -10.0
            elif event[1] == soya.sdlconst.K_q:      soya.IDLER.stop()
            elif event[1] == soya.sdlconst.K_ESCAPE: soya.IDLER.stop()
          if event[0] == soya.sdlconst.KEYUP:
            if   event[1] == soya.sdlconst.K_UP:     self.speed.z = 0.0
            elif event[1] == soya.sdlconst.K_DOWN:   self.speed.z = 0.0
            elif event[1] == soya.sdlconst.K_j:      self.speed.y = 0.0
            elif event[1] == soya.sdlconst.K_k:      self.speed.y = 0.0
            elif event[1] == soya.sdlconst.K_LEFT:   self.rotation_lateral_speed = 0.0
            elif event[1] == soya.sdlconst.K_RIGHT:  self.rotation_lateral_speed = 0.0

      def advance_time(self, proportion):
        self.add_mul_vector(proportion, self.speed)
        self.turn_lateral (self.rotation_lateral_speed  * proportion)
        self.turn_vertical(self.rotation_vertical_speed * proportion)


    camera = MovableCamera(scene)
    camera.set_xyz(16.0, 6.0, 0.0)
    camera.look_at(soya.Point(scene, 16.0, 6.0, 10.0))
    soya.set_root_widget(camera)

    soya.Idler(scene).idle()
    soya.quit()



