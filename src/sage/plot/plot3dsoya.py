##############################################################################
#       Copyright (C) 2006 Josh Kantor and William Stein <wstein@ucsd.edu>
#
#  Distributed under the terms of the GNU General Public License (GPL)
#
#                  http://www.gnu.org/licenses/
##############################################################################


import math, sys, os, os.path, time

from sage.misc.misc import SAGE_DATA, DOT_SAGE

from sage.server.support import EMBEDDED_MODE

try:
    import soya
    import soya.sdlconst
    import soya.pudding as pudding
    import soya.laser
except ImportError, msg:
    msg = "%s\n\nYou must install soya3d first, e.g., install an optional SAGE"%msg + '\n' +\
          'package named something like soya-*.'
    raise RuntimeError, msg



class ProgressMeter(object):
    ESC = chr(27)
    def __init__(self, **kw):
        # What time do we start tracking our progress from?
        self.timestamp = kw.get('timestamp', time.time())
        # What kind of unit are we tracking?
        self.unit = str(kw.get('unit', ''))
        # Number of units to process
        self.total = int(kw.get('total', 100))
        # Number of units already processed
        self.count = int(kw.get('count', 0))
        # Refresh rate in seconds
        self.rate_refresh = float(kw.get('rate_refresh', .5))
        # Number of ticks in meter
        self.meter_ticks = int(kw.get('ticks', 60))
        self.meter_division = float(self.total) / self.meter_ticks
        self.meter_value = int(self.count / self.meter_division)
        self.last_update = None
        self.rate_history_idx = 0
        self.rate_history_len = 10
        self.rate_history = [None] * self.rate_history_len
        self.rate_current = 0.0
        self.last_refresh = 0
        self._cursor = False
        if not EMBEDDED_MODE:
            self.reset_cursor()

    def reset_cursor(self, first=False):
        if self._cursor:
            sys.stdout.write(self.ESC + '[u')
        self._cursor = True
        sys.stdout.write(self.ESC + '[s')

    def update(self, count, **kw):
        now = time.time()
        # Caclulate rate of progress
        rate = 0.0
        # Add count to Total
        self.count += count
        self.count = min(self.count, self.total)
        if self.last_update:
            delta = now - float(self.last_update)
            if delta:
                rate = count / delta
            else:
                rate = count
            self.rate_history[self.rate_history_idx] = rate
            self.rate_history_idx += 1
            self.rate_history_idx %= self.rate_history_len
            cnt = 0
            total = 0.0
            # Average rate history
            for rate in self.rate_history:
                if rate == None:
                    continue
                cnt += 1
                total += rate
            rate = total / cnt
        self.rate_current = rate
        self.last_update = now
        # Device Total by meter division
        value = int(self.count / self.meter_division)
        if value > self.meter_value:
            self.meter_value = value
        if self.last_refresh:
            if (now - self.last_refresh) > self.rate_refresh or \
                (self.count >= self.total):
                    self.refresh()
        else:
            self.refresh()

    def get_meter(self, **kw):
        bar = '-' * self.meter_value
        pad = ' ' * (self.meter_ticks - self.meter_value)
        perc = (float(self.count) / self.total) * 100
        return '[%s>%s] %d%%' % (bar, pad, perc)

    def refresh(self, **kw):
        import sage.server.support
        # Clear line1
        if not EMBEDDED_MODE:
            sys.stdout.write(self.ESC + '[2K')
            self.reset_cursor()
        sys.stdout.write(self.get_meter(**kw))
        # Are we finished?
        if self.count >= self.total:
            sys.stdout.write('\n')
        sys.stdout.flush()
        # Timestamp
        self.last_refresh = time.time()






class coordLabel(pudding.control.SimpleLabel):
   def __init__(self,camera,plot,parent=None):
       pudding.control.SimpleLabel.__init__(self,parent,autosize=True)
       self.camera=camera
       self.plot = plot

   def begin_round(self):
       self.label = "(%s,%s,%s)" % (self.camera.x*(2*self.plot.side)/self.plot.land_size-self.plot.side+self.plot.p[0],self.camera.z*(2*self.plot.side)/self.plot.land_size-self.plot.side+self.plot.p[1],self.camera.y*(2*self.plot.side)/self.plot.land_size)
       self.update()
       self.on_resize()




class MovableCamera(soya.Camera,soya.laser.Laser):
    def __init__(self, parent,plot, light,marker,step=.1,color = (1.0, 0.0, 0.0, 1.0), reflect = 0,pointer=False):
        soya.Camera.__init__(self, parent)
        soya.laser.Laser.__init__(self, parent,color,reflect)
        self.light=light
        self.marker=marker
        self.speed = soya.Vector(self)
        self.rotation_lateral_speed  = 0.0
        self.rotation_vertical_speed = 0.0
        self.rotation_incline_speed = 0.0
        self.impact="none"
        self.impact_x = 0
        self.impact_y = 0
        self.impact_z = 0
        self.direc=soya.Vector(self)
        self.direc.x=0.0
        self.direc.y=0.0
        self.direc.z=-1.0
        self.step=step
        self.plot=plot
        self.pointer=pointer

    def begin_round(self):
        soya.Camera.begin_round(self)
        raypicker = self.get_root()
        self.impact = raypicker.raypick(self.position(), self.direc)

        for event in soya.process_event():
        #  if event[0] == soya.sdlconst.MOUSEMOTION:

        # For mouse motion event, rotate the laser toward the mouse.
        # The formulas are empirical; see soya.cursor for a better algorithm
        # if you want to translate mouse positions into 3D coordinates.

            #	mouse = soya.Point(
             # 		scene,
              #		(float(event[1]) / camera.get_screen_width () - 0.5) *  4.0,
              #		(float(event[2]) / camera.get_screen_height() - 0.5) * -4.0,
              #		0.0,
              #		)
               #     self.look_at(mouse)
            if event[0] == soya.sdlconst.KEYDOWN:
                if   event[1] == soya.sdlconst.K_UP:     self.speed.z = -1.0
                elif event[1] == soya.sdlconst.K_DOWN:   self.speed.z =  1.0
                elif event[1] == soya.sdlconst.K_j:      self.speed.y = -1.0
                elif event[1] == soya.sdlconst.K_k:      self.speed.y =  1.0
                elif event[1] == soya.sdlconst.K_LEFT:   self.rotation_lateral_speed =  3.0
                elif event[1] == soya.sdlconst.K_RIGHT:  self.rotation_lateral_speed = -3.0
                elif event[1] == soya.sdlconst.K_q:      soya.IDLER.stop()
                elif event[1] == soya.sdlconst.K_ESCAPE: soya.IDLER.stop()
                elif event[1] == soya.sdlconst.K_w:       self.rotation_vertical_speed =  3.0
                elif event[1] == soya.sdlconst.K_s:       self.rotation_vertical_speed = -3.0
                elif event[1] == soya.sdlconst.K_f:      soya.toggle_wireframe()
                elif event[1] == soya.sdlconst.K_a:      self.rotation_incline_speed  = 3.0
                elif event[1] == soya.sdlconst.K_d:       self.rotation_incline_speed  = -3.0
                elif event[1] == soya.sdlconst.K_c:
                    if self.impact:
                       self.impact_x = self.impact[0].x
                       self.impact_y = self.impact[0].y
                       self.impact_z = self.impact[0].z
                       if self.pointer:
                           self.marker.x=self.impact[0].x
                           self.marker.y=self.impact[0].y
                           self.marker.z=self.impact[0].z
                           self.light.set_xyz(self.impact[0].x,self.impact[0].y+self.plot.m*self.plot.land_size/(2*self.plot.side),self.impact[0].z)
            if event[0] == soya.sdlconst.KEYUP:
                if   event[1] == soya.sdlconst.K_UP:     self.speed.z = 0.0
                elif event[1] == soya.sdlconst.K_DOWN:   self.speed.z = 0.0
                elif event[1] == soya.sdlconst.K_j:      self.speed.y = 0.0
                elif event[1] == soya.sdlconst.K_k:      self.speed.y = 0.0
                elif event[1] == soya.sdlconst.K_LEFT:   self.rotation_lateral_speed = 0.0
                elif event[1] == soya.sdlconst.K_RIGHT:  self.rotation_lateral_speed = 0.0
                elif event[1] == soya.sdlconst.K_w:       self.rotation_vertical_speed = 0.0
                elif event[1] == soya.sdlconst.K_s:       self.rotation_vertical_speed = 0.0
                elif event[1] == soya.sdlconst.K_a:      self.rotation_incline_speed = 0.0
                elif event[1] == soya.sdlconst.K_d:      self.rotation_incline_speed = 0.0


    def advance_time(self, proportion):
        self.add_mul_vector(proportion*self.step*self.plot.land_size, self.speed)
        self.turn_lateral (self.rotation_lateral_speed  * proportion)
        self.turn_vertical(self.rotation_vertical_speed * proportion)
        self.turn_incline(self.rotation_incline_speed*proportion)


class targetLabel(pudding.control.SimpleLabel):
   def __init__(self,camera,plot,parent=None):
       pudding.control.SimpleLabel.__init__(self,parent,autosize=True)
       self.camera=camera
       self.plot=plot
   def begin_round(self):
       self.label = "(%s,%s,%s)" % ((self.camera.impact_x)*(2*self.plot.side)/self.plot.land_size-self.plot.side+self.plot.p[0],self.camera.impact_z*(2*self.plot.side)/self.plot.land_size-self.plot.side+self.plot.p[1],(self.camera.impact_y)*(2*self.plot.side)/self.plot.land_size)
#       self.label=str(self.camera.impact)
       self.update()
       self.on_resize()





class plot3dsoya:
    """
    A 3d plot object.

    Use the show method to view it.
    """
    def __init__(self,f,p,side,N=10,fineness=25, ignore_bad_values=True):
        self.sampled = False
        self.p = p
        self.f = f
        self.side = float(side)
        self.N = float(N)
        self.fineness = float(fineness)
        #self.land_size = 2**(int(math.ceil(math.log(self.N*2*self.side)/math.log(2))))  + 1
        self.land_size = 2**(int(math.ceil(math.log(self.N)/math.log(2)))) + 1
        self.ignore_bad_values = ignore_bad_values
        from sage.matrix.all import MatrixSpace
        from sage.rings.all import RR
        M = MatrixSpace(RR,self.land_size)
        self.A = M(0)

        try:
            reload(soya)
            reload(soya.sdlconst)

        except ImportError, msg:
            raise RuntimeError, "%s\nYou must install soya3d first."%msg
        SOYA_DATA = '%s/soya3d/data'%SAGE_DATA
        soya.path.append(SOYA_DATA)

    def sample(self):
        """
        Call this to sample the function at the points on the grid.
        Once this has been called it does not have to be called again.

        The show function calls this implicitly if it hasn't already
        been called.
        """
        Total = self.land_size*self.land_size
        p = ProgressMeter(total=Total)

        self.m = abs(self.f(self.p[0]-self.side,self.p[1]-self.side))

        for i in range(self.land_size):

          for j in range(self.land_size):
              try:
                  self.A[i,j] = float(self.f(self.p[0]-self.side + i*(2*self.side)/self.land_size,
                                             self.p[1]-self.side+j*(2*self.side)/self.land_size))
              except (TypeError, ValueError), msg:
                  if self.ignore_bad_values:
                      self.A[i,j] = 0.0
                  else:
                      raise ValueError, msg

              self.m = max(self.m, abs(self.f(self.p[0]-self.side + i*(2*self.side)/self.land_size,
                                              self.p[1]-self.side +j*(2*self.side)/self.land_size)))
              p.update(1)
        self.sampled = True


    def show(self, step=.1,pointer=False):
        """
        Display a plot of this 3d function using soya3d.

        INPUT:
            step -- (default: .1), float, determines how fast you move when flying
                    around the plot of the function.  This can also be set using
                    the show command.
            pointer -- (default: False), bool, whether or not hitting the c key
                    displays a cursor at the pointer you're facing on the surface,
                    along with its coordinates.
        """
        if not self.sampled:
            self.sample()

        try:
            step = float(step)
        except TypeError:
            raise TypeError, "step (=%s) must be coercible to a float"%step

        try:
            pointer = bool(pointer)
        except TypeError:
            raise TypeError, "pointer (=%s) must be coercible to a boolean"%pointer

        soya.init()
        pudding.init()
        scene = soya.World()
        soya.set_root_widget(pudding.core.RootWidget())
        land = soya.Land(scene, self.land_size, self.land_size)

        scale_height=self.m*self.land_size/(2*self.side)

        for i in range(self.land_size):
            for j in range(self.land_size):
                land.set_height(i,j,self.A[i,j])

        land.multiply_height(self.land_size/(2*self.side))

        color={}
        n=0
        quot=(self.fineness*3)

        if(False):
            ground = soya.Material.get("ground")
            grass  = soya.Material.get("grass")
            land.set_material_layer(grass, -self.m*self.land_size/(2*side), m*self.land_size/(2*self.side))
        else:
            for g in range(int(self.fineness)):
                color[(0,(g+1),0)]= soya.Material()
                color[(0,(g+1),0)].diffuse=(float(.1*(g+1)/self.fineness),float(.6*(g+1)/self.fineness),float(0.7),1)
                color[(0,(g+1), 0)].shininess=0


                land.set_material_layer(color[(0,(g+1),0)],-scale_height+2*scale_height*n/quot,-scale_height+2*scale_height*(n+1)/quot)

                n=n+1



            for r in range(int(self.fineness)):
                color[((r+1),0,0)]= soya.Material()
                color[((r+1),0,0)].diffuse=(float(.1+.6*(r+1)/self.fineness),.6,float(.6*(1-(r+1)/self.fineness)),1)
                color[((r+1),0,0)].shininess=0


                land.set_material_layer(color[((r+1),0,0)],-scale_height+2*scale_height*n/quot,-scale_height+2*scale_height*(n+1)/quot)

                n=n+1



            for b in range(int(self.fineness)):
                color[(0,0,(b+1))]= soya.Material()
                color[(0,0,(b+1))].diffuse=(.7,float(.6*(1-(b+1)/self.fineness)),float(.1*(b+1)/self.fineness),1)
                color[(0,0,(b+1))].shininess=0

                land.set_material_layer(color[(0,0,(b+1))],-scale_height+2*scale_height*n/quot,-scale_height+2*scale_height*(n+1)/quot)

                n=n+1

        light1 = soya.Light(scene)
        light1.set_xyz((self.side)*self.land_size/(2*self.side), 2*self.m*self.land_size/(2*self.side), (self.side)*self.land_size/(2*self.side))
        light6  = soya.Light(scene)
        light6.set_xyz(0,0,0)
        cube = None
        if pointer:
            cube = soya.Volume(scene, soya.Shape.get("cube"))
        camera = MovableCamera(scene,self,light6, cube, step=step,pointer=pointer)
        camera.set_xyz((self.side)*self.land_size/(2*self.side), 1.1*self.m*self.land_size/(2*self.side), (self.side)*self.land_size/(2*self.side))
        camera.look_at(soya.Point(scene, (self.side)*self.land_size/(2*self.side), self.m*self.land_size/(2*self.side), (self.side)*self.land_size/(2*self.side)))

        camera.front=.001*2*self.land_size/self.side
        camera.back=100*2*self.land_size/self.side

        soya.root_widget.add_child(camera)
        text1 = coordLabel(camera,self,soya.root_widget)
        if pointer:
            text2 = targetLabel(camera,self,soya.root_widget)
            text2.set_pos_bottom_right(bottom = 10)


        soya.Idler().stop()

        soya.Idler(scene).idle()
        soya.quit()

