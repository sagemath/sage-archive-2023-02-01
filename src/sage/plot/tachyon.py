r"""
Interface to the Tachyon Ray Tracer

AUTHOR:
    -- John E. Stone (johns@megapixel.com) -- wrote tachyon ray tracer
    -- William Stein -- write interface
    -- Joshua Kantor -- 3d function plotting
    -- Tom Boothby -- stuff

TODO:
   -- currently only spheres, lights, and textures are wrapped.  need to add triangles, etc.
"""

from colorsys import hsv_to_rgb

from sage.interfaces.tachyon import tachyon_rt

from sage.ext.sage_object import SageObject

#from sage.ext import fast_tachyon_routines

import os

from math import modf,fabs

class Tachyon(SageObject):
    """
    A scene the can be rendered using the Tachyon ray tracer.

    Note that the coordinates are by default such that z is up,
    positive y is to the *left* and x is toward you.  This is
    not oriented according to the right hand rule.

    EXAMPLES:
    Three spheres on the coordinate axes:

        sage: t = Tachyon(xres=500,yres=500, camera_center=(2,0,0))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
        sage: t.texture('t3', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
        sage: t.texture('t4', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
        sage: t.sphere((0,0.5,0), 0.2, 't2')
        sage: t.sphere((0.5,0,0), 0.2, 't3')
        sage: t.sphere((0,0,0.5), 0.2, 't4')
        sage: t.save()

    Sphere's along the twisted cubic.
        sage: t = Tachyon(xres=512,yres=512, camera_center=(3,0.3,0))
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2,diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in srange(-1,1,0.05):
        ...    k += 1
        ...    t.sphere((i,i^2-0.5,i^3), 0.1, 't%s'%(k%3))
        ...
        sage: t.save()

    Many random spheres:
        sage: t = Tachyon(xres=512,yres=512, camera_center=(2,0.5,0.5), look_at=(0.5,0.5,0.5), raydepth=4)
        sage: t.light((4,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1.0,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.3, opacity=1.0, color=(0,1.0,0))
        sage: t.texture('t2', ambient=0.2, diffuse=0.7, specular=0.5, opacity=0.7, color=(0,0,1.0))
        sage: k=0
        sage: for i in range(100):
        ...    k += 1
        ...    t.sphere((random(),random(), random()), random()/10, 't%s'%(k%3))
        ...
        sage: t.save()         # long (several seconds)


    Points on an elliptic curve:
        sage: t = Tachyon(camera_center=(5,2,2), look_at=(0,1,0))
        sage: t.light((10,3,2), 0.2, (1,1,1))
        sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(1,0,0))
        sage: t.texture('t1', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,1,0))
        sage: t.texture('t2', ambient=0.1, diffuse=0.9, specular=0.5, opacity=1.0, color=(0,0,1))
        sage: E = EllipticCurve('37a')
        sage: P = E([0,0])
        sage: Q = P
        sage: n = 100
        sage: for i in range(n):   # increase 20 for a better plot
        ...    Q = Q + P
        ...    t.sphere((Q[1], Q[0], ZZ(i)/n), 0.1, 't%s'%(i%3))
        ...
        sage: t.save()
    """
    def __init__(self,
                 xres=350, yres=350,
                 zoom = 1.0,
                 antialiasing = False,
                 aspectratio = 1.0,
                 raydepth = 12,
                 camera_center = (-3, 0, 0),
                 updir = (0, 0, 1),
                 look_at = (0,0,0),
                 viewdir = None,
                 projection = 'PERSPECTIVE'):
        self._xres = xres
        self._yres = yres
        self._zoom = zoom
        self._aspectratio = aspectratio
        self._antialiasing = antialiasing
        self._raydepth = raydepth
        self._camera_center = camera_center
        self._updir = updir
        self._projection = projection
        self._objects = []
        if viewdir is None:
            self._viewdir = [look_at[i] - camera_center[i] for i in range(3)]
        else:
            self._viewdir = viewdir



    def __repr__(self):
        return self.str()

    def save(self, filename='sage.png', verbose=0, block=True, extra_opts=''):
        """
            filename -- (default: 'sage.png')
                       output filename; the extension of
                       the filename determines the type.
                       Supported types include:
                         tga -- 24-bit (uncompressed)
                         bmp -- 24-bit Windows BMP (uncompressed)
                         ppm -- 24-bit PPM (uncompressed)
                         rgb -- 24-bit SGI RGB (uncompressed)
                         png -- 24-bit PNG (compressed, lossless)
            verbose -- integer; (default: 0)
                       0 -- silent
                       1 -- some output
                       2 -- very verbose output

            block -- bool (default: True); if False, run the rendering
                     command in the background.

            extra_opts -- passed directly to tachyon command line.
                     Use tachyon_rt.usage() to see some of the possibilities.
        """
        tachyon_rt(self.str(), filename, verbose, block, extra_opts)

    def show(self, verbose=0, extra_opts=''):
        import sage.server.support
        if sage.server.support.EMBEDDED_MODE:
            i = 0
            while os.path.exists('sage%s.png'%i):
                i += 1
            filename = 'sage%s.png'%i
            self.save(filename, verbose=verbose, extra_opts=extra_opts)
        else:
            raise NotImplementedError


    def _res(self):
        return '\nresolution %s %s\n'%(self._xres, self._yres)

    def _camera(self):
        return """
           camera
              zoom %s
              aspectratio %s
              antialiasing %s
              raydepth %s
              center %s
              viewdir %s
              updir %s
           end_camera
        """%(float(self._zoom), float(self._aspectratio),
             int(self._antialiasing),
             int(self._raydepth),
             tostr(self._camera_center),
             tostr(self._viewdir),
             tostr(self._updir))

    def str(self):
        return """
        begin_scene
        %s
        %s
        %s
        end_scene"""%(
            self._res(),
            self._camera(),
            '\n'.join([x.str() for x in self._objects])
            )

    def light(self, center, radius, color):
        self._objects.append(Light(center, radius, color))

    def texture(self, name, ambient=0.2, diffuse=0.8,
                specular=0.0, opacity=1.0,
                color=(1.0,0.0, 0.5), texfunc=0, phong=0, phongsize=.5, phongtype="PLASTIC"):
        self._objects.append(Texture(name, ambient, diffuse,
                                     specular, opacity, color, texfunc,
                                     phong,phongsize,phongtype))

    def sphere(self, center, radius, texture):
        self._objects.append(Sphere(center, radius, texture))

    def cylinder(self, center, axis, radius, texture):
        self._objects.append(Cylinder(center, axis, radius, texture))

    def fcylinder(self, base, apex, radius, texture):
        self._objects.append(FCylinder(base, apex, radius, texture))

    def triangle(self, vertex_1, vertex_2, vertex_3, texture):
	 self._objects.append(Triangle(vertex_1,vertex_2,vertex_3,texture))

    def plot(self,f,(xmin,xmax),(ymin,ymax),texture,max_var=.1,max_depth=5,initial_depth=3,num_colors=False):
        """
        General caveat, this code needs work, its was my proof of concept of a
        way to construct a mesh with nonstandard sampling density that was a
        true mesh without gaps. At high resolutions it becomes very
        slow. pyrex may help but as me and Tom have discussed algorithmically
        the complexity is bad. At some level of resolution it will be cheaper
        to sample uniformly at a high resolution.

        EXAMPLE:
            sage: t = Tachyon(xres=512,yres=512, camera_center=(4,-4,3),viewdir=(-4,4,-3), raydepth=4)
            sage: t.light((4.4,-4.4,4.4), 0.2, (1,1,1))
            sage: def f(x,y): return(float(math.sin(x*y)))
            sage: t.texture('t0', ambient=0.1, diffuse=0.9, specular=0.1,  opacity=1.0, color=(1.0,0,0))
            sage: t.plot(f,(-4,4),(-4,4),"t0",max_depth=5,initial_depth=3,num_colors=60)  # change initial_dept to 5 for better picture
            sage.: t.show()
        """
        Tachyplot(self,f,(xmin,ymin),xmax-xmin,ymax-ymin,texture,max_var = max_var, max_depth=max_depth, initial_depth=initial_depth,num_colors=num_colors)

    def collect(self, objects):
        """Add a set of objects to the scene from a collection"""
        self._objects.extend(objects)




class Tachyplot:
    def __init__(self,tachyon,f,corner,x_side,y_side,texture,max_var=.1,max_depth=5,initial_depth=3,num_colors=False):
         self.tachyon = tachyon
         self.plot_texture=texture
         self.max_depth=max_depth
         self.max_var=.1
         self.initial_depth=initial_depth
         self.x_side=x_side
         self.y_side=y_side
         self.corner=corner
         from sage.rings.rational_field import RationalField
         self.f=f
         self.mesh_point_list={}
         self.trilist=[]
         self.max=f(corner[0],corner[1])
         self.min=self.max
         self.trilist.extend(self.initialrefine(0,0,0)[1])
         if num_colors==False:
             tachyon.collect(self.trilist)
         else:
             for i in range(num_colors):
                 tachyon.texture('tri_col%d'%i,ambient=.1,diffuse = .9,specular = .3, opacity=1.0,color=hue(float(i)/float(num_colors)))

             i = 0
             #print('constructing color map')

             z_range = self.max - self.min
             if z_range != 0:
                 for t in self.trilist:
                     ave = (.33333)*(t._vertex_1[2]+t._vertex_2[2]+t._vertex_3[2])
                     t._texture = 'tri_col%d'%int(num_colors*(ave-self.min)/z_range)

             tachyon.collect(self.trilist)


#    def naive_non_recursive_refine(self,linear_density)
#        l={}
#        for i in range(linear_density):
#            for j in range(linear_density):
#
#l[(i,j)]=self.mesh_refine_2(self.corner[0]+i*x_side*1.0/linear_density,self.corner[1]+j*x_side*1.0/linear_density,x_side*1.0/linear_density,y_side*1.0/linear_density,depth)
#
#        for i in range(linear_densityr):
#            for j in range(linear_densityr):


    def initialrefine(self,x,y,depth):

        d = .5**depth
        d1 = .5**(depth+1)


        center=(float(self.corner[0]+float(x+d1)*self.x_side),\
                float(self.corner[1]+float(y+d1)*self.y_side))


        p1 = (float(self.corner[0]+float(x*self.x_side)),float(self.corner[1]+float(y*self.y_side)))

        p2 = (float(self.corner[0]+float(x+d)*self.x_side),float(self.corner[1]+float(y*self.y_side)))

        p3 = (float(self.corner[0]+x*self.x_side),float( self.corner[1]+float(y+d)*self.y_side))


        p4 = (float(self.corner[0]+float(x+d)*self.x_side) ,\
              float(self.corner[1]+float(y+d)*self.y_side))



        if depth<self.initial_depth:
            tri_list=[]
            point_list=[]
            list1 = self.initialrefine(x,y,depth+1)
            list2 = self.initialrefine(x+d1,y,depth+1)
            list3 = self.initialrefine(x,y+d1,depth+1)
            list4 = self.initialrefine(x+d1,y+d1,depth+1)
            tri_list.extend(list1[1])
            tri_list.extend(list2[1])
            tri_list.extend(list3[1])
            tri_list.extend(list4[1])
            point_list.extend(list1[0])
            point_list.extend(list2[0])
            point_list.extend(list3[0])
            point_list.extend(list4[0])

            self.refine(point_list,tri_list)
            fast_cull_triangles_python(tri_list,self.trilist,p1[0],p2[0],p1[1],p3[1])
            fast_cull_points_python(point_list,p1[0],p2[0],p1[1],p3[1])

            return[point_list,tri_list]

        else:
            return self.meshrefine(x,y,depth)

    def meshrefine_2(self,x,y,x_side,y_side,depth):
        d =.5**depth
        d1 = .5**(depth+1)

        center=(x+.5*x_side,y+.5*y_side,float(self.f(x+.5*x_side,y+.5*y_side)))

        p1 = (x,y, float(self.f(x,y)))

        p2 = (x+x_side,y,self.f(x+x_side,y))

        p3 = (x,y+y_side,self.f(x,y+y_side))

        p4 = (x+x_side,y+y_side,self.f(x+x_side,y+y_side))

        self.max = max(p1[2],p2[2],p3[2],p4[2],self.max)
        self.min = min(p1[2],p2[2],p3[2],p4[2],self.min)

        t1 = (.5*p1[0]+.5*center[0],.5*p1[1]+.5*center[1],.5*p1[2]+.5*center[2])
        t2 = (.5*p2[0]+.5*center[0],.5*p2[1]+.5*center[1],.5*p2[2]+.5*center[2])
        t3 = (.5*p3[0]+.5*center[0],.5*p3[1]+.5*center[1],.5*p3[2]+.5*center[2])
        t4 = (.5*p4[0]+.5*center[0],.5*p4[1]+.5*center[1],.5*p4[2]+.5*center[2])

#        t1 = self.Q('1/2')*(p4[2]-p1[2])
#        t2 = self.Q('1/2')*(p3[2]-p2[2])
        b1 = -self.max_var<t1[2]-self.f(t1[0],t1[1])<self.max_var
        b2 = -self.max_var<t2[2]-self.f(t2[0],t2[1])<self.max_var
        b3 = -self.max_var<t3[2]-self.f(t3[0],t3[1])<self.max_var
        b4 = -self.max_var<t4[2]-self.f(t4[0],t4[1])<self.max_var

#        if max(p1[2]-center[2],p2[2]-center[2],p3[2]-center[2],p4[2]-center[2])> self.max_var and depth<self.max_depth:
        if b1 and b2 and b3 and b4 and depth<self.max_depth:
            tri_list=[]
            point_list=[]
            list1 = self.meshrefine(x,y,depth+1)
            list2 = self.meshrefine(x+d1,y,depth+1)
            list3 = self.meshrefine(x,y+d1,depth+1)
            list4 = self.meshrefine(x+d1,y+d1,depth+1)
            tri_list.extend(list1[1])
            tri_list.extend(list2[1])
            tri_list.extend(list3[1])
            tri_list.extend(list4[1])
            point_list.extend(list1[0])
            point_list.extend(list2[0])
            point_list.extend(list3[0])
            point_list.extend(list4[0])

            self.refine(point_list,tri_list)
            fast_cull_points_python(point_list,p1[0],p2[0],p1[1],p3[1])
            fast_cull_triangles_python(tri_list,self.trilist,p1[0],p2[0],p1[1],p3[1])
            tri_list1=[]
            tri_list2=[]
            tri_list3=[]
            tri_list4=[]

            for t in tri_list:
                if (t._vertex_1[0] == p1[0] and t._vertex_2[0]==p1[0]) or (t._vertex_1[0]==p1[0] and t._vertex_3[0]==p1[0]) or (t._vertex_2[0]==p1[0] and t._vertex_3[0]==p1[0]):
                    tri_list1.append(t)
                elif (t._vertex_1[0] == p2[0] and t._vertex_2[0]==p2[0]) or (t._vertex_1[0]==p2[0] and t._vertex_3[0]==p2[0]) or (t._vertex_2[0]==p2[0] and t._vertex_3[0]==p2[0]):
                    tri_list2.append(t)

                elif (t._vertex_1[1] == p2[1] and t._vertex_2[1]==p2[1]) or (t._vertex_1[1]==p2[1] and t._vertex_3[1]==p2[1]) or (t._vertex_2[1]==p2[1] and t._vertex_3[1]==p2[1]):
                    tri_list3.append(t)

                elif (t._vertex_1[1] == p3[1] and t._vertex_2[1]==p3[1]) or (t._vertex_1[1]==p3[1] and t._vertex_3[1]==p3[1]) or (t._vertex_2[1]==p3[1] and t._vertex_3[1]==p3[1]):
                    tri_list4.append(t)


            return[point_list,[tri_list1,tri_list2,tri_list3,tri_list4]]


        else:
            point_list=[]
            trilist=[]
            mesh_list=[point_list,trilist]
            mesh_list[1].append(Triangle(p1,center,p2,self.plot_texture))
            mesh_list[1].append(Triangle(p1,center,p3,self.plot_texture))
            mesh_list[1].append(Triangle(p2,center,p4,self.plot_texture))
            mesh_list[1].append(Triangle(p3,center,p4,self.plot_texture))
#            mesh_list[1].append(Triangle(p1,p2,p4,self.plot_texture))
#            mesh_list[1].append(Triangle(p1,p3,p4,self.plot_texture))
            mesh_list[0].append(p1)
            mesh_list[0].append(p2)
            mesh_list[0].append(p3)
            mesh_list[0].append(p4)
            return(mesh_list)


    def meshrefine(self,x,y,depth):
        d =.5**depth
        d1 = .5**(depth+1)

        center=(float(self.corner[0]+(x+d1)*self.x_side),\
                float(self.corner[1]+(y+d1)*self.y_side),\
                float(self.f(self.corner[0]+(x+d1)*self.x_side,\
                self.corner[1]+(y+d1)*self.y_side)))

        p1 = (float(self.corner[0]+float(x*self.x_side)),float(self.corner[1]+float(y*self.y_side)), \
            float(self.f(self.corner[0]+float(x*self.x_side),self.corner[1]+float(y*self.y_side))))

        p2 = (float(self.corner[0]+float(x+d)*self.x_side),float(self.corner[1]+float(y*self.y_side)),\
              float(self.f(self.corner[0]+float(x+d)*self.x_side,self.corner[1]+float(y*self.y_side))))

        p3 = (float(self.corner[0]+x*self.x_side),float( self.corner[1]+float(y+d)*self.y_side), \
              float(self.f(self.corner[0]+x*self.x_side, self.corner[1]+float(y+d)*self.y_side)))


        p4 = (float(self.corner[0]+float(x+d)*self.x_side) ,\
              float(self.corner[1]+float(y+d)*self.y_side), \
              float(self.f( self.corner[0]+float(x+d)*self.x_side ,\
              self.corner[1]+float(y+d)*self.y_side )))


        self.max = max(p1[2],p2[2],p3[2],p4[2],self.max)
        self.min = min(p1[2],p2[2],p3[2],p4[2],self.min)

        t1 = (.5*p1[0]+.5*center[0],.5*p1[1]+.5*center[1],.5*p1[2]+.5*center[2])
        t2 = (.5*p2[0]+.5*center[0],.5*p2[1]+.5*center[1],.5*p2[2]+.5*center[2])
        t3 = (.5*p3[0]+.5*center[0],.5*p3[1]+.5*center[1],.5*p3[2]+.5*center[2])
        t4 = (.5*p4[0]+.5*center[0],.5*p4[1]+.5*center[1],.5*p4[2]+.5*center[2])

#        t1 = self.Q('1/2')*(p4[2]-p1[2])
#        t2 = self.Q('1/2')*(p3[2]-p2[2])
        b1 = -self.max_var<t1[2]-self.f(t1[0],t1[1])<self.max_var
        b2 = -self.max_var<t2[2]-self.f(t2[0],t2[1])<self.max_var
        b3 = -self.max_var<t3[2]-self.f(t3[0],t3[1])<self.max_var
        b4 = -self.max_var<t4[2]-self.f(t4[0],t4[1])<self.max_var

#        if max(p1[2]-center[2],p2[2]-center[2],p3[2]-center[2],p4[2]-center[2])> self.max_var and depth<self.max_depth:
        if b1 and b2 and b3 and b4 and depth<self.max_depth:
            tri_list=[]
            point_list=[]
            list1 = self.meshrefine(x,y,depth+1)
            list2 = self.meshrefine(x+d1,y,depth+1)
            list3 = self.meshrefine(x,y+d1,depth+1)
            list4 = self.meshrefine(x+d1,y+d1,depth+1)
            tri_list.extend(list1[1])
            tri_list.extend(list2[1])
            tri_list.extend(list3[1])
            tri_list.extend(list4[1])
            point_list.extend(list1[0])
            point_list.extend(list2[0])
            point_list.extend(list3[0])
            point_list.extend(list4[0])

            self.refine(point_list,tri_list)
            fast_cull_points_python(point_list,p1[0],p2[0],p1[1],p3[1])
            fast_cull_triangles_python(tri_list,self.trilist,p1[0],p2[0],p1[1],p3[1])

            return[point_list,tri_list]

        else:
            point_list=[]
            trilist=[]
            mesh_list=[point_list,trilist]
            mesh_list[1].append(Triangle(p1,center,p2,self.plot_texture))
            mesh_list[1].append(Triangle(p1,center,p3,self.plot_texture))
            mesh_list[1].append(Triangle(p2,center,p4,self.plot_texture))
            mesh_list[1].append(Triangle(p3,center,p4,self.plot_texture))
#            mesh_list[1].append(Triangle(p1,p2,p4,self.plot_texture))
#            mesh_list[1].append(Triangle(p1,p3,p4,self.plot_texture))
            mesh_list[0].append(p1)
            mesh_list[0].append(p2)
            mesh_list[0].append(p3)
            mesh_list[0].append(p4)
            return(mesh_list)

    def refine(self,points,trilist):
        points.sort(lexi_sort)
        done = []
        ##        print('number of triangles in mesh: %d '%len(self.trilist))
        for p in points:
            if p in done:
                pass
            else:
                for t in trilist:

                    if p[0]==t._vertex_1[0] and p[0]==t._vertex_2[0] and (t._vertex_1[1]<p[1] and p[1]<t._vertex_2[1]):
                        trilist.remove(t)
                        done.append(p)
                        trilist.append(Triangle(t._vertex_1,p,t._vertex_3,self.plot_texture))
                        self.trilist.append(Triangle(p,t._vertex_2,t._vertex_3,self.plot_texture))


                        break


                    elif p[0]==t._vertex_1[0] and p[0]==t._vertex_3[0] and (t._vertex_1[1]<p[1] and p[1]<t._vertex_3[1]):
                        trilist.remove(t)
                        done.append(p)
                        trilist.append(Triangle(t._vertex_1,p,t._vertex_2,self.plot_texture))
                        self.trilist.append(Triangle(p,t._vertex_3,t._vertex_2,self.plot_texture))

                        break

                    elif p[0]==t._vertex_2[0] and p[0]==t._vertex_3[0] and (t._vertex_2[1]<p[1] and p[1]<t._vertex_3[1]):
                        trilist.remove(t)
                        done.append(p)
                        trilist.append(Triangle(t._vertex_2,p,t._vertex_1,self.plot_texture))
                        self.trilist.append(Triangle(p,t._vertex_3,t._vertex_1,self.plot_texture))

                        break

                    elif p[1]==t._vertex_1[1] and p[1]==t._vertex_2[1] and (t._vertex_1[0]<p[0] and p[0]<t._vertex_2[0]):
                        trilist.remove(t)

                        done.append(p)
                        self.trilist.append(Triangle(t._vertex_1,p,t._vertex_3,self.plot_texture))
                        trilist.append(Triangle(p,t._vertex_2,t._vertex_3,self.plot_texture))

                        break

                    elif p[1]==t._vertex_1[1] and p[1]==t._vertex_3[1] and (t._vertex_1[0]<p[0] and p[0]<t._vertex_3[0]):
                        trilist.remove(t)
                        done.append(p)

                        self.trilist.append(Triangle(t._vertex_1,p,t._vertex_2,self.plot_texture))
                        trilist.append(Triangle(p,t._vertex_3,t._vertex_2,self.plot_texture))
                        break

                    elif p[1]==t._vertex_2[1] and p[1]==t._vertex_3[1] and (t._vertex_2[0]<p[0] and p[0]<t._vertex_3[0]):
                        trilist.remove(t)
                        done.append(p)
                        self.trilist.append(Triangle(t._vertex_2,p,t._vertex_1,self.plot_texture))
                        trilist.append(Triangle(p,t._vertex_3,t._vertex_1,self.plot_texture))

                        break

                    elif p[0]==t._vertex_1[0] and p[0]==t._vertex_2[0] and (t._vertex_1[1]>p[1] and p[1]>t._vertex_2[1]):
                        trilist.remove(t)
                        done.append(p)
                        self.trilist.append(Triangle(t._vertex_1,p,t._vertex_3,self.plot_texture))
                        trilist.append(Triangle(p,t._vertex_2,t._vertex_3,self.plot_texture))

                        break

                    elif p[0]==t._vertex_1[0] and p[0]==t._vertex_3[0] and (t._vertex_1[1]>p[1] and p[1]>t._vertex_3[1]):
                        trilist.remove(t)
                        done.append(p)
                        self.trilist.append(Triangle(t._vertex_1,p,t._vertex_2,self.plot_texture))
                        trilist.append(Triangle(p,t._vertex_3,t._vertex_2,self.plot_texture))

                        break

                    elif p[0]==t._vertex_2[0] and p[0]==t._vertex_3[0] and (t._vertex_2[1]>p[1] and p[1]>t._vertex_3[1]):
                        trilist.remove(t)
                        done.append(p)
                        self.trilist.append(Triangle(t._vertex_2,p,t._vertex_1,self.plot_texture))
                        trilist.append(Triangle(p,t._vertex_3,t._vertex_1,self.plot_texture))

                        break

                    elif p[1]==t._vertex_1[1] and p[1]==t._vertex_2[1] and (t._vertex_1[0]>p[0] and p[0]>t._vertex_2[0]):
                        trilist.remove(t)
                        done.append(p)

                        trilist.append(Triangle(t._vertex_1,p,t._vertex_3,self.plot_texture))
                        self.trilist.append(Triangle(p,t._vertex_2,t._vertex_3,self.plot_texture))

                        break

                    elif p[1]==t._vertex_1[1] and p[1]==t._vertex_3[1] and (t._vertex_1[0]>p[0] and p[0]>t._vertex_3[0]):
                        trilist.remove(t)
                        done.append(p)

                        trilist.append(Triangle(t._vertex_1,p,t._vertex_2,self.plot_texture))
                        self.trilist.append(Triangle(p,t._vertex_3,t._vertex_2,self.plot_texture))
                        break

                    elif p[1]==t._vertex_2[1] and p[1]==t._vertex_3[1] and (t._vertex_2[0]>p[0] and p[0]>t._vertex_3[0]):
                        trilist.remove(t)
                        done.append(p)

                        trilist.append(Triangle(t._vertex_2,p,t._vertex_1,self.plot_texture))
                        self.trilist.append(Triangle(p,t._vertex_3,t._vertex_1,self.plot_texture))

                        break



class Light:
    def __init__(self, center, radius, color):
        self._center = center
        self._radius = radius
        self._color = color

    def str(self):
        return """
        light center %s
              rad %s
              color %s
        """%(tostr(self._center), float(self._radius),
             tostr(self._color))

class Texture:
    def __init__(self, name, ambient=0.2, diffuse=0.8,
                 specular=0.0, opacity=1.0,
                 color=(1.0,0.0, 0.5), texfunc=0, phong=0, phongsize=0, phongtype="PLASTIC"):
        self._name = name
        self._ambient = ambient
        self._diffuse = diffuse
        self._specular = specular
        self._opacity = opacity
        self._color = color
        self._texfunc = texfunc
        self._phong = phong
        self._phongsize = phongsize
        self._phongtype = phongtype

    def str(self):
        return """
        texdef %s ambient %s diffuse %s specular %s opacity %s
        phong %s %s phong_size %s
        color %s texfunc %s
        """%(self._name,
             self._ambient,
             self._diffuse,
             self._specular,
             self._opacity,
             self._phongtype,
             self._phong,
             self._phongsize,
             tostr(self._color),
             self._texfunc)

class Sphere:
    def __init__(self, center, radius, texture):
        self._center = center
        self._radius = radius
        self._texture = texture

    def str(self):
        return """
        sphere center %s rad %s %s
        """%(tostr(self._center), float(self._radius), self._texture)


class Cylinder:
    def __init__(self, center, axis, radius, texture):
        self._center = center
        self._axis = axis
        self._radius = radius
        self._texture = texture

    def str(self):
        return """
        cylinder center %s axis %s rad %s %s
        """%(tostr(self._center), tostr(self._axis), float(self._radius), self._texture)


class FCylinder:
    def __init__(self, base, apex, radius, texture):
        self._center = base
        self._axis = apex
        self._radius = radius
        self._texture = texture

    def str(self):
        return """
        fcylinder base %s apex %s rad %s %s
        """%(tostr(self._center), tostr(self._axis), float(self._radius), self._texture)

class Triangle:
      def __init__(self,vertex_1,vertex_2,vertex_3,texture):
	  self._vertex_1 = vertex_1
	  self._vertex_2 = vertex_2
	  self._vertex_3 = vertex_3
	  self._texture = texture


      def str(self):
	  return """
	  TRI
	  V0 %s
	  V1 %s
	  V2 %s
	  %s
	  """%(tostr(self._vertex_1), tostr(self._vertex_2),tostr(self._vertex_3), self._texture)

def tostr(s):
    if isinstance(s, str):
        return s
    return ' %s %s %s '%(float(s[0]), float(s[1]), float(s[2]))



def hue(h, s=1, v=1):
    """
      hue(h,s=1,v=1) where 'h' stands for hue,
      's' stands for saturation, 'v' stands for value.
      hue returns a list of rgb intensities (r, g, b)
      All values are in range 0 to 1.

      INPUT:
         h, s, v -- real numbers between 0 and 1.  Note that
                    if any are not in this range they are automatically
                    normalized to be in this range by reducing them
                    modulo 1.
      OUTPUT:
         A valid RGB tuple.

      EXAMPLES:
        sage: hue(0.6)
        (0.0, 0.40000000000000036, 1.0)

        hue is an easy way of getting a broader
        range of colors for graphics

        sage: p = plot(sin, -2, 2, rgbcolor=hue(0.6))

    """
    h = float(h); s = float(s); v = float(v)
    if h != 1:
        h = modf(h)[0]
        if h < 0:
            h += 1
    if s != 1:
        s = modf(s)[0]
        if s < 0:
            s += 1
    if v != 1:
        v = modf(v)[0]
        if v < 0:
            v += 1
    c = hsv_to_rgb(h, s, v)
    return (float(c[0]), float(c[1]), float(c[2]))

def lexi_sort(p,q):
    if p[1]>q[1]:
        return -1
    elif p[1]<q[1]:
        return 1
    elif p[0]<q[0]:
        return -1
    elif p[0]>q[0]:
        return 1
    else:
        return 0






def fast_cull_triangles_python(trilist,final_list,x_min,x_max, y_min, y_max):

	c_len=len(trilist)

	list = []

        for j in range(c_len):
		vertex_1_x = trilist[j]._vertex_1[0]
		vertex_1_y = trilist[j]._vertex_1[1]
		vertex_2_x = trilist[j]._vertex_2[0]
		vertex_2_y = trilist[j]._vertex_2[1]
		vertex_3_x = trilist[j]._vertex_3[0]
		vertex_3_y = trilist[j]._vertex_3[1]
		Bool1 = x_min < vertex_1_x < x_max
		Bool2 = y_min < vertex_1_y < y_max
		Bool3 = x_min < vertex_2_x < x_max
		Bool4 = y_min < vertex_2_y < y_max
		Bool5 = x_min < vertex_3_x < x_max
		Bool6 = y_min < vertex_3_y < y_max

		if Bool1 and Bool2 and Bool3 and Bool4 and Bool5 and Bool6:
			list.append(trilist[j])
			final_list.append(trilist[j])

	c_len = len(list)
        for j in range(c_len):
		trilist.remove(list[j])



def fast_cull_points_python(pointlist,x_min, x_max, y_min, y_max):

	c_len=len(pointlist)
	list = []

        for j in range(c_len):
		vertex_1_x = pointlist[j][0]
		vertex_1_y = pointlist[j][1]
		Bool1 = x_min < vertex_1_x < x_max
		Bool2 = y_min < vertex_1_y < y_max

		if Bool1 and Bool2:
			list.append(pointlist[j])

	c_len = len(list)
        for j in range(c_len):
		pointlist.remove(list[j])




