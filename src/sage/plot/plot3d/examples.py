"""
EXAMPLES of 3d Plots:

First we draw a spiral of spheres:
   sage: S = sum(Sphere(0.5,color=((i-4)/8.0, 0.5, (4-i)/8.0)).translate((sin(i),cos(i),i)) for i in [-4,-3.5,..4])

Next we draw a framed translucent box, and show it with the spheres together:
   sage: B = Box([1,2,4], color=(0.5,0.5,1), opacity=0.5) + frame3d([1,2,4])
   sage: show(B + S)

"""



