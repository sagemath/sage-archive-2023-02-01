import mplot3d
import pylab as p
delta = 0.025
x = y = p.arange(-3.0, 3.0, delta)
X, Y = p.meshgrid(x,y)

Z1 = p.bivariate_normal(X, Y, 1.0, 1.0, 0.0, 0.0)
Z2 = p.bivariate_normal(X, Y, 1.5, 0.5, 1, 1)
Z = Z2-Z1

X = X * 10
Y = Y * 10
Z = Z * 500

fig = p.figure()
ax = mplot3d.Axes3D(fig, axisbg='k')
ax.plot_surface(X,Y,Z, div=40)

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

fig.add_axes(ax)
p.savefig("surface.png", dpi=100)
p.show()
