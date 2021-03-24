import numpy
import math
from pylab import *
from sympy import *
from scipy.optimize import minimize_scalar
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import axes3d, Axes3D

# Function to be minimized
z_str = '3 * a[0] ** 2 + a[1] ** 2 - a[0] * a[1] - 4 * a[0]'

# Formula from string
exec('z = lambda a: ' + z_str)

z_str = z_str.replace('a[0]', 'x')
z_str = z_str.replace('a[1]', 'y')
step = 0

# GRAD(z(x,y) in point (a[0], a[1]))
x = Symbol('x')
y = Symbol('y')
z_d = eval(z_str)
z_prime_by_y = z_d.diff(y)
z_prime_by_x = z_d.diff(x)
def z_grad(a):
    dif_y = str(z_prime_by_y).replace('y', str(a[1]))
    dif_y = dif_y.replace('x', str(a[0]))

    # Partial derivative
    # yprime = z_d.diff(x)
    dif_x = str(z_prime_by_x).replace('y', str(a[1]))
    dif_x = dif_x.replace('x', str(a[0]))

    return numpy.array([eval(dif_y), eval(dif_x)])


def mininize(a):
    # .x returns the solution of optimization minimization problem
    print("На вход r поступают начальные приближения alpha и beta:", a)
    print("На вход r поступает градиент :", z_grad(a))
    global step
    if step < 2:
        step += 1
        l_min = minimize_scalar(lambda l: math.fabs(z(a - l * z_grad(a)))).x
        return a - l_min * z_grad(a)
    direction = a[-1]-a[0]
    l_min = minimize_scalar(lambda l: math.fabs(z(a - l * direction))).x
    return a - l_min * direction


def norm(a):
    return math.sqrt(a[0] ** 2 + a[1] ** 2)


def grad_step(dot):
    return mininize(dot)


dot = [numpy.array([-150.0, 150.0])]
dot.append(grad_step(dot[0]))

eps = 0.0001

while norm(dot[-2] - dot[-1]) > eps:
    dot.append(grad_step(dot[-1]))


def experimentalFunc(_x, _y):
    return 3 * _x ** 2 + _y ** 2 - _x * _y - 4 * _x


def makeData ():
    x = numpy.arange(-200, 200, 1.0)
    y = numpy.arange(-200, 200, 1.0)
    # GRID production from two arrays
    xgrid, ygrid = numpy.meshgrid(x, y)
    print("xgrid[0]", xgrid[0])
    # zgrid = z([xgrid, ygrid])
    zgrid = experimentalFunc(xgrid, ygrid)
    print("Z_GRID:", zgrid)
    return xgrid, ygrid, zgrid


xt, yt, zt = makeData()

# Create a new figure.
fig = plt.figure()
# The Axes contains most of the figure elements: Axis, Tick, Line2D, Text, Polygon, etc., and sets the coordinate system.
ax = plt.axes(projection='3d')
# Create a surface plot.
ax.plot_surface(xt, yt, zt, cmap=cm.cool)
ax.plot([x[0] for x in dot], [x[1] for x in dot], [z(x) for x in dot], color='#000000')

print("amount of steps:", len(dot)-1)

plt.show()