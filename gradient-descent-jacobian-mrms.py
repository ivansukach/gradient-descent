import re as regular
from pylab import *
import sympy as sp
import time

start_time = time.time()

size_of_y = 5
p_step_bdf = 3
tau = 0.01

y = [[1, 1, 1, 1, 1], [2, 2, 2, 2, 2], [3, 3, 3, 3, 3]]
t0 = 0
t = [t0, t0+tau, t0+2*tau, t0+3*tau]
c = [-1/3, 3/2, -3, 11/6]

def rhp0(_t, _y):
    return _t*(_y[0]*_y[1] + _y[2]/_y[3])


def rhp1(_t, _y):
    return _t*(_y[0]*_y[2] + _y[1]/_y[3])


def rhp2(_t, _y):
    return _t*(_y[0]*_y[3] + _y[2]/_y[1])


def rhp3(_t, _y):
    return _t*(_y[1]*_y[2] + _y[0]/_y[3])


def rhp4(_t, _y):
    return _t*(_y[1]*_y[3] + _y[0]/_y[2])


def f(_t, _y):
    return [rhp0(_t, _y), rhp1(_t, _y), rhp2(_t, _y), rhp3(_t, _y), rhp4(_t, _y)]


def rhp0_str(_t, _y):
    return _t+"*(("+_y+")[0]*("+_y+")[1] + ("+_y+")[2]/("+_y+")[3])"


def rhp1_str(_t, _y):
    return _t+"*(("+_y+")[0]*("+_y+")[2] + ("+_y+")[1]/("+_y+")[3])"


def rhp2_str(_t, _y):
    return _t+"*(("+_y+")[0]*("+_y+")[3] + ("+_y+")[2]/("+_y+")[1])"


def rhp3_str(_t, _y):
    return _t+"*(("+_y+")[1]*("+_y+")[2] + ("+_y+")[0]/("+_y+")[3])"


def rhp4_str(_t, _y):
    return _t+"*(("+_y+")[1]*("+_y+")[3] + ("+_y+")[0]/("+_y+")[2])"


def f_str(_t, _y):
    return [rhp0_str(_t, _y), rhp1_str(_t, _y), rhp2_str(_t, _y), rhp3_str(_t, _y), rhp4_str(_t, _y)]


def x_formula():
    x_template = ""
    for i in range(0, p_step_bdf):
        x_template += "tau*beta"+str(i)+"*f(t["+str(i)+"], y["+str(i)+"]) - alpha"+str(i)+"*y["+str(i)+"] "
        if i != p_step_bdf-1:
            x_template += '+ '
    return x_template


# Function to be minimized
def r_formula():
    r_template = "tau*f(t[" + str(p_step_bdf) + "], "+x_f+") - c[" + str(p_step_bdf) + "]*"+x_f
    i = p_step_bdf-1
    while i >= 0:
        r_template += "- c["+str(i)+"]*y["+str(i)+"]"
        i -= 1
    return r_template


def x_s():
    x_template = ""
    for i in range(0, p_step_bdf):
        x_template += str(tau)+"*np.array("+str(f(t[i], y[i]))+")*beta"+str(i)+" - alpha"+str(i)+"*np.array("+str(y[i])+")"
        if i != p_step_bdf-1:
            x_template += " + "
    return x_template


# Function to be minimized
def r_s():
    f_s = f_str(str(t[p_step_bdf]), x)
    r_template = str(tau)+"*np.array("+str(f_s)+") - " + str(c[p_step_bdf])+"*"+x
    i = p_step_bdf-1
    while i >= 0:
        r_template += "-" + str(c[i])+"*np.array("+str(y[i])+")"
        i -= 1
    return r_template.replace("'", "")


x_f = x_formula()
r_f = r_formula()
print("x formula:", x_f)
print("r formula:", r_f)
x = x_s()
r = r_s()
print("x:", x)
print("r:", r)
for i in range(0, p_step_bdf):
    exec("alpha" + str(i) + " = sp.Symbol('alpha" + str(i) + "')")
    exec("beta" + str(i) + " = sp.Symbol('beta" + str(i) + "')")
r = eval(r)
jacobian = [[None]*2*p_step_bdf] * size_of_y
for i in range(0, size_of_y):
    for j in range(0, p_step_bdf):
        jacobian[i][j] = eval('sp.diff(r[' + str(i) + '], alpha' + str(j) + ')')
    for j in range(0, p_step_bdf):
        jacobian[i][p_step_bdf+j] = eval('sp.diff(r[' + str(i) + '], beta' + str(j) + ')')
transposed_jacobian = np.array(jacobian).transpose()
print("r after simplify and replacement: ", r)
print("Jacobian: ", jacobian)
gradient = transposed_jacobian.dot(r)
print("Gradient: ", gradient)
print("--- %s seconds ---" % (time.time() - start_time))

