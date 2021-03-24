import re as regular
from pylab import *
import sympy as sp
import time

start_time = time.time()
size_of_y = 5
p_step_bdf = 3
tau = 0.01


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


# Function to be minimized
def squared_norm_of_r():
    r_norm_template = ""
    for i in range(0, size_of_y):
        r_norm_template += "(tau*f["+str(i)+"](t[3], tau*beta0*f(t[0], y[0]) - alpha0*y[0] + tau*beta1*f(t[1], y[1]) - alpha1*y[1] + " \
                           + "tau*beta2*f(t[2], y[2]) - alpha2*y[2]) - c[3]*(tau*beta0*f["+str(i)+"](t[0], y[0]) - alpha0*y[0]["+str(i)+"] + " \
                            + "tau*beta1*f["+str(i)+"](t[1], y[1]) - alpha1*y[1]["+str(i)+"] + tau*beta2*f["+str(i)+"](t[2], y[2]) - alpha2*y[2]["+str(i)+"]) - " \
                              + "c[2]*y[2]["+str(i)+"] - c[1]*y[1]["+str(i)+"] - c[0]*y[0]["+str(i)+"])**2"
        if i != size_of_y-1:
            r_norm_template += '+'
    return r_norm_template


def unique_strings(prev):
    result = [prev[0]]
    for i in range(1, len(prev)):
        unique = 1
        for j in range(1, len(result)):
            if prev[i] == result[j]:
                unique = 0
        if unique != 0:
            result.append(prev[i])
    return result


normOfR = squared_norm_of_r()
print("Squared norm of r:", normOfR)
t0 = 0
t = [t0, t0+tau, t0+2*tau, t0+3*tau]
c = [-1/3, 3/2, -3, 11/6]
y = [[1, 1, 1, 1, 1], [2, 2, 2, 2, 2], [3, 3, 3, 3, 3]]
normOfR = normOfR.replace('tau', str(tau))
for i in range(0, len(t)):
    normOfR = normOfR.replace('t['+str(i)+']', str(t[i]))
for i in range(0, p_step_bdf):
    normOfR = normOfR.replace('y[' + str(i) + ']', str(y[i]))
f_vectors = regular.findall(r'f\([^,]+,\s?\[[^]]+]\)', normOfR)
print("f vectors:", f_vectors)
for f_vector in f_vectors:
    valueOfFVector = eval(f_vector)
    normOfR = normOfR.replace(f_vector, str(valueOfFVector))
f_components = regular.findall(r'f\[\d]\([^,]+,\s?[^\)]+\)', normOfR)
print("f components:", f_components)
for f_component in f_components:
    index = regular.findall(r'f\[(\d)]\([^,]+,\s?[^\)]+\)', normOfR)[0]
    argument1 = regular.findall(r'f\[\d]\(([^,]+),\s?[^\)]+\)', normOfR)[0]
    argument2 = regular.findall(r'f\[\d]\([^,]+,(\s?[^\)]+)\)', normOfR)[0]
    # normOfR = normOfR.replace(f_component, rhp0_str(argument1, argument2))
    rhpStr = eval('rhp'+index+'_str(argument1, argument2)')
    normOfR = normOfR.replace(f_component, eval('rhp' + index + '_str(argument1, argument2)'))
vectors = regular.findall(r'\[[^],]+,[^]]+]', normOfR)
vectors = unique_strings(vectors)
for vector in vectors:
    normOfR = normOfR.replace(vector, 'np.array(' + vector + ')')
print("Squared norm of r after replacement:", normOfR)
for i in range(0, p_step_bdf):
    exec("alpha" + str(i) + " = sp.Symbol('alpha" + str(i) + "')")
    exec("beta" + str(i) + " = sp.Symbol('beta" + str(i) + "')")
normOfR = str(eval(normOfR))
print("Squared norm of r after simplify:", normOfR)


def gradient():
    dif_alpha = sp.diff(normOfR, alpha0)
    dif_beta = sp.diff(normOfR, beta0)
    alphas = [dif_alpha]
    betas = [dif_beta]
    for i in range(1, p_step_bdf):
        dif_alpha = eval('sp.diff(normOfR, alpha'+str(i)+')')
        dif_beta = eval('sp.diff(normOfR, beta'+str(i)+')')
        alphas.append(dif_alpha)
        betas.append(dif_beta)
    return np.asarray(alphas + betas)


grad = gradient()
print("gradient: ", grad)
print("--- %s seconds ---" % (time.time() - start_time))
