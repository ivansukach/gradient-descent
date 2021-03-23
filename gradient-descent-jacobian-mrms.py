import re as regular
from pylab import *
import sympy as sp

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


# def unique_strings(prev):
#     result = [prev[0]]
#     for i in range(1, len(prev)):
#         unique = 1
#         for j in range(1, len(result)):
#             if prev[i] == result[j]:
#                 unique = 0
#         if unique != 0:
#             result.append(prev[i])
#     return result

x_f = x_formula()
r_f = r_formula()
print("x formula:", x_f)
print("r formula:", r_f)
x = x_s()
r = r_s()
print("x:", x)
print("r:", r)




# normOfR = squared_norm_of_r()
# print("Squared norm of r:", normOfR)
# t0 = 0
# t = [t0, t0+tau, t0+2*tau, t0+3*tau]
# c = [-1/3, 3/2, -3, 11/6]
# r_s = r_s.replace('tau', str(tau))
# for i in range(0, len(t)):
#     normOfR = normOfR.replace('t['+str(i)+']', str(t[i]))
# for i in range(0, p_step_bdf):
#     normOfR = normOfR.replace('y[' + str(i) + ']', str(y[i]))
# f_vectors = regular.findall(r'f\([^,]+,\s?\[[^]]+]\)', normOfR)
# print("f vectors:", f_vectors)
# for f_vector in f_vectors:
#     valueOfFVector = eval(f_vector)
#     normOfR = normOfR.replace(f_vector, str(valueOfFVector))
# f_components = regular.findall(r'f\[\d]\([^,]+,\s?[^\)]+\)', normOfR)
# print("f components:", f_components)
# for f_component in f_components:
#     index = regular.findall(r'f\[(\d)]\([^,]+,\s?[^\)]+\)', normOfR)[0]
#     argument1 = regular.findall(r'f\[\d]\(([^,]+),\s?[^\)]+\)', normOfR)[0]
#     argument2 = regular.findall(r'f\[\d]\([^,]+,(\s?[^\)]+)\)', normOfR)[0]
#     # normOfR = normOfR.replace(f_component, rhp0_str(argument1, argument2))
#     rhpStr = eval('rhp'+index+'_str(argument1, argument2)')
#     normOfR = normOfR.replace(f_component, eval('rhp' + index + '_str(argument1, argument2)'))
# vectors = regular.findall(r'\[[^],]+,[^]]+]', normOfR)
# vectors = unique_strings(vectors)
# for vector in vectors:
#     normOfR = normOfR.replace(vector, 'np.array(' + vector + ')')
# print("Squared norm of r after replacement:", normOfR)



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
        # if j == 0:
        #     print("DIFF beta0: ", eval('sp.diff(r[' + str(i) + '], beta' + str(j) + ')'))
        #     print("R[i]: ", eval('r[' + str(i) + ']'))
transposed_jacobian = np.array(jacobian).transpose()
# r_ab = r
# for i in range(0, p_step_bdf):
#     r_ab = r_ab.replace("alpha"+str(i), "ab["+str(i)+"]").replace("beta"+str(i), "ab["+str(p_step_bdf+i)+"]")
# r_fun = eval(r)

print("r after simplify and replacement: ", r)
print("Jacobian: ", jacobian)
gradient = transposed_jacobian.dot(r)
print("Gradient: ", gradient)
matrix_A = np.array([[1, 2, 3], [4, 5, 6]])
vector_B = np.array([7, 8, 9])
print("matrix_A*vector_B: ", matrix_A.dot(vector_B))


# normOfR = str(eval(normOfR))
# print("Squared norm of r after simplify:", normOfR)
#
#
# def gradient():
#     alpha0 = sp.Symbol('alpha0')
#     dif_alpha = sp.diff(normOfR, alpha0)
#     beta0 = sp.Symbol('beta0')
#     dif_beta = sp.diff(normOfR, beta0)
#     alphas = [dif_alpha]
#     betas = [dif_beta]
#     for i in range(1, p_step_bdf):
#         exec("alpha" + str(i) + " = sp.Symbol('alpha" + str(i) + "')")
#         exec("beta" + str(i) + " = sp.Symbol('beta" + str(i) + "')")
#         dif_alpha = eval('sp.diff(normOfR, alpha'+str(i)+')')
#         dif_beta = eval('sp.diff(normOfR, beta'+str(i)+')')
#         alphas.append(dif_alpha)
#         betas.append(dif_beta)
#     return np.asarray(alphas + betas)
#
#
# grad = gradient()
# print("gradient: ", grad)
