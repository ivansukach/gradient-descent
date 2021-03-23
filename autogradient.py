import autograd.numpy as np
from autograd import grad, jacobian

x = np.array([5,3], dtype=float)

def cost(x):
    return x[0]**2 / x[1] - np.log(x[1])

gradient_cost = grad(cost)
jacobian_cost = jacobian(cost)

print("grad r:", gradient_cost(x))
print("jacobian r:", jacobian_cost(np.array([x, x, x])))

gradient_cost(x)
jacobian_cost(np.array([x,x,x]))