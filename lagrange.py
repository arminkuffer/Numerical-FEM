import numpy as np
import matplotlib.pyplot as plt

def lagrange_basis(x, order, i, x_nodes):
    basis = 1
    for k in range(order + 1):
        if k != i:
            basis *= (x - x_nodes[k]) / (x_nodes[i] - x_nodes[k])
    return basis

def lagrange_polynomial(x, order, x_nodes, f_nodes):
    polynomial = 0
    for i in range(order + 1):
        polynomial += np.multiply(f_nodes[i], lagrange_basis(x, order, i, x_nodes))
    return polynomial

def forward_difference(h):
    polynomial_h = lagrange_polynomial(x + h, o, x_nodes, f_nodes)
    polynomial = lagrange_polynomial(x, o, x_nodes, f_nodes)
    derivative = (polynomial_h - polynomial) / h
    return derivative

def lagrange_basis_derivative(x, order, i, x_nodes):
    derivative = 0
    for m in range(order + 1):
        product = 1
        for k in range(order + 1):
            if k != i and k != m:
                product *= (x - x_nodes[k]) / (x_nodes[i] - x_nodes[k])
        if m != i:
            derivative += (1 / (x_nodes[i] - x_nodes[m])) * product
    return derivative

def lagrange_polynomial_derivative(x, order, x_nodes, f_nodes):
    polynomial_derivative = 0
    for i in range(order + 1):
        polynomial_derivative += np.multiply(f_nodes[i], lagrange_basis_derivative(x, order, i, x_nodes))
    return polynomial_derivative

def plot_lagrange(polynomial, derivative):
    plt.plot(x, polynomial, color="yellow", linewidth=1)
    plt.plot(x_nodes, f_nodes, 'o', color='black')
    plt.plot(x, derivative, color="green", linewidth="1")
    plt.plot(x, fx, color="black", linewidth=1)
    plt.title('Lagrange Polynomial of Order 4')
    plt.xlabel('x')
    plt.ylabel('f(x)')
    plt.grid(linewidth=.05)
    plt.show()

##############################################
h = 10**(-6)                                 # Step size
x_nodes = [0, 1, 2, 3, 4]                    # Data points
f_nodes = (x_nodes / (np.add(1, x_nodes)))**5  # Data points
x = np.linspace(0, 4, 500)                   # Function input
fx = (x / (np.add(1, x)))**5                 # Exact output
o = 4                                        # Order
##############################################
print(lagrange_polynomial_derivative(0.6, 4, x_nodes, f_nodes))
plot_lagrange(lagrange_polynomial(x, o, x_nodes, f_nodes), lagrange_polynomial_derivative(x, o, x_nodes, f_nodes))
forward_difference(h)

