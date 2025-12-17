import numpy as np
from itertools import product
from lagrange2D import linear_quadrilateral_shape_functions as linquadref, linear_quadrilateral_shape_function_derivatives as linquadderivref

def midpoint_rule(a, b):
    return (b - a) * f((a + b) / 2)

def trapezoidal_rule(a, b):
    return 1 / 2 * (b - a) * (f(a) + f(b))

def get_gauss_points_1d(n):
    return points[n - 1]

def get_gauss_weights_1d(n):
    return weights[n - 1]

def f(x):
    return (x**5) / (1 + x)**5

def gaussian_quadrature_1d(x, w, a, b):
    total = 0
    for i in range(len(x)):
        xi = ((b - a) / 2) * x[i] + (a + b) / 2
        total += w[i] * f(xi)
    return (b - a) / 2 * total

def get_gauss_points_2d_reference(n):
    if n == 1:
        return np.array([0, 0])
    else:
        positions = list(product(points[n - 1], repeat=2))
        return np.array(positions)

def get_gauss_weights_2d_reference(n):
    if n == 1:
        return [4.0]
    else:
        weights_2d = list(product(weights[n - 1], repeat=2))
    for i in range(len(weights_2d)):
        weights_2d[i] = weights_2d[i][0] * weights_2d[i][1]
    return weights_2d

def get_position(nodes, xi, eta):
    position = 0
    for i in range(4):  # 4 because there are 4 shape functions for 2D!
        position += linquadref(xi, eta)[i] * nodes[i]
    return position

def compute_jacobian(nodes, xi, eta):
    jacobian_values = np.zeros(3, dtype="object")
    jacobian = np.zeros((2, 2))
    for i in range(4):
        jacobian[0][0] += linquadderivref(xi, eta)[i][0] * nodes[i][0]
        jacobian[0][1] += linquadderivref(xi, eta)[i][1] * nodes[i][0]

        jacobian[1][0] += linquadderivref(xi, eta)[i][0] * nodes[i][1]
        jacobian[1][1] += linquadderivref(xi, eta)[i][1] * nodes[i][1]
    jacobian_values[0] = jacobian
    jacobian_values[1] = np.linalg.det(jacobian)
    if jacobian_values[1] == 0:
        jacobian_values[2] = np.zeros((2, 2))
    else:
        jacobian_values[2] = np.linalg.inv(jacobian)
    return jacobian_values

points = [
    [0],
    [-1 / np.sqrt(3), 1 / np.sqrt(3)],
    [float(-np.sqrt(3 / 5)), 0, float(np.sqrt(3 / 5))],
]
weights = [[2], [1, 1], [5 / 9, 8 / 9, 5 / 9]]

# Example usage
# print(gaussian_quadrature_1d(get_gauss_points_1d(1), get_gauss_weights_1d(1), 0, 4))
# print(gaussian_quadrature_1d(get_gauss_points_1d(2), get_gauss_weights_1d(2), 0, 4))
# print(gaussian_quadrature_1d(get_gauss_points_1d(3), get_gauss_weights_1d(3), 0, 4))
# print(get_gauss_points_2d_reference(3))
# print(get_gauss_weights_2d_reference(3))
# print(get_position(np.array([[2, 1], [4, 1], [4, 3], [2, 2]]), 0.577, -0.577))
# print(compute_jacobian(np.array([[2, 1], [4, 1], [4, 3], [2, 2]]), 0.577, -0.577))
