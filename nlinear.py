import numpy as np
from lagrange2D import linear_quadrilateral_shape_function_derivatives as linquadderivref, linear_quadrilateral_shape_functions as linquadref
from triagplot import plot_quadrilateral_mesh as quadplot
from integration import get_gauss_points_2d_reference as gx2dref, get_gauss_weights_2d_reference as gw2dref, compute_jacobian as getJacobian

# Constants
THERMAL_CONDUCTIVITY = 48
TEMP_BOUNDARY_1 = 600
TEMP_BOUNDARY_2 = 300
RADIUS = 0.1  # Minimum radius for T15 < 450 (max value for y = HEIGHT)
WIDTH = 0.3
HEIGHT = 0.3
EXP_COEFF_1 = 10**6 
EXP_COEFF_2 = 10**3

# Nodes and elements
define_nodes = np.array([[0, 0],
                  [WIDTH / 3, 0],
                  [(2 / 3) * WIDTH, 0],
                  [WIDTH, 0],
                   [0, HEIGHT / 3],
                   [WIDTH / 3, HEIGHT / 3],
                   [(2 / 3) * WIDTH, HEIGHT / 3],
                   [WIDTH, HEIGHT / 3],
                   [0, (2 * HEIGHT) / 3],
                   [WIDTH / 3, (2 * HEIGHT) / 3],
                   [(2 / 3) * WIDTH, (2 * HEIGHT) / 3],
                   [WIDTH - RADIUS * np.sin(np.pi / 6), HEIGHT - RADIUS * np.cos(np.pi / 6)],
                   [WIDTH, HEIGHT - RADIUS],
                   [WIDTH - RADIUS * np.cos(np.pi / 6), HEIGHT - RADIUS * np.sin(np.pi / 6)],
                   [0, HEIGHT],
                   [WIDTH / 3, HEIGHT],
                   [WIDTH / 2, HEIGHT],
                   [WIDTH - RADIUS, HEIGHT]])

define_elements = np.array([[1, 2, 6, 5],
                     [2, 3, 7, 6],
                     [3, 4, 8, 7],
                     [5, 6, 10, 9],
                      [6, 7, 11, 10],
                      [11, 7, 12, 14],
                      [7, 8, 13, 12],
                      [9, 10, 16, 15],
                      [10, 11, 17, 16],
                      [11, 14, 18, 17]])

boundary_conditions = np.array([[1, TEMP_BOUNDARY_1],
                [2, TEMP_BOUNDARY_1],
                [3, TEMP_BOUNDARY_1],
                [4, TEMP_BOUNDARY_1],
                [12, TEMP_BOUNDARY_2],
                [13, TEMP_BOUNDARY_2],     
                [14, TEMP_BOUNDARY_2],
                [18, TEMP_BOUNDARY_2]])

initial_temperature = temperature = np.array([600] * 4 + [250] * 14)

# Functions

def evaluate_stationary_nonlinear(element_nodes, gauss_points, gauss_weights, temperature):
    element_matrix = np.zeros((4, 4))
    element_vector = np.zeros((4,))
    jacobian_temperature = np.zeros((4, 4))
    for k in range(len(gauss_points)):
        jacobian, det_jacobian, inv_jacobian = getJacobian(element_nodes, gauss_points[k][0], gauss_points[k][1])
        shape_functions = linquadref(gauss_points[k][0], gauss_points[k][1])
        for i in range(element_matrix.shape[0]):
            temperature_sum = 0
            for j in range(element_matrix.shape[1]):
                temperature_sum += shape_functions[j] * temperature[j]
            for j in range(element_matrix.shape[1]):
                grad_shape_i = np.dot(linquadderivref(gauss_points[k][0], gauss_points[k][1])[i], inv_jacobian)
                grad_shape_j = np.dot(linquadderivref(gauss_points[k][0], gauss_points[k][1])[j], inv_jacobian)
                element_matrix[i][j] += THERMAL_CONDUCTIVITY * (grad_shape_i @ grad_shape_j) * det_jacobian * gauss_weights[k]
                jacobian_temperature[i][j] += THERMAL_CONDUCTIVITY * (grad_shape_i @ grad_shape_j) * det_jacobian * gauss_weights[k] - shape_functions[i] * shape_functions[j] * (
                    EXP_COEFF_1 * EXP_COEFF_2
                ) / (temperature_sum**2) * np.exp(-(EXP_COEFF_2) / (temperature_sum)) * det_jacobian * gauss_weights[k]
            element_vector[i] += shape_functions[i] * EXP_COEFF_1 * np.exp(-(EXP_COEFF_2 / temperature_sum)) * det_jacobian * gauss_weights[k]
    residual_force = -(np.dot(element_matrix, temperature) - element_vector)
    return jacobian_temperature, residual_force

def assemble_global_system(element_matrix, element_vector, global_matrix, global_rhs, element):
    for i in range(len(element)):
        global_rhs[int(element[i] - 1)] += element_vector[i]
        for j in range(len(element)):
            global_matrix[int(element[i] - 1)][int(element[j] - 1)] += element_matrix[i][j]
    return global_matrix, global_rhs

def apply_boundary_conditions(global_matrix, global_rhs, boundary_conditions):
    for i in range(boundary_conditions.shape[0]):
        global_rhs[int(boundary_conditions[i][0] - 1)] = 0
        for k in range(global_matrix.shape[1]):
            if k == int(boundary_conditions[i][0] - 1):
                global_matrix[int(boundary_conditions[i][0] - 1)][k] = 1
            else:
                global_matrix[int(boundary_conditions[i][0] - 1)][k] = 0
    return global_matrix, global_rhs

def solve_nonlinear_system(nodes, elements, boundary_conditions, temperature, tolerance, max_iterations):
    residual = np.ones((18,))
    iteration = 0
    while np.linalg.norm(residual) > tolerance and iteration < max_iterations:
        element_evaluations = []
        global_assembly = []
        for element in elements:
            element_evaluations = evaluate_stationary_nonlinear(
                nodes[element - 1], gx2dref(3), gw2dref(3), temperature[element - 1]
            )
            if global_assembly == []:
                global_assembly = assemble_global_system(
                    element_evaluations[0],
                    element_evaluations[1],
                    np.zeros((nodes.shape[0], nodes.shape[0])),
                    np.zeros((nodes.shape[0],)),
                    element,
                )
            else:
                global_assembly = assemble_global_system(
                    element_evaluations[0],
                    element_evaluations[1],
                    global_assembly[0],
                    global_assembly[1],
                    element,
                )
        global_system = apply_boundary_conditions(global_assembly[0], global_assembly[1], boundary_conditions)
        delta_temperature = np.linalg.solve(global_system[0], global_system[1])
        temperature = temperature + delta_temperature
        residual = global_system[1]
        iteration += 1
    return temperature

# Example usage
print(
    evaluate_stationary_nonlinear(
        [[0.0, 0.0], [0.0, 1.0], [1.0, 1.0], [1.0, 0.0]],
        [[0.5, 0.5]],
        [1.0],
        [1.0, 2.0, 2.0, 1.0],
    )
)
solution = solve_nonlinear_system(
    define_nodes, define_elements, boundary_conditions, initial_temperature, 10**-8, 10
)
print(solution)
quadplot(define_nodes, define_elements - 1, solution)