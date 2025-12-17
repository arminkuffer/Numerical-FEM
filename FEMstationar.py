import numpy as np
from lagrange2D import linear_quadrilateral_shape_functions as linquadref, linear_quadrilateral_shape_function_derivatives as linquadderivref
from triagplot import plot_quadrilateral_mesh as quadplot
from integration import get_gauss_points_2d_reference as gx2dref, get_gauss_weights_2d_reference as gw2dref, compute_jacobian as getJacobian
import ode as ode
import time as Time
import linear_solver as linsolve

# Variables
thermal_conductivity = 48
high_temperature = 600
low_temperature = 300
radius = 0.02 + 0.06  # Minimum radius for T15 < 450 (max value for y = height)
width = 0.3
height = 0.3
nodes = np.array([[0, 0],
    [width / 3, 0],
    [(2 / 3) * width, 0],
    [width, 0],
    [0, height / 3],
    [width / 3, height / 3],
    [(2 / 3) * width, height / 3],
    [width, height / 3],
    [0, (2 * height) / 3],
    [width / 3, (2 * height) / 3],
    [(2 / 3) * width, (2 * height) / 3],
    [width - radius * np.sin(np.pi / 6), height - radius * np.cos(np.pi / 6)],
    [width, height - radius],
    [width - radius * np.cos(np.pi / 6), height - radius * np.sin(np.pi / 6)],
    [0, height],
    [width / 3, height],
    [width / 2, height],
    [width - radius, height]])

elements = np.array([[1, 2, 6, 5],
    [2, 3, 7, 6],
    [3, 4, 8, 7],
    [5, 6, 10, 9],
    [6, 7, 11, 10],
    [11, 7, 12, 14],
    [7, 8, 13, 12],
    [9, 10, 16, 15],
    [10, 11, 17, 16],
    [11, 14, 18, 17]])

dirichlet_boundary_conditions = np.array([[1, high_temperature],
    [2, high_temperature],
    [3, high_temperature],
    [4, high_temperature],
    [12, low_temperature],
    [13, low_temperature],
    [14, low_temperature],
    [18, low_temperature]])

# Functions:
def compute_element_matrix(element_nodes, gauss_points, gauss_weights):
    element_matrix = np.zeros((4, 4))
    for k in range(len(gauss_points)):
        jacobian, det_jacobian, inv_jacobian = getJacobian(element_nodes, gauss_points[k][0], gauss_points[k][1])
        for i in range(element_matrix.shape[0]):
            for j in range(element_matrix.shape[1]):
                grad_shape_i = np.matmul(linquadderivref(gauss_points[k][0], gauss_points[k][1])[i], inv_jacobian)
                grad_shape_j = np.matmul(linquadderivref(gauss_points[k][0], gauss_points[k][1])[j], inv_jacobian)
                element_matrix[i][j] += thermal_conductivity * (grad_shape_i @ grad_shape_j) * det_jacobian * gauss_weights[k]
    return element_matrix, np.zeros(4)

def assemble_global_system(element_matrix, element_vector, global_matrix, global_rhs, element_nodes):
    for i in range(4):
        global_rhs[int(element_nodes[i] - 1)] += element_vector[i]
        for j in range(4):
            global_matrix[int(element_nodes[i] - 1)][int(element_nodes[j] - 1)] += element_matrix[i][j]
    return global_matrix, global_rhs

def apply_dirichlet_boundary_conditions(global_matrix, global_rhs, boundary_conditions):
    for i in range(boundary_conditions.shape[0]):
        global_rhs[int(boundary_conditions[i][0] - 1)] = boundary_conditions[i][1]
        for k in range(global_matrix.shape[1]):
            if k == int(boundary_conditions[i][0] - 1):
                global_matrix[int(boundary_conditions[i][0] - 1)][k] = 1
            else:
                global_matrix[int(boundary_conditions[i][0] - 1)][k] = 0
    return global_matrix, global_rhs

def solve_system(node_coordinates, element_indices, boundary_conditions):
    element_evaluation = []
    assembled_system = []
    for element in element_indices:
        element_evaluation = compute_element_matrix(node_coordinates[element - 1], gx2dref(3), gw2dref(3))
        if assembled_system == []:
            assembled_system = assemble_global_system(
                element_evaluation[0],
                element_evaluation[1],
                np.zeros((node_coordinates.shape[0], node_coordinates.shape[0])),
                np.zeros((node_coordinates.shape[0],)),
                element,
            )
        else:
            assembled_system = assemble_global_system(
                element_evaluation[0],
                element_evaluation[1],
                assembled_system[0],
                assembled_system[1],
                element,
            )
    assigned_system = apply_dirichlet_boundary_conditions(
        assembled_system[0], assembled_system[1], boundary_conditions
    )
    solution = linsolve.solve_gaussian_elimination(assigned_system[0], assigned_system[1])
    solution = np.array(solution).flatten()
    return solution

start_time = Time.time()
solution = solve_system(nodes, elements, dirichlet_boundary_conditions)
end_time = Time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
print(solution)
quadplot(nodes, elements - 1, solution)


