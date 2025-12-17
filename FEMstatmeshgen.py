import numpy as np
from lagrange2D import linear_quadrilateral_shape_functions as linquadref, linear_quadrilateral_shape_function_derivatives as linquadderivref
from triagplot import plot_quadrilateral_mesh as quadplot
from integration import get_gauss_points_2d_reference as gx2dref, get_gauss_weights_2d_reference as gw2dref, compute_jacobian as getJacobian
from meshgen import meshGen
import ode as ode
import time as Time
import linear_solver as linsolve

# Variables
thermal_conductivity = 48
boundary_temp_high = 600
boundary_temp_low = 300
min_radius = 0.02 + 0.06  # Smallest radius for which T15 < 450 (Max value for y = h)
domain_width = 0.3
domain_height = 0.3

nodes, elements, edges = meshGen(domain_width, domain_height, min_radius, 4)
dirichlet_bc = np.empty((len(edges[2]) + len(edges[4]),), dtype=object)
for i in range(dirichlet_bc.shape[0]):
    if i < len(edges[2]):
        dirichlet_bc[i] = [edges[2][i], boundary_temp_high]
    else:
        dirichlet_bc[i] = [edges[4][i - len(edges[2])], boundary_temp_low]

# Functions:
def evaluate_stat(element_nodes, gpx, gpw):
    element_stiffness_matrix = np.zeros((4, 4))
    for k in range(len(gpx)):
        for i in range(element_stiffness_matrix.shape[0]):
            for j in range(element_stiffness_matrix.shape[1]):
                J, detJ, invJ = getJacobian(element_nodes, gpx[k][0], gpx[k][1])
                gradN_i = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[i], invJ)
                gradN_j = np.matmul(linquadderivref(gpx[k][0], gpx[k][1])[j], invJ)
                element_stiffness_matrix[i][j] += thermal_conductivity * (gradN_i @ gradN_j) * detJ * gpw[k]
    return element_stiffness_matrix, np.zeros(4)

def assemble(element_stiffness_matrix, element_load_vector, global_system_matrix, global_rhs_vector, element_node_indices):
    for i in range(4):
        global_rhs_vector[int(element_node_indices[i] - 1)] += element_load_vector[i]
        for j in range(4):
            global_system_matrix[int(element_node_indices[i] - 1)][int(element_node_indices[j] - 1)] += element_stiffness_matrix[i][j]
    return global_system_matrix, global_rhs_vector

def assignDBC(global_system_matrix, global_rhs_vector, dirichlet_bc):
    for i in range(dirichlet_bc.shape[0]):
        global_rhs_vector[int(dirichlet_bc[i][0] - 1)] = dirichlet_bc[i][1]
        for k in range(global_system_matrix.shape[1]):
            if k == int(dirichlet_bc[i][0] - 1):
                global_system_matrix[int(dirichlet_bc[i][0] - 1)][k] = 1
            else:
                global_system_matrix[int(dirichlet_bc[i][0] - 1)][k] = 0
    return global_system_matrix, global_rhs_vector

def solve(nodes, elements, dirichlet_bc):
    element_evaluation = []
    assembled_system = []
    for element_node_indices in elements:
        element_evaluation = evaluate_stat(nodes[element_node_indices - 1], gx2dref(3), gw2dref(3))
        if assembled_system == []:
            assembled_system = assemble(
                element_evaluation[0],
                element_evaluation[1],
                np.zeros((nodes.shape[0], nodes.shape[0])),
                np.zeros((nodes.shape[0], 1)),
                element_node_indices,
            )
        else:
            assembled_system = assemble(
                element_evaluation[0],
                element_evaluation[1],
                assembled_system[0],
                assembled_system[1],
                element_node_indices,
            )
    assigned_system = assignDBC(assembled_system[0], assembled_system[1], dirichlet_bc)
    solution_vector = np.linalg.solve(assigned_system[0], assigned_system[1])
    solution_vector = np.array(solution_vector).flatten()
    return solution_vector

start_time = Time.time()
sol = solve(nodes, elements, dirichlet_bc)
end_time = Time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
print(sol)
quadplot(nodes, elements - 1, sol)