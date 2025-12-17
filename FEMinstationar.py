import numpy as np
from lagrange2D import linear_quadrilateral_shape_functions as linquadref, linear_quadrilateral_shape_function_derivatives as linquadderivref
from triagplot import plot_quadrilateral_mesh as quadplot
from integration import get_gauss_points_2d_reference as gx2dref, get_gauss_weights_2d_reference as gw2dref, compute_jacobian as getJacobian
import ode as ode
import time as Time
import linear_solver as linsolve


# Variables
material_density = 7800
specific_heat_capacity = 452
thermal_conductivity = 48
initial_temperature = np.array([300, 300, 300, 300])
high_boundary_temperature = 600
low_boundary_temperature = 300
radius = 0.02
domain_width = 0.3
domain_height = 0.3

nodes = np.array([
    [0, 0],
    [domain_width / 3, 0],
    [(2 / 3) * domain_width, 0],
    [domain_width, 0],
    [0, domain_height / 3],
    [domain_width / 3, domain_height / 3],
    [(2 / 3) * domain_width, domain_height / 3],
    [domain_width, domain_height / 3],
    [0, (2 * domain_height) / 3],
    [domain_width / 3, (2 * domain_height) / 3],
    [(2 / 3) * domain_width, (2 * domain_height) / 3],
    [domain_width - radius * np.sin(np.pi / 6), domain_height - radius * np.cos(np.pi / 6)],
    [domain_width, domain_height - radius],
    [domain_width - radius * np.cos(np.pi / 6), domain_height - radius * np.sin(np.pi / 6)],
    [0, domain_height],
    [domain_width / 3, domain_height],
    [domain_width / 2, domain_height],
    [domain_width - radius, domain_height],
])

elements = np.array(
    [
        [1, 2, 6, 5],
        [2, 3, 7, 6],
        [3, 4, 8, 7],
        [5, 6, 10, 9],
        [6, 7, 11, 10],
        [11, 7, 12, 14],
        [7, 8, 13, 12],
        [9, 10, 16, 15],
        [10, 11, 17, 16],
        [11, 14, 18, 17],
    ]
)

dirichlet_boundary_conditions = np.array(
    [
        [1, high_boundary_temperature],
        [2, high_boundary_temperature],
        [3, high_boundary_temperature],
        [4, high_boundary_temperature],
        [12, low_boundary_temperature],
        [13, low_boundary_temperature],
        [14, low_boundary_temperature],
        [18, low_boundary_temperature],
    ]
)

#Functions:
def evaluate_transient(element_nodes, gauss_points, gauss_weights, current_solution, previous_solution, time_integration_method, time_step, theta, is_first_step):
    element_stiffness_matrix = np.zeros((4, 4))
    element_load_vector = np.zeros((4,), dtype=object)
    if is_first_step == 1:
        time_integration_method = 1
    for k in range(len(gauss_points)):
        jacobian, determinant_jacobian, inverse_jacobian = getJacobian(element_nodes, gauss_points[k][0], gauss_points[k][1])
        for i in range(element_stiffness_matrix.shape[0]):
            for j in range(element_stiffness_matrix.shape[1]):
                shape_function_i = linquadref(gauss_points[k][0], gauss_points[k][1])[i]
                shape_function_j = linquadref(gauss_points[k][0], gauss_points[k][1])[j]
                gradient_shape_function_i = np.matmul(linquadderivref(gauss_points[k][0], gauss_points[k][1])[i], inverse_jacobian)
                gradient_shape_function_j = np.matmul(linquadderivref(gauss_points[k][0], gauss_points[k][1])[j], inverse_jacobian)
                mass_matrix = material_density * specific_heat_capacity * (shape_function_i * shape_function_j) * determinant_jacobian * gauss_weights[k]
                damping_matrix = -(thermal_conductivity * (gradient_shape_function_i @ gradient_shape_function_j) * determinant_jacobian * gauss_weights[k])
                stiffness_matrix = 0
                if time_integration_method == 1:
                    result = ode.OST(theta, time_step, np.array(mass_matrix), np.array([damping_matrix, damping_matrix]), np.array([0, 0]), current_solution[j])
                    element_stiffness_matrix[i][j] += result[0]
                    element_load_vector[i] += result[1]
                elif time_integration_method == 2:
                    result = ode.AB2(time_step, np.array(mass_matrix), np.array([damping_matrix, damping_matrix]), np.array([stiffness_matrix, stiffness_matrix]), [current_solution[j], previous_solution[j]])
                    element_stiffness_matrix[i][j] += result[0]
                    element_load_vector[i] += result[1]
                elif time_integration_method == 3:
                    result = ode.AM3(time_step, np.array(mass_matrix), np.array([damping_matrix, damping_matrix, damping_matrix]), np.array([stiffness_matrix, stiffness_matrix, stiffness_matrix]), [current_solution[j], previous_solution[j]])
                    element_stiffness_matrix[i][j] += result[0]
                    element_load_vector[i] += result[1]
                elif time_integration_method == 4:
                    result = ode.BDF2(time_step, mass_matrix, damping_matrix, 0, [current_solution[j], previous_solution[j]])
                    element_stiffness_matrix[i][j] += result[0]
                    element_load_vector[i] += result[1]
    return element_stiffness_matrix, element_load_vector

def assemble(element_stiffness_matrix, element_load_vector, global_system_matrix, global_rhs_vector, element_node_indices):
    for i in range(4):
        global_rhs_vector[int(element_node_indices[i] - 1)] += element_load_vector[i]
        for j in range(4):
            global_system_matrix[int(element_node_indices[i] - 1)][int(element_node_indices[j] - 1)] += element_stiffness_matrix[i][j]
    return global_system_matrix, global_rhs_vector

def assign_dirichlet_boundary_conditions(global_system_matrix, global_rhs_vector, dirichlet_boundary_conditions):
    for i in range(dirichlet_boundary_conditions.shape[0]):
        global_rhs_vector[int(dirichlet_boundary_conditions[i][0] - 1)] = dirichlet_boundary_conditions[i][1]
        for k in range(global_system_matrix.shape[1]):
            if k == int(dirichlet_boundary_conditions[i][0] - 1):
                global_system_matrix[int(dirichlet_boundary_conditions[i][0] - 1)][k] = 1
            else:
                global_system_matrix[int(dirichlet_boundary_conditions[i][0] - 1)][k] = 0
    return global_system_matrix, global_rhs_vector

def solve(node_coordinates, element_indices, dirichlet_boundary_conditions):
    current_solution = np.full((18, 1), 300)
    previous_solution = np.zeros((18, 1))
    time_integration_method = 1
    time_step = 500
    total_time_steps = int(5000 / time_step)
    theta = 0.5
    is_first_step = 0
    for t in range(total_time_steps):
        if t == 0:
            is_first_step = 1
        else:
            is_first_step = 0
        element_evaluation = []
        assembled_system = []
        for element_node_indices in element_indices:
            element_evaluation = evaluate_transient(
                node_coordinates[element_node_indices - 1],
                gx2dref(2),
                gw2dref(2),
                current_solution[element_node_indices - 1],
                previous_solution[element_node_indices - 1],
                time_integration_method,
                time_step,
                theta,
                is_first_step,
            )
            if assembled_system == []:
                assembled_system = assemble(
                    element_evaluation[0],
                    element_evaluation[1],
                    np.zeros((node_coordinates.shape[0], node_coordinates.shape[0])),
                    np.zeros((node_coordinates.shape[0], 1)),
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
        assigned_system = assign_dirichlet_boundary_conditions(
            assembled_system[0], assembled_system[1], dirichlet_boundary_conditions
        )
        temperature_solution = np.linalg.solve(assigned_system[0], assigned_system[1])
        temperature_solution = np.array(temperature_solution).flatten()
        previous_solution = current_solution
        current_solution = temperature_solution
        if temperature_solution[14] > 450:
            print((t + 1) * (time_step))
            return temperature_solution
    return temperature_solution

start_time = Time.time()
sol = solve(nodes, elements, dirichlet_boundary_conditions)
end_time = Time.time()
elapsed_time = end_time - start_time
print(elapsed_time)
print(sol)
print(
    evaluate_transient(
        np.array([[0, 0], [1, 0], [1, 2], [0, 2]]),
        gx2dref(3),
        gw2dref(3),
        np.array([1, 2, 3, 4]),
        np.array([0, 0, 0, 0]),
        1,
        1000,
        0.66,
        1,
    )
)
quadplot(nodes, elements - 1, sol)
