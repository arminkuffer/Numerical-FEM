import numpy as np
import time as t 

def solve_gaussian_elimination(matrix, rhs):
    num_equations = len(rhs)
    solution = np.zeros((num_equations,))
    for i in range(num_equations - 1):
        for j in range(i + 1, num_equations):
            multiplier = matrix[j][i] / matrix[i][i]
            matrix[j] = matrix[j] - multiplier * matrix[i]
            rhs[j] = rhs[j] - multiplier * rhs[i]
    for i in range(num_equations - 1, -1, -1):
        solution[i] = rhs[i] / matrix[i][i] 
        for j in range(i + 1, num_equations):
            solution[i] -= (matrix[i][j] * solution[j]) / matrix[i][i]

    return solution

def solve_gradient_descent(matrix, rhs, initial_guess, tolerance, max_iterations):
    residual = rhs - np.matmul(matrix, initial_guess)
    iteration = 0
    solution = initial_guess
    while np.linalg.norm(residual, ord=2) > tolerance and iteration < max_iterations:
        direction = np.matmul(matrix, residual)
        step_size = np.dot(residual.T, residual) / np.dot(residual.T, direction)
        solution += step_size * residual
        residual -= step_size * direction
        iteration += 1
    print(iteration)
    return solution

def solve_conjugate_gradient(matrix, rhs, initial_guess, tolerance, max_iterations):
    residual = rhs - np.matmul(matrix, initial_guess)
    direction = residual
    iteration = 0
    solution = initial_guess
    new_residual = np.zeros((len(rhs),))
    while np.linalg.norm(residual, ord=2) > tolerance and iteration < max_iterations:
        direction_matrix_product = np.matmul(matrix, direction)
        step_size = np.dot(residual.T, residual) / np.dot(direction.T, direction_matrix_product)
        solution += step_size * direction
        new_residual = residual - step_size * direction_matrix_product
        conjugate_coefficient = np.dot(new_residual.T, new_residual) / np.dot(residual.T, residual)
        direction = new_residual + conjugate_coefficient * direction
        residual = new_residual
        iteration += 1
    print(iteration)
    return solution

"""n = 300
phi = 5.01
rhs = np.ones((n,))
matrix = np.diag(np.full((n), phi)) - 2 * np.diag(np.ones((n - 1)), k=1) - 2 * np.diag(np.ones((n - 1)), k=-1)
time_start = t.time()
solution = solve_gaussian_elimination(matrix, rhs)
time_stop = t.time()
print(solution)
print(time_stop - time_start)
time_start = t.time()
solution = solve_gradient_descent(matrix, rhs, np.zeros((n,)), 0.0000001, 1000)
time_stop = t.time()
print(solution)
print(time_stop - time_start)"""

"""print(solve_conjugate_gradient(np.array([[10.0, 2.0, 10.0],[ 2.0, 40.0, 8.0],[ 10.0, 8.0, 60.0]]),np.array([1.0, 1.0, 2.0]),np.array([0.0, 0.0, 0.0]), 10**(-7), 1000))
print(solve_gradient_descent(np.array([[10.0, 2.0, 10.0],[ 2.0, 40.0, 8.0],[ 10.0, 8.0, 60.0]]),np.array([1.0, 1.0, 2.0]),np.array([0.0, 0.0, 0.0]), 10**(-7), 1000))"""
print(solve_gaussian_elimination(np.array([[10.0,2.0,2.0],[3.0,4.0,4.0],[1.0,8.0,4.0]]),np.array([1.0,1.0,2.0])))
