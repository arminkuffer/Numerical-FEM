# Function I
import numpy as np
def linear_quadrilateral_shape_functions(xi, eta): 
    shape_functions = np.zeros(4)  
    for i in range(len(shape_functions)):
        if i == 0: 
            shape_functions[i] = (1/4)*(1-xi)*(1-eta)
        elif i == 1:
            shape_functions[i] = (1/4)*(1+xi)*(1-eta)
        elif i == 2:
            shape_functions[i] = (1/4)*(1+xi)*(1+eta)
        elif i == 3:
            shape_functions[i] = (1/4)*(1-xi)*(1+eta)
    return shape_functions   

def lagrange_polynomial(shape_functions, function_values):
    polynomial = 0
    for i in range(len(shape_functions)):
        polynomial += shape_functions[i] * function_values[i]
    return polynomial

# Function II
def linear_quadrilateral_shape_function_derivatives(xi, eta):
    derivatives = np.zeros((4, 2)) 
    for i in range(len(derivatives)):
        if i == 0:
            derivatives[i][0] = -(1/4)*(1-eta)  # dN1/dxi
            derivatives[i][1] = -(1/4)*(1-xi)   # dN1/deta
        elif i == 1:
            derivatives[i][0] = (1/4)*(1-eta)  # dN2/dxi
            derivatives[i][1] = -(1/4)*(1+xi)  # dN2/deta
        elif i == 2:
            derivatives[i][0] = (1/4)*(1+eta)  # dN3/dxi
            derivatives[i][1] = (1/4)*(1+xi)   # dN3/deta
        elif i == 3:
            derivatives[i][0] = -(1/4)*(1+eta) # dN4/dxi
            derivatives[i][1] = (1/4)*(1-xi)   # dN4/deta
    return derivatives

def lagrange_polynomial_derivative(derivatives, function_values):
    polynomial_x = 0
    polynomial_y = 0
    for i in range(len(derivatives)):
        polynomial_x += (derivatives[i][0]) * function_values[i]
        polynomial_y += (derivatives[i][1]) * function_values[i]
    return polynomial_x, polynomial_y

xi1 = 0
eta1 = 0
xi2 = 0.577
eta2 = -0.577
function_values_at_nodes = [0, 1, 3, 1]
#print(lagrange_polynomial(linear_quadrilateral_shape_functions(xi1, eta1), function_values_at_nodes))
#print(lagrange_polynomial(linear_quadrilateral_shape_functions(xi2, eta2), function_values_at_nodes))
#print(lagrange_polynomial_derivative(linear_quadrilateral_shape_function_derivatives(xi1, eta1), function_values_at_nodes))
#print(lagrange_polynomial_derivative(linear_quadrilateral_shape_function_derivatives(xi2, eta2), function_values_at_nodes))



