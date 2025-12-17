import numpy as np
import matplotlib.pyplot as plt

def forward_euler(time_steps):
    phi_values = np.zeros(len(time_steps))  # Initialize solution array
    phi_values[0] = 0
    n = 0
    while n <= (len(time_steps) - 2):
        phi_values[n + 1] = phi_values[n] + time_steps[1] * forward_function(phi_values[n], time_steps[n])
        n += 1
    return phi_values

def backward_euler(time_steps):
    phi_values = np.zeros(len(time_steps))  # Initialize solution array
    phi_values[0] = 0
    n = 0
    while n <= (len(time_steps) - 2):
        phi_values[n + 1] = (phi_values[n] + time_steps[1] * backward_function(time_steps[n + 1])) / (6 * time_steps[1] + 1)
        n += 1
    return phi_values

def trapezoidal_rule(time_steps):
    phi_values = np.zeros(len(time_steps))  # Initialize solution array
    phi_values[0] = 0
    n = 0
    while n <= (len(time_steps) - 2):
        phi_values[n + 1] = (
            phi_values[n]
            + time_steps[1] * trapezoidal_function(time_steps[n + 1])
            + 1 / 2 * time_steps[1] * forward_function(phi_values[n], time_steps[n])
        ) / (4 * time_steps[1] + 1)
        n += 1
    return phi_values

def forward_function(phi, time):
    return (time**2) * np.exp(-5 * time) - 6 * phi

def backward_function(time):
    return (time**2) * np.exp(-5 * time)

def trapezoidal_function(time):
    return (1 / 2) * (time**2) * np.exp(-5 * time)

time_step_size = 0.2
t = np.arange(0, 2 + time_step_size, time_step_size)
x_ex = np.linspace(0, 2, 100)
y_ex = np.exp(-5 * x_ex) * ((x_ex**2) - 2 * x_ex + 2) - 2 * np.exp(-6 * x_ex)
"""plt.plot(x_ex,y_ex)
plt.plot(t,forward_euler(t))
plt.plot(t,backward_euler(t), color="red")
plt.plot(t,trapezoidal_rule(t),color="green")
plt.show()"""

#FKT IX:
def OST(theta, time_step, mass_matrix, damping_matrix, stiffness_matrix, solution_vector):
    lhs_rhs = np.zeros((2,), dtype=object)
    lhs_rhs[0] = mass_matrix - theta * time_step * damping_matrix[0]
    lhs_rhs[1] = (
        np.dot((mass_matrix + (1 - theta) * time_step * damping_matrix[1]), solution_vector)
        + time_step * (theta * stiffness_matrix[0] + (1 - theta) * stiffness_matrix[1])
    )
    return lhs_rhs


#FKT X:
def AB2(time_step, mass_matrix, damping_matrix, stiffness_matrix, solution_vectors):
    lhs_rhs = np.zeros((2,), dtype=object)
    lhs_rhs[0] = mass_matrix
    lhs_rhs[1] = (
        np.matmul(solution_vectors[0], mass_matrix)
        + time_step
        * (
            (3 / 2) * (np.matmul(solution_vectors[0], damping_matrix[0]) + stiffness_matrix[0])
            - (1 / 2) * (np.dot(solution_vectors[1], damping_matrix[1]) + stiffness_matrix[1])
        )
    )
    return lhs_rhs


#FKT XI:
def AM3(time_step, mass_matrix, damping_matrix, stiffness_matrix, solution_vectors):
    lhs_rhs = np.zeros((2,), dtype=object)
    lhs_rhs[0] = mass_matrix - time_step * (1 / 12) * 5 * damping_matrix[0]
    lhs_rhs[1] = (
        np.dot(mass_matrix, solution_vectors[0])
        + time_step
        * (
            (1 / 12)
            * (
                5 * stiffness_matrix[0]
                + 8 * np.dot(damping_matrix[1], solution_vectors[0])
                + 8 * stiffness_matrix[1]
                - np.dot(damping_matrix[2], solution_vectors[1])
                - stiffness_matrix[2]
            )
        )
    )
    return lhs_rhs


#FKT XII:
def BDF2(time_step, mass_matrix, damping_matrix, stiffness_matrix, solution_vectors):
    lhs_rhs = np.zeros((2,), dtype=object)
    lhs_rhs[0] = 1.5 * mass_matrix - time_step * damping_matrix
    lhs_rhs[1] = (
        2 * np.dot(mass_matrix, solution_vectors[0])
        - (1 / 2) * np.dot(mass_matrix, solution_vectors[1])
        + time_step * stiffness_matrix
    )
    return lhs_rhs



#Anwendung auf DGL
def DGLOST(phi0,timestep,theta):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        phi[n] = (OST(theta,timestep,[1.0],[-6.0,-6,0],[backward_function((n)*timestep),backward_function((n-1)*timestep)],[phi[n-1]])[1])/(OST(theta,timestep,[1.0],[-6.0,-6,0],[backward_function((n)*timestep),backward_function((n-1)*timestep)],[phi[n-1]])[0])
    return phi

def DGLAB2(phi0,timestep):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        if n == 1: 
            phi[n] = (OST(0.5,timestep,[1.0],[-6.0,-6.0],[backward_function((n-1)*timestep),backward_function((n-2)*timestep)],[phi[n-1]])[1])/(OST(0.5,timestep,[1.0],[-6.0,-6,0],[backward_function((n)*timestep),backward_function((n-1)*timestep)],[phi[n]])[0])
        else: 
            phi[n] = (AB2(timestep,[1.0],[-6.0,-6.0],[backward_function((n-1)*timestep),backward_function((n-2)*timestep)],[phi[n-1],phi[n-2]])[1])/(AB2(timestep,[1.0],[-6.0,-6.0],[backward_function((n-1)*timestep),backward_function((n-2)*timestep)],[phi[n-1],phi[n-2]])[0])
    return phi

def DGLAM3(phi0,timestep):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        if n == 1: 
            phi[n] = (OST(0.5,timestep,[1.0],[-6.0,-6.0],[backward_function((n)*timestep),backward_function((n-1)*timestep)],[phi[n-1]])[1])/(OST(0.5,timestep,[1.0],[-6.0,0],[backward_function((n)*timestep),backward_function((n-1)*timestep)],[phi[n-1]])[0])
        else: 
            phi[n] = (AM3(timestep,[1.0],[-6.0,-6.0,-6.0],[backward_function((n)*timestep),backward_function((n-1)*timestep),backward_function((n-2)*timestep)],[phi[n-1],phi[n-2]])[1])/(AM3(timestep,[1.0],[-6.0,-6.0,-6.0],[backward_function((n)*timestep),backward_function((n-1)*timestep),backward_function((n-2)*timestep)],[phi[n-1],phi[n-2]])[0])
    return phi

def DGLBDF2(phi0,timestep):
    phi = np.zeros(int(2/timestep+1))
    phi[0] = phi0
    for n in range(1,len(phi)):
        if n == 1: 
            phi[n] = (OST(0.5,timestep,[1.0],[-6.0,-6.0],[backward_function((n)*timestep),backward_function((n-1)*timestep)],[phi[n-1]])[1])/(OST(0.5,timestep,[1.0],[-6.0,0],[backward_function((n)*timestep),backward_function((n-1)*timestep)],[phi[n-1]])[0])
        else: 
            phi[n] = (BDF2(timestep,[1.0],[-6.0],[backward_function((n)*timestep)],[phi[n-1],phi[n-2]])[1])/(BDF2(timestep,[1.0],[-6.0],[backward_function((n)*timestep)],[phi[n-1],phi[n-2]])[0])
    return phi

"""plt.plot(t,forward_euler(t),color="orange")
plt.plot(t,DGLOST(0,dt,0.5),color="red")
plt.plot(t,DGLAB2(0,dt),color="blue")
plt.plot(t,DGLAM3(0,dt),color="green")
plt.plot(t,DGLBDF2(0,dt),color="yellow")
plt.plot(x_ex,y_ex,color="black")
plt.show()"""
#print(OST(0.5, 0.2, np.array([1.1]), np.array([1.4, 1.5]), np.array([1.7, 1.8]), np.array([2.0])))
# print(AB2(0.2, np.array([1.1]), np.array([1.5, 1.6]), np.array([1.8, 1.9]), np.array([2.0, 2.1])))
"""print(AM3(0.2, np.array([1.1]), np.array([1.4, 1.5, 1.6]), np.array([1.7, 1.8, 1.9]), np.array([2.0, 2.1])))
print(BDF2(0.2, np.array([1.1]), np.array([1.4]), np.array([1.7]), np.array([2.0, 2.1])))"""
M = np.array([[1.1, 1.2],[1.2, 1.1]],dtype=float)
B = np.array([[1.4, 1.5],[1.5 ,1.4]],dtype=float)
C = np.array([1.7,1.8],dtype=float)
sol = np.array([2.0,3.0],dtype=float)
print(OST(0.5, 0.2, M, np.array([B, B]), np.array([C, C]), sol))
print(AB2(0.2, M, [B, B], [C, C], [sol,sol]))
print(AM3(0.2, M, [B, B, B], [C, C, C], [sol,sol]))
print(BDF2(0.2, M, B, C, [sol,sol]))
