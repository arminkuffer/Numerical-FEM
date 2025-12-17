import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from scipy.ndimage import gaussian_filter
"""def plot2D(x):
    y = np.sin(x)
    plt.plot(x,y,label="f(x) = sin(x)")
    plt.title("Plot von sin(x)")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.grid(True)
    plt.legend()
    plt.show()
plot2D(np.linspace(-np.pi,np.pi,500))"""
node_coordinates = np.array([[-1, -1], [0, -1], [1, -1], [-1, 0], [0, 0], [1, 0], [-1, 1], [0, 1], [1, 1]])
element_indices = np.array([[0, 1, 4, 3], [1, 2, 5, 4], [3, 4, 7, 6], [4, 5, 8, 7]])
temperature_solution = node_coordinates[:, 0]**2 + node_coordinates[:, 1]**2

# Function to plot quadrilateral elements
def plot_quadrilateral_mesh(node_coordinates, element_indices, solution):
    triangles = []
    for element in element_indices:
        triangles.append([element[0],element[1],element[2]])
        triangles.append([element[2],element[3],element[0]])
    triangles = np.array(triangles)
    x = node_coordinates[:,0]
    y = node_coordinates[:,1]
    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')
    triang = Triangulation(x,y,triangles)
    trisurf = ax.plot_trisurf(triang,solution,cmap='hot')
    color_bar = fig.colorbar(trisurf, ax=ax, shrink=0.5, aspect=5)
    color_bar.set_label('Temperature')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('T(x,y)')
    ax.set_title('3D Plot')
    plt.show()
"""quadplot(node_coordinates, element_indices, temperature_solution)"""