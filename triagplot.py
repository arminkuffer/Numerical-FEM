import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.tri import Triangulation
#commented solutions not important for FEM programm
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
nodes = np.array([[-1,-1],[0,-1],[1,-1],[-1,0],[0,0],[1,0],[-1,1],[0,1],[1,1]])
elements = np.array([[0,1,4,3],[1,2,5,4],[3,4,7,6],[4,5,8,7]])
sol = nodes[:,0]**2+nodes[:,1]**2

#FKT. 0
def quadplot(nodes,elements,sol):
    triangles = []
    for element in elements:
        triangles.append([element[0],element[1],element[2]])
        triangles.append([element[2],element[3],element[0]])
    triangles = np.array(triangles)
    x = nodes[:,0]
    y = nodes[:,1]

    fig = plt.figure()
    ax = fig.add_subplot(111,projection='3d')

    triang = Triangulation(x,y,triangles)

    ax.plot_trisurf(triang,sol,cmap='viridis')

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.set_title('3D Plot')
    plt.show()
quadplot(nodes,elements,sol)