import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
r = 0.02+0.06 #r = 0.08 kleinster Radius für den gilt T15 < 450 (Maximalwert für y = h)
b = 0.3
h = 0.3
n = 4

def meshGen(b,h,r,n):
    n_rhq = int(np.round(np.sqrt(1/2)*n))
    n_r = 2*n_rhq+1
    center = np.array([b,h])
    s_x = b/(n_rhq)
    s_y = h/(n_rhq)
    P_base = np.zeros((n_r,2))
    for i in range(n_r):
        if i < n_rhq:
            P_base[i] = np.array([0,b-s_y*i])
        elif i == n_rhq:
            P_base[i] = np.array([0,0])
        else:  
            P_base[i] = np.array([s_x*(i-n_rhq),0])

    V_base = np.zeros((n_r,2))
    for i in range(n_r):
        V_base[i] = center - P_base[i] 
        V_base[i] = V_base[i]*1/n * (1-r/np.linalg.norm(V_base[i]))

    P_nodes = np.zeros((n_r*(n+1),2))
    for i in range(n+1):
        for j in range(n_r):
            P_nodes[n_r*i+j] = P_base[j] + V_base[j]*i

    E = np.zeros((n*(n_r-1),4))
    for i in range(n):
        for j in range(n_r-1):
            E[(n_r-1)*i+j] = np.array([j, j+1, j+n_r+1, j+n_r]) + i*n_r+1
    E = np.int_(E)

    P_edges = np.zeros((5,),dtype=object)
    P_edges[0] = np.arange(1, n_r * (n + 1), n_r)  #north
    P_edges[1] = np.arange(n_r, n_r * (n + 1) + 1, n_r)  #east
    P_edges[2] = np.arange(n_rhq+1, n_r + 1, 1)  #south
    P_edges[3] = np.arange(1, n_rhq + 2, 1)  #west
    P_edges[4] = np.arange(1, n_r + 1, 1) + n_r * n #circle  
    return P_nodes, E , P_edges

#print mesh
"""nodes,elements,edges = meshGen(b,h,r,n)
polygons = [Polygon(nodes[element - 1], closed=True) for element in elements]

# Create the plot
fig, ax = plt.subplots(figsize=(6, 6))

# Plot elements (polygons)
patch_collection = PatchCollection(polygons, edgecolor='blue', facecolor='cyan', alpha=0.5)
ax.add_collection(patch_collection)

# Plot edges (lines)
segments = []
for edge in edges:
    for i in range(len(edge) - 1):
        segments.append([nodes[edge[i] - 1], nodes[edge[i + 1] - 1]])
line_collection = LineCollection(segments, colors='black', linewidths=0.75)
ax.add_collection(line_collection)

# Plot nodes
ax.plot(nodes[:, 0], nodes[:, 1], 'o', color='red')

ax.set_aspect('equal')
plt.show()"""

    