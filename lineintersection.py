import numpy as np
import matplotlib.pyplot as plt
#Blatt 2:
def lineintersection(p2):
    y = 2
    x = [None]*500
    ex = [1]*500
    for i in range(len(p2[0])):
        x[i] = (p2[0][i])/((p2[1][i])-1)
    fig, ax = plt.subplots(figsize=(8,6))
    ax.loglog(delta,np.abs(np.subtract(ex,x)))
    plt.show()
delta = np.logspace(-20,5, 500)
x2 = [delta,1+delta]
lineintersection(x2)