import numpy as np
import matplotlib.pyplot as plt
#Blatt 2:
def lineintersection(p1,p2):
    y = 2
    x = [None]*500
    for i in range(len(p2[0])):
        x[i] = np.abs(1-(y-1)/((p2[1][i]-p1[1])/(p2[0][i]-p1[0])))
    print(x)
    fig, ax = plt.subplots(figsize=(8,6))
    ax.plot(np.log(delta),x)
    plt.show()
delta = np.logspace(-20,5, 500)
x1 = [0,1]
x2 = [delta,1+delta]
print(x2)
lineintersection(x1,x2)