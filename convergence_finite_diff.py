import numpy as np
import matplotlib.pyplot as plt
import finitedifference as fdf

def deriv(x):
    return 5*(x/(1+x))**4 * 1/((1+x)**2)

def plotfivepointCDF(x0,h):
    plt.figure(4)
    plt.loglog(h,np.abs(np.subtract(deriv(x0),fdf.fivepointCDF(x0,h)))) #Plot abs err fivepointCDF
    plt.ylabel(r"$|f'({}) - f'_{err}$({})|")
    plt.xlabel('h')
    plt.title('abs err of different approximation methods for derivatives')
    plt.plot(h,np.power(h,4))
    plt.show()

def plotthreepointCDF(x0,h):
    plt.figure(3)
    plt.loglog(h,np.abs(np.subtract(deriv(x0),fdf.threepointCDF(x0,h)))) #Plot abs err fivepointCDF
    plt.ylabel(r"$|f'({}) - f'_{err}$({})|")
    plt.xlabel('h')
    plt.title('abs err of different approximation methods for derivatives')
    plt.plot(h,np.power(h,2))
    plt.show()

def plotthreepointFDF(x0,h):
    plt.figure(2)
    plt.loglog(h,np.abs(np.subtract(deriv(x0),fdf.threepointFDF(x0,h)))) #Plot abs err threepointFDF
    plt.ylabel(r"$|f'(x) - f'_{err}$(x)|")
    plt.xlabel('h')
    plt.title('abs err of different approximation methods for derivatives')
    plt.plot(h,np.power(h,2))
    plt.show()

def plotFDF(x0,h):
    plt.figure(1)
    plt.loglog(h,np.abs(np.subtract(deriv(x0),fdf.FDF(x0,h)))) #Plot abs err FDF
    plt.ylabel(r"$|f'(x) - f'_{err}$(x)|")
    plt.xlabel('h')
    plt.title('abs err of different approximation methods for derivatives')
    plt.plot(h,np.power(h,1))
    plt.show()

x1 = 0.6
x2 = 2
h = np.logspace(-5,0,500)

plotFDF(x1,h)
plotthreepointFDF(x1,h)
plotthreepointCDF(x1,h)
plotfivepointCDF(x1,h)
plotFDF(x2,h)
plotthreepointFDF(x2,h)
plotthreepointCDF(x2,h)
plotfivepointCDF(x2,h)


