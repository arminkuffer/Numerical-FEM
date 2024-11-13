def FDF(x0,h):
    deriv = (f((x0))-f((x0-h))) /h
    return deriv     

def threepointFDF(x0,h):
    deriv = 1/(2*h)*(-3*f(x0)+4*f(x0+h)-f(x0+2*h))
    return deriv
def threepointCDF(x0,h):
    deriv = 1/(2*h)*(f(x0+h)-f(x0-h))
    return deriv
def fivepointCDF(x0,h):
    deriv = 1/(12*h)*(f(x0-2*h)-8*f(x0-h)+8*f(x0+h)-f(x0+2*h))
    return deriv
def f(x):
    return (x/(1+x))**5

x_0 = 0.6
h = 0.1
print(FDF(x_0,h))
print(threepointFDF(x_0,h))
print(threepointCDF(x_0,h))
print(fivepointCDF(x_0,h))