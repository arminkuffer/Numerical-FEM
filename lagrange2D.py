#FKT. I
import numpy as np
def linquadref(xi,eta): 
    N_i = np.array([0]*4)
    print(len(N_i))
    for i in range(len(N_i)):
        if i == 0: 
            N_i[i] = (1/4)*(1-xi)*(1-eta)
            print(N_i[i])
        elif i == 1:
            N_i[i] = (1/4)*(1+xi)*(1-eta)
        elif i == 2:
            N_i[i] = (1/4)*(1+xi)*(1+eta)      
        elif i == 3:
                N_i[i] = (1/4)*(1-xi)*(1+eta)   
    return N_i   

def LagrangePolynom(N_i,fx):
    poly = 0
    for i in range(len(N_i)):
        poly += N_i[i]*fx[i]
    return poly

#FKT. II
def linquadderivref(xi,eta):
    deriv = np.array([[0,0],[0,0],[0,0],[0,0]])
    for i in range(len(deriv)):
        if i == 0:
            deriv[i][0] = -(1/4)*(1-eta) #dxi
            deriv[i][1] = -(1/4)*(1-xi)  #deta
        elif i == 1:
            deriv[i][0] = (1/4)*(1-eta) #dxi
            deriv[i][1] = -(1/4)*(1+xi)#deta
        elif i == 2:
            deriv[i][0] = (1/4)*(1+eta) #dxi
            deriv[i][1] = (1/4)*(1+xi)#deta
        elif i == 3:
            deriv[i][0] = -(1/4)*(1+eta) #dxi
            deriv[i][1] = (1/4)*(1-xi)#deta
    return deriv
def LagrangePolyDeriv(deriv,fx):
    poly1 = 0
    poly2 = 0
    for i in range(len(deriv)):
        poly1 += (deriv[i][0])*fx[i]
        poly2 += (deriv[i][1])*fx[i]
    return poly1,poly2
xi1 = 0
eta1 = 0
xi2 = 0.577
eta2 = -0.577
f_node = [0,1,3,1]
#print(LagrangePolynom(linquadref(xi1,eta1),f_node))
#print(LagrangePolynom(linquadref(xi2,eta2),f_node))
#print(LagrangePolyDeriv(linquadderivref(xi1,eta1),f_node))
#print(LagrangePolyDeriv(linquadderivref(xi2,eta2),f_node))
print(linquadref(0,0))
print(linquadref(0.577,-0.577))
print(linquadderivref(0,0))
print(linquadderivref(0.577,-0.577))


