def linquadref(xi,eta):
    N_i = [0]*4
    for i in range(len(N_i)):
        match i:
            case 0: 
                N_i[i] = (1/4)*(1-xi)*(1-eta)
            case 1:
                N_i[i] = (1/4)*(1+xi)*(1-eta)
            case 2:
                N_i[i] = (1/4)*(1+xi)*(1+eta)      
            case 3:
                N_i[i] = (1/4)*(1-xi)*(1+eta)   
    return N_i   

def LagrangePolynom(N_i,fx):
    poly = 0
    for i in range(len(N_i)):
        poly += N_i[i]*fx[i]
    return poly


def linquadderivref(xi,eta):
    deriv = [[0,0],[0,0],[0,0],[0,0]]
    for i in range(len(deriv)):
        match i:
            case 0:
                deriv[i][0] = -(1/4)*(1-eta)
                deriv[i][1] = -(1/4)*(1-xi)
            case 1:
                deriv[i][0] = (1/4)*(1-eta)
                deriv[i][1] = -(1/4)*(1+xi)
            case 2:
                deriv[i][0] = (1/4)*(1+eta)
                deriv[i][1] = (1/4)*(1+xi)
            case 3:
                deriv[i][0] = -(1/4)*(1+eta)
                deriv[i][1] = (1/4)*(1-xi)
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
print(LagrangePolynom(linquadref(xi1,eta1),f_node))
print(LagrangePolynom(linquadref(xi2,eta2),f_node))
print(LagrangePolyDeriv(linquadderivref(xi1,eta1),f_node))
print(LagrangePolyDeriv(linquadderivref(xi2,eta2),f_node))


