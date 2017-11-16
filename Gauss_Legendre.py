################################################################################
import xml.etree.ElementTree as et
import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
plt.style.use("ggplot")
################################################################################
ERRORTOL = 1e-6
EPS = 1e-25
ITERATIONCAP = 2000

def XML_Extraction(xmlfile):
    tree = et.parse(xmlfile)
    root = tree.getroot()
    f = lambda x : eval(root[0].text)
    analytical = float(root[1].text)
    return [ f, analytical ]



def Gauss_Legendre_Data(n):
    # find the n roots of L_{n}(x)
    guess = [ ( np.cos((2 * i - 1) * np.pi / (2 * n + 1)) + np.cos(2 * i * np.pi / (2 * n + 1)) ) / 2 for i in range(1, n + 1) ]
    
    # print(guess)
    xi = [ Olver(n, x_0)[0] for x_0 in guess ]
    # compute the weights
    weights = []
    for x in xi:    
        L1_x = Legendre_1(n, x, Legendre_0(n, x))[-1]
        weights.append( 2 / ((1 - x**2) * (L1_x**2)) )
    return [xi, weights]

def Olver(n, x0):
    x = x0        
    xvalues = [x]           
    errors = []
    for i in range(ITERATIONCAP):
        L0_x = Legendre_0(n, x)
        L1_x = Legendre_1(n, x, L0_x)
        L2_x = Legendre_2(n, x, L1_x)
        
        assert(np.abs(L1_x[-1]) > EPS)
        x = x - L0_x[-1] / L1_x[-1] - L2_x[-1] * L0_x[-1]**2 / (2 * L1_x[-1]**3)
        current_error = np.abs(xvalues[-1] - x)
        errors.append(current_error)
        xvalues.append(x)
        if (current_error < ERRORTOL):
            return x, errors
        
    assert(errors[-1] < ERRORTOL)

def Legendre_0(n, x):
    if n == 0:
        return [1]
    if x == -1:
        return [(-1)**i for i in range(n + 1)]
    if x == 1:
        return [1] * (n + 1)
    # result = [ L_0(x), L_1(x), ..., L_n(x) ]
    result = [1, x]
    for i in range (2, n + 1):
        result.append( ((2 * i - 1) * x * result[i - 1] - (i - 1) * result[i - 2]) / i )
    return result

def Legendre_1(n, x, L0):
    if n == 0:
        return [0]
    # result = [ L_0'(x), L_1'(x), ..., L_n'(x) ]
    result = [0, 1]
    for i in range (2, n + 1):
        result.append( ( (2 * i - 1) * (L0[i - 1] + x * result[i - 1]) - (i - 1) * result[i - 2] ) / i )
    return result

def Legendre_2(n, x, L1):
    if n == 0:
        return [0]
    # result = [ L_0''(x), ..., L_n''(x) ]
    result = [0, 0]
    for i in range (2, n + 1):
        result.append( ( (2 * i - 1) * (2 * L1[i - 1] + x * result[i - 1]) - (i - 1) * result[i - 2] ) / i )
    return result
    

def Gauss_Legendre_Quadrature(n, G, f):   
    result = 0
    for i in range(n):
        result += G[1][i] * f(G[0][i])
    return result

def Return_Quadrature(xmlfile, n):
    data = XML_Extraction(xmlfile)
    # n x 2 matrix with nodes and weights
    # with the n + 1 points and weights we can integrate all polynomials of degree <= 2 * n - 1 exactly
    G = Gauss_Legendre_Data(n) 
    numeric = Gauss_Legendre_Quadrature(n, G, data[0])
    analytic = data[1]
    
    return [numeric, analytic]

f = lambda x : 7 * x**12 + 2 * x**3 - 12 * x

G = Gauss_Legendre_Data(5)
print( Gauss_Legendre_Quadrature(5, G, f) )
    