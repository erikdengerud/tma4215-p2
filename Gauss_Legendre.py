################################################################################
import xml.etree.ElementTree as et
import sys
import numpy as np
from math import *
################################################################################
ERRORTOL = 1e-15
EPS = 1e-20
ITERATIONCAP = 100

def XML_Extraction(xmlfile):
    tree = et.parse(xmlfile)
    root = tree.getroot()
    f = lambda x : eval(root[0].text)
    analytical = float((eval(root[1].text)))
    return [f, analytical]

def Gauss_Legendre_Data(n):
    assert(n > 0)
    # find the n roots of L_{n}(x). Since the roots are symmetric about x = 0,
    # we only need to calculate the floor(n / 2) roots in (0, 1]. This also 
    # applies to the corresponding weights
    
    guess = [ ( np.cos((2 * i - 1) * np.pi / (2 * n + 1)) + np.cos(2 * i * np.pi / (2 * n + 1)) ) / 2 for i in range(1, int(np.floor(n / 2)) + 1) ]
    
    xi = []
    weights = []
    for x0 in guess:
        # we keep the already computed derivatives
        x, L1_x = Olver(n, x0)
        xi.append(x)
        weights.append(2 / ((1 - x**2) * L1_x**2))
    # by the symmertri
    xi += [ -x for x in xi ]
    weights += weights
    if n % 2 != 0:
        xi.append(0)
        weights.append( 2 / (Legendre_1(n, 0, Legendre_0(n, 0))**2) )
    
    return [xi, weights]

def Olver(n, x0):
    x = np.float(x0)
    for i in range(ITERATIONCAP):
        L0 = Legendre_0(n, x)
        L0_x = L0[-1]
        L1_x = Legendre_1(n, x, L0)
        L2_x = Legendre_2(n, x, L0_x, L1_x)
        
        # assert(np.abs(L1_x) > EPS)
        delta_x = -L0_x / L1_x - L2_x * L0_x**2 / (2 * L1_x**3)
        x += delta_x
        
        if (np.abs(delta_x) < ERRORTOL):
            return x, L1_x
            
    assert(False)
    # if this assert triggers please decrease
    # ERRORTOL or increase ITERATIONCAP

def Legendre_0(n, x):
    if n == 0:
        return np.array([1])
    if x == -1:
        return np.array([(-1)**i for i in range(n + 1)])
    if x == 1:
        return np.ones(n + 1)
    # result = [ L_0(x), L_1(x), ..., L_n(x) ]
    
    result = np.zeros(n + 1)
    result[0] = 1
    result[1] = x
    
    for i in range (2, n + 1):
        result[i] = ((2 * i - 1) * x * result[i - 1] - (i - 1) * result[i - 2]) / i
        
    return result

def Legendre_1(n, x, L0):
    if n == 0:
        return 0
    if n == 1:
        return 1
    if x == -1:
        return -n * (-1)**n * (n + 1) / 2
    if x == 1:
        return n * (n + 1) / 2
    
    if n % 2 == 0:
        indices = np.arange(1, n, step = 2)
    else:
        indices = np.arange(0, n, step = 2)
        
    return np.sum( (2 * indices + 1) * L0[indices] )

def Legendre_2(n, x, L0_x, L1_x):
    if n <= 1:
        return 0
    if x == -1:
        return (-1)**n * (n - 1) * n * (n + 1) * (n + 2) / 8
    if x == 1:
        return (n - 1) * n * (n + 1) * (n + 2) / 8
    
    return (2 * x * L1_x - n * (n + 1) * L0_x) / (1 - x * x)
    

def Gauss_Legendre_Quadrature(n, G, f):
    assert(n > 0)
    result = 0.
    for i in range(n):
        result += G[1][i] * f(G[0][i]) 
    return result

def Return_Quadrature(xmlfile, n):
    data = XML_Extraction(xmlfile)
    # n x 2 matrix with nodes and weights
    # with the n points and weights we can integrate all polynomials of degree <= 2 * n - 1 exactly
    G = Gauss_Legendre_Data(n) 
    numeric = Gauss_Legendre_Quadrature(n, G, data[0])
    analytic = data[1]
    return [numeric, analytic]
