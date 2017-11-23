################################################################################
import xml.etree.ElementTree as et
import sys
import numpy as np
from math import *
################################################################################
# from decimal import *
# getcontext().prec = 28
ERRORTOL = 1e-14
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
    # we only need to calculate the floor(n / 2) roots in (0, 1]. This
    # also reduces the required calculations for the weights
    
    guess = [ ( np.cos((2 * i - 1) * np.pi / (2 * n + 1)) + np.cos(2 * i * np.pi / (2 * n + 1)) ) / 2 for i in range(1, int(np.floor(n / 2)) + 1) ]
    xi = [ Olver(n, x_0) for x_0 in guess ]
    # the weights
    weights = [ 2 / ((1 - x**2) * Legendre_1(n, x, Legendre_0(n, x))[-1]**2) for x in xi ]
    # the roots to the left of x = 0
    xi_left = [-x for x in xi]
    xi = xi + xi_left
    # the weights
    weights += weights
    if n % 2 != 0:
        xi.append(0)
        weights.append( 2 / (Legendre_1(n, 0, Legendre_0(n, 0))[-1]**2) )
    
    return [xi, weights]

def Olver(n, x0):
    x = x0
    for i in range(ITERATIONCAP):
        L0_x = Legendre_0(n, x)
        L1_x = Legendre_1(n, x, L0_x)
        L2_x = Legendre_2(n, x, L1_x)
        
        L0 = L0_x[-1]
        L1 = L1_x[-1]
        L2 = L2_x[-1]
        
        assert(np.abs(L1) > EPS)
        delta_x = -L0 / L1 - L2 * L0**2 / (2 * L1**3)
        x += delta_x
        
        if (np.abs(delta_x) < ERRORTOL):
            return x
            
    assert(False)

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
    if x == -1:
        return [-i * (-1)**i * (i + 1) / 2 for i in range (n + 1)]
    if x == 1:
        return [i * (i + 1) / 2 for i in range (n + 1)]
    # result = [ L_0'(x), L_1'(x), ..., L_n'(x) ]
    result = [0, 1]
    for i in range (2, n + 1):
        result.append( ( (2 * i - 1) * (L0[i - 1] + x * result[i - 1]) - (i - 1) * result[i - 2] ) / i )
    return result

def Legendre_2(n, x, L1):
    if n == 0:
        return [0]
    # case x == +-1 (?)
    # result = [ L_0''(x), ..., L_n''(x) ]
    result = [0, 0]
    for i in range (2, n + 1):
        result.append( ( (2 * i - 1) * (2 * L1[i - 1] + x * result[i - 1]) - (i - 1) * result[i - 2] ) / i )
    return result
    

def Gauss_Legendre_Quadrature(n, G, f):
    assert(n > 0)
    result = 0
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

"""for i in range (1, 20 + 1):
    guess = [ ( np.cos((2 * j - 1) * np.pi / (2 * i + 1)) + np.cos(2 * j * np.pi / (2 * i + 1)) ) / 2 for j in range(1, int(np.floor(i / 2)) + 1) ]
    for x in guess:
        xi = Olver(i, x)
        assert(Legendre_0(i, xi)[-1] < ERRORTOL)"""


