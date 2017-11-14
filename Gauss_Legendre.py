################################################################################
import xml.etree.ElementTree as et
import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
plt.style.use("ggplot")
################################################################################
ERRORTOL = 1e-15
EPS = 1e-25
ITERATIONCAP = 1000

def XML_Extraction(xmlfile):
    tree = et.parse(xmlfile)
    root = tree.getroot()
    f = lambda x : eval(root[0].text)
    analytical = float(root[1].text)
    return [ f, analytical ]



def Gauss_Legendre_Data(n):
    # find the n roots of L_{n}(x)
    guess = [ np.cos((4 * i - 1) * np.pi / (4 * n + 2)) for i in range(1, n + 1) ]
    xi = [ Olver(n, x_0)[0] for x_0 in guess ]
    # compute the weights
    weights = []
    for x in xi:
        L1_x = Legendre_1(n, x, Legendre_0(n, x))[-1]
        weights.append( 2 / ((1 - x**2) * (L1_x**2)) )
    # weights = [ 2 / (( 1 - x**2 ) * (Legendre_1(n, x, Legendre_0(n, x))[-1])**2) for x in xi ]
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
        result.append( ((2 * n + 1) * x * result[i - 1] - n * result[i - 2] ) / (n + 1) )
    return result

def Legendre_1(n, x, L0):
    if n == 0:
        return [0]
    # result = [ L_0'(x), L_1'(x), ..., L_n'(x) ]
    result = [0, 1]
    for i in range (2, n + 1):
        result.append( ( (2 * n + 1) * (L0[i - 1] + x * result[i - 1]) - n * result[i - 2] ) / (n + 1) )
    return result

def Legendre_2(n, x, L1):
    if n == 0:
        return [0]
    # result = [ L_0''(x), ..., L_n''(x) ]
    result = [0, 0]
    for i in range (2, n + 1):
        result.append( ( (2 * n + 1) * (2 * L1[i - 1] + x * result[i - 1]) - n * result[i - 2] ) / (n + 1) )
    return result
    


def Gauss_Legendre_Quadrature(n, G, f):
    data = Gauss_Legendre_Data(n + 1)
    # with the n + 1 points and weights we can integrate all polynomials of degree <= 2 * n + 1 exactly
    result = 0
    for i in range(len(data[0])):
        result += data[1][i] * f(data[0][i])
    return result

def Return_Quadrature(xmlfile, n):
    G = 1 # ?
    data = XML_Extraction(xmlfile)
    numeric = Gauss_Legendre_Quadrature(n, G, data[0])
    analytic = data[1]
    
    return [numeric, analytic]

"""
n = 10
xvalues = [-1 + 2 * i / 1000 for i in range(1001)]
yvalues = [Legendre_0(n, x)[-1] for x in xvalues]
plt.plot(xvalues, yvalues)

guess = [ np.cos((4 * i - 1) * np.pi / (4 * n + 2)) for i in range(1, n + 1) ]
xi = [ Olver(n, x_0)[0] for x_0 in guess ]
for x in xi:
    plt.plot(x, 0, "ro")
"""

f = lambda x : 2 * x**3 + 3 * x**2 
print(Gauss_Legendre_Quadrature(2, 1, f)) # burde gi 2. Noe er feil med weights (nullpunktene er rett)
    
    