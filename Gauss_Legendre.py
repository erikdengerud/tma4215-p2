################################################################################
import xml.etree.ElementTree as et
import sys
import numpy as np
from math import *
import matplotlib.pyplot as plt
################################################################################


def XML_Extraction(XMLFILE):
    return

def Gauss_Legendre_Data(n):
    return

def Olver(n,x0):
    x = x0        
    xvalues = [x]           
    errors = []
    for i in range(1000):
        L_x = Legendre_0(n, x)[-1]
        L1_x = Legendre_1(n, x)[-1]
        L2_x = Legendre_2(n, x)[-1]
        
        x = x - L_x / L1_x - L2_x * L_x**2 / (2 * L1_x**3)
        current_error = np.abs(xvalues[-1] - x)
        errors.append(current_error)
        xvalues.append(x)
        if (current_error < ERRORTOL):
            return x, errors
        
    assert(errors[-1] < ERRORTOL)

def Legendre_0(n,x):
    if n == 0:
        return [1]
    # result = [ L_0(x), L_1(x), ..., L_n(x) ]
    result = [1, x]
    for i in range (2, n + 1):
        result.append( ((2 * n + 1) * x * result[i - 1] - n * result[i - 2] ) / (n + 1) )
    return result

def Legendre_1(n,x):
    if n == 0:
        return [0]
    # result = [ L_0'(x), L_1'(x), ..., L_n'(x) ]
    result = [0, 1]
    # we need the values of Legendre_0(n, x)
    L0 = Legendre_0(n, x)
    for i in range (2, n + 1):
        result.append( ( (2 * n + 1) * (L0[i - 1] + x * result[i - 1]) - n * result[i - 2] ) / (n + 1) )
    return result

def Legendre_2(n,x,L0,L1):
    if n == 0:
        return [0]
    # result = [ L_0''(x), ..., L_n''(x) ]
    result = [0, 0]
    # we need the values of L1
    L1 = Legendre_1(n, x)
    for i in range (2, n + 1):
        result.append( ( (2 * n + 1) * (2 * L1[i - 1] + x * result[i - 1]) - n * result[i - 2] ) / (n + 1) )
    return result
    
def Gauss_Legendre_Quadrature(n,G,f):
    return

def Return_Quadrature(XMLFILE,n):
    return
    
xvalues = [-1 + 2 * i / 100 for i in range(101)]
yvalues = [ Legendre_0(10, x)[-1] for x in xvalues ]
plt.plot(xvalues, yvalues)

    
    
    
    
    
    
    
    
    
    