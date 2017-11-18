################################################################################
import sys
import Spline_Quadrature as SPQ
from splipy.io import G2
import numpy as np
################################################################################


def Curve_Integral(path):
    


    with G2(path) as file:
        #Make a splipy curve object from file
        curve = file.read()[0]
        
        #Extract knots and order
        Tau = curve.knots(0, True)
        p = curve.order(0) - 1
        
        #Calculate optimal quadrature points and weights
        W, X, itcount = SPQ.Spline_Quadrature(Tau, p)
        
        #Avoid using a for loop at all costs during the evaluation of the 
        # trajectory quadrature, as (unless using numba) this will be far slower
        # than using the c-implemented procedures of numpy. 
        
        #curve.derivative(X) evaluates the derivative in the points X. As this
        # returns a (n*2) array in the format of n*(dx/dt, dy/dt) we use 
        # np.linalg.norm with axis=1 to take the 2-norm s.t we get 
        # ds = sqrt((dx/dt)**2 + (dy/dt)**2). We then dot the weights with the
        # ds-vector. 
        path = np.dot(W, np.linalg.norm(curve.derivative(X), axis = 1))

        return path


print(Curve_Integral('Curve.g2'))