################################################################################
import numpy as np
from splipy.io import G2
import Spline_Quadrature as SQ
################################################################################

def Curve_Integral():
	#Reading in the file
	with G2('Curve.g2') as my_file:
	    curve = my_file.read()[0]
	
	#Gettin weights and nodes from Spline_Quadrature    
	tau = curve.knots(0,True)
	p = curve.order(0) - 1
	W, X, iterations = SQ.Spline_Quadrature(tau, p)

	#Sum over all ds to get the integral
	integral = 0
	for i, x in enumerate(X):
		derivative = curve.derivative(x)
		integral += W[i]*np.linalg.norm(derivative)

	print('Curve_Integral = ', integral)

	return integral

Curve_Integral()