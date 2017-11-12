################################################################################
import numpy as np
import splipy as spl
import matplotlib.pyplot as plt
import sys
################################################################################


def Basis_Plot():
	'''Plots basis splines for different knot vectors'''
	x = np.linspace(0, 4, 100)
	t1 = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
	t2 = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])
	t3 = np.array([0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4])
	t4 = np.array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4])
	knot_vectors = np.array([t1, t2, t3, t4])

	count = 0
	for t in knot_vectors:
		spline = spl.BSplineBasis()
		if count < 2:
			spline.__init__(order=3, knots=t)
			count+=1
		else:
			spline.__init__(order=4, knots=t)								
		matrix = spline.evaluate(x)
		matrix=matrix.T
		#print(spline)
		#print(matrix)
		plt.figure()
		for vec in matrix:
			plt.plot(x, vec.T)
		plt.show()

Basis_Plot()


