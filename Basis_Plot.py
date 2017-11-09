################################################################################
import numpy as np
import splipy as spl
import matplotlib.pyplot as plt
import sys
################################################################################


def Basis_Plot():
	'''Plots basis splines for different knot vectors'''
	x = np.linspace(0,4, 100)
	t1 = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
	t2 = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])
	t3 = np.array([0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4])
	t4 = np.array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4])
	knot_vectors = np.array([t1, t2, t3, t4])

	for t in knot_vectors:
		spline = spl.BSplineBasis()
		spline.__init__(order=3, knots=t)#hva er order her? 3 eller 4 gir best plot. er det forskjell på de ulike vektorene?
										 #er det avhengig av hvor mange like det er på slutten?
		matrix = spline.evaluate(x)
		matrix=matrix.T
		
		plt.figure()
		for vec in matrix:
			plt.plot(x, vec.T)
		plt.show()

Basis_Plot()


