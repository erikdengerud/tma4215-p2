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

	b1 = spl.BSplineBasis(3, knots=t1)
	b2 = spl.BSplineBasis(3, knots=t2)
	b3 = spl.BSplineBasis(4, knots=t3)
	b4 = spl.BSplineBasis(4, knots=t4)


	plt.figure()
	f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
	ax1.plot(x, b1(x))
	ax1.set_title('$\mathbf{\\tau}_1$ ')
	ax2.plot(x, b2(x))
	ax2.set_title('$\mathbf{\\tau}_2$ ')
	ax3.plot(x, b3(x))
	ax3.set_title('$\mathbf{\\tau}_3$ ')
	ax4.plot(x, b4(x))
	ax4.set_title('$\mathbf{\\tau}_4$ ')
	plt.savefig('Basis_Plot.pdf')
	plt.show()


Basis_Plot()


