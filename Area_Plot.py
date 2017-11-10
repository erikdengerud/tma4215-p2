################################################################################
import numpy as np
import splipy as spl
from splipy.io import *
import splipy.surface_factory as spf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
################################################################################


def Area_Plot():

	with G2('Area.g2') as my_surf:
		my_obj = my_surf.read()

	
	x = []
	y = []
	z = []
	for vec in my_obj[0]:
		x.append(vec[0])
		y.append(vec[0])
		z.append(0)

	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.plot_surface(x,y,z)
	plt.show()
	plt.figure()
	plt.plot(y,x)
	plt.show()

Area_Plot()