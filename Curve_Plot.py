################################################################################
import numpy as np
import splipy as spl
from splipy.io import *
import matplotlib.pyplot as plt

import sys
################################################################################
plt.style.use('ggplot')

def Curve_Plot():
	with G2('Curve.g2') as my_file:
		curve = my_file.read()[0] #returns a splinObject
 
		x = curve[:,0]
		y = curve[:,1]

		plt.figure()
		plt.gcf().subplots_adjust(left=0.18)
		plt.title('$[x(t), y(t)]^T$')
		plt.plot(x,y)
		plt.xlabel('$x(t)$')
		plt.ylabel('$y(t)$')
		plt.savefig('Curve_Plot.pdf')
		plt.show()
		
	

Curve_Plot()