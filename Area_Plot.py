################################################################################
import numpy as np
import splipy as spl
from splipy.io import G2
import splipy.surface_factory as spf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
################################################################################


def Area_Plot():
    with G2("Area.g2") as file:
        surface = file.read()[0]

        #2D        
        N = 25
        u = np.linspace(surface.start('u'), surface.end('u'), N)
        v = np.linspace(surface.start('v'), surface.end('v'), N)
        x = surface(u,v)
        plt.plot(x[:,:,0],   x[:,:,1],   'k-')
        plt.plot(x[:,:,0].T, x[:,:,1].T, 'k-')
        plt.axis('equal')
        plt.show()


Area_Plot()