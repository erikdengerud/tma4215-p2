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
        N = 30
        u = np.linspace(surface.start('u'), surface.end('u'), N)
        v = np.linspace(surface.start('v'), surface.end('v'), N)
        X = surface(u, v)
        #Why did we have to find the syntax for this ourselves???????? Jesus Freaking Kennedy
        #Source:
        #https://github.com/sintefmath/Splipy/blob/d2b6dea5055e5c8435b729969aceb0bf170024b2/doc/Tutorial/Factory%20methods.ipynb
        plt.plot(X[:,:,0], X[:,:,1], '-', color="red")
        plt.plot(X[:,:,0], X[:,:,1].T, '-', color="blue")
        plt.axis('equal')
        plt.show()


Area_Plot()