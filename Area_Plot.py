################################################################################
import numpy as np
import splipy as spl
from splipy.io import G2
import splipy.surface_factory as spf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
################################################################################
plt.style.use('ggplot')

def Area_Plot():
    # Reads Area.g2 and stores it in the surface variable
    with G2("Area.g2") as file:
        surface = file.read()[0]

        # Number of points in the plot  
        N = 25
        
        # Creates the area from the surface variable
        u = np.linspace(surface.start('u'), surface.end('u'), N)
        v = np.linspace(surface.start('v'), surface.end('v'), N)
        x = surface(u,v)
        
        # Plots the area
        plt.plot(x[:,:,0],   x[:,:,1],   'k-')
        plt.plot(x[:,:,0].T, x[:,:,1].T, 'k-')
        plt.axis('equal')
        plt.xlabel('x')
        plt.ylabel('y')
        plt.title('Plot of the Area given in Area.g2')
        plt.savefig('Area_Plot.pdf')
        plt.show()


Area_Plot()