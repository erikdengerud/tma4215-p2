# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 15:14:48 2017

@author: William
"""

################################################################################
import numpy as np
import splipy as spl
from splipy.io import G2
import splipy.surface_factory as spf
from Spline_Quadrature import Spline_Quadrature, Assembly, Prepare_Data
################################################################################

with G2("Area.g2") as file:
    surface = file.read()[0]
#2D
N = 25
u = np.linspace(surface.start('u'), surface.end('u'), N)
v = np.linspace(surface.start('v'), surface.end('v'), N)
x = surface(u, v)

tau_u, tau_v = surface.knots(0,True), surface.knots(1,True)

p_u, p_v = surface.order()

p_u += -1
p_v += -1

W_u, X_u, ehh = Spline_Quadrature(tau_u,p_u)
W_v, X_v, ehh = Spline_Quadrature(tau_v,p_v)

#du = surface.derivative(X_u, X_v, d=(1,0))
#dv = surface.derivative(X_u, X_v, d=(0,1))

print('Area calculated from splipy : ',surface.area())

summen = 0

for i in np.arange(len(X_u)):
    for j in np.arange(len(X_v)):
        du = surface.derivative(X_u[i], X_v[j], d=(1,0))
        dv = surface.derivative(X_u[i], X_v[j], d=(0,1))
        summen += np.abs(np.cross(du,dv))
        
print('Area calculated from spline quadrature and Jacoabian integral/sum : ',summen)
