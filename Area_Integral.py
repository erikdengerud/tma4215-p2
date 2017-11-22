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
from Spline_Quadrature import Spline_Quadrature
################################################################################
def Area_Integral():
    # Read Area.g2 file into variable surface
    with G2("Area.g2") as file:
        surface = file.read()[0]
    
    # Fetch knot vectors in first and second direction
    tau_u, tau_v = surface.knots(0,True), surface.knots(1,True)
    
    # Fetch order of the spline for both directions
    p_u, p_v = surface.order()
    
    # Correction to the order
    p_u += -1
    p_v += -1
    
    # Calculate optimal nodes and weights from Spline_Quadrature
    W_u, X_u, dummy = Spline_Quadrature(tau_u,p_u)
    W_v, X_v, dummy = Spline_Quadrature(tau_v,p_v)
    
    # To check if the value is correct
    print('Area calculated by splipy : ',surface.area())
    
    # Define empty variable to add each area component to
    summen = 0
    
    # Main loop, calculates the Jacobian in each node of the surface, multiplies it by its respective weights and adds it to summen
    for i in np.arange(len(X_u)):
        for j in np.arange(len(X_v)):
            du = surface.derivative(X_u[i], X_v[j], d=(1,0))
            dv = surface.derivative(X_u[i], X_v[j], d=(0,1))
            summen += np.abs(W_u[i]*W_v[j]*np.cross(du,dv))
            
    print('Area calculated from spline quadrature and Jacoabian integral/sum : ',summen)

Area_Integral()