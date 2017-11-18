# -*- coding: utf-8 -*-
"""
Created on Sat Nov 18 13:01:30 2017

@author: William
"""

################################################################################
import numpy as np
import splipy as spl
from splipy.io import *
import matplotlib.pyplot as plt
import sys
from Spline_Quadrature import Spline_Quadrature, Assembly, Prepare_Data
################################################################################
import splipy.curve_factory as curves



with G2('Curve.g2') as my_file:
    curve = my_file.read() #returns a splinObject
    
tau = curve[0].knots(0,True)
p = curve[0].order(0) - 1

W, X, ehh = Spline_Quadrature(tau, p)

summen = 0
for i in X:
    derivert = curve[0].derivative(i)
    summen += np.sqrt(np.dot(derivert, derivert))
    
print(summen)