################################################################################
import numpy as np
import scipy as sp
import splipy as spl
import sys
################################################################################


def Prepare_Data(T,p):
	'''Prepare_Data, for creating the initial condition and augmenting the knot
	vector. It must also return n and the constant integral vector in F.
	'''
	#create initial condition

	#augmenting knot vector

	return n, constant_integral_vector_in_F

def Assembly(basis,I,W,X,n):
	'''updating Fn and ∂Fn every time in the Newton iteration.
	Recall that ∂Fn must be permuted and made sparse
	'''
	#update Fn

	#update ∂Fn

	#∂Fn permuted and made sparse

def Spline_Quadrature():
	'''the main program with Newton iteration.'''
	# husk: The arrays obtained after using evaluate must be transposed.
	
	#copied from last project
	#Tolerance
	tol = 1e-11

	#Defining the system as F=(f1, f2, f3)
	F = lambda x,y: np.array([x**4+x-y**3-1, x**2+x*y+y-2])

	#Defining the jacobi matrix of the system
	J = lambda x,y: np.matrix([[4*x**3+1, -3*y**2], [2*x+y, x+1]])


	def Multvariate_Newton(x, y):
		print('Starting at: ', np.array([x,y]))
		vec = np.array([x,y])
		F_value = F(vec[0], vec[1])
		F_norm = np.linalg.norm(F_value, ord=np.inf)
		itcount = 0
		while abs(F_norm)>tol and itcount<100:
			delta = np.linalg.solve(J(vec[0], vec[1]), -F_value)
			vec = vec + delta
			F_value = F(vec[0], vec[1])
			F_norm = np.linalg.norm(F_value, ord=np.inf)

			#update Fn and ∂Fn

			itcount+=1
			print('Iteration ', itcount, '\t xk1=',vec, '\t F_norm=%0.2E' % F_norm)

		if abs(F_norm) > tol:
			itcount = -1 #Does not converge within 100 iterations
		return vec, itcount
