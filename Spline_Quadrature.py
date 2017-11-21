################################################################################
import numpy as np
import scipy as sp
import splipy as spl
import sys
################################################################################

def Prepare_Data(T,p):
	'''Prepare_Data, for creating the initial condition and augmenting the knot
	vector. It must also return n and the constant integral vector in F.
	Done as in the paper: https://ac.els-cdn.com/S004578251630281X/1-s2.0-S004578251630281X-main.pdf?_tid=780e6096-c795-11e7-96c6-00000aacb35e&acdnat=1510483224_a01b14974539b6464a7dce410b080add
	'''

	#Creating the basis
	basis = spl.BSplineBasis(p+1, knots=T)

	#Number of weights of knots. This is also known as 2n.
	n = basis.__len__() -p-1

	#augmenting knot vector: adding to the knot vector if number of basis function is odd. 
	if n%2==1:
		t = (basis.knots[p]+basis.knots[p+1])/2
		basis.insert_knot(t) #(τp+1 + τp+2)/2 between τp+1 and τp+2
		n+=1

	#The constant integral vector in F.
	I = (basis.knots[p+1:]-basis.knots[:-p-1])/(p+1)

	#Creating the initial conditions.
	G = np.array(basis.greville()) 		#greville
	X = np.array(G[1::2]+G[::2])/2 		#nodes
	W = I[1::2]+I[::2]				 	#weights
	
	return basis, I, W, X, n


def Assembly(basis,I,W,X,n):
	'''updating Fn and ∂Fn every time in the Newton iteration.'''

	N = basis.evaluate(X).T
	dN = basis.evaluate(X, d=1).T#sparse=true

	#update Fn
	F = np.zeros(n)
	F[:] = np.dot(N, W) - I
	
	#Permutation of J:
	i_0  = list(range(0, n, 2))
	i_1  = list(range(1, n, 2))
	#updating
	dFde = dN*np.diag(W)
	dFdw = N
	#braiding
	J = np.empty((n,n))
	J[:, i_0] = dFdw
	J[:, i_1] = dFde
	
	#Make J sparse
	J = sp.sparse.csr_matrix(J)
	
	return F, J


def Spline_Quadrature(T, p, printout=False):
	'''the main program with Newton iteration.'''
	
	#Tolerance
	tol = 1e-11

	#Prepare_data
	basis, I, W, X, n = Prepare_Data(T, p)

	#Initial
	norm = 1
	itcount = 0

	while abs(norm)>tol and itcount < 20:

		F, J = Assembly(basis,I,W,X,n)

		delta = sp.sparse.linalg.spsolve(J, F)

		W -= delta[::2]
		X -= delta[1::2]

		norm = np.linalg.norm(delta)
		itcount+=1

		if printout:
			print('Iteration ', itcount, '\t X = ',X ,'\t W = ', W, '\t Norm = %0.2E' % norm)

		if min(X)<T[0] or max(X)>T[-1]:
			print('SINGULAR MATRIX!')
			break


	if abs(norm) > tol:
		itcount = -1 #Does not converge within 20 iterations
	return W, X, itcount


def test_2a():
	t1 = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
	t2 = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])
	t3 = np.array([0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4])
	t4 = np.array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4])

	Spline_Quadrature(t1, 2, printout=True)
	Spline_Quadrature(t2, 2, printout=True)
	Spline_Quadrature(t3, 3, printout=True)
	Spline_Quadrature(t4, 3, printout=True)

#test_2a()
