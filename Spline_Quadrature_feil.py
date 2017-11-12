################################################################################
import numpy as np
import scipy as sp
import splipy as spl
import sys
################################################################################
import time

def Prepare_Data(T,p):
	'''Prepare_Data, for creating the initial condition and augmenting the knot
	vector. It must also return n and the constant integral vector in F.
	Done as in the paper: https://ac.els-cdn.com/S004578251630281X/1-s2.0-S004578251630281X-main.pdf?_tid=780e6096-c795-11e7-96c6-00000aacb35e&acdnat=1510483224_a01b14974539b6464a7dce410b080add
	'''
	T = np.array(T, dtype=float)
	#number of weights of knots
	n = len(T)-p-1 #egt 2n??
	#augmenting knot vector: adding to the knot vector if number of basis function is odd. Dette må først!
	if n%2==1:
		t = (T[p]+T[p+1])/2
		T = np.insert(T, p+1, t) #(τp+1 + τp+2)/2 between τp+1 and τp+2
		n+=1

	#constant integral vector in F.
	I = (T[p+1:]-T[:-p-1])/(p+1)
	#create initial condition
	G = np.zeros(n) #greville innebygd

	for i in range(n):
		G[i] = np.sum(T[i+1:i+1+p])/p #kan gjøres raskere?
	X = (G[1::2]+G[::2])/2 #nodes
	W = I[1::2]+I[::2] #weights
	basis = spl.BSplineBasis().__init__(order=p+1, knots=T)
	
	return basis, I, W, X, n

t1 = np.array([0, 0, 0, 1, 2,  3, 4, 4, 4])


def Assembly(basis,I,W,X,n):
	'''updating Fn and ∂Fn every time in the Newton iteration.
	Recall that ∂Fn must be permuted and made sparse
	Kanskje bare endre der X har endret seg? bedre tid?
	'''
	#make Fn. Trenger kanskje ikke hver gang? Kun oppdatere istedet?
	F = np.zeros(shape=(n,n))
	
	#make ∂Fn. Trenger kanskje ikke hver gang? Kun oppdatere istedet?
	J = np.zeros(shape=(n,n))#.sparse.csr_matrix(n, n, dtype=float), få det til å funke først


	#dspline = spl.SplineObject().__init__(basis, X)

	N = basis.evaluate(X).T
	dN = basis.evaluate(X, d=1).T#sparse=true
	for j in range(n):
		#update Fn
		
		F[j] = np.sum(W*N[j].T) - I[j]
	#for j in range(): fordi J er sparse burde den ha egen looå som går kortere
		#update ∂Fn
	
	temp = dN*sp.sparse.diags(W, 0)
	J = np.concatenate((N, temp), axis=1)#dN*diag(sparse(W))

	#∂Fn permuted and made sparse
	J = sp.sparse.csr_matrix(J)

	return F, J

basis, I, W, X, n = Prepare_Data(t1, 2)

F, J = Assembly(basis, I, W, X, n)
print(F)
print(J.todense())

def Spline_Quadrature(T, p):
	'''the main program with Newton iteration.'''
	# husk: The arrays obtained after using evaluate must be transposed.
	
	#Tolerance
	tol = 1e-11

	#Prepare_data
	basis, I, W, X, n = Prepare_Data(T, p)

	#Assembly
	F, J = Assembly(basis, I, W, X, n)
	delta = sp.sparse.linalg.spsolve(J, F)
	W += delta[0]
	X += delta[1]
	norm = np.linalg.norm(delta)
	itcount += 1
	print('Iteration ', itcount, '\t X = ',X ,'\t W = ', W, '\t F_norm = %0.2E' % norm)

	while abs(norm)>tol and itcount < 20:

		F, J = Assembly(basis,I,W,X,n)
		delta = scipy.sparse.linalg.spsolve(J, F)
		W += delta[0]
		X += delta[1]

		norm = np.linalg.norm(delta)

		itcount+=1
		print('Iteration ', itcount, '\t X = ',X ,'\t W = ', W, '\t F_norm = %0.2E' % norm)

	if abs(norm) > tol:
		itcount = -1 #Does not converge within 100 iterations
	return W, X, itcount
