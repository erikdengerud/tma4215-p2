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
	basis = spl.BSplineBasis()
	basis.__init__(order=p+1, knots=T)
	#number of weights of knots
	n = basis.__len__() -p-1#egt 2n??
	
	#augmenting knot vector: adding to the knot vector if number of basis function is odd. Dette må først!
	if n%2==1:
		t = (basis.knots[p]+basis.knots[p+1])/2
		basis.insert_knot(t) #(τp+1 + τp+2)/2 between τp+1 and τp+2
		n+=1
	#constant integral vector in F.
	I = (basis.knots[p+1:]-basis.knots[:-p-1])/(p+1)
	#create initial condition
	G = np.array(basis.greville()) #greville innebygd
	X = np.array(G[1::2]+G[::2])/2 #nodes
	W = I[1::2]+I[::2] #weights
	
	return basis, I, W, X, n

t1 = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
t2 = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])
t3 = np.array([0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4])
t4 = np.array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4])
knot_vectors = np.array([t1, t2, t3, t4])

def Assembly(basis,I,W,X,n):
	'''updating Fn and ∂Fn every time in the Newton iteration.
	Recall that ∂Fn must be permuted and made sparse
	Kanskje bare endre der X har endret seg? bedre tid?
	'''
	#make Fn. Trenger kanskje ikke hver gang? Kun oppdatere istedet?
	F = np.zeros(n)
	
	#make ∂Fn. Trenger kanskje ikke hver gang? Kun oppdatere istedet?
	J = np.zeros(shape=(n,n))#.sparse.csr_matrix(n, n, dtype=float), få det til å funke først


	#dspline = spl.SplineObject().__init__(basis, X)

	N = basis.evaluate(X).T
	dN = basis.evaluate(X, d=1).T#sparse=true

	#update Fn
	#F = np.squeeze(N.dot(W)-I) # dette skal funke, men gjør det ikke
	for j in range(n):
		F[j] = np.sum(W*N[j].T) - I[j] # F=N*W-I???
	
	#update ∂Fn
	temp = dN*sp.sparse.diags(W, 0)
	J = np.concatenate((N, temp), axis=1)#dN*diag(sparse(W))
	J2 = np.array(N[:,0])
	J2 = np.concatenate((J2, temp[:,0]), axis=1)
	for i in range(1,int(n/2)):
		J2 = np.concatenate((J2, N[:,i]), axis=1)
		J2 = np.concatenate((J2, temp[:,i]), axis=1)
	#print(J2)
	#print(np.asarray(J2))
	#∂Fn permuted and made sparse
	J = sp.sparse.csr_matrix(J)
	#print(J.todense())

	return F, J#, J2

basis, I, W, X, n = Prepare_Data(t1, 2)

#F, J, J2= Assembly(basis, I, W, X, n)
#print(F)
#print(J.todense())
#print(np.linalg.solve(J.todense(), F))
#print(np.linalg.solve(J2, F))


def Spline_Quadrature(T, p):
	'''the main program with Newton iteration.'''
	# husk: The arrays obtained after using evaluate must be transposed.
	
	#Tolerance
	tol = 1e-11

	#Prepare_data
	basis, I, W, X, n = Prepare_Data(T, p)
	print('Starting at: \n W: ', W, '\n X: ', X)
	#Assembly
	F, J = Assembly(basis, I, W, X, n)
	#First iteration
	delta = sp.sparse.linalg.spsolve(J, F)
	W -= delta[:int(n/2)]
	X -= delta[int(n/2):]
	norm = np.linalg.norm(delta)
	itcount = 1
	print('Iteration ', itcount, '\t X = ',X ,'\t W = ', W, '\t Norm = %0.2E' % norm)

	while abs(norm)>tol and itcount < 20:

		F, J = Assembly(basis,I,W,X,n)
		delta = sp.sparse.linalg.spsolve(J, F)
		W -= delta[:int(n/2)]
		X -= delta[int(n/2):]

		norm = np.linalg.norm(delta)

		itcount+=1
		print('Iteration ', itcount, '\t X = ',X ,'\t W = ', W, '\t Norm = %0.2E' % norm)
		if min(X)<T[0]:
			print('SINGULAR MATRIX!')
			break
		if max(X)>T[-1]:
			print('SINGULAR MATRIX!')
			break


	if abs(norm) > tol:
		itcount = -1 #Does not converge within 20 iterations
	return W, X, itcount

Spline_Quadrature(t1, 2)
Spline_Quadrature(t2, 2)
Spline_Quadrature(t3, 3)
Spline_Quadrature(t4, 3)
