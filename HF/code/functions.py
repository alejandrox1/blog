#!/usr/bin/env python

import numpy as np

############################   FUNCTIONS ( coreH.py )   ##############################
def symmetrize(a):
	"""Symmetrize a matrix given a triangular one."""
	return a + a.T - np.diag(a.diagonal())

def eint(a,b,c,d):
	"""Return a unique compound index given four indices.
	We used this function to uniquily determine the given two 
	electron integrals.
	"""
	if a>b:
		ab = a*(a+1)/2 + b
	else:
		ab = b*(b+1)/2 + a
	if c>d:
		cd = c*(c+1)/2 + d
	else:
		cd = d*(d+1)/2 + c
	if ab>cd:
		abcd = ab*(ab+1)/2 + cd
	else:
		abcd = cd*(cd+1)/2 + ab
	return abcd

############################   FUNCTIONS ( hf.py )   ##############################
def tei(twoe, a,b,c,d):
	"""Return the value of the two electron integral
	Example: (12 | 34) = tei(1,2,3,4)
	We use the chemist notation for the scf procedure.

	twoe is a fictionary containing the necessary two electron
	integrals.
	"""
	return twoe.get(eint(a,b,c,d),0.0)

def makefock(Hcore, D, twoe):
	"""Make Fock matrix."""
	F = np.zeros((Hcore.shape[0], Hcore.shape[1]))
	for i in range(Hcore.shape[0]):
		for j in range(Hcore.shape[1]):
			F[i,j] = Hcore[i,j]
			
			for k in range(D.shape[0]):
				for l in range(D.shape[1]):
					F[i,j] += D[k,l]*(tei(twoe,i+1,j+1,k+1,l+1)
						-0.5*tei(twoe,i+1,k+1,j+1,l+1))
	return F

def makedensity(C,D,Nelec):
        """Make Density matrix and store the old one to test for
        convergence.
        """
        OLDD = np.zeros((D.shape[0],D.shape[1]))
        for mu in range(D.shape[0]):
                for nu in range(D.shape[1]):
                        OLDD[mu,nu] = D[mu,nu]
                        D[mu,nu] = 0.0
                        for m in range(Nelec//2):
                                D[mu,nu] += 2*C[mu,m]*C[nu,m]
        return D, OLDD


def fprime(X,F):
        """Put Fock matrix in Orthonormal AO basis."""
        return np.dot(X.T, np.dot(F,X))

def diagonalize(M, S_half):
        """"Diagonalize a matrix and return the Eigenvalues and 
        non-orthogonal Eigenvectors in separate 2D arrays.
        """
        e, Cprime = np.linalg.eigh(M)
        C = np.dot(S_half, Cprime)  
        return e,C

def deltaD(D, OLDD):
	"""Calculate the change in the density matrix."""
	delta = 0.0
	for i in range(D.shape[0]):
		for j in range(D.shape[1]):
			delta += ((D[i,j] - OLDD[i,j])**2.0)
		
	return (np.sqrt(delta/4.0))

def currentenergy(D, Hcore, F):
	"""Calculate energy at each iteration."""
	EN = 0.0
	for mu in range(D.shape[0]):
		for nu in range(D.shape[1]):
			EN += 0.5*D[mu,nu]*(Hcore[mu,nu] + F[mu,nu])
	return EN


