#!/usr/bin/env python

import numpy as np
import functions as func

######## FORM THE CORE HAMILTONIAN
Nelec = 2					# Number of electrons
	
Enuc = np.loadtxt('enuc.dat', dtype=float)	# nuclear repulsion
Sraw = np.loadtxt('s.dat', dtype=float)		# overlap matrix
Traw = np.loadtxt('t.dat', dtype=float)		# kinetic energy matrix
Vraw = np.loadtxt('v.dat', dtype=float)		# potential energy matrix

# Number of basis fuctionss
dim = int((np.sqrt(8*len(Sraw)+1)-1)/2)

# Initialize integrals 
S = np.zeros((dim,dim))
T = np.zeros((dim,dim))
V = np.zeros((dim,dim))

# Store matrices triangularly
for i in range(Sraw.shape[0]):
	row, col    = int(Sraw[i,0]-1), int(Sraw[i,1]-1)
	S[row, col] = Sraw[i,2]
for i in range(Traw.shape[0]):
	row, col    = int(Traw[i,0]-1), int(Traw[i,1]-1)
	T[row, col] = Traw[i,2]
for i in range(Vraw.shape[0]):
	row, col    = int(Vraw[i,0]-1), int(Vraw[i,1]-1)
	V[row, col] = Vraw[i,2]

# For convenience let's fill the whole matrix
S = func.symmetrize(S)  
V = func.symmetrize(V)  
T = func.symmetrize(T)

Hcore = T + V  



######## TWO ELECTRON INTEGRALS
eri = np.loadtxt('eri.dat', dtype=float)	# Electron Repulsion Integrals

twoe = {func.eint(row[0], row[1], row[2], row[3]) : row[4] for row in eri}
