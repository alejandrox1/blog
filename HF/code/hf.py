#!/usr/bin/env python

#~	Hartree-Fock-Roothan Procedure
#~ This script will output the energy (screen) and will produce 'output.out'
#~ 
#~ Code based on : http://joshuagoings.com/2013/04/24/hartree-fock-self-consistent-field-procedure/


from __future__ import division
import sys  
import math  
import csv
import numpy as np  
import functions as func
import coreH as H

# STEP 1: Orthogonalize the basis (we are using the symmetric orthogonalization
#	  which uses the S^(-1/2) as the transformation matrix.
# 	  see Szabo and Ostlund page 143 for more details.
Lval, Lvec   = np.linalg.eig(H.S)
Sval_diag = np.diag(1./np.sqrt(Lval))
S_half    = np.dot(Lvec, np.dot(Sval_diag, Lvec.T))


# This is the main loop. See Szabo and Ostlun page 146.
D = np.zeros((H.dim, H.dim))            	# density matrix
G = np.zeros((H.dim, H.dim))			# used to make the fock matrix

deltaD          = 1.0
deltaE		= 1.0
EN		= 0.0
convergenceD    = 1e-8
convergenceE    = 1e-8
iteration       = 0
OUTPUT = open('output.out', 'w')
OUTPUT.write("Core Hamiltonias:\n {}\n".format(H.Hcore))
while (deltaD>convergenceD) or (deltaE>convergenceE):
	# Header
	iteration += 1
	OUTPUT.write('='*100 + '\n')
	OUTPUT.write('='+' '*48+str(iteration)+' '*48+'='+'\n')
	OUTPUT.write('='*100 + '\n')	

	F = func.makefock(H.Hcore, D, H.twoe)						# Fock matrix
	OUTPUT.write("\nFock Matrix:\n {}\n".format(H.Hcore))	
	
	Fprime = func.fprime(S_half, F) 						# Fock matrix in AO basis
	OUTPUT.write("\nFock Matrix in the AO Representation:\n {}\n".format(Fprime))	

	E, Cprime = np.linalg.eigh(Fprime)						# Find eigenvalues/vectors
	C = np.dot(S_half, Cprime)
	E, C = func.diagonalize(Fprime, S_half)						# change basis
	OUTPUT.write("\nEigenvectors of the Fock Matrix:\n {}\n".format(C))

	D, OLDD = func.makedensity(C, D, H.Nelec)					# new density matrix
	OUTPUT.write("\nDensity Matrix:\n {}\n".format(D))

	oldEN  = EN 
	EN     = func.currentenergy(D, H.Hcore, F)
	deltaD = func.deltaD(D, OLDD)
	deltaE = (EN - oldEN)
	OUTPUT.write("\nChange in the Density Matrix = {}\n".format(deltaD))
	OUTPUT.write("TOTAL E(SCF) = {}\n".format(EN + H.Enuc))
	OUTPUT.write("Change in E(SCF) = {}\n\n\n".format(deltaE))

OUTPUT.close()
print("TOTAL E(SCF) = {}\n".format(EN + H.Enuc))

