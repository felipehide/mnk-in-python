"""
======================================================================
  Any Objective Evolutionary Optimization
  Hernan Aguirre
  Shinshu University
  Nagano, Japan
  ahernan@shinshu-u.ac.jp
  Copyright 1997-2011
  ----------------------------------------------------------------------
  Create multi-objective NKP fitness landscapes (MNKP type)
  with random or near neighbors patterns of epistasis

  Parameters: 
    N  bit string lenght
	K  number of epistatic interactions per bit
	P  number of bits that contribute to fitness (see note 1)
	O  number of objectives
	x  x=R random epistasis, x=N near neighbor epistasis
	L  number of landscapes

  Output:     
    A file containing the definition of the specified landscapes and a 
	random seed per landscape.  
	For each of the N bits there will be a string of 0s and 1s, 
	marking with 1 the bit itself and its K interacting bits (epistasis). 
	The random seeds will be used by the MNK-landscape class to 
	generate 2(K+1) random numbers per bit, which are in turn used by the 
	FitnessFunction of the same class to calculate the contribution to 
	fitness of the bits applying a Reproducible Random Number Generation 
	Method [1].


  Notes: 
    1. In Kauffman's NK-landscapes all bits contribute to fitness, so P=N
	2. If P=N, the MNKP landscapes are known as MNK-landscapes and extend 
	   Kauffman's NK-landscapes to multiple objectives

  Reference:
  [1] M. Shinkai, H. Aguirre, K. Tanaka, "Mutation Strategy Improves GA's 
	Performance on Epistatic Problems", Proc. 2002 Congress on Evolutionary 
	Computation (CEC'02), IEEE Service Center, vol.1, pp.968-973.

  Example:    
    MONKP.exe 4 2 4 2 R 2
    creates landscapes with N=4 bits, K=2 epsistatic interactions
	per bit, P=4 bits that contribute to fitness, O=2 objectives
	per landscape, x=R random patterns of epistatis for all bits, 
	and L=2 landscapes in the file
 
  Output file:  
     L2_N4_K2_P4_O2_RAND (file name)	 
	 ---------- start ----------
		@NKPO N 4 K 2 P 4 O 2 RAND
		@LANDSCAPES 2
		@L 1
		@O 1
		1	1	1	0	
		0	1	1	1	
		1	1	1	0	
		0	1	1	1	

		@O 2
		1	1	0	1	
		1	1	0	1	
		1	0	1	1	
		0	1	1	1	

		@L 2
		@O 1
		1	0	1	1	
		1	1	0	1	
		1	1	1	0	
		1	1	0	1	

		@O 2
		1	1	1	0	
		0	1	1	1	
		0	1	1	1	
		1	0	1	1	

		@SEEDS 2
		2959394231
		2332811499
	 ---------- end ----------
  
  ======================================================================
"""

import numpy as np
import random


def near_epistasis(N,K,P,table):

    right = left = int(K/2)
    if (K%2 == 1):
        right = right+1
    
    for i in range(0,P):
        r = i
        if ((P < N) or i >= N):
            r = randbetween (0,N-1)
        table[i][r] = 1
        for j in range(0,right):
            table[i][(N+r+j+1)% N] = 1
        for j in range(0,left):
            table[i][(N+r+j-1)% N] = 1
    return table

def random_epistasis(N,K,P,table):
    for i in range(0,P):
        r = i
        if (P<N or i >= N):
            r = randbetween(0,N-1)
        table[i][r] = 1
        for j in range(0,K):
            r = randbetween(1,N-1-j)
            z = 0
            x = 0
            for x in range(0,N):
                if (table[i][x] == 0):
                    z += 1
                    if (z == r):
                        break
            table[i][x] = 1
    return table

def randbetween(a,b):
    if (a==b):
        return a
    val = random.randint(0,4294967295)%(b-a+1)+a
    return val

def random_seeds (L,O):
    r_seeds = ""
    for i in range(0,L):
        r_seeds += str(random.randint(0,4294967295)%(4294967295-1+1)+1) + "\n"
    return r_seeds

def main(N,K,P,O,L,pattern):
    print("Insert here confirmation of variables")


    buffer = "@NKPO N " + str(N) + " K " + str(K) + " P " + str(P) + " O " + str(O) + " " + pattern + "\n"
    buffer += "@LANDSCAPES " + str(L) + "\n"

    for count in range(0,L):
        buffer += "@L " + str(count+1) + "\n"
        for obj in range(0,O):
            # Clear NKP landspace memory (not needed in python ???)
            # for i in range(0,P):
            #     for j in range(0,N):
            #         table[i][j] = 0
            buffer += "@O " + str(obj+1) + "\n"
            table = np.zeros([P,N])
            if (pattern == "NEAR"):
                table=near_epistasis(N,K,P,table)
            else:
                table=random_epistasis(N,K,P,table)
            for i in range(0,P):
                for j in range(0,N):
                    buffer += str(int(table[i][j])) + "\t"
                buffer += "\n"
            buffer += "\n"
    
    buffer += "@SEEDS " + str(L) + "\n"

    buffer += random_seeds(L,O)

    with open('output.txt', 'w') as file:
        file.write(buffer)

N=20
K=1
P=20
O=2
L=30
pattern = "RAND" # OR RAND

main(N,K,P,O,L,pattern)
