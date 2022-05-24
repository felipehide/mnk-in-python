"""
======================================================================
  Python Version of the MNK Generator
  --------------------------------------------------------------------
  Create multi-objective NKP fitness landscapes (MNKP type)
  with random or near neighbors patterns of epistasis

  Parameters: 
    N  bit string lenght
	K  number of epistatic interactions per bit
	P  number of bits that contribute to fitness (see note 1)
	O  number of objectives
	x  x=R random epistasis, x=N near neighbor epistasis
	L  number of landscapes
    R  seed for the randomizing functions

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
import sys

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

def main():

    # User inputs
    print(sys.argv)
    if (len(sys.argv) > 1):
        N = int(sys.argv[1])
        K = int(sys.argv[2])
        P = int(sys.argv[3])
        O = int(sys.argv[4])
        pattern = sys.argv[5]
        L = int(sys.argv[6])
        seed = int(sys.argv[7])
    else: # No input given, takes from the keyboard
        N = int(input("N="))
        K = int(input("K="))
        P = int(input("P="))
        O = int(input("O="))
        pattern = input("NEAR(N), RAND(R)=")
        L = int(input("How many landscapes? "))
        seed = int(input("seed="))
    print(pattern)
    if (pattern!="N" and pattern!="R"):
        sys.exit("pattern does not match")

    long_pattern = 'NEAR' if pattern=='N' else 'RAND'
    filename = "L" + str(L) + "_N" + str(N) + "_K" + str(K) + "_P" + str(P) + "_O" + str(O) + "_" + long_pattern 
    print(filename + "( SEED" + str(seed) + ")")

    random.seed(seed)

    buffer = "@NKPO N " + str(N) + " K " + str(K) + " P " + str(P) + " O " + str(O) + " " + long_pattern + "\n"
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
            if (pattern == "N"):
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

    with open(filename, 'w') as file:
        file.write(buffer)

if __name__ == "__main__":
    main()