#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Inputs: 
# path
# landspace
# problem_id

import numpy as np
import random
import sys

class NKLandscape:
    mN = 0
    mK = 0
    mP = 0
    mOb = 0
    EpistasisType = ""
    
    mInitialized = False

    mInitType = ""
    mNkType = ""
    mpEpistasisAddress = 0
    
    mpFitContribution = 0

    mLandscapeId = 0
    mLandscapesInFile = 0

    mpSeeds = []
    mpList = 0              # Stores the seeds to generate the fitness of each bit

    mVerify = False

    def __init__(self,plname, problem_id, verify):

        self.mVerify = verify

        fin = open(plname, 'r')

        self.ReadLandscapeHeader(fin)
        self.GetEpistasisAddress(fin)
        self.GetRandomSeeds(fin)
        fin.close()
        self.GetMemoryVarious()

        self.Initialize(problem_id)

    def ReadLandscapeHeader(self, fin = None):

        line = fin.readline() # Read N x K x P x O x NEAR/RAND
        ch = line.split()
        self.mN  = int(ch[2])
        self.mK  = int(ch[4])
        self.mP  = int(ch[6])
        self.mOb = int(ch[8])
        EpistasisType = ch[9]

        if EpistasisType == "NEAR":
            if (self.mK % 2 == 1 and self.mK != self.mN - 1):
                print("K should be even for NEAR pattern of epistasis")
        elif EpistasisType == "RAND":
            self.mEpistasisType = EpistasisType
        else:
            print("Epistasis type is not correct")
        self.mNkType = "NKP"
        self.mInitType = "PROBLEM_50"

    # Read Epistasis_Address	
    # @LANDSCAPES x
    def GetEpistasisAddress(self, fin):
        line = fin.readline()
        ch = line.split()
        self.mLandscapesInFile = int(ch[1])

        self.mpEpistasisAddress = np.zeros([self.mLandscapesInFile, self.mOb, self.mP, int(self.mK + 1)], dtype=int)

        for x in range(0,self.mLandscapesInFile):
            # Reads @L x
            line = fin.readline()
            ch = line.split()
            self.mL = ch[1]
            for y in range(0,self.mOb):
                # Reads @O x
                line1 = fin.readline()
                ch1 = line1.split()
                self.mOb = int(ch1[1])
                for i in range(0,self.mP):
                    index = 0
                    line2 = fin.readline()
                    ch2 = line2.split()
                    for j in range(0,self.mN):
                        epistasis = int(ch2[j])
                        if (epistasis == 1):
                            self.mpEpistasisAddress[x][y][i][index] = int(j)
                            index = index+1
                line1 = fin.readline() # Gets rid of the blank line beetween @O
                
    def GetRandomSeeds(self, fin):

        self.mpSeeds = np.zeros(self.mLandscapesInFile, dtype = int)

        line = fin.readline()
        ch = line.split()
        if ch[0] != "@SEEDS":
            print("It should be reading @SEEDS")
        for i in range(0,self.mLandscapesInFile):
            line = fin.readline()
            self.mpSeeds[i] = line # Each seed is in a line

    def GetMemoryVarious(self):

        self.mpFitContribution = np.zeros([self.mOb, self.mP])
        self.mpList = np.zeros([self.mOb, self.mP, int(2*(self.mK+1))], dtype=int)

    def Initialize(self, landscape_id):
        if (landscape_id < 1 and landscape_id > self.mLandscapesInFile):
            print("There is no landscape " + landscape_id)
            print("Total number of landscapes in file : " + self.mLandscapesInFile)
            print("Available landscapes are 1,..., " + self.mLandscapesInFile)
        self.mLandscapeId = landscape_id - 1
        
        self.InitList(self.mLandscapeId)

        self.mInitialized = True

    def InitList(self, lid):
        if (self.mVerify == True):
            print ("Seeds[" + str(lid) + "]" + str(self.mpSeeds[lid]))

        random.seed(self.mpSeeds[lid])
        
        for x in range(0,self.mOb):
            for i in range(0, self.mP):
                for j in range(0,int( 2* (self.mK + 1))):
                    self.mpList[x][i][j] = int(random.randint(0,4294967295)) # 0xFFFFFFFF from C++, 32 bits size
                    if (self.mVerify == True):
                        print ("List[" + str(x) + "][" + str(i) + "][" + str(j) + "]=" + str(self.mpList[x][i][j]))

    def FitnessFunction (self, pgenotype):
        pbit_string = pgenotype

        pfitness = np.zeros(self.mOb)

        if self.mInitialized == False:
            sys.exit("Problem must be initialized first")
        for i in range(0,self.mOb):
            pfitness[i] = 0
            for j in range(0,self.mP):
                self.mpFitContribution[i][j] = 0
        
        self.EvalString(pbit_string)

        # Assings the fitness to each objective
        for i in range(0,self.mOb):
            for j in range(0,self.mP):
                pfitness[i] += self.mpFitContribution[i][j]
                if (self.mVerify == True):
                    print("F[" + str(i) + "][" + str(j) + "] = " + str(self.mpFitContribution[i][j]))
            if (self.mVerify == True):
                print("Ob: " + str(i) + " " + str(pfitness[i]))
        for i in range(0,self.mOb):
            pfitness[i] = pfitness[i] / self.mP
            if (self.mVerify == True):
                print("Ob: " + str(i) + " " + str(pfitness[i]))
        return pfitness

    def EvalString(self, pbits):
        for y in range(0,self.mOb):
            for i in range(0,self.mP):
                rseed = 0
                for j in range(0,self.mK+1):
                    bit = pbits[self.mpEpistasisAddress[self.mLandscapeId][y][i][j]]
                    if (bit == 0):
                        z = j
                    elif (bit == 1):
                        z = j + self.mK + 1
                    rseed = rseed ^ self.mpList[y][i][z]
                if (self.mVerify == True):
                    print("rseed: " + str(rseed))
                random.seed(rseed)

                self.mpFitContribution[y][i] += random.random() # A number from 0 to 1

    def EvalStringDebug(self, pbits):
        for y in range(0,self.mOb):
            for i in range(0,self.mP):
                rseed = 0
                for j in range(0,self.mK+1):
                    bit = pbits[self.mpEpistasisAddress[self.mLandscapeId][y][i][j]]
                    print(bit)
                    if (bit == 0):
                        z = j
                    elif (bit == 1):
                        z = j + self.mK + 1
                    rseed = rseed ^ self.mpList[y][i][z]
                random.seed(rseed)

                self.mpFitContribution[y][i] += random.random() # A number from 0 to 1


    def print_for_debug(self):
        with open("debug.txt", "w") as f:
            for i in self.mpList:
                f.write(str(i))
                f.write("\n")