#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Inputs: 
# path
# landspace
# problem_id

from platform import mac_ver
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
    mpList = 0

    mVerify = False


    def __init__(self,pdir, plname, problem_id, verify):

        self.mVerify = verify

        fin = open(pdir + plname, 'r')

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

        for i in range(0,self.mLandscapesInFile):
            # Reads @L x
            line = fin.readline()
            ch = line.split()
            self.mL = ch[1]
            for j in range(0,self.mOb):
                # Reads @O x
                line1 = fin.readline()
                ch1 = line1.split()
                self.mOb = int(ch1[1])
                for k in range(0,self.mP):
                    index = k
                    line2 = fin.readline()
                    ch2 = line2.split()
                    for l in ch2:
                        if l == 1:
                            self.mpEpistasisAddress[i][j][k][index+1] = int(l)
                line1 = fin.readline() # Gets rid of the blank line beetween @O

    def GetRandomSeeds(self, fin):

        self.mpSeeds = np.array(self.mLandscapesInFile)

        line = fin.readline()
        ch = line.split()
        if ch[0] != "@SEEDS":
            print("It should be reading @SEEDS")
        for i in range(0,self.mLandscapesInFile):
            line = fin.readline()
            np.append(self.mpSeeds,line) # Each seed is in a line

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
            print ("Seeds[" + lid + "]" + self.mpSeeds[lid])
        
        for x in range(0,self.mOb):
            for i in range(0, self.mP):
                for j in range(0,self.mK):
                    self.mpList[x][i][j] = int(random.randint(0,1))
                    if (self.mVerify == True):
                        print ("List[" + x + "][" + i + "][" + j + "]=" + self.mpList[x][i][j])
        
        return

    def FitnessFunction (self, pgenotype, num_objectives, pfitness):
        pbit_string = pgenotype

        if self.mInitialized == False:
            print("ble")
            sys.exit("Problem must be initialized first")
        for i in range(0,self.mOb):
            pfitness[i] = 0
            for j in range(0,self.mP):
                self.mpFitContribution[i][j] = 0
        
        self.EvalString(pbit_string)

        for i in range(0,self.mOb):
            for j in range(0,self.mP):
                pfitness[i] += self.mpFitContribution[i][j]
                if (self.mVerify == True):
                    print("F[" + i + "][" + j + "] = " + self.mpFitContribution[i][j])
            if (self.mVerify == True):
                print("Ob: " + i + " " + pfitness[i])
        for i in range(0,self.mOb):
            pfitness[i] = pfitness[i] / self.mP
            if (self.mVerify == True):
                print("Ob: " + i + " " + pfitness[i])
        return pfitness

    def EvalString(self, pbits):

        for y in range(0,self.mOb): # Is this necessary?
            for i in range(0,self.mP):
                rseed = 0
                for j in range(0,self.mK+1):
                    bit = pbits[self.mpEpistasisAddress[self.mLandscapeId][y][i][j]]
                    if (bit==0):
                        z = j
                    elif (bit == 1):
                        z = j + self.mK + 1
                    rseed = rseed ^ self.mpList[y][i][z]
                    rseed += 1
            if (self.mVerify == True):
                print("rseed: " + rseed)
            self.mpFitContribution[y][i] += random.random() # A number from 0 to 1



        return

