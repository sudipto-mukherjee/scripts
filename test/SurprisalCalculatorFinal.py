#!/usr/env/bin python
import sys, os
from scipy.io import mmread
import numpy as np


# Functions


# Classes

class SurprisalCalculator:

  def __init__(self, sparse1, sparse2, Verbose=False):
    self.sparse1 = sparse1
    self.sparse2 = sparse2
    self.Verbose = Verbose

  def prepareSparseMatrices(self, state):
    jStateList = []
    sparse1 = self.sparse1
    sparse2 = self.sparse2 

    index1 = sparse1.row == state
    index2 = sparse2.row == state

    col1 = sparse1.col[index1]
    col2 = sparse2.col[index2]
    counts1 = sparse1.data[index1]
    counts2 = sparse2.data[index2]

    for colIter in range(col1.shape[0]):
      jStateList.append(col1[colIter])

    for colIter2 in range(col2.shape[0]):
      while not jStateList.__contains__(col2[colIter2]):
        jStateList.append(col2[colIter2])

    jStateList.sort()

    for i in range(len(jStateList)):
      if not col1.__contains__(jStateList[i]):
        col1 = np.insert(col1, i, jStateList[i])
        counts1 = np.insert(counts1, i, 0)

      if not col2.__contains__(jStateList[i]):
        col2 = np.insert(col2, i, jStateList[i])
        counts2 = np.insert(counts2, i, 0)

    oldCounts = counts1.astype(np.float64)
    newCounts = counts2.astype(np.float64)
    #print "oldCounts: ", oldCounts
    #print "newCounts: ", newCounts
    return oldCounts, newCounts

  def calculateEntropy(self, counts):
    
    entropySub = 0.0
    tCounts = np.sum(counts)
 
    for i in range(len(counts)):
      if counts[i] != 0:
         entropySub += counts[i]*np.log(counts[i])

    if (tCounts == 0):
      entropy = 0.0

    else:
      entropy = (1.0/tCounts)*(tCounts*np.log(tCounts) - entropySub)
   
    return entropy

  def getCovarianceMatrices(self, sparse):
    """Returns the covariance matrix of a multinomial distribution for p_i ~ Counts."""
    covMatrix = np.zeros( (sparse.shape[0], sparse.shape[1]) )
    
    for state in range(sparse.shape[0]):
      jStateArr = sparse.col[sparse.row == state]
      counts = sparse.data[sparse.row == state].astype(np.float64)
      tCounts = np.sum(counts)
      selfTrans = sparse.data[(sparse.row == state)*(sparse.col == state)][0] 

      for k in range(len(jStateArr)):
        jState = jStateArr[k]    
        covMatrix[state][jState] += -selfTrans*counts[k]/tCounts
        if state == jState:
          covMatrix[state][jState] += selfTrans
    return covMatrix

  def EstimateSurprisalVarianceAnalytical(self, state):
    covMatrix1 = self.getCovarianceMatrices(self.sparse1)
    covMatrix2 = self.getCovarianceMatrices(self.sparse2)
    print 'Covariance Matrix 1', covMatrix1
    #variance = []
    sizeM = self.sparse1.shape[1]
    diagCovMatrix = np.zeros ( (2*sizeM, 2*sizeM), np.float)
    diagCovMatrix[0:sizeM,0:sizeM] = covMatrix1
    diagCovMatrix[sizeM:2*sizeM, sizeM:2*sizeM] = covMatrix2

    #for state in range(self.sparse1.shape[0]):
    sensitivities1 = np.zeros(sizeM)
    sensitivities2 = np.zeros(sizeM)
    jStateArr1 = self.sparse1.col[self.sparse1.row == state]
    jStateArr2 = self.sparse2.col[self.sparse2.row == state]
    print 'jStateArr1', jStateArr1
    print 'jStateArr2', jStateArr2
    allJStates = np.unique(np.append(jStateArr1, jStateArr2))
    counts1, counts2 = self.prepareSparseMatrices(state)
    tCounts1 = np.sum(counts1).astype(np.float64)
    tCounts2 = np.sum(counts2).astype(np.float64)
    print 'total counts', tCounts1+tCounts2
    counts1 = counts1.astype(np.float64)
    counts2 = counts2.astype(np.float64)
    i = 0
    print 'allJstate', allJStates
    print 'counts1', counts1
    print 'counts2', counts2
    for jState in allJStates:
#print "Sensitivities1: tCounts1+tCounts2 jTotCounts1+jTotCounts2 counts1[i]", tCounts1+tCounts2, jTotCounts1+jTotCounts2, counts1[i]
      if counts1[i] != 0:
        sensitivities1[jState] += (1.0/(tCounts1+tCounts2))*(np.log(tCounts1/(tCounts1 + tCounts2))-np.log(counts1[i]/(counts1[i] + counts2[i])))

      if counts2[i] != 0:
        sensitivities2[jState] += (1.0/(tCounts1+tCounts2))*(np.log(tCounts2/(tCounts1 + tCounts2))-np.log(counts2[i]/(counts1[i] + counts2[i])))

      i += 1
      
    sensitivities = np.append(sensitivities1, sensitivities2)
    print 'sensitivities', sensitivities
    cq = np.dot(diagCovMatrix, sensitivities)
    variance = np.dot(np.transpose(sensitivities), cq)

    return variance
    
  def EstimateSurprisalVarianceBootstrap(self, oldCounts, newCounts, nBootstraps=1000):
    """Called for every state!."""
    totalOldCounts = np.sum(oldCounts)
    pOld = np.divide(oldCounts, totalOldCounts)

    totalNewCounts = np.sum(newCounts)
    pNew = np.divide(newCounts, totalNewCounts)

    surprisals = []

    sampledOldCounts = np.random.multinomial(totalOldCounts, pOld, size=nBootstraps)
    sampledNewCounts = np.random.multinomial(totalNewCounts, pNew, size=nBootstraps)

    for trial in range(nBootstraps):
      thisSurprisal = self.calcForState(sampledOldCounts[trial,:], sampledNewCounts[trial,:])
      print "Surprisal: ", thisSurprisal
      surprisals.append(thisSurprisal)
 
    return np.array(surprisals).var()

       

  def calcAll(self):
    surprisal = [] 

    for state in range(self.sparse1.shape[0]):
      if (self.Verbose):
        print("Working on State %d\n"%state)
      oldCounts, newCounts = self.prepareSparseMatrices(state) 
      surprisal.append(self.calcForState(oldCounts, newCounts))

    return surprisal


  def calcForState(self, oldCounts, newCounts):
    combCounts = oldCounts + newCounts
    tCombCounts = np.sum(oldCounts+newCounts)
    tOldCounts = np.sum(oldCounts)
    tNewCounts = np.sum(newCounts)

    combEntropy = self.calculateEntropy(combCounts)
    oldEntropy = self.calculateEntropy(oldCounts)
    newEntropy = self.calculateEntropy(newCounts)

    #print "combEntropy", combEntropy, "oldEntropy", oldEntropy, "newEntropy", newEntropy
   

    return (tCombCounts*combEntropy - tOldCounts*oldEntropy - tNewCounts*newEntropy) 

  def calcWeightedSurprisal(self):
    surprisalWeighted = []
    for state in range(self.sparse1.shape[0]):
      if (self.Verbose):
        print("Working on state %d"%state)

      oldCounts, newCounts = self.prepareSparseMatrices(state)
      combCounts = oldCounts + newCounts
      tCombCounts = np.sum(oldCounts+newCounts)     
      tOldCounts = np.sum(oldCounts)
      tNewCounts = np.sum(newCounts)
      combEntropy = self.calculateEntropy(combCounts)
      oldEntropy = self.calculateEntropy(oldCounts)
      newEntropy = self.calculateEntropy(newCounts)
      surprisalWeighted.append(combEntropy - (tOldCounts/tCombCounts)*oldEntropy - (tNewCounts/tCombCOunts)*newEntropy)
    return surprisalWeighted


# Main

if __name__ == '__main__':
  sparse1 = mmread(sys.argv[1])
  sparse2 = mmread(sys.argv[2])
  surprisalCalc = SurprisalCalculator(sparse1, sparse2)
  print ("Variances analytical: ", surprisalCalc.EstimateSurprisalVarianceAnalytical(15))

  oldCounts, newCounts = surprisalCalc.prepareSparseMatrices(15)
  print "Variance Bootstrapped: ", surprisalCalc.EstimateSurprisalVarianceBootstrap(oldCounts, newCounts) 
  #surprisal = surprisalCalc.calcAll()  
  #data = open('surprisal.dat', 'w')
  #for i in range(len(surprisal)):
   # data.write('%f\n'%(surprisal[i]))

            
        
