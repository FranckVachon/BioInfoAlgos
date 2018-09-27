import numpy as np
import math
from copy import deepcopy

#from main import CostFunctionSimple

alphabet = ["a", "t", "c", "g", "_"]

class NeedlemanWunshAlgo:
    def __init__(self, s1, s2, costFunctionObj):
        # The sequences we want to align
        self.s1 = s1
        self.s2 = s2
        # An object that encapsulate the function cost - so we can easily plug in a different one
        self.cost = costFunctionObj
        #A matrix containing the result - where the dynamic programming takes place
        self.results = []

        self.initMatrix()


    def execNeedlemanWunsh(self):
        # Applies needleman-wunsh algo to s1, s2 to determine optimal alignment(s) using dy prog
        #WARNING: index in python start @ 0, so take care...

        for i,s1i in enumerate(self.s1):
            i+=1
            for j,s2j in enumerate(self.s2):
                j+=1
                dum_costs = []

                #cost of substitution
                dum_costs.append(self.results[i-1,j-1] + self.cost.calculateCost(s1i,s2j))

                #cost of insertion
                dum_costs.append(self.results[i-1,j] + self.cost.calculateCost(s1i,None))

                #cost of deletion
                dum_costs.append(self.results[i,j-1] + self.cost.calculateCost(None,s2j))

                #Now we select the min & put in the corresponding spot in
                # self.results if we're not in the initialized part of the matrix
                self.results[i,j] = min(dum_costs )

                #Now we have the result matrix - need to reconstruct the path... We could store
                # pointers. Or we could just do nearest neighboors, choise mininum values &
                #  rebuild the paths from that.

        solutions = self.reconstructPaths()


        return solutions

    def execNeedlemanWunshSimilarity(self):
        # Applies needleman-wunsh algo to s1, s2 to determine optimal alignment(s) using dy prog
        #WARNING: index in python start @ 0, so take care...

        for i,s1i in enumerate(self.s1):
            i+=1
            for j,s2j in enumerate(self.s2):
                j+=1
                dum_costs = []

                #cost of substitution
                dum_costs.append(self.results[i-1,j-1] + self.cost.calculateCost(s1i,s2j))

                #cost of insertion
                dum_costs.append(self.results[i-1,j] + self.cost.calculateCost(s1i,None))

                #cost of deletion
                dum_costs.append(self.results[i,j-1] + self.cost.calculateCost(None,s2j))

                #Now we select the min & put in the corresponding spot in
                # self.results if we're not in the initialized part of the matrix
                self.results[i,j] = max(dum_costs )

                #Now we have the result matrix - need to reconstruct the path... We could store
                # pointers. Or we could just do nearest neighboors, choise mininum values &
                #  rebuild the paths from that.

        solutions = self.reconstructPaths()


        return solutions

    def initMatrix(self):
        self.results = np.full((len(self.s1)+1, len(self.s2)+1), None)
        self.results[:,0] = [x*self.cost.getInitCost() for x in range(0,self.results.shape[0])]
        self.results[0,:] = [x*self.cost.getInitCost() for x in range(0,self.results.shape[1])]

    def reconstructPaths(self):
        """ This returns a sequence of the kind:
        acgt_acc_gtcu

        TODO: YOU SHOULD POP FROM S1,S2 COPIES BECAUSE THEN YOU ENSURE YOU USE EVERYTHING, DON'T WORRY ABOUT THE INDEXING AT THE EDGE CASES ETC.... AS SEEN IN CLASSE
        :return:
        """

        # This will contain partial sequence objects, so we can just
        # iterate over them & complete them
        partialSequences = []

        # This will contain the completed alignement sequences.
        # There may be a few equivalent one
        completeAlignmentSequences = []
        #Starting bottom right corner.
        # The side sequence is the "reference" sequence
        partialSequences.append(PartialSequence((self.results.shape[0]-1,self.results.shape[1]-1),self.s1[-1]))

        #now we can just iterate over all the partial sequences, send them to the neighbors thing
        seq = partialSequences.pop()
        while seq:
            #TODO: check  if I should compare to s1 or s2...
            while len(seq.sequenceString) < len(self.s1):
                #Because we will update it in place, so if we have a fork it won't work...
                oldSequence = deepcopy(seq)
                neighboursCells = self.getNeighboursMin(seq.currentPos)
                # add the first new neighbor to the current sequence.
                # This will determine what letter to append to the alignment sequence
                nextNeighbour = neighboursCells.pop() if neighboursCells else None
                if nextNeighbour:
                    seq.addLetter(self.nextLetter(seq, nextNeighbour))
                    seq.currentPos = nextNeighbour
                # if there are other neighbours left, then we have a branch. New sequences are
                #  created and added to the pile of partialSequences. We create them & add them
                # to partialSequences. They'll be handled by the while seq. But we really treat
                #  just one sequences at a time, and pile on others for later treatment by the
                # first while
                nextNeighbour = neighboursCells.pop() if neighboursCells else None

                while nextNeighbour:
                    newSequence = deepcopy(oldSequence)
                    newSequence.addLetter(self.nextLetter(newSequence,nextNeighbour))
                    newSequence.currentPos = nextNeighbour
                    partialSequences.append(newSequence)
                    nextNeighbour = neighboursCells.pop() if neighboursCells else None
            # Out of the current sequence while
            # This is a full sequence now - move it over & do the next one
            # But we inverse it before so that the string is in same direction as input
            seq.setSequenceOrder()
            completeAlignmentSequences.append(seq)
            if len(partialSequences)>0:
                seq = partialSequences.pop()
            else:
                seq = None

        print(self.results)
        return completeAlignmentSequences


    def nextLetter(self,seq, nextNeighbourPosition):
        """Returns the next letter (atcg or _) for that sequence from s1 perspective (side sequence)"""
        # gets the letter corresponding to that coordinate, based on corresponding places in S1 Mind the indexes...

        #TODO: test the indexing (0,1 to start?)
        # If (i-1, j-1) or (i, j-1) between seq/next then we return a letter from s1
        # Otherwise it's a delection so _

        if seq.currentPos[0] != nextNeighbourPosition[0] and seq.currentPos[1] != nextNeighbourPosition[1]:
            #same letter, return a letter
            return self.s1[max(nextNeighbourPosition[0]-1,0)]
        elif seq.currentPos[0] == nextNeighbourPosition[0]:
            return "_"

        elif seq.currentPos[1] == nextNeighbourPosition[1]:
            return self.s1[max(nextNeighbourPosition[0]-1,0)]
        else:
            print("Something wrong nextLetter()")


    def getNeighboursMin(self, pos):
        '''
        :param tuple for position (as coherent with self.results dimensions, which is +1 compared to corresponding s1, s2 position
        :return: a list of tuples representing the position of the minimum neighboors. Returns all of them if there is a tie in the value
        '''


        alignOrSubstitution = (pos[0]-1,pos[1]-1)
        delete = (pos[0],pos[1]-1)
        insert = (pos[0]-1,pos[1])
        minNeighbours = []
        scores = [(self.results[insert],insert),(self.results[delete], delete),(self.results[alignOrSubstitution],alignOrSubstitution)]
        scores.sort(key=lambda tup: tup[0])

        minScore = scores[0][0]
        minNeighbours.append(scores.pop(0)[1])

        while len(scores):
            s = scores.pop()
            if s[0]==minScore:
                minNeighbours.append(s[1])
        return minNeighbours


class PartialSequence:
    def __init__(self, currentPosition, sequenceString=""):
        """
        :param currentPosition: the position in self.result. Warning: this position is +1 with regard to (s1, s2)
        :param sequenceString:
        """
        self.currentPos=currentPosition
        self.sequenceString=sequenceString

    def addLetter(self, letter):
        self.sequenceString += letter

    def setSequenceOrder(self):
        # reverses the string order
        self.sequenceString = self.sequenceString[::-1]

# Cost functions

class CostFunctionSimple:
    # the simple cost function where del/insert and substitution all cost the same
    def __init__(self, indel):
        #cost attributed to delete/insert etc.sssss
        self.indel = indel

    def calculateCost(self, s1i, s2j):

        return self.indel if s1i != s2j else 0
    def getInitCost(self):
        """Returns the cost to be used to initialize the first col/row of the cost matrix
        - may vary depending on algo"""
        return self.indel

class CostFunctionWeighted:
    # the simple cost function where del/insert and substitution all cost the same
    def __init__(self, costDelete, costSubstitute):
        #cost attributed to delete/insert etc.sssss
        self.indel = costDelete
        self.costSubstitute = costSubstitute

    def calculateCost(self, s1i, s2j):

        #IF one element is None, then we have a deletion/insert. If both are present it's either = (no cost) or substitution (cost)
        if s1i and s2j:
            return self.costSubstitute if s1i!=s2j else 0
        else:
            return self.indel

        return self.costAll if s1i != s2j else 0

    def getInitCost(self):
        """Returns the cost to be used to initialize the first col/row of the cost matrix
        - may vary depending on algo"""
        return self.indel



class CostFunctionSimilarity:
    # the simple cost function where del/insert and substitution all cost the same
    def __init__(self, match, mismatch, indel):
        #cost attributed to delete/insert etc.sssss
        self.match = match
        self.mismatch = mismatch
        self.indel = indel

    def calculateCost(self, s1i, s2j):

        #IF one element is None, then we have a deletion/insert. If both are present it's either = (no cost) or substitution (cost)
        if s1i and s2j:
            return self.mismatch if s1i!=s2j else self.match
        else:
            return self.indel

        return self.mismatch if s1i != s2j else self.match

    def getInitCost(self):
        """Returns the cost to be used to initialize the first col/row of the cost matrix
        - may vary depending on algo"""
        return self.indel



s1 = ["a","a","c","t","g","g","a","t","g","c","t","t"]
s2 = ["a","c","t","g","g","a","c","t","g","a","g","t"]
#s2 = ["t","c","g","a"]
#s1 = ["c","c","g","a"]


needleSimple = NeedlemanWunshAlgo(s1,s2,CostFunctionSimple(1))
#needleWeighted = NeedlemanWunshAlgo(s1,s2,CostFunctionWeighted(1,2))

print("Simple cost function:")

solutions = needleSimple.execNeedlemanWunsh()
#solutions = needleSimple.execNeedlemanWunshSimilarity()


print("S1: ",s1)
print("S2: ",s2)
for s in solutions:
    print("FI: ",list(s.sequenceString))

print("Weighted cost function:")
#solutions = needleWeighted.execNeedlemanWunsh()
#for s in solutions:
    #print(s.sequenceString)
#    pass