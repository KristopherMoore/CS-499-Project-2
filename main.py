##Program Information######################################################
##
 # @file main.py
 #
 # @version 1.03
 #          Kristopher Moore (8 May 2019) -- Command-Line User Interface
 #          Samantha Muellner (25 April 2019)
 #          Initial Program Build.
##

import sys
import numpy as np
import itertools

#main module responsible for handling dataset and calling functions
def main():

    #initializations
    executionFlag = True
    targetSequenceFileName = ""
    readSequenceFileName = ""

    #greeting message
    print("\nCS:499-Bioinformatics, Local Alignment Using Seeds Tool:\n_____________________________________________________\n\n")

    #request location of TARGET sequence from user (The sequence we want to align to)
    targetSequenceFileName = requestSequence("TARGET")

    #request location of READ sequence from user (The sequence we want to use for alignment)
    readSequenceFileName = requestSequence("READ")

    #program loop, only allow exit when the user has specifically requested to do so
    while executionFlag:
        print("\nSelect a function to Execute\n\n"
        +"\t1: Perform Local Alignment\n"
        +"\t2: Change (Target)Sequence\n"
        +"\t3: Change (Read)Sequence\n"
        +"\t0: Quit\n")

        userInput = input("\tENTER NUMBER:")

        #check userChoices
        if(userInput == "0"):
            print("\nEXITING...")
            executionFlag = False

		#main tool function, initiate all actions
        elif(userInput == "1"):
            print("\nReading Target Sequence, Creating Target dictionary...")
            targetDict = readSequences(targetSequenceFileName)

            print("\nGenerating Seeds from Read Sequence, Creating Seeds dictionary...")
            seedsDict = generateSeeds(readSequenceFileName)

            print("\nMatching Seeds on our Target, Creating Matched dictionary...")
            matchedDict = matchSeedOnReference(targetDict, seedsDict)

            print("\nExtending matches, Creating FinalSeeds dictionary...")
            finalDict = extendMatchedSeeds(targetDict, matchedDict)

            print("\nPerforming Final Local Alignment...")
            localAlignment(targetDict, finalDict)

        #allow users to modify the Target Sequence Location
        elif(userInput == "2"):
            targetSequenceFileName = requestSequence("TARGET")

        #allow users to modify the Read Sequence Location
        elif(userInput == "3"):
            readSequenceFileName = requestSequence("READ")

        #invalid input
        else:
            print("Invalid Selection, Please select from the list")


#---AUSTIN---
#readSequences is responsible for reading in a FASTA formatted file, that may contain
#one or more sequences. If more than one sequence exists the data must be split into
#seperate dictionary structs and returned to the calling method
def readSequences( file_name ):
    dicts = []
    with open(file_name, "r") as fasta:
        line = fasta.readline()
        tmp = {}
        key = ""
        while line:
            if '>' in line:
                key = line.strip('\n')
                tmp[key] = ""
            else:
                tmp[key] = tmp[key] + line.strip('\n')

            line = fasta.readline()

    return tmp




#---KRIS---
#generateSeeds is responsible for taking in a read sequence from a FASTA/FASTQ file
#and utilzing the K-Mer Exact Match Seeding process to develop K-mer seeds to be
#matched upon our reference sequence, return dictionary of seeds to calling method
def generateSeeds( file_name ):

   sequence = readSequences( file_name )
   dict = {}

   #loop through each element, then at each element parse possibility of string of current index to index + k.
   for element in sequence:
      subString = ""
      byteSubString = b""

      if '>' in element:
        continue

      i = int(element)
      while i < element + k:
         if(element + k > len(sequence)):
            break
         if(sequence[i] == "\n"):
            i = element + k  #ensure we escape our while loop, no relevent mers to be found here
            break

         subString += sequence[i]
         i = i+1

      #check if the key already exists if so then add to its count, else make the count 1
      if subString in dict:
         dict[subString] = dict.get(subString, 0) + 1
      elif(len(subString) == k):
         dict[subString] = 1

   return dict



#---Jennie---
#matchSeedOnReference is responsible for searching the reference sequence(s) for
#matching sets of given seeds, once found in sequence, added to a dictionary of
#"matched" seeds to be returned to calling method
def matchSeedOnReference( reference, seeds ):
    matched = {}
    index = 0
    for seed in seeds.values():
        for refString in reference.values():
            if seed in refString and seed not in matched.values():
                matched[index] = seed
                index += 1

    return matched



#---Samantha---
# extendMatchedSeeds is responsible for utilizing the calculated matchedSeeds, and extending
# them along the referenceString, where applicable. This process is detailed further in
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4695832

def extendMatchedSeeds( reference, matchedSeeds ):
    sequence = []

    # get the alignment for every matched seed
    for item in matchedSeeds.values():
        sequence.append(sw(reference, item))

    return sequence



# helper function for SW
def matrix(a, b, match_score=3, gap_cost=2):
    H = np.zeros((len(a) + 1, len(b) + 1), np.int)

    for i, j in itertools.product(range(1, H.shape[0]), range(1, H.shape[1])):
        match = H[i - 1, j - 1] + (match_score if a[i - 1] == b[j - 1] else - match_score)
        delete = H[i - 1, j] - gap_cost
        insert = H[i, j - 1] - gap_cost
        H[i, j] = max(match, delete, insert, 0)
    return H

# helper function for SW
def traceback(H, b, b_='', old_i=0):
    # flip H to get index of **last** occurrence of H.max() with np.argmax()
    H_flip = np.flip(np.flip(H, 0), 1)

    i_, j_ = np.unravel_index(H_flip.argmax(), H_flip.shape)
    i, j = np.subtract(H.shape, (i_ + 1, j_ + 1))  # (i, j) are **last** indexes of H.max()

    if H[i, j] == 0:
        return b_, j

    b_ = b[j - 1] + '-' + b_ if old_i - i > 1 else b[j - 1] + b_

    return traceback(H[0:i, 0:j], b, b_, i)

# helper function for extendMatchedSeeds
def sw(a, b, match_score=3, gap_cost=2):
    a, b = a.upper(), b.upper()

    H = matrix(a, b, match_score, gap_cost)
    b_, pos = traceback(H, b)
    return pos, pos + len(b_)




#---Jake---
#localAlignment is responsible for perfoming the final local sequence alignment of the
#matched and extended seeds to the reference sequence(s), utilizing the Smith-Waterman Algorithm
def localAlignment( reference, finalSeeds ):
	for key in reference:
		refString = str(reference[key])
		finalString = str(finalSeeds)
		print(refString)
		print(finalString)
		begining , ending = sw(refString,finalString);
		print(reference[key][begining:ending]);



#helper function to request a sequence from the user
def requestSequence(fileType):
   UI = input("\nEnter File name for the (" + str(fileType) + ") Sequence (.FASTA format): ")
   return UI


#ensure we can forward declare our methods.
if __name__=="__main__":
   main()
