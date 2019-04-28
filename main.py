##Program Information######################################################
##
 # @file main.py
 #
 # @version 1.01
 #          Kristopher Moore (15 April 2019)
 #          Samantha Muellner (25 April 2019)
 #          Initial Program Build.
##

import sys
import numpy as np
import itertools

#main module responsible for handling dataset and calling functions
def main():
	readSequences("sample.fasta")


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
            
    print(tmp)
    return tmp
	



#---KRIS---
#generateSeeds is responsible for taking in a read sequence from a FASTA/FASTQ file
#and utilzing the K-Mer Exact Match Seeding process to develop K-mer seeds to be
#matched upon our reference sequence, return dictionary of seeds to calling method
def generateSeeds( file_name ):
	print("TODO")



#---Jennie---
#matchSeedOnReference is responsible for searching the reference sequence(s) for
#matching sets of given seeds, once found in sequence, added to a dictionary of
#"matched" seeds to be returned to calling method
def matchSeedOnReference( reference, seeds ):
	print("TODO")



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
	begining , ending = sw(reference,finalSeeds);
	print(reference[begining:ending]);


#ensure we can forward declare our methods.
if __name__=="__main__":
   main()
