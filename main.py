##Program Information######################################################
##
 # @file main.py
 #
 # @version 1.00
 #          Kristopher Moore (15 April 2019)
 #          Initial Program Build.
##

import sys

#main module responsible for handling dataset and calling functions
def main():
	readSequences("TODO")


#---AUSTIN---
#readSequences is responsible for reading in a FASTA formatted file, that may contain
#one or more sequences. If more than one sequence exists the data must be split into 
#seperate dictionary structs and returned to the calling method
def readSequences( file_name ):
	print("TODO")



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
#extendMatchedSeeds is responsible for utilizing the calculated matchedSeeds, and extending
#them along the referenceString, where applicable. This process is detailed further in
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4695832 
def extendMatchedSeeds( reference, matchedSeeds ):
	print("TODO")



#---Jake---
#localAlignment is responsible for perfoming the final local sequence alignment of the 
#matched and extended seeds to the reference sequence(s), utilizing the Smith-Waterman Algorithm
def localAlignment( reference, finalSeeds ):
	print("TODO")


#ensure we can forward declare our methods.
if __name__=="__main__":
   main()