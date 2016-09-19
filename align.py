#!/bin/python

# Leah Briscoe 304-145-856

from math import *
import sys, re

GAP='-'
gap_penalty = -5

# Assumes only two fasta entries
def read_fasta(fasta):
	a=''
	a_name=''
	b=''
	b_name=''
	state=0
	for line in open(fasta, 'r'):
		line = line.rstrip('\n')
		line = line.rstrip('\r')
		if 0 == state:
			if None == re.search('^>', line):
				print('error: not valid fasta format')
			a_name = line[1:]
			state = state + 1
		elif 1 == state:
			if None == re.search('^>', line):
				a = a + line
			else:
				state = state + 1
				b_name = line[1:]
		else:
			b = b + line

	return (a,a_name,b,b_name)

def read_scoring(scoring):
	state = 0
	header = []
	sm = {} #create a dictionary
	for line in open(scoring, 'r'): #for every row
		line = line.rstrip('\n')#remove end of line character
		line = line.rstrip('\r')#remove carriage return
		if None == re.search('^#', line): #any line other than the header is split to create array of scores.
			a = line.split(' ')
			if 0 == state:
				header = a
				state = state + 1
			else:
				if len(a) != 1 + len(header): # number of scores must match header
					print('error: the number of entries did not match the header')
					sys.exit(2)
				for i in range(len(header)):
					sm[a[0] + header[i]] = int(a[i+1]) #e.g. AA,4 and AR,-1
   
	return sm

#Function to perform global alignment between two sequences, from end to end.
def global_alignment(a,b,sm):

	#initialize the dynamic programming matrix by taking into account the gap penalty
	walk = [[ gap_penalty * i if j == 0 else (gap_penalty * j if i == 0 else 0) for i in xrange(len(b) +1)] for j in xrange(len(a) +1)]

	#initialize two other matrices to store the i and j indices of the optimal path for walk[i][j]
	Path_i = [[ -1 if j == 0 and i == 0 else 0 for i in xrange(len(b) + 1)] for j in xrange(len(a) +1)]

	Path_j = [[ -1 if j == 0 and i == 0 else 0 for i in xrange(len(b) + 1)] for j in xrange(len(a) +1)]


	# for each cell consider 3 possible options, take the max
	#begin to fill the walk matrix
	for i in xrange(1, len(a) +1):
		for j in xrange(1, len(b) +1):
			options = [
				walk[i-1][j] + gap_penalty, #gap in seq a 

				walk[i][j-1] + gap_penalty, # gap in seq b 
				
				walk[i-1][j-1] + sm[a[i-1]+b[j-1]] #match/mismatch
			]
			
			#choose the best option to continue this alignment. Here, gaps are favored over match
			walk[i][j] = max(options)
			
			#store the indices of the best option
			Path_i[i][j], Path_j[i][j] = [(i-1, j), (i,j-1), (i-1, j-1)][options.index(walk[i][j])]
			
			
	#traceback
	a_aln, b_aln = '', ''
	i, j = len(a), len(b)
	# Creating the string backwards. 
	while Path_i[i][j] != -1 and Path_j[i][j] != -1:
		if Path_i[i][j] == i-1 and Path_j[i][j] == j-1:		
			#align a[i-1] to b[j-1]
			a_aln += a[i-1]
			b_aln += b[j-1]
			
		elif Path_i[i][j] == i and Path_j[i][j] == j-1:
			#align gap in a to amino acid in b
			a_aln += GAP # GAP = "-"
			b_aln += b[j-1]
		elif Path_i[i][j] == i-1 and Path_j[i][j] == j:
			#align gap in b to amino acid in a
			a_aln += a[i-1]
			b_aln += GAP
		
		i, j = Path_i[i][j], Path_j[i][j] #The contents of each cell are the indices for the matching letters for the step before
	
	return a_aln[::-1], b_aln[::-1], walk[len(a)][len(b)] #flip the strings backwards and return score also

def main():
	# get arguments
	
	if 3 != len(sys.argv):
		print('error: two input files required')
		sys.exit(2)
	
	#Sequences
	fasta = sys.argv[1]
	
	#Substitution Matrix
	scoring = sys.argv[2]
	


	# read in input files
	(a,a_name,b,b_name)=read_fasta(fasta)
	(sm)=read_scoring(scoring)

	#Perform global alignment
	(ag,bg,sg)=global_alignment(a,b,sm)

	# print output
	print "#GLOBAL SCORE=", sg
	print a_name,ag
	print b_name,bg


if __name__ == '__main__':
	main()

