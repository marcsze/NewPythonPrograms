#! python

# Need to try using a cluster based approach.  Use both consensus and reference OTU sequence distance from a seed sequence.
# From here generate distance of files of interest also based on consensus and reference OTU sequences.
# Track which OTU these test set files lie closest to and then assign new OTUs
	# Need to get the consensus file for the testset.OTUrep.fasta
	# Need to get a seed sequence (although this can be any sequence)
	# Need to rewrite code accordingly

	
#MergebyConsensusShared.py

## Example command line entry:
## python MergeSharedAgain.py testset.OTUrep.fasta testset.count_table clusterSetSeqList.fasta output.txt
# Load the needed modules for the program
import sys, random

# Read in a Command arguments for consensus sequence file and a degapped aligned and screened fasta file
# Input other instructions from here
def commandLine():
	commands = sys.argv
	OTUTestFasta = commands[1]
	countsFile = commands[2]
	OTURefFasta = commands[3]
	outputFile = commands[4]
		
	return OTUTestFasta, countsFile, OTURefFasta, outputFile
	
	
def createArrays(OTUTestFasta, countsFile):
	# Make the first dictionary that will be used matching OTUs
	testOTUs = open(OTUTestFasta, 'r')
	
	x = 1
	wholeSequence = []
	sequenceName = []
	TestOTUsequences = {}
	for line in testOTUs:
		if x%2 == 0:
			goodSequence = line.strip('\n')
			wholeSequence.append(goodSequence)
		else:
			name = line[1:]
			nameStrip = name.strip('\n')
			sequenceName.append(nameStrip)
		x = x + 1

	for i in range(len(sequenceName)):
		
		TestOTUsequences[sequenceName[i]] = wholeSequence[i]

	# Make the second dictionary that will store the total count information
	x = 1
	groupNames = []
	TestSequenceCount = {}
	countTable = open(countsFile, 'r')
	
	for line in countTable:
		if x == 1:
			y = 0
			for i in line.split("\t"):
				if y >= 2:
					name = i.strip('\n')
					groupNames.append(name)
					
				else:
					groupNames = groupNames
				
				y = y + 1			
			x = x + 1
		
		else:
			y = 0
			tempcount = {}
			for i in line.split("\t"):
				if y == 0:
					seqName = i
					
				elif y >= 2:
					name = i.strip('\n')
					name = eval(name)
					tempcount[groupNames[y-2]] = name
				
				y = y + 1
				
				
				TestSequenceCount[seqName] = tempcount
	
	
	return TestOTUsequences, TestSequenceCount
	
	
def getCenteringSequence(TestOTUsequences):
	numberSelection = random.randrange(1, len(TestOTUsequences)+1)
	x = 1
	for i in TestOTUsequences:
		sequence = TestOTUsequences[i]
		if numberSelection == x:
			CenterSequence = sequence
		x = x + 1

	return(CenterSequence)
	

	
def main():
	# Need to create a way to judge accuarcy of calls
	OTUTestFasta, countsFile, OTURefFasta, outputFile = commandLine()
	TestOTUsequences, TestSequenceCount = createArrays(OTUTestFasta, countsFile)
	CenterSequence = getCenteringSequence(TestOTUsequences)
	
	

	
	
	
	
if __name__ == '__main__': main()