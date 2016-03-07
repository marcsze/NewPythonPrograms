#! python

# Need to try using a cluster based approach.  Use both consensus and reference OTU sequence distance from a seed sequence.
# From here generate distance of files of interest also based on consensus and reference OTU sequences.
# Track which OTU these test set files lie closest to and then assign new OTUs
	# Need to get the consensus file for the testset.OTUrep.fasta
	# Need to get a seed sequence (although this can be any sequence)
	# Need to rewrite code accordingly

	
#MergebyConsensusShared.py

## Example command line entry:
## python MergeSharedAgain.py testset.OTUrep.fasta testset.count_table clusterSetSeqList.fasta consensus.fasta output.txt

# Load the needed modules for the program
import sys, random

# Read in a Command arguments for consensus sequence file and a degapped aligned and screened fasta file
# Input other instructions from here
def commandLine():
	commands = sys.argv
	OTUTestFasta = commands[1]
	countsFile = commands[2]
	OTURefFasta = commands[3]
	OTUdataFile = commands[4]
	outputFile = commands[5]
		
	return OTUTestFasta, countsFile, OTURefFasta, OTUdataFile, outputFile
	
	
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
	
def makeRefFile(OTURefFasta, OTUdataFile):
	# Create a reference with the representative sequence
	OTUIDrefernce = open(OTUdataFile, 'r')
	otuID = []
	x = 2
	for line in OTUIDrefernce:
		if x%2 == 0:
			otu = line[1:]
			otuID.append(otu.strip('\n'))
		x = x + 1
	OTUIDrefernce.close()
	
	
	reference = open(OTURefFasta, 'r')
	x = 1
	refSequence = []
	for line in reference:
		if x%2 == 0:
			refSequence.append(line.strip('\n'))
		else:
			refSequence = refSequence
		x = x + 1
	reference.close()

	OTUSequences = {}
	for i in range(len(otuID)):
		OTUSequences[otuID[i]] = refSequence[i]
	
	return OTUSequences

def getRefOTUSimScores(CenterSequence, OTUSequences):
	RefOTUSim = {}
	for i in OTUSequences:
		sequence = OTUSequences[i]
		total = 0
		match = 0
		matchLength = 1
		for j, nucleotide in enumerate(sequence):
			if nucleotide == CenterSequence[j]:
				match = match + 1
				total = total + 1
			else:
				total = total + 1
			
			if matchLength == 300:
				break
		
			matchLength = matchLength + 1
		
		RefOTUSim[i] = 100 - ((match/total)*100)
	
	return RefOTUSim
		
	
def findTotalDifference(TestOTUSim, TestOTUSim2, RefOTUSim, RefOTUSim2):
	tempDictionary = {}
	minDifference = 100
	for i in TestOTUSim:
		
		TestScore = TestOTUSim[i]
		for j in RefOTUSim:
	
			RefScore = RefOTUSim[j]
			tempDifference = abs(TestScore-RefScore)
			TestScore2 = TestOTUSim2[i]
			RefScore2 = RefOTUSim2[j]
			tempDifference2 = abs(TestScore2-RefScore2)
			TotalDifference = tempDifference + tempDifference2
			if TotalDifference < minDifference:
				tempDictionary[i] = j
				minDifference = TotalDifference
			else:
				tempDictionary = tempDictionary
			
	return tempDictionary		
	
	
def main():
	# Need to create a way to judge accuarcy of calls
	OTUTestFasta, countsFile, OTURefFasta, OTUdataFile, outputFile = commandLine()
	TestOTUsequences, TestSequenceCount = createArrays(OTUTestFasta, countsFile)
	OTUSequences = makeRefFile(OTURefFasta, OTUdataFile)
	# Generate Information for first axis Reference
	CenterSequence = getCenteringSequence(TestOTUsequences)
	RefOTUSim = getRefOTUSimScores(CenterSequence, OTUSequences)
	# Generate Information for second axis Reference
	CenterSequence2 = getCenteringSequence(TestOTUsequences)
	RefOTUSim2 = getRefOTUSimScores(CenterSequence2, OTUSequences)
	# Generate Information for first axis Test Set
	TestOTUSim = getRefOTUSimScores(CenterSequence, TestOTUsequences)
	# Generate Information for second axis Test Set
	TestOTUSim2 = getRefOTUSimScores(CenterSequence2, TestOTUsequences)
	
	# Find the cumilative difference between the two axes
	tempDictionary = findTotalDifference(TestOTUSim, TestOTUSim2, RefOTUSim, RefOTUSim2)
	
	print(tempDictionary)
	

	
	
	
	
if __name__ == '__main__': main()