#! python

#MergebyConsensusShared.py

## Example command line entry:
## python MergebyConsesusShared.py testset.OTUrep.fasta testset.count_table consensus.fasta clusterSetSeqList.fasta output.txt

# Load the needed modules for the program
import sys

# Read in a Command arguments for consensus sequence file and a degapped aligned and screened fasta file
# Input other instructions from here
def commandLine():
	commands = sys.argv
	OTUTestFasta = commands[1]
	countsFile = commands[2]
	OTURefFasta = commands[3]
	OTURef2Fasta = commands[4]
	outputFile = commands[5]
		
	return OTUTestFasta, countsFile, OTURefFasta, OTURef2Fasta, outputFile

	
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


def makeRefFile(OTURefFasta, OTURef2Fasta):
	# Create a reference with the consensus sequences
	reference = open(OTURefFasta, 'r')
	x = 1
	refSequence = []
	otuID = []
	for line in reference:
		if x%2 == 0:
			refSequence.append(line)
		else:
			otuID.append(line[1:])			
			refSequence = refSequence
		x = x + 1
	reference.close()

	OTUSequences = {}
	for i in range(len(otuID)):
		name = otuID[i]
		name = name.strip('\n')
		OTUSequences[name] = refSequence[i]
	
	# Create a reference with the representative sequence
	reference2 = open(OTURef2Fasta, 'r')
	x = 1
	refSequences = []
	for line in reference2:
		if x%2 == 0:
			refSequences.append(line)
		
		x = x + 1
		
	reference2.close()
	
	OTUSequences2 = {}
	for i in range(len(otuID)):
		name = otuID[i]
		name = name.strip('\n')
		OTUSequences2[name] = refSequences[i]
	
	
	
	return OTUSequences, OTUSequences2
	
# Seperate consensus sequence by N and remove '' from the list
def seperateReference(OTUSequences):
	SeperatedDatabase = {}
	
	for i in OTUSequences:
		tempStorage = []
		sequence = OTUSequences[i]
		seqLength = len(sequence)
		for j in sequence.split("N"):
			name = j.strip('\n')
			if name == '':
				tempStorage = tempStorage
			else:
				tempStorage.append(name)
			
		SeperatedDatabase[i] = tempStorage
	
		
	return SeperatedDatabase

# Generate a total length of the consensus sequences without N's	
def getOTUconsensusTotal(SeperatedDatabase):
	OTUConsenTotal = {}
	for i in SeperatedDatabase:
		consensus = SeperatedDatabase[i]
		total = 0
		for j in consensus:
			length = len(j)
			total = total + length
		OTUConsenTotal[i] = total
	return OTUConsenTotal

# Generate first round of matching using the consensus sequences instead of the rep OTU
def getMatch(TestOTUsequences, SeperatedDatabase, OTUConsenTotal):
	finalizedCompletedMatch = {}
	finalizedUnMatched = {}
	OTUlength = len(SeperatedDatabase)
	
	print("Trying to Match OTUs...")
	for i in TestOTUsequences:
		#print(i)
		sequence = TestOTUsequences[i]
		counter = 0
		for otu in SeperatedDatabase:
			counter = counter + 1
			OTUConsensus = SeperatedDatabase[otu]
			total = OTUConsenTotal[otu]
			previous = 0
			PrevchunkLength = 0
			prevLocation = 0
			x = 0
			for chunk in OTUConsensus:
				length = len(OTUConsensus)
				
				if length == 1:
					if chunk in sequence:
						x = x + len(chunk)
					
				else:
					# Using <test>.find() to get match and where the match is
					if chunk in sequence:
						location = sequence.find(chunk, (prevLocation))
						x = x + len(chunk)
						prevLocation = location + 1
						
						# Break out of functions if preLocation is at all behind the current match
						if prevLocation <= 0:
							break
					if prevLocation <= 0 :
						break
				if prevLocation <= 0:
					break
					
			if length == 1 and (x/total) == 1:
				finalizedCompletedMatch[i] = otu
				break
			elif (x/total) == 1:
				finalizedCompletedMatch[i] = otu
				break
			elif counter == OTUlength:
				finalizedUnMatched[i] = sequence
	
	return finalizedCompletedMatch, finalizedUnMatched
	
# Generate second round of matching using rep OTU
def getMatch2(finalizedCompletedMatch, finalizedUnMatched, OTUSequences2):
	print("Trying to Match OTUs Round 2...")
	for i in finalizedUnMatched:
		testSequence = finalizedUnMatched[i]
		
		for j in OTUSequences2:
			x = 0
			refSequence = OTUSequences2[j]
			
			for k, nucleotide in enumerate(testSequence[0:306]):
				if nucleotide == refSequence[k]:
					x = x + 1
				else:
					x = x
				
				if k == 163 and x/(k+1) <= 0.96:
					break
				elif k == 50 and x/(k+1) <= 0.6:
					break
				elif k == 305 and x/(k+1) >= 0.99:
					finalizedCompletedMatch[i] = j
					break
		
			if k == 305 and x/(k+1) >= 0.99:
				break
	
	# Removing matched OTUs from the unMatched dictionary
	for i in finalizedCompletedMatch:
		if i in finalizedUnMatched:
			del finalizedUnMatched[i]

				
	return finalizedCompletedMatch, finalizedUnMatched

# Adding any unmatched sequences to a new OTU	
def UnmatchedSeqToOTU(SeperatedDatabase, finalizedCompletedMatch, finalizedUnMatched):
	totalOTUs = len(SeperatedDatabase)
	for i in finalizedUnMatched:
		totalOTUs = totalOTUs + 1
		name = "Otu{0}".format(totalOTUs)
		finalizedCompletedMatch[i] = name

	return finalizedCompletedMatch
	
# Combine OTUs that are duplicated
def combineDuplicates(finalizedCompletedMatch, TestSequenceCount):
	OTUDict = {}
	for i in finalizedCompletedMatch:
		otuID = finalizedCompletedMatch[i]
		counts = TestSequenceCount[i]
		tempStorage = {}
		if otuID in OTUDict:
			storedCounts = OTUDict[otuID]
			for j in storedCounts:
				tempStorage[j] = counts[j] + storedCounts[j]
			OTUDict[otuID] = tempStorage
		else:
			OTUDict[otuID] = counts
			
	return OTUDict
		

		
# Making a tab delimited table with the necessary information		
def makeSharedFile(OTUDict, outputFile):
	outfile = open(outputFile, 'w')
	x = 1
	for i in OTUDict:
		series = OTUDict[i]
		total = len(series)
		if x == 1:
			y = 0
			print("{0}".format("OTU"), end ='\t', file = outfile)
			for j in series:
				if y == len(series) - 1:
					print("{0}".format(j), end ='\n', file = outfile)
				else:
					print("{0}".format(j), end ='\t', file = outfile)
					y = y + 1
			y = 0
			for k in series:
				OTUValue = series[k]
				if y == 0:
					print("{0}".format(i), end ='\t', file = outfile)
					print("{0}".format(OTUValue), end ='\t', file = outfile)
				elif y == total-1:
					print("{0}".format(OTUValue), end ='\n', file = outfile)
				else:
					print("{0}".format(OTUValue), end ='\t', file = outfile)
				y = y + 1
		
		else:
			y = 0
			for k in series:
				OTUValue = series[k]
				if y == 0:
					print("{0}".format(i), end ='\t', file = outfile)
					print("{0}".format(OTUValue), end ='\t', file = outfile)
				elif y == total-1:
					print("{0}".format(OTUValue), end ='\n', file = outfile)
				else:
					print("{0}".format(OTUValue), end ='\t', file = outfile)
				y = y + 1
		x = x + 1
		
	outfile.close()
		
def main():
	# Need to create a way to judge accuarcy of calls
	OTUTestFasta, countsFile, OTURefFasta, OTURef2Fasta, outputFile = commandLine()
	TestOTUsequences, TestSequenceCount = createArrays(OTUTestFasta, countsFile)
	OTUSequences, OTUSequences2 = makeRefFile(OTURefFasta, OTURef2Fasta)
	SeperatedDatabase = seperateReference(OTUSequences)
	OTUConsenTotal = getOTUconsensusTotal(SeperatedDatabase)
	finalizedCompletedMatch, finalizedUnMatched = getMatch(TestOTUsequences, SeperatedDatabase, OTUConsenTotal)
	finalizedCompletedMatch, finalizedUnMatched = getMatch2(finalizedCompletedMatch, finalizedUnMatched, OTUSequences2)
	finalizedCompletedMatch = UnmatchedSeqToOTU(SeperatedDatabase, finalizedCompletedMatch, finalizedUnMatched)
	OTUDict = combineDuplicates(finalizedCompletedMatch, TestSequenceCount)
	makeSharedFile(OTUDict, outputFile)
	

	
	
	
	
if __name__ == '__main__': main()





















