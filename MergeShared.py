#! python

#MergeShared.py

## Example command line entry:
## python MergeShared.py testset.OTUrep.fasta testset.count_table clusterSetSeqList.fasta cluster.taxonomy


# Load the needed modules for the program
import sys


# Read in a Command arguments for consensus sequence file and a degapped aligned and screened fasta file
# Input other instructions from here
def commandLine():
	commands = sys.argv
	OTUTestFasta = commands[1]
	countsFile = commands[2]
	OTURefFasta = commands[3]
	TaxFile = commands[4]
		
	return OTUTestFasta, countsFile, OTURefFasta, TaxFile
	
	
	
# Read in necessary data and make two dictionaries
# One will store sequence information and the other will store counts for every sample
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
	
	
# Create a dictonary for each OTU and consensus sequence
def makeRefFile(OTURefFasta, TaxFile):
	reference = open(OTURefFasta, 'r')
	x = 1
	refSequence = []
	for line in reference:
		if x%2 == 0:
			refSequence.append(line)
		else:
			refSequence = refSequence
		x = x + 1
	reference.close()
	
	taxonomy = open(TaxFile, 'r')
	otuID = []
	x = 0
	for line in taxonomy:
		if x > 0:
			otu, total, taxID = line.split("\t")
			otuID.append(otu)
		x = x + 1
	taxonomy.close()
	
	referenceData = {}
	for i in range(len(otuID)):
		referenceData[otuID[i]] = refSequence[i]
		
	return otuID, referenceData
	

# Seperate reference sequences into smaller 25bp parts number of 25bp parts based on size
def seperateReference(otuID, referenceData):
	SeperatedDatabase = {}
	
	for i in otuID:
		tempStorage = []
		sequence = referenceData[i]
		seqLength = len(sequence)
		DivNum = (seqLength // 25) + 1
		x = 0
		
		for j in range(DivNum):
			if j == (DivNum-1):
				final = sequence[(0+x):]
				EndPart = final.strip('\n')
				tempStorage.append(EndPart)
			else:
				tempStorage.append(sequence[(0+x):(25+x)])
				
			x = x + 25
		
		SeperatedDatabase[i] = tempStorage
		
		
	return SeperatedDatabase
	
	
# Level 3 > 98.5% match of all sequence matches
def getThirdMatch(TestOTUsequences, TestSequenceCount, OTUSeperatedDatabase, otuID):
	SeperatedDatabase = {}
	finalizedCompletedMatch = {}
	finalizedUnMatched = {}
	
	# Make a 25 bp key to make searching quicker similar to function seperateReference 
	for i in TestOTUsequences:
		sequence = TestOTUsequences[i]
		seqLength = len(sequence)
		DivNum = (seqLength // 25) + 1
		x = 0
		tempStorage = []
		for k in range(DivNum):
			if k == (DivNum-1):
				final = sequence[(0+x):]
				EndPart = final.strip('\n')
				tempStorage.append(EndPart)
			else:
				tempStorage.append(sequence[(0+x):(25+x)])
				
			x = x + 25
		
		SeperatedDatabase[i] = tempStorage

	
	# Search each 25bp chunk with the 25bp chunk and log % matching if less than 90% after 7 move on to next sequence
	print("Trying to match OTUs...")
	
	for i in SeperatedDatabase:
		sequence = SeperatedDatabase[i]
		
		for k in range(len(otuID)):
			tempRefKey = OTUSeperatedDatabase[otuID[k]]
			x = 0
			total = 0
			try:
				for l in range(13):
					tempkeySequence = tempRefKey[l]
					try:
						if tempkeySequence in sequence[l]:
							x = x + 25
						#Modify this code to account for changes in the middle of the sequence
						else:
							testSequence = sequence[l]
							y = 0
							for nucleotide in range(len(tempkeySequence)):
								if tempkeySequence[nucleotide] == testSequence[nucleotide]:
									y = y + 1
								else:
									y = y
								
								if nucleotide == (len(tempkeySequence)-1):
									x = x + y
																
					except IndexError:
						x = x
						
					total = (x/((l+1)*25))*100
					
					if l == 1 and total <= 80:
						break
					elif l == 6 and total <= 90:
						break
					elif l == 12 and total >= 99:
						finalizedCompletedMatch[i] = otuID[k]
						break
						
			except IndexError:
				total = (x/((l+1)*25))*100
				if l == 1 and total <= 80:
					break
				elif l == 6 and total <= 90:
					break
				elif l == 12 and total >= 99:
					finalizedCompletedMatch[i] = otuID[k]
					break
							
			if total >= 99: break
			
			if k == (len(otuID)-1) and (i not in finalizedCompletedMatch):
				finalizedUnMatched[i] = sequence

	return finalizedCompletedMatch, finalizedUnMatched
	
	
def main():
	# Need to create a way to judge accuarcy of calls
	OTUTestFasta, countsFile, OTURefFasta, TaxFile = commandLine()
	TestOTUsequences, TestSequenceCount = createArrays(OTUTestFasta, countsFile)
	otuID, referenceData = makeRefFile(OTURefFasta, TaxFile)
	SeperatedDatabase = seperateReference(otuID, referenceData)
	finalizedCompletedMatch, finalizedUnMatched = getThirdMatch(TestOTUsequences, TestSequenceCount, SeperatedDatabase, otuID)
	counts = {}
	for i in finalizedCompletedMatch:
		otu = finalizedCompletedMatch[i]
		counts[otu] = counts.get(otu,0) + 1
		#print(i, otu)
	print(counts)
	retrievedOTUList = list(counts.items())
	retrievedOTUList.sort()
	#print(retrievedOTUList)
	print(len(finalizedCompletedMatch), len(finalizedUnMatched))
	

if __name__ == '__main__': main()








