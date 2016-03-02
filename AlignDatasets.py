#! python

#AlignDatasets.py

## Example command line entry:
## python AlignDatasets.py SmallclusterPart.fasta SmallclusterPart.groups clusterSetSeqList.fasta cluster.taxonomy


# Load the needed modules for the program
import sys, re


# Read in a Command arguments for consensus sequence file and a degapped aligned and screened fasta file
# Input other instructions from here
def commandLine():
	commands = sys.argv
	FastaFile = commands[1]
	GroupsFile = commands[2]
	ConsensusFile = commands[3]
	TaxFile = commands[4]
		
	return FastaFile, GroupsFile, ConsensusFile, TaxFile


# Create a dictionary for each sample and their sequences from the fasta file
def makeSampleArray(FastaFile, GroupsFile):
	Groups = open(GroupsFile, 'r')
	DataFileTemp = {}
	storage = []
	#This is not working correctly it is adding to the existing sets instead only to specific ones.
	for line in Groups:
		seqs, sample = line.split("\t")
		sampleStrip = sample.strip('\n')
		DataFileTemp.setdefault(sampleStrip, []).append(seqs)
		
	Groups.close()	
	Fasta = open(FastaFile, 'r')
	
	x = 1
	wholeSequence = []
	sequenceName = []
	sequenceInfo = {}
	for line in Fasta:
		if x%2 == 0:
			goodSequence = line.strip('\n')
			wholeSequence.append(goodSequence)
		else:
			name = line[1:]
			nameStrip = name.strip('\n')
			sequenceName.append(nameStrip)
		x = x + 1
	
	Fasta.close()
	
	for i in range(len(sequenceName)):
		
		sequenceInfo[sequenceName[i]] = wholeSequence[i]

	OverallData = {}
	
	
	for j in DataFileTemp:
		sequenceTable = {}
		tempSequenceList = DataFileTemp[j]

		for k in tempSequenceList:
			sequenceTable[k] = sequenceInfo[k]
		
		OverallData[j] = sequenceTable
		
	return OverallData


# Create a dictonary for each OTU and consensus sequence
def makeRefFile(ConsensusFile, TaxFile):
	reference = open(ConsensusFile, 'r')
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

	
# Create a 10bp key in the middle of each ~25 bp sequence to use for aligning
def makeTenBP(otuID, SeperatedDatabase):
	tenBPDatabase = {}
	
	for i in otuID:
		sequenceList = SeperatedDatabase[i]
		tempStorage = []
		
		for j in range(len(sequenceList)):
			if j == (len(sequenceList)-1):
				total = len(sequenceList[j]) // 2
				keySequence = sequenceList[j]
				tempStorage.append(keySequence[(total-5):(total+5)])
				
			
			else:
				keySequence = sequenceList[j]
				tempStorage.append(keySequence[8:18])
			
		tenBPDatabase[i] = tempStorage
		
	return tenBPDatabase

# Use the OTU-Regex dictionary to parse the fasta file dictionary.
# Level 1 match all 10bp keys, output to be matched file, unmatched file, 
# matched sequence list, unmatched sequence list
def getFirstMatch(tenBPDatabase, groupData, otuID):
	finalizedUnMatched = {}
	finalizedCompletedMatch = {}
	print("Level 1 of Alignment...")
	for i in groupData:
		completedMatch = {}
		unMatchedSeq = {}
		sequenceList = groupData[i]
		print(i)
		for j in sequenceList:
			x = 0
			tempSequence = sequenceList[j]
			for k in otuID:
				tempMatchList = tenBPDatabase[k]
				total = 0
				n = 0
				
				for l in tempMatchList:
								
					if l in tempSequence[(0+n):(25+n)]:
						total = total + 1
						n = n + 25
					else:
						total = total
						break
						
				if total == len(tempMatchList):
					completedMatch[j] = k
					break
				
				else:
					completedMatch = completedMatch
					if x == len(otuID)-1:
						unMatchedSeq[j] = sequenceList[j]
					else:
						unMatchedSeq = unMatchedSeq

				x = x + 1
		
		finalizedUnMatched[i] = unMatchedSeq
		finalizedCompletedMatch[i] = completedMatch
		
	return 	finalizedCompletedMatch, finalizedUnMatched

# Level 2 == 98.5% match of all 10bp keys
def getSecondMatch(finalizedCompletedMatch, finalizedUnMatched, tenBPDatabase, otuID):
	storedOverallData = {}
	# Make a 25 bp key to make searching quicker similar to function seperateReference 
	for i in finalizedUnMatched:
		SeperatedDatabase = {}
		sequenceList = finalizedUnMatched[i]
		for j in sequenceList:
			sequence = sequenceList[j]
			tempStorage = []
			seqLength = len(sequence)
			DivNum = (seqLength // 25) + 1
			x = 0
			for k in range(DivNum):
				if k == (DivNum-1):
					final = sequence[(0+x):]
					EndPart = final.strip('\n')
					tempStorage.append(EndPart)
				else:
					tempStorage.append(sequence[(0+x):(25+x)])
				
				x = x + 25
		
			SeperatedDatabase[j] = tempStorage
		storedOverallData[i] = SeperatedDatabase		
	
	# Search each 25bp chunk with the 10bp chunk and log % matching if less than 90% after 7 move on to next sequence
	print("Level 2 of Alignment...")
	for i in storedOverallData:
		print(i)
		sequenceList = storedOverallData[i]
		for j in sequenceList:
			individualKey = sequenceList[j]

			for k in range(len(otuID)):

				tempRefKey = tenBPDatabase[otuID[k]]
				x = 0
				total = 0
				try:
					for l in range(len(individualKey)):

						tempkeySequence = tempRefKey[l]
						try:
							if tempkeySequence in individualKey[l]:
								x = x + 10
							elif tempkeySequence[0:9] in individualKey[l]:
								x = x + 9
							elif tempkeySequence[0:8] in individualKey[l]:
								x = x + 8
							else:
								x = x + 7	
						except IndexError:
							x = x
						
						total = (x/((l+1)*10))*100
					
						if l == 1 and total <= 80:
							break
						elif l == ((len(individualKey)//2)-1) and total <= 90:
							break
						elif l == (len(individualKey)-1) and total >= 98.5:
						
							temp = finalizedCompletedMatch[i]
							temp[j] = otuID[k]
							finalizedCompletedMatch[i] = temp
						
							updateUnMatched =  finalizedUnMatched[i]
							del updateUnMatched[j]
							finalizedUnMatched[i] = updateUnMatched
							break
				
				except IndexError:
					total = (x/((l+1)*10))*100
					if l == 1 and total <= 80:
						break
					elif l == ((len(individualKey)//2)-1) and total <= 90:
						break
					elif l == (len(individualKey)-1) and total >= 98.5:
						temp = finalizedCompletedMatch[i]
						temp[j] = otuID[k]
						finalizedCompletedMatch[i] = temp
					
						updateUnMatched =  finalizedUnMatched[i]
						del updateUnMatched[j]
						finalizedUnMatched[i] = updateUnMatched
						break
							
				if total >= 98.5: break
						
	return finalizedCompletedMatch, finalizedUnMatched
					
# Level 3 < 98.5% match of all sequence matches
def getThirdMatch(finalizedCompletedMatch, finalizedUnMatched, OTUSeperatedDatabase, otuID):
	storedOverallData = {}
	# Make a 25 bp key to make searching quicker similar to function seperateReference 
	for i in finalizedUnMatched:
		SeperatedDatabase = {}
		sequenceList = finalizedUnMatched[i]
		for j in sequenceList:
			sequence = sequenceList[j]
			tempStorage = []
			seqLength = len(sequence)
			DivNum = (seqLength // 25) + 1
			x = 0
			for k in range(DivNum):
				if k == (DivNum-1):
					final = sequence[(0+x):]
					EndPart = final.strip('\n')
					tempStorage.append(EndPart)
				else:
					tempStorage.append(sequence[(0+x):(25+x)])
				
				x = x + 25
		
			SeperatedDatabase[j] = tempStorage
		storedOverallData[i] = SeperatedDatabase
	# Search each 25bp chunk with the 25bp chunk and log % matching if less than 90% after 7 move on to next sequence
	print("Level 3 of Alignment...")
	for i in storedOverallData:
		print(i)
		sequenceList = storedOverallData[i]
		for j in sequenceList:
			individualKey = sequenceList[j]

			for k in range(len(otuID)):

				tempRefKey = OTUSeperatedDatabase[otuID[k]]
				x = 0
				total = 0
				try:
					for l in range(len(individualKey)):

						tempkeySequence = tempRefKey[l]
						try:
							if tempkeySequence in individualKey[l]:
								x = x + 25
							elif tempkeySequence[0:24] in individualKey[l]:
								x = x + 24
							elif tempkeySequence[0:23] in individualKey[l]:
								x = x + 23
							elif tempkeySequence[0:22] in individualKey[l]:
								x = x + 22
							elif tempkeySequence[0:21] in individualKey[l]:
								x = x + 21
							elif tempkeySequence[0:20] in individualKey[l]:
								x = x + 20
							else:
								x = x + 17
								
						except IndexError:
							x = x
						
						total = (x/((l+1)*25))*100
					
						if l == 1 and total <= 80:
							break
						elif l == ((len(individualKey)//2)-1) and total <= 90:
							break
						elif l == (len(individualKey)-1) and total >= 98.5:
						
							temp = finalizedCompletedMatch[i]
							temp[j] = otuID[k]
							finalizedCompletedMatch[i] = temp
						
							updateUnMatched =  finalizedUnMatched[i]
							del updateUnMatched[j]
							finalizedUnMatched[i] = updateUnMatched
							break
				
				except IndexError:
					total = (x/((l+1)*10))*100
					if l == 1 and total <= 80:
						break
					elif l == ((len(individualKey)//2)-1) and total <= 90:
						break
					elif l == (len(individualKey)-1) and total >= 98.5:
						temp = finalizedCompletedMatch[i]
						temp[j] = otuID[k]
						finalizedCompletedMatch[i] = temp
					
						updateUnMatched =  finalizedUnMatched[i]
						del updateUnMatched[j]
						finalizedUnMatched[i] = updateUnMatched
						break
							
				if total >= 98.5: break
						
	return finalizedCompletedMatch, finalizedUnMatched
		
# If there is 97% or greater matching assign into that OTU and move on
# Create a new dictionary with the sample and OTU assignment
# Can use the degaped fasta file after everything has been completed



# Need to store sequences that don't fall into any category for later use



#Might need to cluster remaining sequences between each other and assign a unique OTU to each one
# Add this to the existing sample and OTU dictionary






# Count all the OTUs up and create a new dictionary with these values in it.










# Create a new shared file of the groups in the fasta file based on this information.





# Write this file as a tab text deliminated file.







# Run the program
def main():

	FastaFile, GroupsFile, ConsensusFile, TaxFile = commandLine()
	groupData = makeSampleArray(FastaFile, GroupsFile)
	otuID, referenceData = makeRefFile(ConsensusFile, TaxFile)
	OTUSeperatedDatabase = seperateReference(otuID, referenceData)
	tenBPDatabase = makeTenBP(otuID, OTUSeperatedDatabase)
	
	finalizedCompletedMatch, finalizedUnMatched = getFirstMatch(tenBPDatabase, groupData, otuID)
	
	finalizedCompletedMatch, finalizedUnMatched = getSecondMatch(finalizedCompletedMatch, finalizedUnMatched, tenBPDatabase, otuID)
	
	finalizedCompletedMatch, finalizedUnMatched = getThirdMatch(finalizedCompletedMatch, finalizedUnMatched, OTUSeperatedDatabase, otuID)
	
	test = finalizedCompletedMatch["SRR327639"]
	test2 = finalizedUnMatched["SRR327639"]
	print(len(test))
	print(len(test2))

if __name__ == '__main__': main()

