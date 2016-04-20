#! python

# First attempt at using a Kmer based approach to identify contamination. 

	
# Load the needed modules for the program
import sys, re


# Read in a Command arguments for files to match
# Input other instructions from here 
def commandLine():
	commands = sys.argv
	fastafile = commands[1]
	length = int(commands[2])
	
	return fastafile, length

# Read in data and create dictionary
def makeDataArray(inputfile, type):
	inputfile = open(inputfile, 'r')
	if type == "fasta":
		# Create specific environment variables
		x = 1
		names = []
		sequence = []
		DataDict = {}
		# Read data in and create a dictionary
		for line in inputfile:
			if x%2 != 0:
				newLine = re.sub('>', '', line)
				names.append(newLine.strip('\t\n'))
			else:
				sequence.append(line.strip('\n'))
			
			x = x + 1
			
		inputfile.close()
		
		for i, seqName in enumerate(names):
			DataDict[seqName] = sequence[i]
	# deal with data that is a group file
	elif type == "group":
		DataDict = {}
		for line in inputfile:
			seqName, group = line.split('\t')
			DataDict[seqName] = group.strip('\n')
	# deal with data that is a map file
	else:
		DataDict = {}
		for line in inputfile:
			number, seqname = line.split('\t')
			DataDict[number] = seqname.strip('\n')
		
	return DataDict	

# Function to generate total Kmers
def generateKmersList(sequence, length):
	tempKmerList = []
	total = len(sequence)
	for j,nucleotide in enumerate(sequence):
		if j+length > total-1: break
		else:
			kmerCheck = sequence[j:j+length]
			tempKmerList.append(kmerCheck)
		
	return tempKmerList

# Function to remove duplicates from website
# http://www.dotnetperls.com/duplicates-python
def removeDups(tempKmerList):
	output = []
	seen = set()
	for value in tempKmerList:
        # If value has not been encountered yet,
        # ... add it to both list and set.
		if value not in seen:
			output.append(value)
			seen.add(value)
	return output

#Function to get total unique kmers in a given dictionary
def getTotalKmers(fastaDict, length):
	kmerNumDict = {}
	for i in fastaDict:
		sequence = fastaDict[i]
		tempKmerList = generateKmersList(sequence, length)
		output = removeDups(tempKmerList)
		kmerNumDict[i] = len(output)
	return kmerNumDict

# Look at what Kmers are highly represented across sequences	
# Store the Kmers 
# Store the number of times each kmer comes up
	# Probably looking for something that shows up in all sequences
	
	
# Run the main program
def main():
	fastafile, length = commandLine()
	fastaDict = makeDataArray(fastafile, "fasta")
	kmerNumDict = getTotalKmers(fastaDict, length)
	print(kmerNumDict)

if __name__ == '__main__': main()