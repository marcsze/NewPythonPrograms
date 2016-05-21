#! python

# New program to compare two different files and output a number of columns on the 
# data infromation of samples that are the same 
# Uses the certain convertToList.py programs to edit the existing text file to make 
# it more amenable for downstream use

# Read in necessary programs and functions
import sys
from convertToList import readData, identifySamples

# Read in a Command arguments for the program
# Input other instructions from here
def commandLine():
	commands = sys.argv
	textToConvert = commands[1]
	pattern = commands[2]
	mappingFile = commands[3]
	return textToConvert, pattern, mappingFile

# Read in test data as a dictionary to modify file as needed
def readTestFile(mappingFile, pattern):
	data = open(mappingFile, 'r')
	testDict = {}
	x = 0
	for line in data:
		if x != 0:
			sampleID, barcode, linker, group = line.split("\t")
			if pattern in group:
				testDict[group.strip('\n')] = sampleID

		x = x + 1
	return testDict

# Make the comparison between the data and test data
def compareData(dataList, testDict):
	noMatchList = []
	YESMatchDict = {}
	for i, group in enumerate(dataList):
		try:
			sample = testDict[group]
			YESMatchDict[group] = sample
		except KeyError:
			noMatchList.append(group)
	return YESMatchDict, noMatchList
	
	
	

def main():
	# Need to change this so that existing commandLine is what is read.
	textToConvert, pattern, mappingFile  = commandLine()
	dataList = readData(textToConvert)
	goodData = identifySamples(dataList, pattern)
	testDict = readTestFile(mappingFile, pattern)
	YESMatchDict, noMatchList = compareData(dataList, testDict)
	
	print(YESMatchDict)
if __name__ == '__main__': main()
