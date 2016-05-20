#! python

## A program to convert a group of samples into a single column list
## based on user defined inputs of patterns to match.  
## Should be able to find patterns in any general text file
## Used to complement the compareDoc.py program
## The sample ID input is case sensitive
## An example of how to use is below:
	# python convertToList.py samplesToMatch.txt "FD" outputtest.txt


# Load the needed modules for the program
import sys

# Read in a Command arguments for the program
# Input other instructions from here
def commandLine():
	commands = sys.argv
	textToConvert = commands[1]
	pattern = commands[2]
	outputfile = commands[3]
	return textToConvert, pattern, outputfile
	
# Read in data as a list to modify with another function to identify pattern of interest
def readData(textToConvert):
	data = open(textToConvert, 'r')
	dataList = []
	for line in data:
		for sample in line.split(" "):
			sampleStrip = sample.strip('\n')
			dataList.append(sampleStrip)
	return dataList
	
# This function finds samples with the given user specificed pattern to keep
def identifySamples(dataList, pattern):
	goodData = []
	for i in dataList:
		if pattern in i:
			goodData.append(i)
	return goodData
	
# This function will output a single column txt file of the kept samples
def saveDataAsColumn(goodData, outputfile):
	outfile = open(outputfile, 'w')
	for i, sampleID in enumerate(goodData):
		if i == len(goodData) - 1:
			print("{0}".format(sampleID), end ='', file = outfile)
		else:
			print("{0}".format(sampleID), end ='\n', file = outfile)
	outfile.close()

# This runs the all the necessary functions of the program 
def main():
	textToConvert, pattern, outputfile = commandLine()
	dataList = readData(textToConvert)
	goodData = identifySamples(dataList, pattern)
	saveDataAsColumn(goodData, outputfile)
	
if __name__ == '__main__': main()










