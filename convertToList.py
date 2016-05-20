#! python


## A program to convert a group of samples into a single column list
## based on user defined inputs of patterns to match.  
## Should be able to find patterns in any general text file
## Used to complement the compareDoc.py program



# Load the needed modules for the program
import sys

# Read in a Command arguments for the program
# Input other instructions from here
def commandLine():
	commands = sys.argv
	textToConvert = commands[1]
	pattern = commands[2]
		
	return textToConvert, pattern
	
# Read in data as a list to modify with another function to identify pattern of interest
def readData(textToConvert):
	data = open(textToConvert, 'r')
	dataList = []
	for line in data:
		for sample in line.split(" "):
			sampleStrip = sample.strip('\n')
			dataList.append(sampleStrip)
		
	return dataList
	

def identifySamples(dataList):
	
	
def main():
	textToConvert, pattern = commandLine()
	dataList = readData(textToConvert)
	print(dataList)
	
	
if __name__ == '__main__': main()










