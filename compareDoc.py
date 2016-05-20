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
	return textToConvert, pattern




def main():
	# Need to change this so that existing commandLine is what is read.
	textToConvert, pattern  = commandLine()
	dataList = readData(textToConvert)
	goodData = identifySamples(dataList, pattern)
	print(goodData)
	
if __name__ == '__main__': main()
