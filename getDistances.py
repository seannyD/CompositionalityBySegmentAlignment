# This script must be run with python3

# The mantel test code comes from Jon Carr:
# http://www.jonwcarr.net/blog/2014/9/19/a-guide-to-the-mantel-test-for-linguists
# https://github.com/jwcarr/MantelTest


from lingpy import *
import codecs
import Mantel


def getAlignmentWeights(filename):
	""" Read the alignment weights file"""
	
	# Open file
	o = codecs.open(filename,'r','utf-8')
	d = o.read()
	o.close()
	d = d.replace("\r","")
	# get the data
	dat = [x.split(",")[1:] for x in d.split("\n") if len(x)>1]
	# get column names
	colNames = dat[0]
	#(exclude the first line, which should be the column names)
	dat = dat[1:]
	
	# Build a dictionary of transition weights between segments
	# (read the matrix)
	ret = {}
	for i in range(len(dat)):
		row = dat[i]
		for j in range(len(row)):
			ret[(colNames[i],colNames[j])] = float(dat[i][j])
	
	return ret
	
def readWordFile(wordFile):
	""" Read the language file- one word per line, 
	with the meaning dimensions divided by tabs"""
	
	# Read the file text
	wordFileLink = codecs.open(wordFile,'r','utf-8')
	wordFileText = wordFileLink.read()
	wordFileLink.close()
	wordFileText = wordFileText.replace("\r","")
	
	# Get the first column (words)
	words = [line.split(',')[0] for line in wordFileText.split('\n')[1:] if len(line)>0]
	
	# Get the meaning columns (a list of lists)
	meanings = [line.split(",")[1:] for line in wordFileText.split('\n')[1:] if len(line)>0]
	
	return words,meanings

def getMeaningDistances(meanings):
	# matrix to store 
	meaningDistances = []
	for m1 in range(len(meanings)):
		for m2 in range(len(meanings)):
			meaning1 = meanings[m1]
			meaning2 = meanings[m2]
			# the Mantel test only requires the upper triangle of the matrix
			if m1 > m2:
				d = sum([meaning1[i]==meaning2[i] for i in range(len(meaning1))]) / float(len(meaning1))
				meaningDistances.append(d)
	return meaningDistances
			
def getWordDistances(words, alignmentWeights, gap_weight):
	""" Get word distances from optimal string alignment"""
	
	# This makes a csv file with a matrix of distances between each word pair
	out = "Words,"+",".join(words)
	
	# to store the distances for the mantel test:
	distanceList = []
		
	# For every word
	for n1 in range(len(words)):
		out += "\n"+words[n1]+","
		# for every word
		for n2 in range(len(words)): 
				wordA = words[n1]
				wordB = words[n2]
				# get the alignment (Needleman-Wunsch algorithm nw_align from LingPy)
				# similarity score
				almA, almB, simAB = nw_align(wordA, wordB, scorer=alignmentWeights, gap=gap_weight)
				almAX, almBX, simAA = nw_align(wordA, wordA, scorer=alignmentWeights, gap=gap_weight)
				almAY, almBY, simBB = nw_align(wordB, wordB, scorer=alignmentWeights, gap=gap_weight)
				
				#  Divided by length of longest string?
				#sim = sim / float(max([len(x) for x in almA+almB]))
				# No: divide by sum of self-self distance:
				dist = 1 - ((2 * simAB) / (simAA + simBB))
				# TODO: check that we're calculating distance, and not similarity. 
				#  i.e. is "1 -" correct?
				
				# show the alignments
				print(' '.join(almA)+"\n" + ' '.join(almB) + "     distance={0}".format(dist))
				print("----")
				out += str(dist)+","
				# the Mantel test only requires the upper triangle of the matrix
				if n1 > n2:
					distanceList.append(dist)
	return distanceList, out

def getDistances(wordFile, weightsFile, gap_weight, outFile):
	""" Get word and meaning distances """
	
	# Read the language file
	words,meanings = readWordFile(wordFile)
	
	# get distances between meanings
	meaningDistances = getMeaningDistances(meanings)
	
	# Read the alignment weights file
	alignmentWeights = getAlignmentWeights(weightsFile)
	
	# Check that the alignment file has all the necessary letters
	alignmentLetters = [x[0] for x in alignmentWeights.keys()]
	for word in words:
		for letter in word:
			if not letter in alignmentLetters:
				print("Error: Letter "+letter+" not found in alignment matrix for word "+word)


	# get the distances
	wordDistances, wordDistancesString = getWordDistances(words,alignmentWeights,gap_weight) 
	
	# Write the word distances to a file
	o = codecs.open(outFile,"w",'utf-8')
	o.write(wordDistancesString)
	o.close()
	return wordDistances,meaningDistances
	

def getCompositionality(languageFile,AlignmentWeights,gap_weight,outFile, mantel_test_method='pearson',mantel_test_tail='upper',mantel_test_perms=10000):
	""" Given a langauge file and an alignment file, 
	Work out the word and meaning distances, then compute a 
	Mantel test between them """
	wordDistances,meaningDistances = getDistances(languageFile,AlignmentWeights, gap_weight,outFile)

	# Do mantel test using Jon Carr's code
	# 
	mantelResult = Mantel.test(wordDistances,meaningDistances,perms=mantel_test_perms, method=mantel_test_method, tail=mantel_test_tail)

	# Print results
	print("Mantel test result")
	print("r=",mantelResult[0],"p=",mantelResult[1], "z=",mantelResult[2])
	return mantelResult

