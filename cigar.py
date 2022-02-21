import sys
from Bio.Seq import Seq 
from Bio import pairwise2
from Bio.pairwise2 import align, format_alignment


templates = [] #list to olds files containing sequences from the template
NEATs = [] #list to hold files containing sequences from the NEAT reads
cigars = []

def alginSeqs():
    """
    This function uses pairwise2 to locally align the NEAT sequence to the template sequence.
    Calls the makeCigar function after each alignment.
    """
    for tmplt, neat in zip(templates, NEATs):
        test = pairwise2.align.localms(tmplt, neat, match=2, mismatch=-1, open=-.5, extend=-.1)
        alignment = format_alignment(*test[0]).split()
        makeCigar(alignment)

def makeCigar(alignment):
    """
    Iterate through the alignment and counts the number of matches(including mismatches), deletions, and insertions. Formatted into a cigar string.
    Calls inserCigar to input the cigar string to its respective location in the samfile.
    param: alignment: the alignment of two sequences, index: keeps track of which sequences are being aligned
    """
    templateSeq = alignment[1]
    NEATSeq = alignment[-2]
    cigCount = ''
    currChar = ''
    cigString = ''
    for char in range(0, len(NEATSeq)):
        if templateSeq[char] == '-': #insertion
            if currChar == 'I': #more insertions
                cigCount = cigCount + 1
            else:# new insertion    
                cigString = cigString + str(cigCount) + currChar
                currChar = 'I'
                cigCount = 1

        elif NEATSeq[char] == '-': #deletion
            if currChar == 'D': #more deletions
                cigCount = cigCount + 1
            else:# new deletion
                cigString = cigString + str(cigCount) + currChar
                currChar = 'D'
                cigCount = 1

        else: #match
            if currChar == 'M': #more matches
                cigCount = cigCount + 1
            else:# new match
                cigString = cigString + str(cigCount) + currChar
                currChar = 'M'
                cigCount = 1

    cigString = cigString + str(cigCount) + currChar    
    cigars.append(cigString)


def insertCigar():#made for a file input
    """
    Replaces the template sequence with a cigar string in a new file called "output.tsam"
    """
    i = 0
    try:
        with open("output.tsam", 'w') as output:
            with open (filename , 'r') as samFile:
                for line in samFile:
                    templateSeq = line.split('\t')[5]
                    output.write(line.replace(templateSeq, cigars[i], 1))
                    i = i+1
    except (FileNotFoundError, IOError):
        print("Wrong file or file path")
    finally:
        samFile.close()


filename = sys.argv[1] #Takes file as command line input (testing purposes)
#This can be changed to take a file from gen_reads.py
try:
    with open (filename , 'r') as samFile:
        for line in samFile:
            templateSeq = line.split('\t')[5]
            templates.append(Seq(templateSeq))

            NEATSeq = line.split('\t')[9]
            NEATs.append(Seq(NEATSeq))
    alginSeqs()        
except (FileNotFoundError, IOError):
    print("Wrong file or file path")
finally:
    samFile.close()
insertCigar()
