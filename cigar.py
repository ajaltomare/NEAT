"""
Input:File (SAM file with template)
Output:File (SAM file with Cigar string instead of template)

Functions:
    file handeling
        Parse Template string and NEAT string
    ALign using BLAST
    create CIGAR string using alignment, and insert in SAM file
"""
#import pysam
#import io
import sys
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbimakeblastdbCommandline
#import subprocess
import io
import pathlib

from Bio.Seq import Seq 
from Bio import pairwise2


# WorkDirectory = pathlib.Path().resolve() #where the sequence files will be saved
# #Can be made to take in path from gen_reads

templates = [] #list to olds files containing sequences from the template
NEATs = [] #list to hold files containing sequences from the NEAT reads
alignments = []#list to holds alignment results
# fileCntr = 1

# def makeFastqFiles(textInput, name, counter):
#     """
#     Writes sequence into FASTSQ file and adds it to a list
#     param: textInput: sequnece
#     retun: fastq file
#     """
#     # f = io.StringIO(">seq" "\n" + textInput + "\n")
#     f = open(name + str(counter) + ".fastq", 'w+')
#     f.write(">" + name + str(counter) + "\n" + textInput + "\n")
#     f.seek(0)
#     return f

def BLAST(counter):
    """
    This will take in two files, each containg a single sequence and run them through the blastn command and then insert the alignment into alignments[]
    param: counter, keeps track of which pair of template and NEAT sequence to BLAST
    """

    db = NcbimakeblastdbCommandline(dbtype="nucl", input_file = str(WorkDirectory) + "/template" + str(counter) + ".fastq")
    db()
    
    alignmentCommand = NcbiblastnCommandline(query = str(WorkDirectory) + "/template" + str(counter) + ".fastq", 
                    subject = str(WorkDirectory) + "/NEAT" + str(counter) + ".fastq", outfmt = 17)
    #print(alignmentCommand) --> blastn -outfmt 1 -query seq1 -subject seq2
    # subprocess.run(alignmentCommand, shell=True)
    alignments.append(alignmentCommand())

def alginSeqs():
    """
    
    """
    NeatIndex = 0
    for seq in templates:
        alignments.append(pairwise2.align.localms(seq, NEATs[NeatIndex], match=2, mismatch=-1, open=-.5, extend=-.1))
        

def insertCigar(f, alignments):
    """
    Replaces the template sequence with a cigar string
    param: f: sam file, alignments: list with the alignments from 
    return: sam file
    """


#make into 2 functions
filename = sys.argv[1] #Takes file as command line input (testing purposes)
#This can be changed to take a file from gen_reads.py
try:
    with open (filename , 'r') as samFile:
        for line in samFile:
            templateSeq = line.split('\t')[5]
            #templates.append(makeFastqFiles(templateSeq, "template",fileCntr))
            templates.append(Seq(templateSeq))

            NEATSeq = line.split('\t')[9]
            #NEATs.append(makeFastqFiles(NEATSeq, "NEAT", fileCntr))
            NEATs.append(Seq(NEATSeq))

            # BLAST(fileCntr)
            # fileCntr += 1

except (FileNotFoundError, IOError):
    print("Wrong file or file path")
finally:
    samFile.close()    


#path to current working directory

#f = pysam.AlignmentFile(filename, "rb")
    #ValueError: file does not contain alignment data (temp_1.sam)


# Open the sam file
#     f = open(filename, 'rb')
#     i = 0
#     for line in f:
#         line = line.split()
#         i+=1
#         print(i)
#         if i == 5:
#             templates.append(line)
#         elif i==9:
#             SEQs.append(line)