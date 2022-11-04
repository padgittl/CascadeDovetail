import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

def writeFasta(fastaFile):
    recordDict = {}
    for record in SeqIO.parse(fastaFile,'fasta'):
        # >Scaffold_1;HRSCAF=21
        if record.id not in recordDict:
            recordDict[record.id] = record
            SeqIO.write(recordDict[record.id], record.id + ".fasta", "fasta")
            
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <fasta file> "
if len(sys.argv) != 2:
    print usage
    sys.exit()

fastaFile = sys.argv[1]

writeFasta(fastaFile)


