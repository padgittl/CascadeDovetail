import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq


###############
# SUBROUTINES #
###############



def readFileList(fileList):
    frameshiftDict = {}
    uniqDict = {}
    frameshiftPairs = open('pairs_with_frameshift.txt','w')
    with open(fileList,'r') as F:
        for line in F:
            fileName = line.strip()
            # print fileName
            frameshiftDict = readFasta(fileName,frameshiftDict)
            for pairID in frameshiftDict:
                if pairID not in uniqDict:
                    uniqDict[pairID] = 1
                    frameshiftPairs.write("%s\n" % (pairID))



def readFasta(fastaFile,frameshiftDict):
    fileID = fastaFile
    # print fileID
    for record in SeqIO.parse(fastaFile,'fasta'):
        record.seq = record.seq.upper()
        if '!' in record.seq or 'N' in record.seq:
            # if '!' in record.seq:
            #print record.id
            if fileID not in frameshiftDict:
                frameshiftDict[fileID] = 1
    return(frameshiftDict)

    
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <fasta file list>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

fileList = sys.argv[1]

readFileList(fileList)


