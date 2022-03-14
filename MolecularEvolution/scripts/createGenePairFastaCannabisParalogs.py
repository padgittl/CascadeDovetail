import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np


###############
# SUBROUTINES #
###############


'''
## Alignment 3: score=284.0 e_value=9.3e-11 N=7 csNC_044370.1&hl1531 plus
  3-  0:        XP_030500325.1  HUMLU_CAS0006395.t1           0
'''
def parseCollinearityFile(collinearityFile):
    genePairDict = {}
    genePairList = []
    with open(collinearityFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                numberBlock,geneID1,geneID2,eValue = line.strip().split('\t')
                # IDs = sorted([ID1,ID2], key=lambda x:contigLengths[x])
                # only collinear paralogs from cannabis
                if "XP_" in geneID1 and "XP_" in geneID2:
                    IDs = sorted([geneID1,geneID2])
                    if (IDs[0],IDs[1]) not in genePairDict:
                        genePairDict[(IDs[0],IDs[1])] = True
                        genePairList.append((IDs[0],IDs[1]))
                        #print geneID1,geneID2
                        #print IDs[0],IDs[1]
                    #else:
                        #print IDs[0],IDs[1]
                        #print geneID1,geneID2
    return(genePairList)


def readFasta(fastaFile):
    recordDict = {}
    recordList = []
    for record in SeqIO.parse(fastaFile,"fasta"):
        record.name = record.name.replace(record.name,'')
        record.description = record.description.replace(record.description,'')
        if record.id not in recordDict:
            recordDict[record.id] = record
            recordList.append(record)
    #print len(recordList)
    return(recordDict,recordList)


def createFasta(genePairList,recordDict):
    #print len(genePairList)
    for geneID1,geneID2 in genePairList:
        pairRecordList = []
        if geneID1 in recordDict and geneID2 in recordDict:
            pairRecordList.append((recordDict[geneID1]))
            pairRecordList.append((recordDict[geneID2]))
            #print geneID1,geneID2
            SeqIO.write(pairRecordList,geneID1 + "_vs_" + geneID2 + ".fasta", "fasta")
        else:
            print "problem"
            sys.exit()

 

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <collinearity file> <fasta file>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

collinearityFile = sys.argv[1]
fastaFile = sys.argv[2]

genePairList = parseCollinearityFile(collinearityFile)
recordDict,recordList = readFasta(fastaFile)
createFasta(genePairList,recordDict)
