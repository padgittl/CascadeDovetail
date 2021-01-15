import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq


###############
# SUBROUTINES #
###############


def readESTFastaFile(estFastaFile):
    estDict = {}
    for record in SeqIO.parse(estFastaFile,"fasta"):
        accessionID = record.id
        accessionDesc = record.description
        accessionSeq = record.seq
        estDict[accessionID] = record
        #print record
    return estDict


def readBlastnOutputFileList(blastnOutputFileList,estDict):
    with open(blastnOutputFileList,'r') as BFL:
        for line in BFL:
            # each line is the blast output file for a different scaffold
            fullPath = line.strip() 
            fileName = os.path.basename(fullPath)  
            baseName,fileExt = os.path.splitext(fileName)
            getScaffoldID = re.search('(.+)_vs',baseName)
            scaffoldID = getScaffoldID.group(1)
            hitList = getESTHits(fullPath,estFastaFile,estDict)
            #print hitList
            recordList = []
            for contigName in hitList:
                estIDs = hitList[contigName]
                for estID in estIDs:
                    if estID in estDict:
                        recordList.append(estDict[estID])
                        SeqIO.write(recordList, contigName + "_vs_trichomeESTs.fasta", "fasta")
                    else:
                        print "error!"
                        sys.exit()


def getESTHits(pathToFile,estFastaFile,estDict):
    hitDict = {}
    for line in open(pathToFile,'r'):
        contigID,estID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore = line.strip().split("\t")
        if contigID not in hitDict:
            hitDict[contigID] = []
        if estID not in hitDict[contigID]:
            hitDict[contigID].append(estID)
    return hitDict

            
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <est fasta file> <blastn output file> "
if len(sys.argv) != 3:
    print usage
    sys.exit()

estFastaFile = sys.argv[1]
blastnOutputFileList = sys.argv[2]

estDict = readESTFastaFile(estFastaFile)
readBlastnOutputFileList(blastnOutputFileList,estDict)



