import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############


def readUniprotFastaFile(uniprotFastaFile):
    uniprotDict = {}
    for record in SeqIO.parse(uniprotFastaFile,"fasta"):
        # sp|F4HVA6|TAF6B_ARATH
        getAccessionID = re.search('.+\|(.+)',record.id)
        accessionID = getAccessionID.group(1)
        record.id = record.id.replace(record.id,accessionID)
        record.name = record.name.replace(record.name,'')
        record.description = record.description.replace(record.description,'')
        accessionID = record.id
        #print record.id,record.name,record.description
        uniprotDict[accessionID] = record
        #print record.id
    return(uniprotDict)


def readBlastxOutputFileList(blastxOutputFileList,uniprotDict):
    with open(blastxOutputFileList,'r') as BFL:
        for line in BFL:
            # each line in the file list represents a different scaffold
            fullPath = line.strip() 
            fileName = os.path.basename(fullPath)  
            baseName,fileExt = os.path.splitext(fileName)
            getScaffoldID = re.search('(.+)_vs',baseName)
            scaffoldID = getScaffoldID.group(1)
            hitList = getUniprotHits(fullPath,uniprotDict)
            #print hitList
            # new recordList initialized for each scaffold
            recordList = []
            for contigName in hitList:
                uniprotIDs = hitList[contigName]
                for uniprotID in uniprotIDs:
                    if uniprotID in uniprotDict:
                        recordList.append(uniprotDict[uniprotID])
                        SeqIO.write(recordList, contigName + "_vs_refSeqCSativa.fasta", "fasta")
                    else:
                        #print "error!"
                        print contigName,uniprotID
                        sys.exit()


def getUniprotHits(pathToFile,uniprotDict):
    hitDict = {}
    for line in open(pathToFile,'r'):
        contigID,uniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore = line.strip().split("\t")
        if contigID not in hitDict:
            hitDict[contigID] = []
        if uniprotID not in hitDict[contigID]:
            hitDict[contigID].append(uniprotID)
            #print contigID,uniprotID
    return(hitDict)

            
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <fasta file> <blastx output file> "
if len(sys.argv) != 3:
    print usage
    sys.exit()

uniprotFastaFile = sys.argv[1]
blastxOutputFileList = sys.argv[2]

uniprotDict = readUniprotFastaFile(uniprotFastaFile)
readBlastxOutputFileList(blastxOutputFileList,uniprotDict)




