import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq


###############
# SUBROUTINES #
###############


def readFileList(fileList):
    with open(fileList,'r') as F:
        for line in F:
            fileName = line.strip()
            modifyFastaDefline(fileName)


def modifyFastaDefline(fastaFile):
    recordDict = {}
    recordList = []
    fileID = fastaFile
    getOGG = re.search('(OG.+)_aligned',fileID)
    oggID = getOGG.group(1)
    #print oggID
    for record in SeqIO.parse(fastaFile,'fasta'):
        #print record.id
        idStuff = record.id.split('_')
        if len(idStuff) == 4:
            genus,species,geneIDPrefix,geneIDSuffix = idStuff
            speciesID = genus + "_" + species
            #print speciesID
        else:
            genus,species,geneID = idStuff
            speciesID = genus + "_" + species
            #print speciesID
        record.id = record.id.replace(record.id,speciesID)
        record.name = record.name.replace(record.name,'')
        record.description = record.description.replace(record.description,'')
        #print record.id
        if record.id not in recordDict:
            recordDict[record.id] = record
            recordList.append(record)
    SeqIO.write(recordList, oggID + ".fasta", "fasta")
        
    
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <fasta file list> "
if len(sys.argv) != 2:
    print usage
    sys.exit()

fileList = sys.argv[1]

readFileList(fileList)


