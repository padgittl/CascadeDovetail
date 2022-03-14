import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq


###############
# SUBROUTINES #
###############

def readProteinFasta(proteinFasta):
    seqDict = {}
    for record in SeqIO.parse(proteinFasta,"fasta"):    
        if record.id not in seqDict:
            seqDict[record.id] = record
    return(seqDict)
        

# OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffoldID,geneID,newGeneID,mRNAStart,mRNAStop,cdsStart,cdsStop))
def readMapFile(mapFile,seqDict):
    mapList = []
    with open(mapFile,'r') as F:
        for line in F:
            scaffoldID,oldGeneID,newGeneID,mRNAStart,mRNAStop,cdsStart,cdsStop = line.strip().split('\t')
            if newGeneID in seqDict:
                mapList.append((scaffoldID,newGeneID,int(cdsStart),int(cdsStop)))
    return(mapList)


def createFilesForMCScanX(mapList):
    OUT = open('hopGenes.gff','w')
    # hl1001
    for scaffoldID,geneID,cdsStart,cdsStop in mapList:
        # sp# gene starting_position ending_position
        # sp#, sp is the two-letter short name for the species; # is the chromosome number.
        # species is hl 
        getScaffoldNumber = re.search('Scaffold_(.+)',scaffoldID)
        scaffoldNumber = getScaffoldNumber.group(1)
        sp = 'hl' + scaffoldNumber
        OUT.write("%s\t%s\t%s\t%s\n" % (sp,geneID,cdsStart,cdsStop))


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <gene map file> <protein fasta>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

mapFile = sys.argv[1]
proteinFasta = sys.argv[2]

seqDict = readProteinFasta(proteinFasta)
mapList = readMapFile(mapFile,seqDict)
createFilesForMCScanX(mapList)
