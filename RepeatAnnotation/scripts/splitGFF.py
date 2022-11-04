import sys, re, os

################
# SUBROUTINES ##
################

def readGFF(gffFile):
    scaffoldDict = {}
    with open(gffFile,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                scaffoldID,source,feature,start,end,score,strand,frame,attribute  = line.strip().split("\t")
                if scaffoldID not in scaffoldDict:
                    scaffoldDict[scaffoldID] = []
                scaffoldDict[scaffoldID].append((source,feature,start,end,score,strand,frame,attribute))
    return(scaffoldDict)

def splitGFF(scaffoldDict):
    for scaffoldID in scaffoldDict:
        OUT = open(scaffoldID + ".gff",'w')
        for source,feature,start,end,score,strand,frame,attribute in scaffoldDict[scaffoldID]:
            OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffoldID,source,feature,start,end,score,strand,frame,attribute))

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <gff file>\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

gffFile = sys.argv[1]

scaffoldDict = readGFF(gffFile)
splitGFF(scaffoldDict)
