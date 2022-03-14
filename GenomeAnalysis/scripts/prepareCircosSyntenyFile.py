import sys, re, os, math


def readGFF(gffFile):
    coordDict = {}
    with open(gffFile,'r') as GFF:
        for line in GFF:
            # hl1001  HUMLU_CAS0065845.t1     30674   31037
            scaffoldID,geneID,geneStart,geneStop = line.strip().split('\t')
            if 'hl' in scaffoldID:
                scaffoldNumber = scaffoldID.replace('hl','').strip()
                scaffoldID = 'Scaffold_' + scaffoldNumber
                #print scaffoldID
                if geneID not in coordDict:
                    coordDict[geneID] = (scaffoldID,int(geneStart),int(geneStop))
                    # print scaffoldID
            else:
                print "problem parsing gff file"
                sys.exit()
    return(coordDict)


def readScaffoldLengthsFile(scaffoldLengthsFile):
    scaffoldNumber = 10
    scaffoldLengthDict = {}
    scaffoldLenList = []
    with open(scaffoldLengthsFile,'r') as LEN:
        for line in LEN:
            scaffoldID,scaffoldLen = line.strip().split()
            scaffoldLenList.append((scaffoldID,int(scaffoldLen)))
    scaffoldLenList.sort(key=lambda x:x[1], reverse=True)
    longestScaffoldList = scaffoldLenList[0:scaffoldNumber]
    for scaffoldID,scaffoldLen in longestScaffoldList:
        if scaffoldID not in scaffoldLengthDict:
            scaffoldLengthDict[scaffoldID] = scaffoldLen
    return(scaffoldLengthDict,longestScaffoldList)


def readGeneMapFile(geneMapFile,scaffoldLengthDict):
    geneMapDict = {}
    with open(geneMapFile,'r') as F:
        for line in F:
            scaffoldID,oldGeneID,newGeneID,mRNAStart,mRNAStop,cdsStart,cdsStop = line.strip().split('\t')
            if scaffoldID in scaffoldLengthDict:
                if newGeneID not in geneMapDict:
                    geneMapDict[newGeneID] = scaffoldID
    return(geneMapDict)


def parseCollinearityFile(collinearityFile,coordDict,geneMapDict):
    blockDict1 = {}
    blockDict2 = {}
    strandDict = {}
    with open(collinearityFile,'r') as F:
        for line in F:
            # Alignment
            if not line.startswith('#'):
                # print line
                # 0-  0:        XP_030485978.1  XP_030486287.1        0
                fullNumberBlock,geneID1,geneID2,eValue = line.strip().split('\t')
                numberBlock,extra = fullNumberBlock.split('-')
                numberBlock = int(numberBlock)
                IDs = sorted([geneID1,geneID2])
                fileID = IDs[0] + "_vs_" + IDs[1]
                if geneID1 in geneMapDict and geneID2 in geneMapDict:
                    scaffoldID1,geneStart1,geneStop1 = coordDict[IDs[0]]
                    scaffoldID2,geneStart2,geneStop2 = coordDict[IDs[1]]
                    #print scaffoldID1,scaffoldID2
                    if 'Scaffold' in scaffoldID1 and 'Scaffold' in scaffoldID2:
                        if numberBlock not in blockDict1:
                            blockDict1[numberBlock] = {}
                        if numberBlock not in blockDict2:
                            blockDict2[numberBlock] = {}
                        if scaffoldID1 not in blockDict1[numberBlock]:
                            blockDict1[numberBlock][scaffoldID1] = []
                        if scaffoldID2 not in blockDict2[numberBlock]:
                            blockDict2[numberBlock][scaffoldID2] = []
                        blockDict1[numberBlock][scaffoldID1].append((geneStart1,geneStop1))
                        blockDict2[numberBlock][scaffoldID2].append((geneStart2,geneStop2))
                    else:
                        print "problem parsing scaffold ID"
                        sys.exit()
            if line.startswith('## Alignment'):
                # ## Alignment 0: score=437.0 e_value=7e-20 N=10 csNC_044370.1&csNC_044370.1 plus
                line = line.strip().split()
                strand = line[7]
                blockNumber = line[2].replace(':','').strip()
                blockNumber = int(blockNumber)
                if blockNumber not in strandDict:
                    strandDict[blockNumber] = strand
                    # print blockNumber,strand
    blockCoords1 = {}
    blockCoords2 = {}
    for blockID1 in blockDict1:
        if blockID1 not in blockCoords1:
            blockCoords1[blockID1] = {}
        for scaffoldID1 in blockDict1[blockID1]:
            #print blockID,scaffoldID
            blockDict1[blockID1][scaffoldID1].sort(key=lambda x: x[0], reverse=False)
            blockStart1 = blockDict1[blockID1][scaffoldID1][0][0] 
            blockDict1[blockID1][scaffoldID1].sort(key=lambda x: x[1], reverse=True)
            blockStop1 = blockDict1[blockID1][scaffoldID1][0][1]
            if scaffoldID1 not in blockCoords1[blockID1]:
                blockCoords1[blockID1][scaffoldID1] = (blockStart1,blockStop1)

    for blockID2 in blockDict2:
        if blockID2 not in blockCoords2:
            blockCoords2[blockID2] = {}
        for scaffoldID2 in blockDict2[blockID2]:
            blockDict2[blockID2][scaffoldID2].sort(key=lambda x: x[0], reverse=False)
            blockStart2 = blockDict2[blockID2][scaffoldID2][0][0]
            blockDict2[blockID2][scaffoldID2].sort(key=lambda x: x[1], reverse=True)
            blockStop2 = blockDict2[blockID2][scaffoldID2][0][1]
            if scaffoldID2 not in blockCoords2[blockID2]:
                blockCoords2[blockID2][scaffoldID2] = blockStart2,blockStop2
    return(blockCoords1,blockCoords2,strandDict)

'''
segdup01284 hs1 102024400 102025440
segdup01284 hs3 111883743 111884767
segdup02286 hs1 152617218 152618252
segdup02286 hs3 111883745 111884756
'''
def createFile(blockCoords1,blockCoords2,outName):
    OUT = open(outName,'w')
    for blockID in blockCoords1:
        for scaffoldID1 in blockCoords1[blockID]:
            for scaffoldID2 in blockCoords2[blockID]:
                blockStart1,blockStop1 = blockCoords1[blockID][scaffoldID1]
                blockStart2,blockStop2 = blockCoords2[blockID][scaffoldID2]
                newBlockID = 'block' + str(blockID)
                OUT.write("%s\t%s\t%s\t%s\n" % (newBlockID,scaffoldID1,blockStart1,blockStop1))
                OUT.write("%s\t%s\t%s\t%s\n" % (newBlockID,scaffoldID2,blockStart2,blockStop2))

############
### MAIN ###
############


usage = "Usage: " + sys.argv[0] + " <scaffold lengths file> <gene-scaffold map file> <collinearity file> <gff file> <out name>\n"
if len(sys.argv) != 6:
    print usage
    sys.exit()

scaffoldLengthsFile = sys.argv[1]
geneMapFile = sys.argv[2]
collinearityFile = sys.argv[3]
gffFile = sys.argv[4]
outName = sys.argv[5]

coordDict = readGFF(gffFile)
scaffoldLengthDict,longestScaffoldList = readScaffoldLengthsFile(scaffoldLengthsFile)
geneMapDict = readGeneMapFile(geneMapFile,scaffoldLengthDict)
blockCoords1,blockCoords2,strandDict = parseCollinearityFile(collinearityFile,coordDict,geneMapDict)

createFile(blockCoords1,blockCoords2,outName)
