import sys, re, os

################
# SUBROUTINES ##
################

'''
Scaffold_9      94574   
type    count   percentage
DNA     0       0.0     
LTR     10918   11.544  
LINE    0       0.0     
SINE    0       0.0     
Satellite       0       0.0     
Simple  2496    2.639   
Mobile Element  116     0.123   
rRNA    0       0.0     
Other   0       0.0     
Total Repeat Content    13530   14.306  
Non-repeat      81044   85.694 
'''

def readFileList(fileList,assemblyLen):

    dnaCount = 0
    ltrCount = 0
    lineCount = 0
    sineCount = 0
    satelliteCount = 0
    simpleCount = 0
    mobileElementCount = 0
    rRNACount = 0
    otherCount = 0
    totalRepeatCount = 0
    nonRepeatCount = 0

    with open(fileList,'r') as FL:
        for line in FL:
            fileName = line.strip()
            dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount = readRepeatFile(fileName,dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount)

    percDNA = round(float(dnaCount) / assemblyLen * 100,3)
    percLTR = round(float(ltrCount) / assemblyLen * 100,3)
    percLINE = round(float(lineCount) / assemblyLen * 100,3)
    percSINE = round(float(sineCount) / assemblyLen * 100,3)
    percSatellite = round(float(satelliteCount) / assemblyLen * 100,3)
    percSimple = round(float(simpleCount) / assemblyLen * 100,3)
    percMobileElement = round(float(mobileElementCount) / assemblyLen * 100,3)
    percrRNA = round(float(rRNACount) / assemblyLen * 100,3)
    percOTher = round(float(otherCount) / assemblyLen * 100,3)

    totalRepeatPercentage = percDNA + percLTR + percLINE + percSINE + percSatellite + percSimple + percMobileElement + percrRNA + percOTher
    nonRepeatPercentage = 100 - totalRepeatPercentage

    print("Assembly Length:\t%s\t" % (assemblyLen))
    print("type\tcount\tpercentage")
    print("DNA\t%s\t%s\t" % (dnaCount,percDNA))
    print("LTR\t%s\t%s\t" % (ltrCount,percLTR))
    print("LINE\t%s\t%s\t" % (lineCount,percLINE))
    print("SINE\t%s\t%s\t" % (sineCount,percSINE))
    print("Satellite\t%s\t%s\t" % (satelliteCount,percSatellite))
    print("Simple\t%s\t%s\t" % (simpleCount,percSimple))
    print("Mobile Element\t%s\t%s\t" % (mobileElementCount,percMobileElement))
    print("rRNA\t%s\t%s\t" % (rRNACount,percrRNA))
    print("Other\t%s\t%s\t" % (otherCount,percOTher))
    print("Total Repeat Content\t%s\t%s\t" % (totalRepeatCount,totalRepeatPercentage))
    print("Non-Repeat\t%s\t%s\t" % (nonRepeatCount,nonRepeatPercentage))
    

def readRepeatFile(repeatFile,dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount):
    with open(repeatFile,'r') as RF:
        for line in RF:
            if 'Scaffold' in line:
                scaffoldID,scaffoldLength = line.strip().split('\t')
            if 'DNA' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                dnaCount += repCount
            if 'LTR' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                ltrCount += repCount
            if 'LINE' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                lineCount += repCount
            if 'SINE' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                sineCount += repCount
            if 'Satellite' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                satelliteCount += repCount
            if 'Simple' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                simpleCount += repCount
            if 'Mobile' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                mobileElementCount += repCount
            if 'rRNA' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                rRNACount += repCount
            if 'Other' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                otherCount += repCount
            if 'Total Repeat Content' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                totalRepeatCount += repCount
            if 'Non-repeat' in line:
                repType,repCount,repPercentage = line.strip().split('\t')
                repCount = int(repCount)
                nonRepeatCount += repCount
    return(dnaCount,ltrCount,lineCount,sineCount,satelliteCount,simpleCount,mobileElementCount,rRNACount,otherCount,totalRepeatCount,nonRepeatCount)


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <file list> <assembly length>\n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

fileList = sys.argv[1]
assemblyLen = sys.argv[2]

assemblyLen = int(assemblyLen)

readFileList(fileList,assemblyLen)
