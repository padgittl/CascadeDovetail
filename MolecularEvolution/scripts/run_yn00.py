import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

#fileName = os.path.basename(fullPath)
#filePrefix,fileExt = os.path.splitext(fileName)
#baseName,extra = os.path.splitext(filePrefix)


def readFileList(fileList):
    with open(fileList,'r') as FL:
        for line in FL:
            fileName = line.strip()
            #print fileName
            commands(fileName)


def commands(ctlFile):
    fullPath = ctlFile.strip()
    fileName = os.path.basename(fullPath)
    #print fullPath
    # HUMLU_CAS0000020.t1_vs_HUMLU_CAS0000407.t1_yn00.ctl
    fileNamePrefix,fileNameSuffix = fileName.split('.yn00.ctl')
    outputFile = fileName.replace('ctl','txt')
    errFileName = fileNamePrefix + ".err"
    os.system("mv " + fullPath + " yn00.ctl")
    # ~/bin/paml4.9j/bin/yn00 > err
    os.system("~/bin/paml4.9j/bin/yn00 > " + errFileName)
    os.system("mv " + outputFile + " " + fileNamePrefix + "/.")
    os.system("mv 2YN.dN " + fileNamePrefix + "/.")
    os.system("mv 2YN.dS " + fileNamePrefix + "/.")
    os.system("mv 2YN.t " + fileNamePrefix + "/.")
    os.system("mv " + errFileName + " " + fileNamePrefix + "/.")
    os.system("mv rst " + fileNamePrefix + "/.")
    os.system("mv rst1 " + fileNamePrefix + "/.")
    os.system("mv rub " + fileNamePrefix + "/.")
    os.system("mv yn00.ctl " + fileNamePrefix + "/.")

         
############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <file list> "
if len(sys.argv) != 2:
    print usage
    sys.exit()

fileList = sys.argv[1] 

readFileList(fileList)



