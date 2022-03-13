import sys, re, os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import numpy as np

###############
# SUBROUTINES #
###############


def readOrthogroupIDFile(orthogroupIDFile):
    selectedOrthogroupIDDict = {}
    with open(orthogroupIDFile,'r') as O:
        for line in O:
            orthogroupID = line.strip()
            selectedOrthogroupIDDict[orthogroupID] = 1
    return(selectedOrthogroupIDDict)


def readTSV(orthogroupTSV,selectedOrthogroupIDDict):
    orthogroupDict = {}
    with open(orthogroupTSV,'r') as F:
        for line in F:
            if not line.startswith('Orthogroup'):
                orthoGroupInfo = line.strip().split('\t')
                orthogroupID = orthoGroupInfo[0]
                if orthogroupID in selectedOrthogroupIDDict:
                    # cannabis
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    cannabisGeneID = orthoGroupInfo[1:][0]
                    fullCannabisGeneID = "cannabis_sativa_" + cannabisGeneID
                    orthogroupDict[orthogroupID].append(fullCannabisGeneID)

                    # hops
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    hopGeneID = orthoGroupInfo[1:][1]
                    fullHopGeneID = "humulus_lupulus_" + hopGeneID
                    orthogroupDict[orthogroupID].append(fullHopGeneID)
                        
                    # mulberry
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    mulberryGeneID = orthoGroupInfo[1:][2]
                    fullMulberryGeneID = "morus_notabilis_" + mulberryGeneID
                    orthogroupDict[orthogroupID].append(fullMulberryGeneID)

                    # Parasponia andersonii
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    parasponiaGeneID = orthoGroupInfo[1:][3]
                    fullParasponiaGeneID = "parasponia_andersonii_" + parasponiaGeneID
                    orthogroupDict[orthogroupID].append(fullParasponiaGeneID)

                    # prunus persica
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    prunusGeneID = orthoGroupInfo[1:][4]
                    fullPrunusGeneID = "prunus_persica_" + prunusGeneID
                    orthogroupDict[orthogroupID].append(fullPrunusGeneID)

                    # Trema orientale
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    tremaGeneID = orthoGroupInfo[1:][5]
                    fullTremaGeneID = "trema_orientale_" + tremaGeneID
                    orthogroupDict[orthogroupID].append(fullTremaGeneID)
                        
                    # Vitis vinifera
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    vitisGeneID = orthoGroupInfo[1:][6]
                    fullVitisGeneID = "vitis_vinifera_" + vitisGeneID
                    orthogroupDict[orthogroupID].append(fullVitisGeneID)

                    # ziziphus jujuba
                    if orthogroupID not in orthogroupDict:
                        orthogroupDict[orthogroupID] = []
                    ziziphusGeneID = orthoGroupInfo[1:][7]
                    fullZiziphusGeneID = "ziziphus_jujuba_" + ziziphusGeneID
                    orthogroupDict[orthogroupID].append(fullZiziphusGeneID)
    return(orthogroupDict)



def readFasta(fastaFile,speciesID):
    recordDict = {}
    recordList = []
    for record in SeqIO.parse(fastaFile,"fasta"):
        fullGeneID = speciesID + "_" + record.id
        record.id = record.id.replace(record.id,fullGeneID)
        record.name = record.name.replace(record.name,'')
        record.description = record.description.replace(record.description,'')
        if record.id not in recordDict:
            recordDict[record.id] = record
    return(recordDict)


def createOrthogroupFastas(orthogroupDict,cannabisRecordDict,humulusRecordDict,morusRecordDict,parasponiaRecordDict,prunusRecordDict,tremaRecordDict,vitisRecordDict,ziziphusRecordDict):
    # orthogroupDict[orthogroupID].append(fullZiziphusGeneID)
    for orthogroupID in orthogroupDict:
        # cannabis,hop,morus,parasponia,prunus,trema,vitis,ziziphus,total
        fullOrthogroupRecordDict = {}
        fullOrthogroupRecordList = []
        outName = orthogroupID + ".fasta"
        for geneID in orthogroupDict[orthogroupID]:
            #print speciesID,orthogroupID,geneID
            if geneID in cannabisRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(cannabisRecordDict[geneID])
                
            if geneID in humulusRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(humulusRecordDict[geneID])

            if geneID in morusRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(morusRecordDict[geneID])

            if geneID in parasponiaRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(parasponiaRecordDict[geneID])

            if geneID in prunusRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(prunusRecordDict[geneID])

            if geneID in tremaRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(tremaRecordDict[geneID])

            if geneID in vitisRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(vitisRecordDict[geneID])

            if geneID in ziziphusRecordDict:
                if geneID not in fullOrthogroupRecordDict:
                    fullOrthogroupRecordDict[geneID] = 1
                    fullOrthogroupRecordList.append(ziziphusRecordDict[geneID])
                    
        SeqIO.write(fullOrthogroupRecordList, outName, "fasta")



############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <orthogroup tsv file> <orthogroup ID file> <cannabis fasta> <hop fasta> <morus fasta> <parasponia fasta> <prunus fasta> <trema fasta> <vitis fasta> <ziziphus fasta>"
if len(sys.argv) != 11:
    print usage
    sys.exit()

orthogroupTSV = sys.argv[1]
orthogroupIDFile = sys.argv[2]
cannabisFasta = sys.argv[3]
humulusFasta = sys.argv[4]
morusFasta = sys.argv[5]
parasponiaFasta = sys.argv[6]
prunusFasta = sys.argv[7]
tremaFasta = sys.argv[8]
vitisFasta = sys.argv[9]
ziziphusFasta = sys.argv[10]

selectedOrthogroupIDDict = readOrthogroupIDFile(orthogroupIDFile)
orthogroupDict = readTSV(orthogroupTSV,selectedOrthogroupIDDict)

cannabisRecordDict = readFasta(cannabisFasta,"cannabis_sativa")
humulusRecordDict = readFasta(humulusFasta,"humulus_lupulus")
morusRecordDict = readFasta(morusFasta,"morus_notabilis")
parasponiaRecordDict = readFasta(parasponiaFasta,"parasponia_andersonii")
prunusRecordDict = readFasta(prunusFasta,"prunus_persica")
tremaRecordDict = readFasta(tremaFasta,"trema_orientale")
vitisRecordDict = readFasta(vitisFasta,"vitis_vinifera")
ziziphusRecordDict = readFasta(ziziphusFasta,"ziziphus_jujuba")

createOrthogroupFastas(orthogroupDict,cannabisRecordDict,humulusRecordDict,morusRecordDict,parasponiaRecordDict,prunusRecordDict,tremaRecordDict,vitisRecordDict,ziziphusRecordDict)

