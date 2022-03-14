import sys, re, os
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import numpy as np


###############
# SUBROUTINES #
###############


# Orthogroup      cannabis_sativa humulus_lupulus morus_notabilis parasponia_andersonii   prunus_per
# OG0000000       0       477     0       0       0       0       0       0       477


def readTSV(orthogroupGeneCountFile):
    orthogroupDict = {}
    with open(orthogroupGeneCountFile,'r') as F:
        for line in F:
            if not line.startswith('Orthogroup'):
                # Orthogroup      cannabis_sativa humulus_lupulus morus_notabilis parasponia_andersonii   prunus_persica  trema_orientale vitis_vinifera  ziziphus_jujuba Total
                orthogroupID,cannabis_sativa,humulus_lupulus,morus_notabilis,parasponia_andersonii,prunus_persica,trema_orientale,vitis_vinifera,ziziphus_jujuba,total = line.strip().split('\t')
                if orthogroupID not in orthogroupDict:
                    orthogroupDict[orthogroupID] = (int(cannabis_sativa),int(humulus_lupulus),int(morus_notabilis),int(parasponia_andersonii),int(prunus_persica),int(trema_orientale),int(vitis_vinifera),int(ziziphus_jujuba),int(total))
    return(orthogroupDict)


def writeNewOrthogroupFile(orthogroupDict,maxGeneNumber):
    OUT = open('orthogroupCountFile.cafe.unfiltered.tsv','w')
    OUT.write("Description\tID\tcannabisSativa\thumulusLupulus\tmorusNotabilis\tparasponiaAndersonii\tprunusPersica\ttremaOrientale\tvitisVinifera\tziziphusJujuba\n")
    for orthogroupID in orthogroupDict:
        cannabis_sativa,humulus_lupulus,morus_notabilis,parasponia_andersonii,prunus_persica,trema_orientale,vitis_vinifera,ziziphus_jujuba,total = orthogroupDict[orthogroupID]
        # orthogroup should have at least one gene from every species and fewer than 100 genes for each species
        #if cannabis_sativa > 0 and humulus_lupulus > 0 and morus_notabilis > 0 and parasponia_andersonii > 0 and prunus_persica > 0 and trema_orientale > 0 and vitis_vinifera > 0 and ziziphus_jujuba > 0:
        #if cannabis_sativa < maxGeneNumber and humulus_lupulus < maxGeneNumber and morus_notabilis < maxGeneNumber and parasponia_andersonii < maxGeneNumber and prunus_persica < maxGeneNumber and trema_orientale < maxGeneNumber and vitis_vinifera < maxGeneNumber and ziziphus_jujuba < maxGeneNumber:
        OUT.write("\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (orthogroupID,cannabis_sativa,humulus_lupulus,morus_notabilis,parasponia_andersonii,prunus_persica,trema_orientale,vitis_vinifera,ziziphus_jujuba))


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <orthogroup tsv file> "
if len(sys.argv) != 2:
    print usage
    sys.exit()

orthogroupTSV = sys.argv[1]

# A reference genome for pea provides insight into legume genome evolution 
# https://www.nature.com/articles/s41588-019-0480-1
# Supplementary info, setting max group size to 20 or more
# "The larger the orthogroup size the more likely the orthogroup includes spurious homologs with Ks values spreading throughout the Ks range (eg. groups with 20 or more genes)."
# maxOrthogroupSize = 100
maxGeneNumber = 100

orthogroupDict = readTSV(orthogroupTSV)

writeNewOrthogroupFile(orthogroupDict,maxGeneNumber)


