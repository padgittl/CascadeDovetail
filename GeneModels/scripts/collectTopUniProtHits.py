import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import math

###############
# SUBROUTINES #
###############


def readRepeatFileFromUniProtPlants(repeatFileUniProtPlants):
    repeatsFromUniProtPlants = {}
    with open(repeatFileUniProtPlants,'r') as F:
        for line in F:
            # sp|Q9XEC1|Y4443_ARATH   sp|Q9XEC1|Y4443_ARATH Probable transposase-like protein At4g04430 OS=Arabidopsis thaliana OX=3702 GN=At4g04430 PE=3 SV=1 
            fullUniprotID,fullUniprotDescription = line.strip().split('\t')
            getUniprotID = re.search('\|(.+)\|',fullUniprotID)
            uniprotID = getUniprotID.group(1)
            IDStuff,uniprotDesc = fullUniprotDescription.split(fullUniprotID)
            uniprotDescription = uniprotDesc.strip()
            #print uniprotID,uniprotDescription
            if uniprotID not in repeatsFromUniProtPlants:
                repeatsFromUniProtPlants[uniprotID] = uniprotDescription
    return(repeatsFromUniProtPlants)


def filterTwilightThresholdGenes(alignmentLength):
    # https://www.rostlab.org/papers/1999_twilight/paper.html#ref5
    # "changed new formula"
    # 480 * x**(-0.32*(1+exp(-x/1000)))
    lengthThreshold = 450
    L = alignmentLength
    if L <= lengthThreshold:
        percentIdentityThreshold = 480 * L**(-0.32*(1+math.exp(-float(L)/1000)))
    if L > lengthThreshold:
        percentIdentityThreshold = 19.5
    return(percentIdentityThreshold)


def readBlastpOutputFile_gene_vs_uniprot(pepBlastpOutputFile1,eValueThreshold,queryCoverageThreshold):
    bestHits_gene_vs_uniprot = {}
    bestScore = {}
    with open(pepBlastpOutputFile1,'r') as BF1:
        for line in BF1:
            # 001620F.g33.t1  sp|O80438|MAK3_ARATH    74.444  180     46      0       17      196     11      190     6.62e-100       287     88
            geneID,fullUniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            # storing all gene IDs with a hit in this direction
            pIdent = float(pIdent)
            length = int(length)
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            extractUniprotID = re.search('sp\|(.+)\|',fullUniprotID)
            uniprotID = extractUniprotID.group(1)
            percentIdentityThreshold = filterTwilightThresholdGenes(length)
            if eValue < eValueThreshold:
                if queryCov >= queryCoverageThreshold:
                    if pIdent >= percentIdentityThreshold:
                        #print geneID,uniprotID,length,pIdent,percentIdentityThreshold
                        if geneID in bestHits_gene_vs_uniprot:
                            if bestScore[geneID] < bitScore:
                                bestScore[geneID] = bitScore
                                bestHits_gene_vs_uniprot[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
                        else:
                            bestScore[geneID] = bitScore
                            bestHits_gene_vs_uniprot[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
                    # else:
                    #    print geneID,uniprotID,length,pIdent,percentIdentityThreshold
    return(bestHits_gene_vs_uniprot)



def readBlastpOutputFile_uniprot_vs_gene(pepBlastpOutputFile2,eValueThreshold,queryCoverageThreshold):
    bestHits_uniprot_vs_gene = {}
    bestScore = {}
    with open(pepBlastpOutputFile2,'r') as BF2:
        for line in BF2:
            fullUniprotID,geneID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            pIdent = float(pIdent)
            length = int(length)
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            extractUniprotID = re.search('sp\|(.+)\|',fullUniprotID)
            uniprotID = extractUniprotID.group(1)
            percentIdentityThreshold = filterTwilightThresholdGenes(length)
            if eValue < eValueThreshold:
                if queryCov >= queryCoverageThreshold:
                    if pIdent >= percentIdentityThreshold:
                        if geneID in bestHits_uniprot_vs_gene:
                            if bestScore[geneID] < bitScore:
                                bestScore[geneID] = bitScore
                                bestHits_uniprot_vs_gene[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
                        else:
                            bestScore[geneID] = bitScore
                            bestHits_uniprot_vs_gene[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
    return(bestHits_uniprot_vs_gene)


def getTopHitsBothDirections(bestHits_gene_vs_uniprot,bestHits_uniprot_vs_gene):
    bestHitsBothDirections = {}
    bestScoreBothDirections = {}
    not_in_gene_to_uniprot = open('not_in_gene_to_uniprot.txt','w')
    not_in_uniprot_to_gene = open('not_in_uniprot_to_gene.txt','w')
    for geneID1 in bestHits_gene_vs_uniprot:
        uniprotID1,percIdentity1,eValue1,bitScore1,queryCov1 = bestHits_gene_vs_uniprot[geneID1]
        #print uniprotID1
        if geneID not in bestHits_uniprot_vs_gene:
            not_in_uniprot_to_gene.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID1,uniprotID1,percIdentity1,eValue1,bitScore1,queryCov1))
        if geneID1 in bestHitsBothDirections:
            if bestScoreBothDirections[geneID1] < bitScore1:
                bestScoreBothDirections[geneID1] = bitScore1
                bestHitsBothDirections[geneID1] = (uniprotID1,percIdentity1,eValue1,bitScore1,queryCov1)
        else:
            bestScoreBothDirections[geneID1] = bitScore1
            bestHitsBothDirections[geneID1] = (uniprotID1,percIdentity1,eValue1,bitScore1,queryCov1)
    for geneID2 in bestHits_uniprot_vs_gene:
        uniprotID2,percIdentity2,eValue2,bitScore2,queryCov2 = bestHits_uniprot_vs_gene[geneID2]
        if geneID2 not in bestHits_gene_vs_uniprot:
            not_in_gene_to_uniprot.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID2,uniprotID2,percIdentity2,eValue2,bitScore2,queryCov2))
        #print geneID2,uniprotID2,percIdentity2,eValue2,bitScore2,queryCov2
        if geneID2 in bestHitsBothDirections:
            if bestScoreBothDirections[geneID2] < bitScore2:
                bestScoreBothDirections[geneID2] = bitScore2
                bestHitsBothDirections[geneID2] = (uniprotID2,percIdentity2,eValue2,bitScore2,queryCov2)
        else:
            bestScoreBothDirections[geneID2] = bitScore2
            bestHitsBothDirections[geneID2] = (uniprotID2,percIdentity2,eValue2,bitScore2,queryCov2)
    return(bestHitsBothDirections)


def getDescriptionFromFasta(fastaFile):
    uniprotDescriptionDict = {}
    for record in SeqIO.parse(fastaFile,"fasta"):
        # sp|F4HVA6|TAF6B_ARATH
        getUniprotID = re.search('sp\|(.+)\|',record.id)
        uniprotID = getUniprotID.group(1)
        #print uniprotID
        if uniprotID not in uniprotDescriptionDict:
            uniprotDescriptionDict[uniprotID] = record.description
            #print geneID,record.description
    return(uniprotDescriptionDict)


########
# MAIN #
########


usage = "Usage: " + sys.argv[0] + " <blastp output file 1> <blastp output file 2> <uniprot plant fasta file> <uniprot keyword TE fasta file> <repeat file from uniprot plants> <query coverage threshold>\n"
if len(sys.argv) != 7:
    print usage
    sys.exit()

blastpOutputFile1 = sys.argv[1]
blastpOutputFile2 = sys.argv[2]
uniprotPlantFastaFile = sys.argv[3]
# key word TE = KWTE
uniprotKWTEFastaFile = sys.argv[4]
repeatFileUniProtPlants = sys.argv[5]
queryCovThreshold = sys.argv[6]

eValueThreshold = 1e-3
queryCoverageThreshold = int(queryCovThreshold)
# queryCoverageThreshold = 20

repeatsFromUniProtPlants = readRepeatFileFromUniProtPlants(repeatFileUniProtPlants)
bestHits_gene_vs_uniprot = readBlastpOutputFile_gene_vs_uniprot(blastpOutputFile1,eValueThreshold,queryCoverageThreshold)
OUT_bothDirections = open('bestHitsBothDirections.txt','w')
for geneID in bestHits_gene_vs_uniprot:
    uniprotID,pIdent,eValue,bitScore,queryCov = bestHits_gene_vs_uniprot[geneID]
    OUT_bothDirections.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID,uniprotID,pIdent,eValue,bitScore,queryCov))

bestHits_uniprot_vs_gene = readBlastpOutputFile_uniprot_vs_gene(blastpOutputFile2,eValueThreshold,queryCoverageThreshold)
for geneID in bestHits_uniprot_vs_gene:
    uniprotID,pIdent,eValue,bitScore,queryCov = bestHits_uniprot_vs_gene[geneID]
    OUT_bothDirections.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID,uniprotID,pIdent,eValue,bitScore,queryCov))


uniprotPlantDescriptionDict = getDescriptionFromFasta(uniprotPlantFastaFile)
# these are TE genes from uniprot -- these will potentially include plant TEs that are harder to identify by basic keywords like 'retrotransposon'
uniprotKWTEDescriptionDict = getDescriptionFromFasta(uniprotKWTEFastaFile)

bestHitsBothDirections = getTopHitsBothDirections(bestHits_gene_vs_uniprot,bestHits_uniprot_vs_gene)


OUT_NOREP = open('topHitsRepeatsRemovedQueryCov' + str(queryCoverageThreshold) + '.txt','w')
OUT_NOREP.write("geneID\tuniprotID\tpercentIdentity\teValue\tbitScore\tqueryCoverage\tuniprotDescription\n")
for geneID in bestHitsBothDirections:
    uniprotID,percIdentity,eValue,bitScore,queryCov = bestHitsBothDirections[geneID]
    # require that uniprot hits not be repeat-associated 
    if uniprotID not in repeatsFromUniProtPlants:
        if uniprotID not in uniprotKWTEDescriptionDict:
            if uniprotID in uniprotPlantDescriptionDict:
                uniprotDescription = uniprotPlantDescriptionDict[uniprotID]
                OUT_NOREP.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID,uniprotID,percIdentity,eValue,bitScore,queryCov,uniprotDescription))
            #else:
            #print "problem with uniprot dict"
            #sys.exit()


OUT_REP = open('topHitsRepeatsQueryCov' + str(queryCoverageThreshold) + '.txt','w')
OUT_REP.write("geneID\tuniprotID\tpercentIdentity\teValue\tbitScore\tqueryCoverage\tuniprotDescription\n")
for geneID in bestHitsBothDirections:
    uniprotID,percIdentity,eValue,bitScore,queryCov = bestHitsBothDirections[geneID]
    if uniprotID in repeatsFromUniProtPlants or uniprotID in uniprotKWTEDescriptionDict:
        if uniprotID in uniprotPlantDescriptionDict:
            uniprotDescription = uniprotPlantDescriptionDict[uniprotID]
            OUT_REP.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID,uniprotID,percIdentity,eValue,bitScore,queryCov,uniprotDescription))


