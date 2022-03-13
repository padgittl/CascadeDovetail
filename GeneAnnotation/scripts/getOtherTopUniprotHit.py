import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import math

###############
# SUBROUTINES #
###############

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
            # HUMLU_CAS0065976.t1     Q89703 
            geneID,uniprotID,pIdent,length,mismatch,gapOpen,qStart,qEnd,sStart,sEnd,eValue,bitScore,queryCov = line.strip().split("\t")
            pIdent = float(pIdent)
            length = int(length)
            eValue = float(eValue)
            bitScore = float(bitScore)
            queryCov = int(queryCov)
            percentIdentityThreshold = filterTwilightThresholdGenes(length)
            if eValue < eValueThreshold:
                if queryCov >= queryCoverageThreshold:
                    if pIdent >= percentIdentityThreshold:
                        if geneID in bestHits_gene_vs_uniprot:
                            if bestScore[geneID] < bitScore:
                                bestScore[geneID] = bitScore
                                bestHits_gene_vs_uniprot[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
                        else:
                            bestScore[geneID] = bitScore
                            bestHits_gene_vs_uniprot[geneID] = (uniprotID,pIdent,eValue,bitScore,queryCov)
    return(bestHits_gene_vs_uniprot)


def readBlastpOutputFile_uniprot_vs_gene(pepBlastpOutputFile2,eValueThreshold,queryCoverageThreshold):
    bestHits_uniprot_vs_gene = {}
    bestScore = {}
    with open(pepBlastpOutputFile2,'r') as BF2:
        for line in BF2:
            # sp|O55756|VF157_IIV6    HUMLU_CAS0023436.t1
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
    for geneID1 in bestHits_gene_vs_uniprot:
        uniprotID1,percIdentity1,eValue1,bitScore1,queryCov1 = bestHits_gene_vs_uniprot[geneID1]
        #print uniprotID1
        if geneID1 in bestHitsBothDirections:
            if bestScoreBothDirections[geneID1] < bitScore1:
                bestScoreBothDirections[geneID1] = bitScore1
                bestHitsBothDirections[geneID1] = (uniprotID1,percIdentity1,eValue1,bitScore1,queryCov1)
        else:
            bestScoreBothDirections[geneID1] = bitScore1
            bestHitsBothDirections[geneID1] = (uniprotID1,percIdentity1,eValue1,bitScore1,queryCov1)

    for geneID2 in bestHits_uniprot_vs_gene:
        uniprotID2,percIdentity2,eValue2,bitScore2,queryCov2 = bestHits_uniprot_vs_gene[geneID2]
        #print uniprotID2
        #if geneID2 not in bestHits_gene_vs_uniprot:
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


usage = "Usage: " + sys.argv[0] + " <blastp output file 1> <blastp output file 2> <uniprot fasta file> <uniprot set, e.g. Virus or Bacteria, or UniProtTEs> <query coverage threshold, e.g. 90> <species, e.g. hop or can>\n"
if len(sys.argv) != 7:
    print usage
    sys.exit()


blastpOutputFile1 = sys.argv[1]
blastpOutputFile2 = sys.argv[2]
fastaFile = sys.argv[3]
uniprotSet = sys.argv[4]
queryCoverageThresholdValue = sys.argv[5]
species = sys.argv[6]


eValueThreshold = 1e-3
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3052304/
queryCoverageThreshold = int(queryCoverageThresholdValue)

bestHits_gene_vs_uniprot = readBlastpOutputFile_gene_vs_uniprot(blastpOutputFile1,eValueThreshold,queryCoverageThreshold)
bestHits_uniprot_vs_gene = readBlastpOutputFile_uniprot_vs_gene(blastpOutputFile2,eValueThreshold,queryCoverageThreshold)

uniprotDescriptionDict = getDescriptionFromFasta(fastaFile)

bestHitsBothDirections = getTopHitsBothDirections(bestHits_gene_vs_uniprot,bestHits_uniprot_vs_gene)


OUT = open(species + 'Top' + uniprotSet + 'Hits_queryCov' + queryCoverageThresholdValue + '.txt','w')
OUT.write("geneID\tuniprotID\tpercentIdentity\teValue\tbitScore\tqueryCoverage\tuniprotDescription\n")
for geneID in bestHitsBothDirections:
    uniprotID,percIdentity,eValue,bitScore,queryCov = bestHitsBothDirections[geneID]
    if uniprotID in uniprotDescriptionDict:
        uniprotDescription = uniprotDescriptionDict[uniprotID]
        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (geneID,uniprotID,percIdentity,eValue,bitScore,queryCov,uniprotDescription))

        
    



