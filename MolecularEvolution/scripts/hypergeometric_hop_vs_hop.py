import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
from scipy.stats import hypergeom
import numpy as np
from statsmodels.stats.multitest import fdrcorrection

###############
# SUBROUTINES #
###############

def readTopHitsFile(topHitsFile):
    topHits = {}
    with open(topHitsFile,'r') as HIT:
        for line in HIT:
            # XP_030489797.1  uniprotPlantNonRepeat   Q1PE49  sp|Q1PE49|4ON1_ARATH Protein At-4/1 OS=Arabidopsis thaliana OX=3702 GN=At4g26020 PE=1 SV=1
            if 'uniprotPlantNonRepeat' in line:
                geneID,source,uniprotID,uniprotGeneName = line.strip().split('\t')
                topHits[geneID] = (uniprotID,uniprotGeneName)
    return(topHits)

def parseCollinearityFile(collinearityFile):
    collinearGeneDict = {}
    collinearGeneList = []
    with open(collinearityFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                numberBlock,geneID1,geneID2,eValue = line.strip().split('\t')
                #print numberBlock,geneID1,geneID2,eValue
                if 'HUMLU_CAS' in geneID1 and 'HUMLU_CAS' in geneID2:
                    #print geneID1
                    if geneID1 not in collinearGeneDict:
                        collinearGeneDict[geneID1] = 1
                        collinearGeneList.append(geneID1)
                    if geneID2 not in collinearGeneDict:
                        collinearGeneDict[geneID2] = 1
                        collinearGeneList.append(geneID2)
    #print len(collinearGeneDict)
    #print len(collinearGeneList)
    return(collinearGeneDict,collinearGeneList)


# Entry   Entry name      Protein names   Organism        Gene ontology IDs       Gene ontology (GO)
# all go terms
def readGODescriptionFile(allGOTermFile,categorySpecificGOTermDict):
    goDescriptionDict = {}
    uniprotGODict = {}
    with open(allGOTermFile,'r') as GO:
        for line in GO:
            if 'Entry' not in line:
                line = line.strip().split('\t')
                uniprotID = line[0]
                #print uniprotID
                uniprotName = line[1]
                uniprotDescription = line[2]
                species = line[3]
                # requiring a length > 0 means that the uniprot id has an associated go term
                if len(line[4:]) > 0:
                    #print line[5]
                    if ';' in line[5]:
                        goTermStuff = line[5].split(';')
                        for item in goTermStuff:
                            goTermDesc,goTerm = item.split('[GO:')
                            goTermDesc = goTermDesc.strip()
                            goTerm,bracket = goTerm.split(']')
                            goTerm = "GO:"  + goTerm
                            #print goTerm
                            # this if-statement requires that the go term be in the category-specific dict
                            if goTerm in categorySpecificGOTermDict:
                                if goTerm not in goDescriptionDict:
                                    goDescriptionDict[goTerm] = goTermDesc
                                if uniprotID not in uniprotGODict:
                                    uniprotGODict[uniprotID] = []
                                uniprotGODict[uniprotID].append(goTerm)
                                # otherwise, there's no ';' - meaning, there's only one go term for that uniprot id
                    else:
                        # calcium ion binding [GO:0005509]
                        #print line[5]
                        getGOTerm = re.search('(.+)\[(GO:.+)\]',line[5])
                        goTermDesc = getGOTerm.group(1)
                        goTermDesc = goTermDesc.strip()
                        goTerm = getGOTerm.group(2)
                        if goTerm in categorySpecificGOTermDict:
                            if goTerm not in goDescriptionDict:
                                goDescriptionDict[goTerm] = goTermDesc
                            if uniprotID not in uniprotGODict:
                                uniprotGODict[uniprotID] = []
                            uniprotGODict[uniprotID].append(goTerm)
                            #print goTerm,goTermDesc
                # else, no go term
                #else:
                #    #print line[3:]
                #    continue
    # uniprotGODict[uniprotID].append(goTerm)
    #for uniprotID in uniprotGODict:
    #    for goTerm in uniprotGODict[uniprotID]:
    #        print uniprotID,goTerm
    return(goDescriptionDict,uniprotGODict)


# Entry   Entry name      Protein names   Organism        Gene ontology IDs       Gene ontology (biological process)
# category specific go term file
def readCategorySpecificGOTermFile(categorySpecificGOFileWithDesc):
    categorySpecificGOTermDict = {}
    with open(categorySpecificGOFileWithDesc,'r') as GO:
        for line in GO:
            if 'Entry' not in line:
                line = line.strip().split('\t')
                uniprotID = line[0]
                uniprotName = line[1]
                uniprotDescription = line[2]
                species = line[3]
                #print line
                # integral component of membrane [GO:0016021]
                # anchored component of membrane [GO:0031225]; plasma membrane [GO:0005886]
                # re the if statement below: the brackets enclose the go term and the formatting is from the column corresponding to the category-specific go terms...
                # the last column will always have go terms but if they don't have the brackets and a description, they're not category specific
                # so in selecting the last column, if it only includes go terms and no description, those are indicative of the full set of go terms (all categories) for a given uniprot ID
                if '[' in line[-1] and ']' in line[-1]:
                    if ';' in line[-1]:
                        # category-specific go terms
                        allGOTerms = line[-1].split(';')
                        #print goTermStuff
                        #print allGOTerms
                        for goTermStuff in allGOTerms:
                            goTermDesc,goTerm = goTermStuff.strip().split('[GO:')
                            goTermDesc = goTermDesc.strip()
                            goTerm,bracket = goTerm.strip().split(']')
                            goTerm = "GO:" + goTerm
                            #print uniprotID,goTermDesc,goTerm
                            if goTerm not in categorySpecificGOTermDict:
                                categorySpecificGOTermDict[goTerm] = goTermDesc
                    else:
                        goTermStuff = line[-1].strip()
                        #print goTerm
                        goTermDesc,goTerm = goTermStuff.strip().split('[GO:')
                        goTermDesc = goTermDesc.strip()
                        goTerm,bracket = goTerm.strip().split(']')
                        goTerm = "GO:" + goTerm
                        #print uniprotID,goTerm,goTermDesc
                        if goTerm not in categorySpecificGOTermDict:
                            categorySpecificGOTermDict[goTerm] = goTermDesc
    return(categorySpecificGOTermDict)


# find associations between uniprot IDs and go terms
def getGOAssociations(topHits,uniprotGODict,goDescriptionDict):
    fullGOTermList = []
    for geneID in topHits:
        uniprotID,uniprotGeneName = topHits[geneID]
        if uniprotID in uniprotGODict:
            for goTerm in uniprotGODict[uniprotID]:
                #print uniprotID,goTerm
                if goTerm in goDescriptionDict:
                    #print uniprotID,goTerm 
                    goDesc = goDescriptionDict[goTerm]
                    fullGOTermList.append((geneID,uniprotID,uniprotGeneName,goTerm,goDesc))
    return(fullGOTermList)


# get go term counts, in preparation for hypergeometric test
def getGOCountsForHypergeometric(fullGOTermList,collinearGeneDict,goTermCategory):
    OUT = open('genesWithGOTerms_' + goTermCategory + '.txt','w')
    pathwayDesc = {}      
    pathwayToGenes = {}

    populationCount = {}
    populationTotal = 0   
    sampleCount = {}
    sampleTotal = 0
    uniqueGeneIDDict = {}
    for hopGeneID,uniprotID,uniprotDesc,goTerm,goDesc in fullGOTermList:
        #print hopGeneID,uniprotID,uniprotDesc,goTerm,goDesc
        # sample total is just the total number of collinear genes
        populationTotal += 1
        pathwayDesc[goTerm] = goDesc
        if goTerm not in populationCount:
            populationCount[goTerm] = 0
        populationCount[goTerm] += 1
        if hopGeneID in collinearGeneDict:
            if hopGeneID not in uniqueGeneIDDict:
                uniqueGeneIDDict[hopGeneID] = 1
                OUT.write("%s\t%s\t%s\n" % (hopGeneID,uniprotID,uniprotDesc))
            OUT.write("\t%s\t%s\n" % (goTerm,goDesc))
            sampleTotal += 1
            if goTerm not in sampleCount:
                sampleCount[goTerm] = 0
            sampleCount[goTerm] += 1
            if goTerm not in pathwayToGenes:
                pathwayToGenes[goTerm] = []
            pathwayToGenes[goTerm].append((uniprotID,uniprotDesc,hopGeneID))
    return(pathwayDesc,populationCount,populationTotal,sampleTotal,sampleCount,pathwayToGenes)


def computeSignificance(populationCount,populationTotal,sampleCount,sampleTotal,pathwayToGenes,sampleCountThreshold,pathwayDesc):
    pathwaysUnsorted = []
    uncorrected_pValues = []
    i = 0
    for goTermID in sampleCount:
        #if sampleCount[goTermID] >= sampleCountThreshold:
        N = sampleTotal # sample size, n in wikipedia
        k = sampleCount[goTermID] # sample successes, k in wikipedia
        M = populationTotal       # population size, N in wikipedia
        n = populationCount[goTermID] # population successes, K in wikipedia
        #n = populationCount            # population successes, K in wikipedia
        #p = hypergeom.sf(k, M, n, N)   # p value
        expK = round(float(N*n)/M, 5) # expected successes
        i += 1
        #if k > expK:
        #p = hypergeom.sf(k, M, n, N)
        p = 2*min(hypergeom.cdf(k, M, n, N),hypergeom.sf(k, M, n, N))
        pathwaysUnsorted.append((goTermID,p,k,expK))
        uncorrected_pValues.append(p)

    pathwayStats = []
    rejected,corrected_pValues = fdrcorrection(uncorrected_pValues, alpha=0.05, method='indep', is_sorted=False)
    fullDataList = zip(pathwaysUnsorted,corrected_pValues)
    # sort by p-value 
    fullDataList.sort(key = lambda x:x[0][1])
    for i in range(len(fullDataList)):
        data,FDR = fullDataList[i]
        goTermID,p,k,expK = data
        foldChange = float(k)/expK
        goTermIDDescription = pathwayDesc[goTermID]
        #if FDR < fdrThreshold:
        pathwayStats.append((goTermID,k,expK,p,FDR,foldChange))
    # sort pathway list with stats, ascending by q value
    pathwayStats = sorted(pathwayStats, key=lambda q: q[4])
    return(pathwayStats,pathwayToGenes)


############
# MAIN  ####
############


usage = "Usage: " + sys.argv[0] + " <top hit file> <collinearity file> <go term file> <category-specific go term file> <go term category>\n"
if len(sys.argv) != 6:
    print usage
    sys.exit()

topHitsFile = sys.argv[1]
collinearityFile = sys.argv[2]
allGOTermFile = sys.argv[3]
categorySpecificGOTermFile = sys.argv[4]
goTermCategory = sys.argv[5]

# sampleCountThreshold = 2
sampleCountThreshold = 3
fdrThreshold = 0.05

# functional annotations
topHits = readTopHitsFile(topHitsFile)

# collinearity file from mcscanx
collinearGeneDict,collinearGeneList = parseCollinearityFile(collinearityFile)

# go term stuff
categorySpecificGOTermDict = readCategorySpecificGOTermFile(categorySpecificGOTermFile)
goDescriptionDict,uniprotGODict = readGODescriptionFile(allGOTermFile,categorySpecificGOTermDict)
fullGOTermList = getGOAssociations(topHits,uniprotGODict,goDescriptionDict)

# hypergeometric
pathwayDesc,populationCount,populationTotal,sampleTotal,sampleCount,pathwayToGenes = getGOCountsForHypergeometric(fullGOTermList,collinearGeneDict,goTermCategory)

pathwayStats,pathwayToGenes = computeSignificance(populationCount,populationTotal,sampleCount,sampleTotal,pathwayToGenes,sampleCountThreshold,pathwayDesc)

outputFile1 = "hypergeometric_" + goTermCategory + "_geneIDsIncluded.txt"
GO1 = open(outputFile1,'w')
GO1.write("hopGeneID\tuniprotID\tgoTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tgoTermDescription\n")

outputFile2 = "hypergeometric_" + goTermCategory + ".txt"
GO2 = open(outputFile2,'w')
GO2.write("goTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tgoTermDescription\n")

outputFile3 = "hypergeometric_" + goTermCategory + "_geneIDsIncluded_allFDRValues.txt"
GO3 = open(outputFile3,'w')
GO3.write("hopGeneID\tuniprotID\tgoTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tgoTermDescription\n")

outputFile4 = "hypergeometric_" + goTermCategory + "_allFDRValues.txt"
GO4 = open(outputFile4,'w')
GO4.write("goTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tgoTermDescription\n")

for item in pathwayStats:
    goTermID,obsK,expK,p,qValue,foldChange = item
    #print qValue
    #p = round(p, 5)
    #qValue = round(qValue, 5)
    if qValue <= fdrThreshold:
        #GO1.write("%s\t%s\t%s\t%s\t%s\n" % (goTermID,qValue,obsK,expK,pathwayDesc[goTermID])) 
        #GO1.write("%s\tFDR=%.3f\t%d genes, %.2e expected\t%s\n" % (goTermID,qValue,k,expK,pathwayDesc[goTermID]))
        GO2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (goTermID,p,qValue,obsK,expK,foldChange,pathwayDesc[goTermID])) 
        for uniprotID,uniprotDesc,hopGeneID in pathwayToGenes[goTermID]:
            GO1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (hopGeneID,uniprotID,goTermID,p,qValue,obsK,expK,foldChange,pathwayDesc[goTermID]))
    GO4.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (goTermID,p,qValue,obsK,expK,foldChange,pathwayDesc[goTermID]))
    for uniprotID,uniprotDesc,hopGeneID in pathwayToGenes[goTermID]:
        GO3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (hopGeneID,uniprotID,goTermID,p,qValue,obsK,expK,foldChange,pathwayDesc[goTermID]))

