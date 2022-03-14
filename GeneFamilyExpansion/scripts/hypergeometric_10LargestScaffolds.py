import sys, re, os
import numpy as np
from scipy.stats import hypergeom
from statsmodels.stats.multitest import fdrcorrection
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.switch_backend('agg')
from matplotlib.lines import Line2D
import pandas as pd
import seaborn as sns


###############
# SUBROUTINES #
###############


def readLengthsfile(lengthsFile):
    scaffoldLengthDict = {}
    with open(lengthsFile,'r') as F:
        # Scaffold_1531 476495644
        for line in F:
            scaffoldID,scaffoldLen = line.strip().split(' ')
            scaffoldLengthDict[scaffoldID] = int(scaffoldLen)
            #print scaffoldID,scaffoldLen
    return(scaffoldLengthDict)


# Scaffold_1531   MAKER0000001.t1 HUMLU_CAS0000001.t1.p1  13449   23487   13449   23487
def readGeneMapFile(geneMapFile,scaffoldLengthDict):
    filteredGeneMapDict = {}
    with open(geneMapFile,'r') as F:
        for line in F:
            # OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffoldID,transcriptID,newGeneID,geneStart,geneStop,cdsStart,cdsStop))
            scaffoldID,oldGeneID,newGeneID,geneStart,geneStop,cdsStart,cdsStop = line.strip().split('\t')
            if scaffoldID in scaffoldLengthDict:
                scaffoldLen = scaffoldLengthDict[scaffoldID]
                filteredGeneMapDict[newGeneID] = scaffoldID
    return(filteredGeneMapDict)


# expanded and contracted gene families in hop
def readFamilyFile(familyFile):
    expansions = {}
    contractions = {}
    expansionCount = 0
    contractionCount = 0
    asteriskGroups = {}
    with open(familyFile,'r') as F:
        for line in F:
            if 'The labeled CAFE tree:' not in line and 'Overall' not in line:
                if 'humulusLupulus' in line:
                    line = line.strip().split(',')
                    stuff,firstOrthogroupID = line[0].split('humulusLupulus<10>:')
                    firstOrthogroupID = firstOrthogroupID.strip()
                    firstOrthogroupID,firstFamilyState = firstOrthogroupID.split('[')
                    firstFamilyState,extra = firstFamilyState.split(']')
                    #print firstOrthogroupID,firstFamilyState
                    if '*' in firstFamilyState:
                        asteriskGroups[orthogroupID] = 1
                    # remove * character to get + or -
                    firstFamilyState = firstFamilyState[0]
                    if firstFamilyState == '+':
                        if firstOrthogroupID not in expansions:
                            expansions[firstOrthogroupID] = 1
                    elif firstFamilyState == '-':
                        if firstOrthogroupID not in contractions:
                            contractions[firstOrthogroupID] = 1
                    else:
                        print "something wrong with orthogroupID"
                        sys.exit()
                        
                    for orthogroupID in line[1:]:
                        orthogroupID,familyState = orthogroupID.split('[')
                        familyState,extra = familyState.split(']')
                        #print familyState
                        if '*' in familyState:
                            asteriskGroups[orthogroupID] = 1
                        familyState = familyState[0]
                        if familyState == '+':
                            if orthogroupID not in expansions:
                                expansions[orthogroupID] = 1
                        elif familyState == '-':
                            if orthogroupID not in contractions:
                                contractions[orthogroupID] = 1
                        else:
                            print "something wrong with orthogroupID"
                            sys.exit()

    print "expansions: " + str(len(expansions))
    print "contractions: " + str(len(contractions))
    # groups with asterisks are rapidly evolving
    expansion_asterisk_count = 0
    contraction_asterisk_count = 0
    for orthogroupID in contractions:
        if orthogroupID in asteriskGroups:
            contraction_asterisk_count += 1
    for orthogroupID in expansions:
        if orthogroupID in asteriskGroups:
            expansion_asterisk_count += 1

    print "contraction asterisk count: " + str(contraction_asterisk_count)
    print "expansion asterisk count: " + str(expansion_asterisk_count)
    return(expansions,contractions)



# includes all orthogroups
def readGeneCountTSV(orthogroupGeneCountFile):
    orthogroupGeneCountDict = {}
    with open(orthogroupGeneCountFile,'r') as F:
        for line in F:
            if 'Description' not in line:
                orthogroupID,cannabis_sativa,humulus_lupulus,morus_notabilis,parasponia_andersonii,prunus_persica,trema_orientale,vitis_vinifera,ziziphus_jujuba = line.strip().split('\t')
                #print orthogroupID
                if orthogroupID not in orthogroupGeneCountDict:
                    orthogroupGeneCountDict[orthogroupID] = (int(cannabis_sativa),int(humulus_lupulus),int(morus_notabilis),int(parasponia_andersonii),int(prunus_persica),int(trema_orientale),int(vitis_vinifera),int(ziziphus_jujuba))
    return(orthogroupGeneCountDict)


# all orthogroups
def getSpeciesSpecificOrthogroupCounts(orthogroupGeneCountDict):
    speciesSpecificOrthogroupDict = {}
    for orthogroupID in orthogroupGeneCountDict:
        cannabis_sativa,humulus_lupulus,morus_notabilis,parasponia_andersonii,prunus_persica,trema_orientale,vitis_vinifera,ziziphus_jujuba = orthogroupGeneCountDict[orthogroupID]
        if orthogroupID not in speciesSpecificOrthogroupDict:
            speciesSpecificOrthogroupDict[orthogroupID] = {}

        if cannabis_sativa > 0:
            if 'cannabis_sativa' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['cannabis_sativa'] = cannabis_sativa

        if humulus_lupulus > 0:
            if 'humulus_lupulus' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['humulus_lupulus'] = humulus_lupulus
            
        if morus_notabilis > 0:
            if 'morus_notabilis' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['morus_notabilis'] = morus_notabilis

        if parasponia_andersonii > 0:
            if 'parasponia_andersonii' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['parasponia_andersonii'] = parasponia_andersonii
        
        if prunus_persica > 0:
            if 'prunus_persica' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['prunus_persica'] = prunus_persica
            
        if trema_orientale > 0:
            if 'trema_orientale' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['trema_orientale'] = trema_orientale

        if vitis_vinifera > 0:
            if 'vitis_vinifera' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['vitis_vinifera'] = vitis_vinifera

        if ziziphus_jujuba > 0:
            if 'ziziphus_jujuba' not in speciesSpecificOrthogroupDict[orthogroupID]:
                speciesSpecificOrthogroupDict[orthogroupID]['ziziphus_jujuba'] = ziziphus_jujuba
    # print speciesSpecificOrthogroupDict
    return(speciesSpecificOrthogroupDict)


def readOrthogroupTSV(orthogroupTSV,speciesSpecificOrthogroupDict):
    orthogroupDict = {}
    with open(orthogroupTSV,'r') as F:
        for line in F:
            if not line.startswith('Orthogroup'):
                orthoGroupInfo = line.strip().split('\t')
                orthogroupID = orthoGroupInfo[0]
                # print line
                if orthogroupID in speciesSpecificOrthogroupDict:
                    # cannabis
                    if 'cannabis_sativa' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'cannabis_sativa' not in orthogroupDict:
                            orthogroupDict['cannabis_sativa'] = {}
                        if orthogroupID not in orthogroupDict['cannabis_sativa']:
                            orthogroupDict['cannabis_sativa'][orthogroupID] = []
                        cannabisGeneID = orthoGroupInfo[1:][0]
                        #fullCannabisGeneID = "cannabis_sativa_" + cannabisGeneID
                        orthogroupDict['cannabis_sativa'][orthogroupID].append(cannabisGeneID)
                        # print cannabisGeneID
                        # hops
                    if 'humulus_lupulus' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'humulus_lupulus' not in orthogroupDict:
                            orthogroupDict['humulus_lupulus'] = {}
                        if orthogroupID not in orthogroupDict['humulus_lupulus']:
                            orthogroupDict['humulus_lupulus'][orthogroupID] = []
                        hopGeneID = orthoGroupInfo[1:][1]
                        #fullHopGeneID = "humulus_lupulus_" + hopGeneID
                        orthogroupDict['humulus_lupulus'][orthogroupID].append(hopGeneID)
                        #print hopGeneID
                        
                    # mulberry
                    if 'morus_notabilis' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'morus_notabilis' not in orthogroupDict:
                            orthogroupDict['morus_notabilis'] = {}
                        if orthogroupID not in orthogroupDict['morus_notabilis']:
                            orthogroupDict['morus_notabilis'][orthogroupID] = []
                        mulberryGeneID = orthoGroupInfo[1:][2]
                        #fullMulberryGeneID = "morus_notabilis_" + mulberryGeneID
                        orthogroupDict['morus_notabilis'][orthogroupID].append(mulberryGeneID)

                    # parasponia andersonii
                    if 'parasponia_andersonii' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'parasponia_andersonii' not in orthogroupDict:
                            orthogroupDict['parasponia_andersonii'] = {}
                        if orthogroupID not in orthogroupDict['parasponia_andersonii']:
                            orthogroupDict['parasponia_andersonii'][orthogroupID] = []
                        parasponiaGeneID = orthoGroupInfo[1:][3]
                        #print orthogroupID,parasponiaGeneID
                        #fullParasponiaGeneID = "parasponia_andersonii_" + parasponiaGeneID
                        orthogroupDict['parasponia_andersonii'][orthogroupID].append(parasponiaGeneID)
                    
                    # prunus persica
                    if 'prunus_persica' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'prunus_persica' not in orthogroupDict:
                            orthogroupDict['prunus_persica'] = {}
                        if orthogroupID not in orthogroupDict['prunus_persica']:
                            orthogroupDict['prunus_persica'][orthogroupID] = []
                        prunusGeneID = orthoGroupInfo[1:][4]
                        #fullPrunusGeneID = "prunus_persica_" + prunusGeneID
                        orthogroupDict['prunus_persica'][orthogroupID].append(prunusGeneID)

                    # trema orientale
                    if 'trema_orientale' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'trema_orientale' not in orthogroupDict:
                            orthogroupDict['trema_orientale'] = {}
                        if orthogroupID not in orthogroupDict['trema_orientale']:
                            orthogroupDict['trema_orientale'][orthogroupID] = []
                        tremaGeneID = orthoGroupInfo[1:][5]
                        #fullTremaGeneID = "trema_orientale_" + tremaGeneID
                        orthogroupDict['trema_orientale'][orthogroupID].append(tremaGeneID)
                    
                    # vitis vinifera
                    if 'vitis_vinifera' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'vitis_vinifera' not in orthogroupDict:
                            orthogroupDict['vitis_vinifera'] = {}
                        if orthogroupID not in orthogroupDict['vitis_vinifera']:
                            orthogroupDict['vitis_vinifera'][orthogroupID] = []
                        vitisGeneID = orthoGroupInfo[1:][6]
                        #fullVitisGeneID = "vitis_vinifera_" + vitisGeneID
                        orthogroupDict['vitis_vinifera'][orthogroupID].append(vitisGeneID)
                    
                    # ziziphus jujuba
                    if 'ziziphus_jujuba' in speciesSpecificOrthogroupDict[orthogroupID]:
                        if 'ziziphus_jujuba' not in orthogroupDict:
                            orthogroupDict['ziziphus_jujuba'] = {}
                        if orthogroupID not in orthogroupDict['ziziphus_jujuba']:
                            orthogroupDict['ziziphus_jujuba'][orthogroupID] = []
                        ziziphusGeneID = orthoGroupInfo[1:][7]
                        #fullZiziphusGeneID = "ziziphus_jujuba_" + ziziphusGeneID
                        orthogroupDict['ziziphus_jujuba'][orthogroupID].append(ziziphusGeneID)
                        # print orthogroupDict
    return(orthogroupDict)


# uniprot hits for all species
def readAnnotationFile(annotationFile,hitDict):
    with open(annotationFile,'r') as A:
        for line in A:
            if 'uniprotPlantNonRepeat' in line:
                geneID,source,uniprotID,uniprotDescription = line.strip().split('\t')
                #print geneID,source,uniprotID,uniprotDescription
                if geneID not in hitDict:
                    hitDict[geneID] = (uniprotID,uniprotDescription)
    return(hitDict)


def getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,speciesID):
    # orthogroupDict['humulus_lupulus'][orthogroupID].append(hopGeneID)
    # uniprotID,pIdent,eValue,bitScore,queryCov,uniprotDescription = topHits[geneID]
    # retrieve all genes for a species and orthogroup
    # print orthogroupDict
    for orthogroupID in orthogroupDict[speciesID]:
        for geneIDStuff in orthogroupDict[speciesID][orthogroupID]:
            if ',' in geneIDStuff:
                geneIDStuff = geneIDStuff.split(',')
                for geneID in geneIDStuff:
                    geneID = geneID.strip()
                    #print geneID
            else:
                geneID = geneIDStuff.strip()
                # print geneID
            if geneID in hitDict:
                # print geneID
                uniprotID,uniprotDescription = hitDict[geneID]
                allGenesInOrthogroups[geneID] = orthogroupID
                if orthogroupID in contractions:
                    if geneID not in allContractedGeneIDs:
                        allContractedGeneIDs[geneID] = (orthogroupID,uniprotID,uniprotDescription)
                if orthogroupID in expansions:
                    if geneID not in allExpandedGeneIDs:
                        allExpandedGeneIDs[geneID] = (orthogroupID,uniprotID,uniprotDescription)
    return(allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups)


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
                        for goTermStuff in allGOTerms:
                            goTermDesc,goTerm = goTermStuff.strip().split('[GO:')
                            goTermDesc = goTermDesc.strip()
                            goTerm,bracket = goTerm.strip().split(']')
                            goTerm = "GO:" + goTerm
                            #print goTermDesc,goTerm
                            if goTerm not in categorySpecificGOTermDict:
                                categorySpecificGOTermDict[goTerm] = goTermDesc
                    else:
                        goTermStuff = line[-1].strip()
                        #print goTerm
                        goTermDesc,goTerm = goTermStuff.strip().split('[GO:')
                        goTermDesc = goTermDesc.strip()
                        goTerm,bracket = goTerm.strip().split(']')
                        goTerm = "GO:" + goTerm
                        #print goTerm,goTermDesc
                        if goTerm not in categorySpecificGOTermDict:
                            categorySpecificGOTermDict[goTerm] = goTermDesc
    return(categorySpecificGOTermDict)


# find associations between uniprot IDs and go terms
def getGOAssociations(hitDict,uniprotGODict,goDescriptionDict,allGenesInOrthogroups,filteredGeneMapDict):
    fullGOTermList = []
    # hitDict[geneID] = (uniprotID,uniprotDescription)
    # this list only includes genes that occur in orthogroups
    for geneID in hitDict:
        if geneID in allGenesInOrthogroups:
            orthogroupID = allGenesInOrthogroups[geneID]
            uniprotID,uniprotDescription = hitDict[geneID]
            if uniprotID in uniprotGODict:
                for goTerm in uniprotGODict[uniprotID]:
                    if goTerm in goDescriptionDict:
                        goDesc = goDescriptionDict[goTerm]
                        if 'HUMLU_CAS' in geneID:
                            if geneID in filteredGeneMapDict:
                                fullGOTermList.append((orthogroupID,geneID,uniprotID,uniprotDescription,goTerm,goDesc))
                        else:
                            fullGOTermList.append((orthogroupID,geneID,uniprotID,uniprotDescription,goTerm,goDesc))
    return(fullGOTermList)


# get go term counts, in preparation for hypergeometric test
def getGOCountsForHypergeometric(fullGOTermList,selectedGeneIDs,goTermCategory,familyCategory):
    OUT = open(familyCategory + '_genesWithGOTerms_' + goTermCategory + '.txt','w')
    pathwayDesc = {}  
    pathwayToGenes = {}
    populationCount = {}
    populationTotal = 0
    sampleCount = {}
    sampleTotal = 0
    uniqueGeneIDDict = {}
    for orthogroupID,geneID,uniprotID,uniprotDesc,goTerm,goDesc in fullGOTermList:
        # 'fullGOTermList' contains genes from all species in an orthogroup that have a uniprot and GO term - but here, only interested in functional enrichment of hop
        populationTotal += 1
        pathwayDesc[goTerm] = goDesc
        if goTerm not in populationCount:
            populationCount[goTerm] = 0
        populationCount[goTerm] += 1
        if geneID in selectedGeneIDs:
            #if 'HUMLU_CAS' in geneID:
            if geneID not in uniqueGeneIDDict:
                uniqueGeneIDDict[geneID] = 1
                OUT.write("%s\t%s\t%s\n" % (geneID,uniprotID,uniprotDesc))
            OUT.write("\t%s\t%s\n" % (goTerm,goDesc))
            sampleTotal += 1
            if goTerm not in sampleCount:
                sampleCount[goTerm] = 0
            sampleCount[goTerm] += 1
            if goTerm not in pathwayToGenes:
                pathwayToGenes[goTerm] = []
            pathwayToGenes[goTerm].append((uniprotID,uniprotDesc,geneID))
    return(pathwayDesc,populationCount,populationTotal,sampleTotal,sampleCount,pathwayToGenes)


def computeSignificance(populationCount,populationTotal,sampleCount,sampleTotal,pathwayToGenes,sampleCountThreshold,familyCategory,goTermCategory,pathwayDesc):
    pathwaysUnsorted = []
    uncorrected_pValues = []
    i = 0
    DATA_OUT = open(goTermCategory + '_' + familyCategory + '_data.txt','w')
    DATA_OUT.write("goTermID\tk\texpK\tp\tqValue\tfoldChange\tgoTermIDDescription\n")
    for goTermID in sampleCount:
        N = sampleTotal # sample size, n in wikipedia
        k = sampleCount[goTermID]
        M = populationTotal       # population size, N in wikipedia
        #n = populationCount[goTermID] # population successes, K in wikipedia
        n = populationCount[goTermID] # population successes, K in wikipedia
        expK = round(float(N*n)/M, 4) # expected successes
        #expK = round(sampleTotal[goTermID]*float(n)/M, 3) # expected successes
        #expK = sampleTotal[goTermID]*float(n)/M
        i += 1
        if familyCategory == 'expanded':
            if k > expK:
                p = hypergeom.sf(k, M, n, N)   # p value
                pathwaysUnsorted.append((goTermID,p,k,expK))
                uncorrected_pValues.append(p)
        else:
            if k < expK:
                p = hypergeom.cdf(k, M, n, N)
                pathwaysUnsorted.append((goTermID,p,k,expK))
                uncorrected_pValues.append(p)

    pathwayStats = []
    pathwayStatsForPlot = []
    rejected,corrected_pValues = fdrcorrection(uncorrected_pValues, alpha=0.05, method='indep', is_sorted=False)
    fullDataList = zip(pathwaysUnsorted,corrected_pValues)
    # sort by p-value
    fullDataList.sort(key = lambda x:x[0][1])
    maxNumberInPlot = 10
    fdrThresholdForPlot = 1e-5
    for i in range(len(fullDataList)):
        data,FDR = fullDataList[i]
        goTermID,p,k,expK = data
        foldChange = float(k)/expK
        goTermIDDescription = pathwayDesc[goTermID]
        if FDR < fdrThreshold:
            pathwayStats.append((goTermID,k,expK,p,FDR,foldChange))
            DATA_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (goTermID,k,expK,p,FDR,foldChange,goTermIDDescription))
        if FDR < fdrThresholdForPlot:
            ### pathwayStats.append((goTermID,k,expK,p,FDR,foldChange))
            pathwayStatsForPlot.append((goTermID,k,expK,p,FDR,foldChange))
            #less molecularFunction_expanded_data.txt | sort -k2gr -k5g -k4g -k6gr | less
            #[Linux@anybase hypergeometric_02052022]$ less molecularFunction_expanded_data.txt | sort -k6gr -k5g -k4g -k2gr | less -S
    if familyCategory == 'expanded':
        # in this case, sorting x[1] largest to smallest, these are largest observed cases
        # these are sorted according to largest observed go terms because otherwise sorting on fdr and pvalue first turns up examples where it's 1 observed and 0.2 expected... this goes on for a while and even though it's a tiny bit more statistically significant, biologically it doesnt say much
        # pathwayStats = sorted(pathwayStats, key=lambda x: (x[4], x[3], -x[1], -x[5]))
        pathwayStats = sorted(pathwayStats, key=lambda x: (-x[1], x[4], x[3], -x[5]))
        pathwayStatsForPlot = sorted(pathwayStatsForPlot, key=lambda x: (-x[1], x[4], x[3], -x[5]))
        pathwayStatsForPlot = pathwayStatsForPlot[0:maxNumberInPlot]
    else:
        # contracted
        # here, sorting -x[2], largest to smallest expected
        # pathwayStats = sorted(pathwayStats, key=lambda x: (x[4], x[3], -x[2], -x[5]))
        pathwayStats = sorted(pathwayStats, key=lambda x: (-x[2], x[4], x[3], -x[5]))
        pathwayStatsForPlot = sorted(pathwayStatsForPlot, key=lambda x: (-x[1], x[4], x[3], -x[5]))
        pathwayStatsForPlot = pathwayStatsForPlot[0:maxNumberInPlot]
    return(pathwayStats,pathwayToGenes,pathwayStatsForPlot)


def formatData(dataList,pathwayDesc,familyCategory,goTermCategory):
    # x = ['A', 'B', 'C', 'D']
    # y1 = [10, 20, 10, 30]
    # y2 = [20, 25, 15, 25]
    TOP_OUT = open(goTermCategory + '_top_go_terms_' + familyCategory + '.txt','w')
    TOP_OUT.write("goTermID\tk\texpK\tp\tqValue\n")
    goTermIDs = []
    goTermIDDescriptions = []
    expKList = []
    kList = []
    dataListForFigure = []
    if familyCategory == 'expanded':
        dataList.sort(key=lambda x: x[1], reverse=True)
    else:
        dataList.sort(key=lambda x: x[2], reverse=True)
    for goTermID,k,expK,p,qValue,foldChange in dataList:
        goTermIDDescription = pathwayDesc[goTermID]
        goTermIDDescriptions.append(goTermIDDescription)
        goTermIDs.append(goTermID)
        goTermLabel = goTermID + "\n" + goTermIDDescription
        TOP_OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (goTermID,k,expK,p,qValue,foldChange,goTermIDDescription))
        if familyCategory == 'expanded':
            # smaller number is expected
            expKList.append(expK)
            kList.append(k)
            dataListForFigure.append(('Expected k',goTermLabel,expK))
            dataListForFigure.append(('Observed k',goTermLabel,k))
        else:
            # contracted
            # larger number is expected
            expKList.append(expK)
            kList.append(k)
            dataListForFigure.append(('Expected k',goTermLabel,expK))
            dataListForFigure.append(('Observed k',goTermLabel,k))
    return(goTermIDs,expKList,kList,goTermIDDescriptions,dataListForFigure)


'''
speciesCentricDataDict['ziziphus_jujuba'].append(float(z_jujuba))
        
    dataList = []
    for speciesID in speciesCentricDataDict:
        #print speciesID,speciesCentricDataDict[speciesID][0],speciesCentricDataDict[speciesID][1],
speciesCentricDataDict[speciesID][2]
        # 1 -- 'Percentage of genes in orthogroups'
        # 2 -- 'Percentage of orthogroups containing species'
        # 3 -- 'Percentage of genes in species-specific orthogroups'
        dataList.append(("genes_in_orthogroups",speciesID,speciesCentricDataDict[speciesID][0]))
        dataList.append(("orthogroups_containing_species",speciesID,speciesCentricDataDict[speciesID][1]))
        dataList.append(("species_specific",speciesID,speciesCentricDataDict[speciesID][2]))
    #print dataList
    #dataArray = np.array(dataList)
    return(dataList)

def plot(dataList,outName,plotTitle):
    dataList = sorted(dataList)
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    df = pd.DataFrame(dataList, columns=['Orthogroup category','Species','Percentage'])
    sns.barplot(data=df, x='Species', y='Percentage', hue='Orthogroup category')
    xticklabels = ax1.get_xticklabels()
    ax1.set_xticklabels(xticklabels, rotation=45, ha='right', rotation_mode='anchor')
    plt.tight_layout()
    plt.title(plotTitle)
    plt.savefig(outName + ".png", dpi=600)
    plt.savefig(outName + ".pdf")
    plt.savefig(outName + ".svg")
'''

def plot(goTermIDs,bar1,bar2,goTermIDDescriptions,outName,familyCategory,goTermCategory,dataListForFigure):
    #dataListForFigure = sorted(dataListForFigure)
    #dataArray = np.array(dataListForFigure) 
    #print dataArray
    #x = []
    #for i,j in zip(goTermIDs,goTermIDDescriptions):
    #    label = i + "\n" + j
    #    x.append(label)

    #y1 = bar1
    #y2 = bar2

    #fig = plt.figure(figsize=(15,12))
    #fig = plt.figure(figsize=(12.5,10))
    ###fig = plt.figure(figsize=(11.25,9))
    fig = plt.figure(figsize=(9,9))
    #fig = plt.figure(figsize=(10,8))
    #fig = plt.figure(figsize=(8.75,7))
    #fig = plt.figure(figsize=(7.5,6))
    #fig = plt.figure(figsize=(5,4))
    ax1 = fig.add_subplot(111)
    # dataListForFigure.append((goTermIDs,expK,k))
    # ('Expected k',goTermID,expK))
    df = pd.DataFrame(dataListForFigure, columns=['GOTermCountCategory','GOTermID','GOTermCount'])
    # columns=['Orthogroup category','Species','Percentage']
    #barData = sns.barplot(data=df, x='GOTermID', y='GOTermCount', hue='GOTermCountCategory')
    barData = sns.barplot(data=df, x='GOTermCount', y='GOTermID', hue='GOTermCountCategory')
    yticklabels = barData.get_yticklabels()
    barData.set_yticklabels(yticklabels, rotation=25, ha='right', rotation_mode='anchor')
    #xticklabels = barData.get_xticklabels()
    #barData.set_xticklabels(xticklabels, rotation=45, ha='right', rotation_mode='anchor')

    #left = 0

    #plt.barh(x, y1, color='#42033D')
    #plt.barh(x, y2, left=y1, color='#FF9000')

    #xticklabels = ax1.get_xticklabels()
    #ax1.set_xticklabels(x, rotation=45, ha='right', rotation_mode='anchor')
    #ax1.set_yticks(y_pos, labels=x)
    #plt.tight_layout()
    #custom_lines = [Line2D([0], [0], color='#42033D', lw=4),
    #                Line2D([0], [0], color='#FF9000', lw=4)]
    #if familyCategory == 'expanded':
    #    ax1.legend(custom_lines, ['Expected k','Observed k'])
    #else:
    #    ax1.legend(custom_lines, ['Observed k','Expected k'])
    plt.legend()
    plt.savefig(goTermCategory + outName + ".png", dpi=600)
    plt.savefig(goTermCategory + outName + ".pdf")
    plt.savefig(goTermCategory + outName + ".svg")


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <orthogroup gene count file> <orthogroup tsv file> <cafe family file> <hop> <cannabis> <mulberry> <parasponia> <prunus> <trema> <vitis> <ziziphus> <scaffold lengths file> <gene map file> <all GO term file> <category-specific GO term file> <GO term category>"
if len(sys.argv) != 17:
    print usage
    sys.exit()

orthogroupGeneCountFile = sys.argv[1]
orthogroupTSV = sys.argv[2]
familyFile = sys.argv[3]
hopAnnotationFile = sys.argv[4]
cannabisAnnotationFile = sys.argv[5]
mulberryAnnotationFile = sys.argv[6]
parasponiaAnnotationFile = sys.argv[7]
prunusAnnotationFile = sys.argv[8]
tremaAnnotationFile = sys.argv[9]
vitisAnnotationFile = sys.argv[10]
ziziphusAnnotationFile = sys.argv[11]
lengthsFile = sys.argv[12]
geneMapFile = sys.argv[13]
allGOTermFile = sys.argv[14]
categorySpecificGOTermFile = sys.argv[15]
goTermCategory = sys.argv[16]


scaffoldLengthDict = readLengthsfile(lengthsFile)
filteredGeneMapDict = readGeneMapFile(geneMapFile,scaffoldLengthDict)

expansions,contractions = readFamilyFile(familyFile)
orthogroupGeneCountDict = readGeneCountTSV(orthogroupGeneCountFile)

speciesSpecificOrthogroupDict = getSpeciesSpecificOrthogroupCounts(orthogroupGeneCountDict)
#print speciesSpecificOrthogroupDict

orthogroupDict = readOrthogroupTSV(orthogroupTSV,speciesSpecificOrthogroupDict)
#print orthogroupDict

hitDict = {}
hitDict = readAnnotationFile(cannabisAnnotationFile,hitDict)
hitDict = readAnnotationFile(hopAnnotationFile,hitDict)
hitDict = readAnnotationFile(mulberryAnnotationFile,hitDict)
hitDict = readAnnotationFile(parasponiaAnnotationFile,hitDict)
hitDict = readAnnotationFile(prunusAnnotationFile,hitDict)
hitDict = readAnnotationFile(tremaAnnotationFile,hitDict)
hitDict = readAnnotationFile(vitisAnnotationFile,hitDict)
hitDict = readAnnotationFile(ziziphusAnnotationFile,hitDict)

allContractedGeneIDs = {}
allExpandedGeneIDs = {}
allGenesInOrthogroups = {}

allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'humulus_lupulus')
allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'cannabis_sativa')
allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'morus_notabilis')
allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'parasponia_andersonii')
allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'prunus_persica')
allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'trema_orientale')
allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'vitis_vinifera')
allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups = getAnnotationForOrthogroup(orthogroupDict,hitDict,contractions,expansions,allContractedGeneIDs,allExpandedGeneIDs,allGenesInOrthogroups,'ziziphus_jujuba')

# hypergeometric thresholds
sampleCountThreshold = 2
fdrThreshold = 0.05
#fdrThreshold = 0.01

# go term stuff
categorySpecificGOTermDict = readCategorySpecificGOTermFile(categorySpecificGOTermFile)
goDescriptionDict,uniprotGODict = readGODescriptionFile(allGOTermFile,categorySpecificGOTermDict)
fullGOTermList = getGOAssociations(hitDict,uniprotGODict,goDescriptionDict,allGenesInOrthogroups,filteredGeneMapDict)

# hypergeometric
# contracted
cPathwayDesc,cPopulationCount,cPopulationTotal,cSampleTotal,cSampleCount,cPathwayToGenes = getGOCountsForHypergeometric(fullGOTermList,allContractedGeneIDs,goTermCategory,'contracted')
# expanded
ePathwayDesc,ePopulationCount,ePopulationTotal,eSampleTotal,eSampleCount,ePathwayToGenes = getGOCountsForHypergeometric(fullGOTermList,allExpandedGeneIDs,goTermCategory,'expanded')


# contracted
cPathwayStats,cPathwayToGenes,cPathwayStatsForPlot = computeSignificance(cPopulationCount,cPopulationTotal,cSampleCount,cSampleTotal,cPathwayToGenes,sampleCountThreshold,'contracted',goTermCategory,cPathwayDesc)
# expanded
ePathwayStats,ePathwayToGenes,ePathwayStatsForPlot = computeSignificance(ePopulationCount,ePopulationTotal,eSampleCount,eSampleTotal,ePathwayToGenes,sampleCountThreshold,'expanded',goTermCategory,ePathwayDesc)

# contracted
c_goTermIDs,c_expKList,c_kList,c_goTermIDDescriptions,c_dataListForFigure = formatData(cPathwayStatsForPlot,cPathwayDesc,'contracted',goTermCategory)
plot(c_goTermIDs,c_expKList,c_kList,c_goTermIDDescriptions,'_top_contracted_bar_chart','contracted',goTermCategory,c_dataListForFigure)

# expanded
e_goTermIDs,e_expKList,e_kList,e_goTermIDDescriptions,e_dataListForFigure = formatData(ePathwayStatsForPlot,ePathwayDesc,'expanded',goTermCategory)
plot(e_goTermIDs,e_expKList,e_kList,e_goTermIDDescriptions,'_top_expanded_bar_chart','expanded',goTermCategory,e_dataListForFigure)

# contracted
outputFile1 = "contracted_hypergeometric_" + goTermCategory + "_geneIDsIncluded.txt"
GO1 = open(outputFile1,'w')
GO1.write("geneID\tuniprotID\tgoTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tpopulationCount\tpopulationTotal\tgoTermDescription\n")

outputFile2 = "contracted_hypergeometric_" + goTermCategory + ".txt"
GO2 = open(outputFile2,'w')
GO2.write("goTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tpopulationCount\tpopulationTotal\tgoTermDescription\n")

for cItem in cPathwayStats:    
    # goTermID,p,k,expK = data
    c_goTermID,c_obsK,c_expK,c_p,c_qValue,c_foldChange = cItem
    #print qValue
    if c_qValue <= fdrThreshold:
        #c_p = round(c_p, 5)
        #c_qValue = round(c_qValue, 5)
        GO2.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (c_goTermID,c_p,c_qValue,c_obsK,c_expK,c_foldChange,cPopulationCount,cPopulationTotal,cPathwayDesc[c_goTermID])) 
        for cUniprotID,cUniprotDesc,cGeneID in cPathwayToGenes[c_goTermID]:
            GO1.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (cGeneID,cUniprotID,c_goTermID,c_p,c_qValue,c_obsK,c_expK,c_foldChange,cPopulationCount,cPopulationTotal,cPathwayDesc[c_goTermID]))


# expanded
# ePopulationCount,ePopulationTotal 
outputFile3 = "expanded_hypergeometric_" + goTermCategory + "_geneIDsIncluded.txt"
GO3 = open(outputFile3,'w')
GO3.write("geneID\tuniprotID\tgoTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tpopulationCount\tpopulationTotal\tgoTermDescription\n")

outputFile4 = "expanded_hypergeometric_" + goTermCategory + ".txt"
GO4 = open(outputFile4,'w')
GO4.write("goTermID\tpValue\tqValue\tobsK\texpK\tfoldChange\tpopulationCount\tpopulationTotal\tgoTermDescription\n")

#ePathwayStats.sort(key=lambda x: (x[4], x[3]))
for eItem in ePathwayStats:
    e_goTermID,e_obsK,e_expK,e_p,e_qValue,e_foldChange = eItem
    if e_qValue <= fdrThreshold:
        #e_p = round(e_p, 5)
        #e_qValue = round(e_qValue, 5)
        GO4.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (e_goTermID,e_p,e_qValue,e_obsK,e_expK,e_foldChange,ePopulationCount,ePopulationTotal,ePathwayDesc[e_goTermID]))
        for eUniprotID,eUniprotDesc,eGeneID in ePathwayToGenes[e_goTermID]:
            GO3.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (eGeneID,eUniprotID,e_goTermID,e_p,e_qValue,e_obsK,e_expK,e_foldChange,ePopulationCount,ePopulationTotal,ePathwayDesc[e_goTermID]))


