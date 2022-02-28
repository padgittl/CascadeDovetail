import sys, re, os
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


def readGeneFasta(geneFasta):
    geneIDs = {}
    for record in SeqIO.parse(geneFasta,"fasta"):
        if record.id not in geneIDs:
            geneIDs[record.id] = 1
    return(geneIDs)


# PF13612.7       Transposase DDE domain 
def readRepeatDomainFile(repeatDomainFile):
    repeatDomainDict = {}
    with open(repeatDomainFile,'r') as R:
        for line in R:
            pfamID,pfamDescription = line.strip().split('\t')
            if pfamID not in repeatDomainDict:
                repeatDomainDict[pfamID] = pfamDescription
    return(repeatDomainDict)


def readHMMFile(hmmFile,eValueThreshold,repeatDomainDict):
    OUT_HMM = open('pfamDomainsBelowEvalue' + str(eValueThreshold) + '.txt','w')
    domainCountDict = {}
    domainCount = 0
    repeatDict = {}
    nonRepeatDict = {}
    geneLenDict = {}
    with open(hmmFile,'r') as HMM:
        for line in HMM:
            if not line.startswith('#'):
                if not line.startswith('-'):
                    if not line.isspace():
                        domainInfo = line.strip().split()[0:22]
                        accessionID = domainInfo[1]
                        geneID = domainInfo[3]
                        geneLen = domainInfo[5]
                        eValue = domainInfo[6]
                        score = domainInfo[7]
                        domainAlignmentStart = domainInfo[17]
                        domainAlignmentStop = domainInfo[18]
                        alignmentProbability = domainInfo[21]
                        accessionDesc = line.strip().split()[22:]
                        accessionDesc = ' '.join(accessionDesc)
                        eValue = float(eValue)
                        geneLen = int(geneLen)
                        domainAlignmentStart = int(domainAlignmentStart)
                        domainAlignmentStop = int(domainAlignmentStop)
                        if eValue < eValueThreshold:
                            OUT_HMM.write("%s\t%s\t%s\t%s\n" % (geneID,accessionID,eValue,accessionDesc))
                            if geneID not in domainCountDict:
                                domainCountDict[geneID] = 0
                            domainCountDict[geneID] += 1
                            if geneID not in geneLenDict:
                                geneLenDict[geneID] = geneLen
                            # accessionID is in repeat dict 
                            if accessionID in repeatDomainDict:
                                if geneID not in repeatDict:
                                    repeatDict[geneID] = []
                                repeatDict[geneID].append((accessionID,accessionDesc,domainAlignmentStart,domainAlignmentStop,geneLen))
                                # accessionID is in non-repeat dict
                            else:
                                if geneID not in nonRepeatDict:
                                    nonRepeatDict[geneID] = []
                                nonRepeatDict[geneID].append((accessionID,accessionDesc,domainAlignmentStart,domainAlignmentStop,geneLen))
    #print len(domainCountDict)
    #for geneID in domainCountDict:
    #    print geneID,domainCountDict[geneID]
    return(nonRepeatDict,repeatDict,domainCountDict,geneLenDict)


def organizeGenes(geneIDs,nonRepeatDict,repeatDict):
    genesWithBoth = {}
    genesWithRepeatOnly = {}
    genesWithNonRepeatOnly = {}
    genesWithoutPfam = {}
    for geneID in geneIDs:
        # gene has repeat domain
        if geneID in repeatDict:
            # gene has non-repeat domain
            if geneID in nonRepeatDict:
                if geneID not in genesWithBoth:
                    genesWithBoth[geneID] = []
                for repAccessionID,repAccessionDesc,repAlignStart,repAlignStop,repGeneLen in repeatDict[geneID]:
                    genesWithBoth[geneID].append((repAccessionID,repAccessionDesc,repAlignStart,repAlignStop,repGeneLen))
                for nonRepAccessionID,nonRepAccessionDesc,nonRepAlignStart,nonRepAlignStop,nonRepGeneLen in nonRepeatDict[geneID]:
                    genesWithBoth[geneID].append((nonRepAccessionID,nonRepAccessionDesc,nonRepAlignStart,nonRepAlignStop,nonRepGeneLen))
            # gene only has repeat domain
            else:
                if geneID not in genesWithRepeatOnly:
                    genesWithRepeatOnly[geneID] = []
                for repAccessionID,repAccessionDesc,repAlignStart,repAlignStop,repGeneLen in repeatDict[geneID]:
                    genesWithRepeatOnly[geneID].append((repAccessionID,repAccessionDesc,repAlignStart,repAlignStop,repGeneLen))
        # gene does not have repeat domain
        else:
            # gene has non-repeat domain and no repeat domain
            if geneID in nonRepeatDict:
                if geneID not in genesWithNonRepeatOnly:
                    genesWithNonRepeatOnly[geneID] = []
                for nonRepAccessionID,nonRepAccessionDesc,nonRepAlignStart,nonRepAlignStop,nonRepGeneLen in nonRepeatDict[geneID]:
                    genesWithNonRepeatOnly[geneID].append((nonRepAccessionID,nonRepAccessionDesc,nonRepAlignStart,nonRepAlignStop,nonRepGeneLen))
            # gene has no pfam annotation
            else:
                if geneID not in genesWithoutPfam:
                    genesWithoutPfam[geneID] = 1
    return(genesWithBoth,genesWithRepeatOnly,genesWithNonRepeatOnly,genesWithoutPfam)


def createCountArray(geneLen):
    countArray = [0]*geneLen
    return(countArray)


def updateArray(countArray,start,stop,posMarker):
    for i in range(start,stop):
        countArray[i] = posMarker
    return(countArray)


def assessAndReassignGenesWithBoth(genesWithBoth,nonRepeatDict,repeatDict,geneLenDict,repeatCoverageThreshold,conditionalRepeatGenes):
    OUT_reassignedToRepeat_aboveRepCovThreshold = open('reassignedToRepeat_aboveRepCovThreshold.txt','w')
    reassignedToRepeat = {}
    isConditional = {}
    isNotConditional = {}
    repPosMarker = 1
    for geneID in genesWithBoth:
        # looping through "nonRepeatDict" will retrieve only the domains that are non-repeat
        for nonRepAccessionID,nonRepAccessionDesc,nonRepAlignStart,nonRepAlignStop,nonRepGeneLen in nonRepeatDict[geneID]:
            # check to see if non-repeat domain is in the conditional dict
            if nonRepAccessionID in conditionalRepeatGenes:
                # print geneID,nonRepAccessionID,nonRepAccessionDesc
                if geneID not in isConditional:
                    isConditional[geneID] = 1
            # here, check if non-repeat domain is NOT in conditional dict
            # the idea is that if at least one of the domains is not in the conditional dict, the gene will not be re-categorized as a repeat
            else:
                if geneID not in isNotConditional:
                    isNotConditional[geneID] = 1
    # loop through "genesWithBoth" again...
    # here, the idea is to check the coverage of repeat domains in cases where genes have both repeat- and non-repeat-associated domains
    for geneID in genesWithBoth:
        # if geneID is in conditional dict, may have some repeat assocation
        if geneID in isConditional:
            # but gene ID is NOT in "isNotConditional"
            # this means that none of the genes were conditional
            if geneID not in isNotConditional:
                # then add the geneID to the reassignment dict...
                if geneID not in reassignedToRepeat:
                    reassignedToRepeat[geneID] = 1
            else:
                # if genes have both repeat and non repeat, check the coverage of the repeat domains
                # initialize count array here
                geneLen = geneLenDict[geneID]
                countArray = createCountArray(geneLen)
                # print geneID,geneLen,countArray
                for repAccessionID,repAccessionDesc,repAlignStart,repAlignStop,repGeneLen in repeatDict[geneID]:
                    # print geneID,geneLen,repAccessionID,repAccessionDesc,repAlignStart,repAlignStop,repGeneLen
                    countArray = updateArray(countArray,repAlignStart,repAlignStop,repPosMarker)
                    # print geneID,geneLen,countArray
                totalRepPosCount = countArray.count(repPosMarker)
                # print geneID,geneLen,totalRepPosCount
                repCoverage = float(totalRepPosCount)/geneLen
                if repCoverage > repeatCoverageThreshold:
                    if geneID not in reassignedToRepeat:
                        reassignedToRepeat[geneID] = repCoverage
                    OUT_reassignedToRepeat_aboveRepCovThreshold.write("%s\t%s\n" % (geneID,repCoverage))
                    for repAccessionID,repAccessionDesc,repAlignStart,repAlignStop,repGeneLen in repeatDict[geneID]:
                        OUT_reassignedToRepeat_aboveRepCovThreshold.write("\t%s\t%s\n" % (repAccessionID,repAccessionDesc))
                    for nonRepAccessionID,nonRepAccessionDesc,nonRepAlignStart,nonRepAlignStop,nonRepGeneLen in nonRepeatDict[geneID]:
                        OUT_reassignedToRepeat_aboveRepCovThreshold.write("\t%s\t%s\n" % (nonRepAccessionID,nonRepAccessionDesc))
    return(reassignedToRepeat)


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hmmscan domtblout file> <repeat domain file> <gene fasta file>\n"
if len(sys.argv) != 4:
    print(usage)
    sys.exit()


hmmFile = sys.argv[1]
repeatDomainFile = sys.argv[2]
geneFasta = sys.argv[3]

eValueThreshold = 1e-3
# https://phytozome.jgi.doe.gov/pz/portal.html#!info?alias=Org_Ptrichocarpa
# began with threshold 0.3; found that it could be more stringent without losing any interesting gene examples
#repeatCoverageThreshold = 0.2
repeatCoverageThreshold = 0.3

# original conditionalRepeatGenes dict -- some of these examples were re-assigned to the primary repeat list
# conditionalRepeatGenes = {'PF04434.18':1, 'PF16588.6':1, 'PF02892.16':1, 'PF13696.7':1, 'PF00098.24':1, 'PF13917.7':1, 'PF14392.7':1, 'PF14111.7':1, 'PF14372.7':1, 'PF09331.12':1, 'PF13962.7':1, 'PF00385.25':1, 'PF15136.7':1, 'PF03372.24':1, 'PF02902.20':1, 'PF18907.1':1, 'PF00574.24':1, 'PF12353.9':1, 'PF11604.9':1, 'PF13634.7':1, 'PF00230.21':1, 'PF05193.22':1}

conditionalRepeatGenes = {'PF04434.18':1, 'PF16588.6':1, 'PF02892.16':1, 'PF13917.7':1, 'PF14372.7':1, 'PF09331.12':1, 'PF13962.7':1, 'PF00385.25':1, 'PF02902.20':1, 'PF18907.1':1, 'PF00574.24':1, 'PF11604.9':1, 'PF13634.7':1, 'PF00230.21':1, 'PF05193.22':1}

# these domains co-occur with repeat genes, and it seems like some genes with strong repeat association are included in the "both" set because of the presence of these domains - however, in the "non-repeat" set, they often occur by themselves, where they are the only domain - in this case, their function is not well-defined
# Added 05/11/2021 // PF13634.7 // Nucleoporin FG repeat region // HUMLU_CAS0003036.t1
# Added 05/11/2021 // PF00230.21 // Major intrinsic protein // HUMLU_CAS0048483.t1; HUMLU_CAS0032823.t1
# Added 05/11/2021 // PF05193.22 // Peptidase M16 inactive domain // HUMLU_CAS0006187.t1
# Added 05/11/2021 // PF11604.9 // Copper binding periplasmic protein CusF // HUMLU_CAS0062539.t1; HUMLU_CAS0050529.t1
# Added 05/11/2021 // PF12353.9 // Eukaryotic translation initiation factor 3 subunit G // strong repeat association 
# Added 05/06/2021 // PF15136.7 // Uncharacterised protein family UPF0449 // HUMLU_CAS0006722.t1
# Added 05/06/2021 // PF03372.24 // Endonuclease/Exonuclease/phosphatase family
# Added 05/06/2021 // PF02902.20 // Ulp1 protease family, C-terminal catalytic domain
# Added 05/06/2021 // PF18907.1 // Family of unknown function (DUF5662); (HUMLU_CAS0014859.t1)
# Added 05/06/2021 // PF00574.24 // Clp protease
# SWIM zn finger --> https://pubmed.ncbi.nlm.nih.gov/12151216/
# C2H2 zn finger --> https://www.nature.com/articles/nbt.3128
# BED zn finger  --> https://pfam.xfam.org/family/PF02892.16
# Added 04/30/2021 // PF00385.25      Chromo (CHRromatin Organisation MOdifier) domain
# PF13696.7 --> co-occurs with repeat domains
# PF00098.24 --> http://pfam.xfam.org/family/zf-CCHC "The motifs are mostly from retroviral gag proteins (nucleocapsid). Prototype structure is from HIV. Also contains members involved in eukaryotic gene regulation, such as C. elegans GLH-1."
# PF13917.7 --> http://pfam.xfam.org/family/PF13917.7 "The motifs are mostly from retroviral gag proteins (nucleocapsid). Prototype structure is from HIV. Also contains members involved in eukaryotic gene regulation, such as C. elegans GLH-1."
### Added back to repeat list 04/29/2021 # PF14291.7 Domain of unknown function (DUF4371) --> (seems to occur by itself most of the time and has nr similarity to "zinc finger MYM-type protein 1-like [Cannabis sativa]": HUMLU_CAS0051023.t1, HUMLU_CAS0041690.t1, HUMLU_CAS0011021.t1; occasionally co-occurs with other repeat: HUMLU_CAS0000341.t1)
### Added back to repeat list 04/29/2021 # PF13952.7 Domain of unknown function (DUF4216) --> another DUF that frequently occurs by itself and can't clearly be called repeat-associated (e.g. HUMLU_CAS0046319.t1, HUMLU_CAS0046312.t1, HUMLU_CAS0043760.t1; does co-occur with repeat, e.g. HUMLU_CAS0046358.t1)
### Added back to repeat list 04/29/2021 # PF04937.16 Protein of unknown function (DUF 659) --> http://pfam.xfam.org/family/PF04937.16 "Transposase-like protein with no known function." (occurs by itself: HUMLU_CAS0006586.t1, HUMLU_CAS0007618.t1, HUMLU_CAS0009828.t1, HUMLU_CAS0003849.t1; occurs with other repeat domain: HUMLU_CAS0016894.t1, HUMLU_CAS0045125.t1, HUMLU_CAS0006890.t1)
# PF14111.7 Domain of unknown function (DUF4283) --> "LINEs between Species: Evolutionary Dynamics of LINE-1 Retrotransposons across the Eukaryotic Tree of Life" https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5203782/ (this domain sometimes occurs by itself, e.g. HUMLU_CAS0066255.t1, HUMLU_CAS0066354.t1, HUMLU_CAS0046369.t1, HUMLU_CAS0066615.t1, and sometimes co-occurs with repeat domains, e.g. HUMLU_CAS0046263.t1, HUMLU_CAS0046338.t1, HUMLU_CAS0049810.t1, HUMLU_CAS0043927.t1, HUMLU_CAS0044001.t1) 
### Added back to repeat list 04/29/2021 # PF13960.7 Domain of unknown function (DUF4218) --> (this domain occurs by itself, e.g. HUMLU_CAS0066369.t1, HUMLU_CAS0066238.t1, HUMLU_CAS0066375.t1, HUMLU_CAS0066428.t1 and also co-occurs with repeat domains, e.g. HUMLU_CAS0066402.t1, HUMLU_CAS0066421.t1, HUMLU_CAS0066502.t1)
# PF14372.7 Domain of unknown function (DUF4413) --> (conditional repeat based on association observed with HUMLU_CAS0025411.t1; co-occuring repeat domains)
# PF09331.12 Domain of unknown function (DUF1985) --> co-occuring with repeat domains in HUMLU_CAS0039176.t1 and HUMLU_CAS0056123.t1 
# PF13962.7 Domain of unknown function --> HUMLU_CAS0016430.t1 shows repeat association but normally this domain co-occurs with Ankyrin repeats...
# e-value thresholds --> https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2203978/


geneIDs = readGeneFasta(geneFasta)
repeatDomainDict = readRepeatDomainFile(repeatDomainFile)
nonRepeatDict,repeatDict,domainCountDict,geneLenDict = readHMMFile(hmmFile,eValueThreshold,repeatDomainDict)

genesWithBoth,genesWithRepeatOnly,genesWithNonRepeatOnly,genesWithoutPfam = organizeGenes(geneIDs,nonRepeatDict,repeatDict)

reassignedToRepeat = assessAndReassignGenesWithBoth(genesWithBoth,nonRepeatDict,repeatDict,geneLenDict,repeatCoverageThreshold,conditionalRepeatGenes)

OUT_genesWithBoth1 = open('genesWithBoth.txt','w')
OUT_genesWithBoth2 = open('genesWithBoth.tsv','w')
OUT_genesWithRepeatOnly1 = open('genesWithRepeatOnly.txt','w')
OUT_genesWithRepeatOnly2 = open('genesWithRepeatOnly.tsv','w')
OUT_genesWithNonRepeatOnly1 = open('genesWithNonRepeatOnly.txt','w')
OUT_genesWithNonRepeatOnly2 = open('genesWithNonRepeatOnly.tsv','w')
OUT_reassignedToRepeat1 = open('reassignedToRepeat.txt','w')
OUT_reassignedToRepeat2 = open('reassignedToRepeat.tsv','w')
OUT_genesWithoutPfam = open('genesWithoutPfam.txt','w')

# nonRepAlignStart,nonRepAlignStop,nonRepGeneLen
for geneID1 in genesWithBoth:
    # these genes don't have "conditional" repeats
    if geneID1 not in reassignedToRepeat:
        OUT_genesWithBoth1.write("%s\n" % (geneID1))
        #if "PF04434" in genesWithBoth[geneID1]:
        for accID1,accDesc1,alignStart1,alignStop1,geneLen1 in genesWithBoth[geneID1]:
            #if "PF04434" in accID1:
            OUT_genesWithBoth1.write("\t%s\t%s\n" % (accID1,accDesc1))
            OUT_genesWithBoth2.write("%s\t%s\t%s\n" % (geneID1,accID1,accDesc1))
    # these genes have "conditional" repeats...
    else:
        OUT_reassignedToRepeat1.write("%s\n" % (geneID1))
        for accID1,accDesc1,alignStart1,alignStop1,geneLen1 in genesWithBoth[geneID1]:
            OUT_reassignedToRepeat1.write("\t%s\t%s\n" % (accID1,accDesc1))
            OUT_reassignedToRepeat2.write("%s\t%s\t%s\n" % (geneID1,accID1,accDesc1))


for geneID2 in genesWithRepeatOnly:
    OUT_genesWithRepeatOnly1.write("%s\n" % (geneID2))
    for accID2,accDesc2,alignStart2,alignStop2,geneLen2 in genesWithRepeatOnly[geneID2]:
        OUT_genesWithRepeatOnly1.write("\t%s\t%s\n" % (accID2,accDesc2))
        OUT_genesWithRepeatOnly2.write("%s\t%s\t%s\n" % (geneID2,accID2,accDesc2))


for geneID3 in genesWithNonRepeatOnly:
    OUT_genesWithNonRepeatOnly1.write("%s\n" % (geneID3))
    for accID3,accDesc3,alignStart3,alignStop3,geneLen3 in genesWithNonRepeatOnly[geneID3]:
        OUT_genesWithNonRepeatOnly1.write("\t%s\t%s\n" % (accID3,accDesc3))
        OUT_genesWithNonRepeatOnly2.write("%s\t%s\t%s\n" % (geneID3,accID3,accDesc3))


for geneID4 in genesWithoutPfam:
    OUT_genesWithoutPfam.write("%s\n" % (geneID4))


number_of_transcripts_with_domain = len(domainCountDict)
total_number_of_transcripts = len(geneIDs)
percent_transcripts_with_domains = float(number_of_transcripts_with_domain) / total_number_of_transcripts * 100
print("number_of_transcripts_with_domain\t%s\t" % (number_of_transcripts_with_domain))
print("total_number_of_transcripts\t%s\t" % (total_number_of_transcripts))
print("percent_transcripts_with_domains\t%s\t" % (percent_transcripts_with_domains))
