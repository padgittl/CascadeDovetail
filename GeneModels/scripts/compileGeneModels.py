import os,re,sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


# this subroutine is for the purpose of storing transcript lengths to verify with lengths in transdecoder bed file
def readStringtieTranscriptGFF(stringtieTranscriptGFF):
    stringtieTranscriptLengthDict = {}
    with open(stringtieTranscriptGFF,'r') as S:
        for line in S:
            if not line.startswith('#'):
                if not line.isspace():
                    scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                    transcriptInfo = attribute.split(',')
                    getTranscriptID = re.search('(MSTRG\.\d+\.\d+)',transcriptInfo[2])
                    transcriptID = getTranscriptID.group(1)
                    #print transcriptID
                    start = int(start)
                    end = int(end)
                    matchLen = end - start + 1
                    if transcriptID not in stringtieTranscriptLengthDict:
                        stringtieTranscriptLengthDict[transcriptID] = 0
                    stringtieTranscriptLengthDict[transcriptID] += matchLen
        return(stringtieTranscriptLengthDict)


def readTransdecoderFasta(transdecoderFastaFile):
    transdecoderLengthDict = {}
    for record in SeqIO.parse(transdecoderFastaFile,"fasta"):
        if record.id not in transdecoderLengthDict:
            transdecoderLengthDict[record.id] = len(record.seq)
    return(transdecoderLengthDict)
        

# this subroutine is where the business happens - sorting transcript and orf lengths
def readTransdecoderBedFile(transdecoderBedFile,stringtieTranscriptLengthDict,transdecoderCDSLengthDict,transdecoderProteinLengthDict):
    # D[geneID].append((transcriptID1, transcriptLength1, orfID, orfLength))
    # D[geneID].sort(key = lambda x:x[3], reverse=True)
    geneDict = {}
    transdecoderTranscriptCount = 0
    OUT_orfAndCDSDifferentLengths = open('orfsAndCDSWithDifferentLengths.txt','w')
    OUT_orfAndCDSDifferentLengths.write("transcriptID\torfLen\tcdsLen\tcdsLenDividedByCodonLen\tproteinLen\torfInfo\n")
    with open(transdecoderBedFile,'r') as B:
        for line in B:
            if 'track name' not in line:
                line = line.strip().split('\t')
                transcriptID = line[0]
                getGeneID = re.search('(MSTRG\.\d+)',transcriptID)
                geneID = getGeneID.group(1)
                if geneID not in geneDict:
                    geneDict[geneID] = []
                transcriptLength = line[2]
                transcriptLength = int(transcriptLength)
                attribute = line[3:]
                joinedAttribute = ' '.join(attribute)
                # verify that transcript length is actually what is being collected
                if transcriptID in stringtieTranscriptLengthDict:
                    transcriptLen = stringtieTranscriptLengthDict[transcriptID]
                    if transcriptLen != transcriptLength:
                        print transcriptID,transcriptLength,transcriptLen
                        sys.exit()
                    else:
                        splitAttribute = joinedAttribute.split(';')
                        prefix,orfID = splitAttribute[0].split('ID=')
                        orfInfo = splitAttribute[2].split(',')
                        # ORF_type:complete_len:487_(+)
                        if 'complete_len' in orfInfo[0]:
                            getOrfLen = re.search('ORF_type:.+:(\d+)_',orfInfo[0])
                            orfLen = getOrfLen.group(1)
                            orfLen = int(orfLen)
                            if orfID in transdecoderCDSLengthDict:
                                transdecoderTranscriptCount += 1
                                cdsLen = transdecoderCDSLengthDict[orfID]
                                cdsLenDividedByCodonLen = float(cdsLen) / 3
                                proteinLen = transdecoderProteinLengthDict[orfID]
                                geneDict[geneID].append((transcriptID,transcriptLen,orfID,orfLen,cdsLen,proteinLen))
                                if orfLen != cdsLenDividedByCodonLen:
                                    OUT_orfAndCDSDifferentLengths.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (transcriptID,orfLen,cdsLen,cdsLenDividedByCodonLen,proteinLen,orfInfo[0]))
                                    #if orfLen > cdsLenDividedByCodonLen:
                                    #print transcriptID + ", lengths not equal, cds: " + str(cdsLenDividedByCodonLen) + " orf: " + str(orfLen) + ", " + orfInfo[0]
                                    #print "cds and orf lengths not equal: " + transcriptID,orfLen,cdsLen
                                    #sys.exit()
                            else:
                                print "orfID not in transdecoderCDSLengthDict: "  + orfID
                                sys.exit()
                else:
                    print "transcriptID not in bed file: " + transcriptID
                    sys.exit()

    transdecoderLongestCDSDict = {}
    OUT_examplesWhereLongestTranscriptNotLongestCDS = open('longestTranscriptNotLongestCDS.txt','w')
    transdecoderLongestCDSCount = 0
    # OUT_examplesWhereLongestTranscriptNotLongestCDS.write("# lo = longest orf; lcds = longest cds; lt = longest transcript\n")
    OUT_examplesWhereLongestTranscriptNotLongestCDS.write("# lcds = longest cds; lt = longest transcript\n")
    OUT_examplesWhereLongestTranscriptNotLongestCDS.write("# lcdsTranscriptID\tltTranscriptID\tlcdsTranscriptLen\tltTranscriptLen\tlcdsORFLen\tltORFLen\n")
    for geneID in geneDict:
        if len(geneDict[geneID]) > 0:
            geneDict[geneID].sort(key=lambda x:x[4], reverse=True)
            # longest cds (lcds)
            lcdsTranscriptID,lcdsTranscriptLen,lcdsORFID,lcdsORFLen,lcdsCDSLen,lcdsProteinLen = geneDict[geneID][0]
            if lcdsTranscriptID not in transdecoderLongestCDSDict:
                transdecoderLongestCDSCount += 1
                transdecoderLongestCDSDict[lcdsORFID] = (geneID,lcdsTranscriptID,lcdsTranscriptLen,lcdsORFLen,lcdsCDSLen,lcdsProteinLen)
            # longest transcript (lt)
            geneDict[geneID].sort(key=lambda x:x[1], reverse=True)
            ltTranscriptID,ltTranscriptLen,ltORFID,ltORFLen,ltCDSLen,ltProteinLen = geneDict[geneID][0]
            #print "longest transcript: " + ltTranscriptID,ltTranscriptLen,ltORFID,ltORFLen,ltCDSLen,ltProteinLen
            if lcdsTranscriptID != ltTranscriptID:
                #print loTranscriptID,ltTranscriptID,loTranscriptLen,ltTranscriptLen,loORFLen,ltORFLen
                OUT_examplesWhereLongestTranscriptNotLongestCDS.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (lcdsTranscriptID,ltTranscriptID,lcdsTranscriptLen,ltTranscriptLen,lcdsORFLen,ltORFLen))
    print "full transcript count (complete_len): " + str(transdecoderTranscriptCount)
    print "longest CDS per transcript count (complete_len): " + str(transdecoderLongestCDSCount)
    return(transdecoderLongestCDSDict)


# read in scaffold lengths and also collect largest scaffolds
def readLengthsFile(lengthsFile,maxNumberLongestScaffolds):
    lengthsDict = {}
    lengthsList = []
    topLengthScaffoldDict = {}
    with open(lengthsFile,'r') as L:
        for line in L:
            scaffoldID,scaffoldLength = line.strip().split()
            lengthsDict[scaffoldID] = int(scaffoldLength)
            lengthsList.append((scaffoldID,int(scaffoldLength)))
    # sort by descending length
    lengthsList.sort(key=lambda x:x[1], reverse=True)
    longestLengthScaffolds = lengthsList[:10]
    for topSetScaffoldID,topSetScaffoldLen in longestLengthScaffolds:
        topLengthScaffoldDict[topSetScaffoldID] = topSetScaffoldLen
    return(lengthsDict,lengthsList,topLengthScaffoldDict)


# get IDs with protein length at least 100 aa
def filterProteinsByLength(proteinFastaFile,minLen):
    filterProteinsByLengthDict = {}
    for record in SeqIO.parse(proteinFastaFile,'fasta'):
        if len(record.seq) >= minLen:
            if record.id not in filterProteinsByLengthDict:
                filterProteinsByLengthDict[record.id] = len(record.seq)
    return(filterProteinsByLengthDict)


def getLongestMakerTranscript(gffFile):
    geneDict = {}
    with open(gffFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                if not line.isspace():
                    scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                    start = int(start)
                    end = int(end)
                    length = end - start + 1
                    if feature == 'mRNA':
                        getTranscriptID = re.search('ID=(.+);Parent=',attribute)
                        transcriptID = getTranscriptID.group(1)
                        # MAKER0067600.t1
                        geneID,suffix = transcriptID.split('.t')
                        if geneID not in geneDict:
                            geneDict[geneID] = []
                        geneDict[geneID].append((transcriptID,start,end,length))
    makerLongestTranscriptPerGeneDict = {}
    for geneID in geneDict:
        geneDict[geneID].sort(key=lambda x:x[3], reverse=True)
        transcriptID,start,end,mRNALength = geneDict[geneID][0]
        if transcriptID not in makerLongestTranscriptPerGeneDict:
            makerLongestTranscriptPerGeneDict[transcriptID] = (start,end,mRNALength)
    return(makerLongestTranscriptPerGeneDict)


# read in full gff, containing all gene models for maker and transdecoder individually - no filtering
def readGFF(gffFile):
    gffDict = {}
    mRNACoordDict = {}
    geneIDs = {}
    cdsLengthDict = {}
    with open(gffFile,'r') as F:
        for line in F:
            if not line.startswith('#'):
                if not line.isspace():
                    scaffoldID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                    # adding the strand to the scaffold ID is for overlapping cases on opposite strands - allow this type of overlap to occur
                    newScaffoldID = scaffoldID + ";" + strand
                    if feature == 'gene':
                        getGeneID = re.search('ID=(.+);Name=',attribute)
                        geneID = getGeneID.group(1)
                        if geneID not in geneIDs:
                            geneIDs[geneID] = (int(start),int(end)) 
                    if feature == 'mRNA':
                        getTranscriptID = re.search('ID=(.+);Parent=',attribute)
                        transcriptID = getTranscriptID.group(1)
                        if newScaffoldID not in mRNACoordDict:
                            mRNACoordDict[newScaffoldID] = []
                        mRNACoordDict[newScaffoldID].append((transcriptID,int(start),int(end)))
                    
                    if feature == 'CDS':
                        start = int(start)
                        end = int(end)
                        cdsExonLen = end - start + 1
                        # 09/27/2021 -->
                        # this used to be HUMLU_CAS, but was switched to MAKER to prevent confusion
                        if 'MAKER' in attribute:
                            getTranscriptID = re.search('Parent=(.+);',attribute)
                            transcriptID = getTranscriptID.group(1)
                        elif 'MSTRG' in attribute:
                            getTranscriptID = re.search('Parent=(.+)',attribute)
                            transcriptID = getTranscriptID.group(1)
                        else:
                            print 'gff parsing error'
                            sys.exit()
                        # 2D dictionary
                        if newScaffoldID not in gffDict:
                            gffDict[newScaffoldID] = {}
                        if transcriptID not in gffDict[newScaffoldID]:
                            gffDict[newScaffoldID][transcriptID] = []
                        gffDict[newScaffoldID][transcriptID].append((start,end,cdsExonLen))
                        # get concatenated cds exon lengths
                        if transcriptID not in cdsLengthDict:
                            cdsLengthDict[transcriptID] = 0
                        cdsLengthDict[transcriptID] += cdsExonLen
    
    cdsCoordDict = {}
    sortedCDSDict = {}
    originalScaffoldIDSortedCDSDict = {}
    for newScaffoldID in gffDict:
        # also store this data to dict with original scaffoldID, this is important for re-naming genes in later subroutine 
        oldScaffoldID,strand = newScaffoldID.split(';')
        if oldScaffoldID not in originalScaffoldIDSortedCDSDict:
            originalScaffoldIDSortedCDSDict[oldScaffoldID] = []
        # store newScaffoldIDs to dictionaries
        if newScaffoldID not in cdsCoordDict:
            cdsCoordDict[newScaffoldID] = {}
        if newScaffoldID not in sortedCDSDict:
            sortedCDSDict[newScaffoldID] = []
        # get first and last coding exon positions, sort by starting position in ascending order
        for transcriptID in gffDict[newScaffoldID]:
            gffDict[newScaffoldID][transcriptID].sort(key=lambda x:x[0], reverse=False)
            firstCDSStart = gffDict[newScaffoldID][transcriptID][0][0]
            lastCDSStop = gffDict[newScaffoldID][transcriptID][-1][1]
            cdsLength = cdsLengthDict[transcriptID]
            # this dict is keyed only by scaffold ID and is intended to sort transcript IDs by length
            sortedCDSDict[newScaffoldID].append((transcriptID,firstCDSStart,lastCDSStop,cdsLength))
            originalScaffoldIDSortedCDSDict[oldScaffoldID].append((transcriptID,firstCDSStart,lastCDSStop,cdsLength))
            # this dict is keyed by scaffold ID and transcriptID, to store all coding exons in a list
            if transcriptID not in cdsCoordDict[newScaffoldID]:
                cdsCoordDict[newScaffoldID][transcriptID] = []
            # loop through dict containing all coding exons
            for cdsExonStart,cdsExonEnd,cdsExonLen in gffDict[newScaffoldID][transcriptID]:
                cdsCoordDict[newScaffoldID][transcriptID].append((cdsExonStart,cdsExonEnd))
            cdsCoordDict[newScaffoldID][transcriptID].sort(key=lambda x:x[0], reverse=False)
        # sort by scaffoldID, cds length in descending order
        sortedCDSDict[newScaffoldID].sort(key=lambda x:x[3], reverse=True)
        originalScaffoldIDSortedCDSDict[oldScaffoldID].sort(key=lambda x:x[3], reverse=True)
    return(cdsCoordDict,mRNACoordDict,geneIDs,cdsLengthDict,sortedCDSDict,originalScaffoldIDSortedCDSDict)


def createCountArray(scaffoldLen):
    countArray = [0]*scaffoldLen
    return(countArray)


def updateCountArray(countArray,start,stop):
    # the max value corresponds to the maximum position count per coding exon
    maxValue = 0
    # countArray[stop+1] += 1
    for i in range(start,stop):
        countArray[i] += 1
        if maxValue < countArray[i]:
            maxValue = countArray[i]
    return(countArray,maxValue)


# get maker genes that are not overlapping with transdecoder
def getNonOverlappingMakerGenes(makerCDSDict,transdecoderCDSDict,lengthsDict,sortedTransdecoderCDSDict,sortedMakerCDSDict,transdecoderLongestCDSDict,makerLongestTranscriptPerGeneDict):
    makerGenesToKeep = {}
    transdecoderGenesToKeep = {}
    # create count array, giving priority to transdecoder gene models
    for scaffoldID in sortedTransdecoderCDSDict:
        oldScaffoldID,strand = scaffoldID.split(';')
        scaffoldLen = lengthsDict[oldScaffoldID]
        countArray = createCountArray(scaffoldLen)
        # this dictionary is ordered by the cds start position in ascending order
        for transcriptID,firstCDSStart,lastCDSStop,cdsLength in sortedTransdecoderCDSDict[scaffoldID]:
            if transcriptID in transdecoderLongestCDSDict:
                if transcriptID not in transdecoderGenesToKeep:
                    transdecoderGenesToKeep[transcriptID] = (scaffoldID,firstCDSStart,lastCDSStop,cdsLength)
                # this dict is not the sorted-dict, it's keyed by scaffold and transcriptIDs and contains all coding exons in a list
                for cdsExonStart,cdsExonStop in transdecoderCDSDict[scaffoldID][transcriptID]:
                    countArray,maxCount = updateCountArray(countArray,cdsExonStart,cdsExonStop)
        # check if scaffold ID is in maker set
        if scaffoldID in sortedMakerCDSDict:
            for mTranscriptID,mFirstCDSStart,mLastCDSStop,mCDSLength in sortedMakerCDSDict[scaffoldID]:
                if mTranscriptID in makerLongestTranscriptPerGeneDict:
                    cdsExonMaxCountList = []
                    for mCDSExonStart,mCDSExonStop in makerCDSDict[scaffoldID][mTranscriptID]:
                        #print scaffoldID,mTranscriptID,mCDSExonStart,mCDSExonStop,mCDSLen
                        countArray,maxCount = updateCountArray(countArray,mCDSExonStart,mCDSExonStop)
                        cdsExonMaxCountList.append(maxCount)
                    cdsExonMaxCountSet = set(cdsExonMaxCountList)
                    #print mTranscriptID,cdsExonMaxCountSet 
                    if len(cdsExonMaxCountSet) == 1:
                        if list(cdsExonMaxCountSet)[0] == 1:
                            if mTranscriptID not in makerGenesToKeep:
                                makerGenesToKeep[mTranscriptID] = (scaffoldID,mFirstCDSStart,mLastCDSStop,mCDSLength)

    # also want to check cases where scaffold isn't present among transdecoder genes
    for scaffoldID in sortedMakerCDSDict:
        # check here
        if scaffoldID not in sortedTransdecoderCDSDict:
            oldScaffoldID,strand = scaffoldID.split(';')
            scaffoldLen = lengthsDict[oldScaffoldID]
            countArray = createCountArray(scaffoldLen)
            if scaffoldID not in makerGenesToKeep:
                makerGenesToKeep[scaffoldID] = {}
            for mTranscriptID,mFirstCDSStart,mLastCDSStop,mCDSLength in sortedMakerCDSDict[scaffoldID]: 
                if mTranscriptID in makerLongestTranscriptPerGeneDict:
                    cdsExonMaxCountList = []
                    for mCDSExonStart,mCDSExonStop in makerCDSDict[scaffoldID][mTranscriptID]:
                        #print scaffoldID,mTranscriptID,mCDSExonStart,mCDSExonStop,mCDSLen
                        countArray,maxCount = updateCountArray(countArray,mCDSExonStart,mCDSExonStop)
                        cdsExonMaxCountList.append(maxCount)
                    cdsExonMaxCountSet = set(cdsExonMaxCountList)
                    #print mTranscriptID,cdsExonMaxCountSet
                    if len(cdsExonMaxCountSet) == 1:
                        if list(cdsExonMaxCountSet)[0] == 1:
                            if mTranscriptID not in makerGenesToKeep:
                                makerGenesToKeep[mTranscriptID] = (scaffoldID,mFirstCDSStart,mLastCDSStop,mCDSLength)
    return(makerGenesToKeep,transdecoderGenesToKeep)


# want to create a mapping for every gene model to the updated naming scheme, not just non-overlapping or repeat-filtered examples
def renameGenes(originalScaffoldIDSortedMakerCDSDict,originalScaffoldIDSortedTransdecoderCDSDict,makerGeneIDs,transdecoderGeneIDs,lengthsList):
    compiledGeneDict = {}
    compileGenesForGeneCount = {}
    # add all genes from maker and transdecoder to one dict
    for scaffoldID,scaffoldLength in lengthsList:
        if scaffoldID not in compiledGeneDict:
            compiledGeneDict[scaffoldID] = []
        if scaffoldID not in compileGenesForGeneCount:
            compileGenesForGeneCount[scaffoldID] = []
        if scaffoldID in originalScaffoldIDSortedTransdecoderCDSDict:
            for transcriptID,cdsStart,cdsStop,cdsLength in originalScaffoldIDSortedTransdecoderCDSDict[scaffoldID]:
                geneInfo = transcriptID.split('.')
                geneID = geneInfo[0] + '.' + geneInfo[1]
                geneStart,geneStop = transdecoderGeneIDs[geneID]
                compileGenesForGeneCount[scaffoldID].append((geneID,geneStart))
                compiledGeneDict[scaffoldID].append((transcriptID,cdsStart,cdsStop,geneStart,geneStop))
        if scaffoldID in originalScaffoldIDSortedMakerCDSDict:
            for makerTranscriptID,makerGeneCDSStart,makerGeneCDSStop,makerCDSLen in originalScaffoldIDSortedMakerCDSDict[scaffoldID]:
                makerGeneID,suffix = makerTranscriptID.split('.t')
                maker_geneStart,maker_geneStop = makerGeneIDs[makerGeneID]
                compileGenesForGeneCount[scaffoldID].append((makerGeneID,maker_geneStart))
                compiledGeneDict[scaffoldID].append((makerTranscriptID,makerGeneCDSStart,makerGeneCDSStop,maker_geneStart,maker_geneStop))

    geneCount = 0
    geneCountDict = {}
    for scaffoldID,scaffoldLength in lengthsList:
        if scaffoldID in compileGenesForGeneCount:
            compileGenesForGeneCount[scaffoldID].sort(key=lambda x:x[1], reverse=False)
            for geneID,geneStart in compileGenesForGeneCount[scaffoldID]:
                #print geneID
                if geneID not in geneCountDict:
                    geneCount += 1
                    geneCountDict[geneID] = geneCount

    # do renaming here
    geneIDPrefix = 'HUMLU_CAS'
    offset = 7
    newGeneIDMappingDict = {}
    OUT = open('combinedGeneModels.txt','w')
    for scaffoldID,scaffoldLength in lengthsList:
        if scaffoldID in compiledGeneDict:
            # sort coords based on start pos, in ascending order
            # sorted according to cds start, but the sorting isn't required, because the geneCount is already done in the above subroutine - the suffix is from the existing transcriptID
            compiledGeneDict[scaffoldID].sort(key=lambda x:x[1], reverse=False)
            for transcriptID,cdsStart,cdsStop,geneStart,geneStop in compiledGeneDict[scaffoldID]:
                # add a transcript number identifier, like '.t1'
                # these are transdecoder genes
                if '.p' in transcriptID:
                    geneInfo = transcriptID.split('.')
                    # ['MSTRG', '865', '1', 'p1']
                    geneID = geneInfo[0] + '.' + geneInfo[1]
                    geneCountValue = geneCountDict[geneID]
                    prefix = geneInfo[0]
                    geneNumber = geneInfo[1]
                    transcriptNumber = geneInfo[2]
                    orfNumber = geneInfo[3]
                    newGeneIDSuffix = '.t' + transcriptNumber + '.' + orfNumber
                # these are maker genes
                else:
                    geneID,suffix = transcriptID.split('.t')
                    geneCountValue = geneCountDict[geneID]
                    oldGeneIDPrefix,geneIDSuffix = transcriptID.split('.t')
                    newGeneIDSuffix = '.t' + geneIDSuffix + '.p1'
                geneCountLength = len(str(geneCountValue))
                zeroPadding = offset - geneCountLength
                countWithZeroPad = str(geneCountValue).zfill(geneCountLength+zeroPadding)
                newGeneID = geneIDPrefix + countWithZeroPad + newGeneIDSuffix
                OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffoldID,transcriptID,newGeneID,geneStart,geneStop,cdsStart,cdsStop))
                #print transcriptID
                if transcriptID not in newGeneIDMappingDict:
                    newGeneIDMappingDict[transcriptID] = (newGeneID,scaffoldID)
    return(newGeneIDMappingDict)


        
# create pep and cds fasta files with updated gene IDs
def createFastaFiles(transdecoderProteinFastaFile,transdecoderCDSFastaFile,makerProteinFastaFile,makerCDSFastaFile,makerGenesToKeep,transdecoderGenesToKeep,newGeneIDMappingDict,transdecoderFilterProteinsByLengthDict,makerFilterProteinsByLengthDict,topLengthScaffoldDict):
    pepRecordDict = {}
    pepRecordList = []
    cdsRecordDict = {}
    cdsRecordList = []
    
    # largest 10 scaffolds
    top10ScafPepRecordDict = {}
    top10ScafPepRecordList = []
    top10ScafCDSRecordDict = {}
    top10ScafCDSRecordList = []

    # transdecoderGenesToKeep
    transdecoderOnlyPepRecordDict = {}
    transdecoderOnlyPepRecordList = []
    transdecoderOnlyCDSRecordDict = {}
    transdecoderOnlyCDSRecordList = []

    top10ScafTransdecoderOnlyPepRecordDict = {}
    top10ScafTransdecoderOnlyPepRecordList = []
    top10ScafTransdecoderOnlyCDSRecordDict = {}
    top10ScafTransdecoderOnlyCDSRecordList = []

    # transdecoder
    for record in SeqIO.parse(transdecoderProteinFastaFile,"fasta"):
        if record.id in newGeneIDMappingDict:
            #print record.id
            if record.id in transdecoderGenesToKeep:
                #print record.id
                if record.id in transdecoderFilterProteinsByLengthDict:
                    record.seq = record.seq.upper()
                    newTransdecoderRecord,scaffoldID = newGeneIDMappingDict[record.id]
                    record.id = record.id.replace(record.id,newTransdecoderRecord)
                    record.name = record.name.replace(record.name,'')
                    record.description = record.description.replace(record.description,'')
                    if record.id not in pepRecordDict:
                        pepRecordDict[record.id] = 1
                        pepRecordList.append(record)
                    if record.id not in transdecoderOnlyPepRecordDict:
                        transdecoderOnlyPepRecordDict[record.id] = 1
                        transdecoderOnlyPepRecordList.append(record)
                    # filter for longest scaffold lengths
                    if scaffoldID in topLengthScaffoldDict:
                        if record.id not in top10ScafPepRecordDict:
                            top10ScafPepRecordDict[record.id] = 1
                            top10ScafPepRecordList.append(record)
                        if record.id not in top10ScafTransdecoderOnlyPepRecordDict:
                            top10ScafTransdecoderOnlyPepRecordDict[record.id] = 1
                            top10ScafTransdecoderOnlyPepRecordList.append(record)

        else:
            print record.id + " not in newGeneIDMappingDict"
            sys.exit()

    for record in SeqIO.parse(transdecoderCDSFastaFile,"fasta"):
        if record.id in newGeneIDMappingDict:
            if record.id in transdecoderGenesToKeep:
                if record.id in transdecoderFilterProteinsByLengthDict:
                    record.seq = record.seq.upper()
                    newTransdecoderRecord,scaffoldID = newGeneIDMappingDict[record.id]
                    record.id = record.id.replace(record.id,newTransdecoderRecord)
                    record.name = record.name.replace(record.name,'')
                    record.description = record.description.replace(record.description,'')
                    if record.id not in cdsRecordDict:
                        cdsRecordDict[record.id] = 1
                        cdsRecordList.append(record)
                    if record.id not in transdecoderOnlyCDSRecordDict:
                        transdecoderOnlyCDSRecordDict[record.id] = 1
                        transdecoderOnlyCDSRecordList.append(record)
                    # filter for longest scaffold lengths
                    if scaffoldID in topLengthScaffoldDict:
                        if record.id not in top10ScafCDSRecordDict:
                            top10ScafCDSRecordDict[record.id] = 1
                            top10ScafCDSRecordList.append(record)
                        if record.id not in top10ScafTransdecoderOnlyCDSRecordDict:
                            top10ScafTransdecoderOnlyCDSRecordDict[record.id] = 1
                            top10ScafTransdecoderOnlyCDSRecordList.append(record)
                            
        else:
            print record.id + " not in newGeneIDMappingDict"
            sys.exit()

    # maker
    for record in SeqIO.parse(makerProteinFastaFile,"fasta"):
        if record.id in newGeneIDMappingDict:
            if record.id in makerGenesToKeep:
                if record.id in makerFilterProteinsByLengthDict:
                    record.seq = record.seq.upper()
                    newMakerRecord,scaffoldID = newGeneIDMappingDict[record.id]
                    record.id = record.id.replace(record.id,newMakerRecord)
                    record.name = record.name.replace(record.name,'')
                    record.description = record.description.replace(record.description,'')
                    if record.id not in pepRecordDict:
                        pepRecordDict[record.id] = 1
                        pepRecordList.append(record)
                    # filter for longest scaffold lengths
                    if scaffoldID in topLengthScaffoldDict:
                        if record.id not in top10ScafPepRecordDict:
                            top10ScafPepRecordDict[record.id] = 1
                            top10ScafPepRecordList.append(record)

        else:
            print record.id + " not in newGeneIDMappingDict"
            sys.exit()

    for record in SeqIO.parse(makerCDSFastaFile,"fasta"):
        if record.id in newGeneIDMappingDict:
            if record.id in makerGenesToKeep:
                if record.id in makerFilterProteinsByLengthDict:
                    record.seq = record.seq.upper()
                    newMakerRecord,scaffoldID = newGeneIDMappingDict[record.id]
                    record.id = record.id.replace(record.id,newMakerRecord)
                    record.name = record.name.replace(record.name,'')
                    record.description = record.description.replace(record.description,'')
                    if record.id not in cdsRecordDict:
                        cdsRecordDict[record.id] =1
                        cdsRecordList.append(record)
                    # filter for longest scaffold lengths
                    if scaffoldID in topLengthScaffoldDict:
                        if record.id not in top10ScafCDSRecordDict:
                            top10ScafCDSRecordDict[record.id] = 1
                            top10ScafCDSRecordList.append(record)

        else:
            print record.id + " not in newGeneIDMappingDict"
            sys.exit()

    # transdecoder + maker --> full assembly
    SeqIO.write(pepRecordList,"combinedGeneModels.fullAssembly.pep.fasta", "fasta")
    SeqIO.write(cdsRecordList,"combinedGeneModels.fullAssembly.cds.fasta", "fasta") 

    # transdecoder only --> full assembly
    SeqIO.write(transdecoderOnlyPepRecordList,"transdecoder.fullAssembly.pep.fasta", "fasta")
    SeqIO.write(transdecoderOnlyCDSRecordList,"transdecoder.fullAssembly.cds.fasta", "fasta")

    # transdecoder + maker --> ten scaffolds
    SeqIO.write(top10ScafPepRecordList,"combinedGeneModels.tenScaffolds.pep.fasta", "fasta")
    SeqIO.write(top10ScafCDSRecordList,"combinedGeneModels.tenScaffolds.cds.fasta", "fasta")

    # transdecoder only --> ten scaffolds
    SeqIO.write(top10ScafTransdecoderOnlyPepRecordList,"transdecoder.tenScaffolds.pep.fasta", "fasta")
    SeqIO.write(top10ScafTransdecoderOnlyCDSRecordList,"transdecoder.tenScaffolds.cds.fasta", "fasta")


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <transdecoder protein fasta> <transdecoder cds fasta> <full transdecoder gff> <maker protein fasta> <maker cds fasta> <full maker gff> <scaffold lengths file> <stringtie transcript gff> <transdecoder bed file> <minimum amino acid length>\n"
if len(sys.argv) != 11:
    print usage
    sys.exit()


transdecoderProteinFastaFile = sys.argv[1]
transdecoderCDSFastaFile = sys.argv[2]
fullTransdecoderGFF = sys.argv[3]
makerProteinFastaFile = sys.argv[4]
makerCDSFastaFile = sys.argv[5]
fullMakerGFF = sys.argv[6]
lengthsFile = sys.argv[7]
stringtieTranscriptGFF = sys.argv[8]
transdecoderBedFile = sys.argv[9]
minLen = sys.argv[10]

minLen = int(minLen)
maxNumberLongestScaffolds = 10

lengthsDict,lengthsList,topLengthScaffoldDict = readLengthsFile(lengthsFile,maxNumberLongestScaffolds)

stringtieTranscriptLengthDict = readStringtieTranscriptGFF(stringtieTranscriptGFF)
transdecoderCDSLengthDict = readTransdecoderFasta(transdecoderCDSFastaFile)
transdecoderProteinLengthDict = readTransdecoderFasta(transdecoderProteinFastaFile)

transdecoderLongestCDSDict = readTransdecoderBedFile(transdecoderBedFile,stringtieTranscriptLengthDict,transdecoderCDSLengthDict,transdecoderProteinLengthDict)
makerLongestTranscriptPerGeneDict = getLongestMakerTranscript(fullMakerGFF)

transdecoderFilterProteinsByLengthDict = filterProteinsByLength(transdecoderProteinFastaFile,minLen)
makerFilterProteinsByLengthDict = filterProteinsByLength(makerProteinFastaFile,minLen)

transdecoderCDSDict,transdecoder_mRNACoordDict,transdecoderGeneIDs,transdecoderCDSLengthDict,sortedTransdecoderCDSDict,originalScaffoldIDSortedTransdecoderCDSDict = readGFF(fullTransdecoderGFF)
makerCDSDict,maker_mRNACoordDict,makerGeneIDs,makerCDSLengthDict,sortedMakerCDSDict,originalScaffoldIDSortedMakerCDSDict = readGFF(fullMakerGFF)

makerGenesToKeep,transdecoderGenesToKeep = getNonOverlappingMakerGenes(makerCDSDict,transdecoderCDSDict,lengthsDict,sortedTransdecoderCDSDict,sortedMakerCDSDict,transdecoderLongestCDSDict,makerLongestTranscriptPerGeneDict)

newGeneIDMappingDict = renameGenes(originalScaffoldIDSortedMakerCDSDict,originalScaffoldIDSortedTransdecoderCDSDict,makerGeneIDs,transdecoderGeneIDs,lengthsList)

createFastaFiles(transdecoderProteinFastaFile,transdecoderCDSFastaFile,makerProteinFastaFile,makerCDSFastaFile,makerGenesToKeep,transdecoderGenesToKeep,newGeneIDMappingDict,transdecoderFilterProteinsByLengthDict,makerFilterProteinsByLengthDict,topLengthScaffoldDict)

