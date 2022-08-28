import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO


###############
# SUBROUTINES #
###############


def readFileList(fileList):
    seqDict = {}
    with open(fileList,'r') as F:
        for line in F:
            fileName = line.strip()
            seqDict = createDictKeyedByOrthogroupAndSpecies(fileName,seqDict)
    # for every orthogroup, the species order will be the same
    for orthogroupID in seqDict:
        seqDict[orthogroupID].sort(key=lambda x:x[0])
    orthogroupPosDict = concatenate(seqDict)
    speciesSpecificSeqDict = createSequence(orthogroupPosDict,seqDict)
    for speciesID in speciesSpecificSeqDict:
        speciesSpecificSeqDict[speciesID].sort(key=lambda x:x[0])
    speciesCoordDict = createFastas(speciesSpecificSeqDict)
    createCoordFile(speciesCoordDict)
    createRaxmlPartitionFile(speciesCoordDict)
    createPartitionFinderFile(speciesCoordDict)


def createDictKeyedByOrthogroupAndSpecies(fastaFile,seqDict):
    fileID = fastaFile
    getOGG = re.search('(OG.+)\.fasta',fileID)
    orthogroupID = getOGG.group(1)
    # print orthogroupID
    if orthogroupID not in seqDict:
        seqDict[orthogroupID] = []
    for record in SeqIO.parse(fastaFile,'fasta'):
        record.seq = record.seq.upper()
        # record.id is the speciesID
        seqDict[orthogroupID].append((record.id,record.seq))
        #print orthogroupID,record.id,record.seq
    # else:
    #   print orthogroupID
    return(seqDict)
        

def concatenate(seqDict):
    # https://biopython-tutorial.readthedocs.io/en/latest/notebooks/04%20-%20Sequence%20Annotation%20objects.html
    number_of_species = 8
    orthogroupPosDict = {}
    # orthogroups are not ordered
    for orthogroupID in seqDict:
        posDict = {}
        # speciesIDs are sorted
        for speciesID,sequence in seqDict[orthogroupID]:
            # loop-through sequence in codon-increments
            for i in range(0, len(sequence), 3):
                codon = sequence[i:i+3]
                # check to make sure there's no gap in codon
                if '-' not in codon:
                    # i is the codon position
                    if i not in posDict:
                        posDict[i] = []
                    posDict[i].append(codon)

        # loop through codon positions
        # for every position, there are 8 codons, corresponding to each of the 8 species, per position
        for codonPos in posDict:
            # because there were several conditions imposed above, not all species may be represented in list for that position -- but we are only interested in cases that feature all species
            if len(posDict[codonPos]) == number_of_species:
                #for codon in posDict[codonPos]:
                if orthogroupID not in orthogroupPosDict:
                    orthogroupPosDict[orthogroupID] = {}
                if codonPos not in orthogroupPosDict[orthogroupID]:
                    orthogroupPosDict[orthogroupID][codonPos] = 1
            #else:
            #    print "incorrect number of species"
            #    sys.exit()
    return(orthogroupPosDict)
            

def createSequence(orthogroupPosDict,seqDict):
    speciesSpecificSeqDict = {}
    # contains the sequences and speciesIDs for each orthogroup
    for orthogroupID in seqDict:
        for speciesID,sequence in seqDict[orthogroupID]:
            if speciesID not in speciesSpecificSeqDict:
                speciesSpecificSeqDict[speciesID] = []
            # ortholog sequence
            # initialize empty string to add codons to
            if orthogroupID in orthogroupPosDict:
                catCodon = ''
                # orthogroupPosDict[orthogroupID][codonPos] 
                for codonPos in orthogroupPosDict[orthogroupID]:
                    codon = sequence[codonPos:codonPos+3]
                    codon = str(codon)
                    catCodon += codon
                # add string of codons to list
                speciesSpecificSeqDict[speciesID].append((orthogroupID,catCodon))
    return(speciesSpecificSeqDict)


def createFastas(speciesSpecificSeqDict):
    # seqLenThreshold is related to the number of codons that end up being partitioned for partitionfinder...
    seqLenThreshold = 50
    speciesCoordDict = {}
    # speciesSpecificSeqDict[speciesID].append((orthogroupID,catCodon))
    for speciesID in speciesSpecificSeqDict:
        # this is going to be the full string of codons from all orthogroups
        concatenatedSequence = ''
        concatenatedSequenceList = []
        outName = speciesID  + ".fasta"
        seqStart = 0
        seqStop = 0
        if speciesID not in speciesCoordDict:
            speciesCoordDict[speciesID] = []
        for orthogroupID,concatenatedCodon in speciesSpecificSeqDict[speciesID]:
            # add codons to full string
            seqLen = len(concatenatedCodon)
            if seqLen >= seqLenThreshold:
                concatenatedSequence += concatenatedCodon
                seqStart = seqStop
                seqStop += seqLen
                #seqStart = seqStop + 1
                #seqStop += (seqLen + 1)
                speciesCoordDict[speciesID].append((orthogroupID,seqStart,seqStop,seqLen))
        # creating a new SeqRecord object
        # https://biopython.org/wiki/SeqRecord
        newRecord = SeqRecord(
            Seq(concatenatedSequence),
            id=speciesID,
            name='',
            description='')
        # print newRecord
        concatenatedSequenceList.append(newRecord)
        # print speciesID,concatenatedSequenceList,newRecord
        print("%s\t%s\n" % (speciesID,len(newRecord.seq)))
        SeqIO.write(concatenatedSequenceList, outName, "fasta")
    return(speciesCoordDict)


def createCoordFile(speciesCoordDict):
    for speciesID in speciesCoordDict:
        species_out = open(speciesID + ".coords.txt",'w')
        species_out.write("orthogroupID\tseqStart\tseqStop\tseqLen\n")
        for orthogroupID,seqStart,seqStop,seqLen in speciesCoordDict[speciesID]:
            species_out.write("%s\t%s\t%s\t%s\n" % (orthogroupID,seqStart,seqStop,seqLen))


def createRaxmlPartitionFile(speciesCoordDict):
    # http://cme.h-its.org/exelixis/web/software/raxml/hands_on.html
    '''
    DNA, p1=1-60\3,2-60\3
    DNA, p2=3-60\3

    DNA,Gene1_codonPos1And2=1-2916\3,2-2916\3
    DNA,Gene1_codonPos3=3-2916\3
    DNA,Gene2_codonPos1And2=2917-4023\3,2918-4023\3
    DNA,Gene2_codonPos3=2919-4023\3
    DNA,Gene3_codonPos1And2=4024-4398\3,4025-4398\3
    DNA,Gene3_codonPos3=4026-4398\3
    DNA,Gene4_codonPos1And2=4399-5919\3,4400-5919\3
    DNA,Gene4_codonPos3=4401-5919\3
    '''
    for speciesID in speciesCoordDict:
        geneCount = 0
        PAR = open(speciesID + '_raxml_partition.txt','w')
        for orthogroupID,seqStart,seqStop,seqLen in speciesCoordDict[speciesID]:
            geneCount += 1
            pos1Label = str(seqStart + 1) + "-" + str(seqStop) + "\\3," + str(seqStart + 2) + "-" + str(seqStop) + "\\3"  
            pos2Label = str(seqStart + 3) + "-" + str(seqStop) + "\\3"
            PAR.write("DNA,Gene%s_codonPos1And2=%s\n" % (geneCount,pos1Label))
            PAR.write("DNA,Gene%s_codonPos3=%s\n" % (geneCount,pos2Label))
        PAR.write("\n")


# http://www.robertlanfear.com/partitionfinder/tutorial/
def createPartitionFinderFile(speciesCoordDict):
    for speciesID in speciesCoordDict:
        geneCount = 0
        PAR = open(speciesID + '_partitionfinder.txt','w')
        PAR.write("## ALIGNMENT FILE ##\n")
        PAR.write("alignment = allSpecies.fullCodon.phy;\n")
        PAR.write("\n")
        PAR.write("## BRANCHLENGTHS: linked | unlinked ##\n")
        PAR.write("branchlengths = linked;\n")
        PAR.write("\n")
        PAR.write("## MODELS OF EVOLUTION: all | allx | mrbayes | beast | gamma | gammai | <list> ##\n")
        PAR.write("models = GTR, GTR+G, GTR+I+G;\n")
        PAR.write("\n")
        PAR.write("# MODEL SELECCTION: AIC | AICc | BIC #\n")
        PAR.write("model_selection = aicc;\n")
        PAR.write("\n")
        PAR.write("## DATA BLOCKS: see manual for how to define ##\n")
        PAR.write("[data_blocks]\n")
        for orthogroupID,seqStart,seqStop,seqLen in speciesCoordDict[speciesID]:
            '''
            Gene1_pos1 = 1-789\3;
            Gene1_pos2 = 2-789\3;
            Gene1_pos3 = 3-789\3;
            Gene2_pos1 = 790-1449\3;
            '''
            geneCount += 1
            pos1Label = str(seqStart + 1) + "-" + str(seqStop) + "\\3;"  
            pos2Label = str(seqStart + 2) + "-" + str(seqStop) + "\\3;"
            pos3Label = str(seqStart + 3) + "-" + str(seqStop) + "\\3;" 
            PAR.write("Gene%s_pos1 = %s\n" % (geneCount,pos1Label))
            PAR.write("Gene%s_pos2 = %s\n" % (geneCount,pos2Label))
            PAR.write("Gene%s_pos3 = %s\n" % (geneCount,pos3Label))
        PAR.write("\n")
        PAR.write("## SCHEMES, search: all | user | greedy | rcluster | rclusterf | kmeans ##\n")
        PAR.write("[schemes]\n")
        PAR.write("search = greedy;\n")


'''
Gene1_pos1 = 1-789\3;
Gene1_pos2 = 2-789\3;
Gene1_pos3 = 3-789\3;
Gene2_pos1 = 790-1449\3;
Gene2_pos2 = 791-1449\3;
Gene2_pos3 = 792-1449\3;
Gene3_pos1 = 1450-2208\3;
Gene3_pos2 = 1451-2208\3;
Gene3_pos3 = 1452-2208\3;
'''

############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <fasta file list>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

fileList = sys.argv[1]


readFileList(fileList)



