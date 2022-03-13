import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from StringIO import StringIO


###############
# SUBROUTINES #
###############


def readFileList(fileList,codonTable):
    seqDict = {}
    with open(fileList,'r') as F:
        for line in F:
            fileName = line.strip()
            seqDict = createDictKeyedByOrthogroupAndSpecies(fileName,seqDict)
    # for every orthogroup, the species order will be the same
    for orthogroupID in seqDict:
        seqDict[orthogroupID].sort(key=lambda x:x[0])
    orthogroupPosDict = concatenate(seqDict,codonTable)
    speciesSpecificSeqDict = createFDTvSequence(orthogroupPosDict,seqDict)
    for speciesID in speciesSpecificSeqDict:
        speciesSpecificSeqDict[speciesID].sort(key=lambda x:x[0])
    speciesCoordDict = createFDTvFastas(speciesSpecificSeqDict)
    createCoordFile(speciesCoordDict)


def createDictKeyedByOrthogroupAndSpecies(fastaFile,seqDict):
    fileID = fastaFile
    getOGG = re.search('(OG.+)\.fasta',fileID)
    orthogroupID = getOGG.group(1)
    # print orthogroupID
    if orthogroupID not in seqDict:
        seqDict[orthogroupID] = []
    for record in SeqIO.parse(fastaFile,'fasta'):
        record.seq = record.seq.upper()
        seqDict[orthogroupID].append((record.id,record.seq))
    # else:
    #   print orthogroupID
    return(seqDict)
        

def concatenate(seqDict,codonTable):
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
                # check to make sure there's no gap in codon, and that codon is 4D
                if '-' not in codon:
                    if codon in codonTable:
                        # only if two conditions are satisfied, is codon appended to list for that codon position
                        if i not in posDict:
                            posDict[i] = []
                        posDict[i].append(codon)

        for codonPos in posDict:
            # because there were several conditions imposed above, not all species may be represented in list for that position -- but we are only interested in cases that feature all species
            if len(posDict[codonPos]) == number_of_species:
                fdAminoAcidList = []
                for codon in posDict[codonPos]:
                    aminoAcid = codonTable[codon]
                    #print codonPos,codon,aminoAcid
                    fdAminoAcidList.append(aminoAcid)
                # convert list to set, to get uniq amino acids
                fdAminoAcidSet = set(fdAminoAcidList)
                # only interested in cases where different codons encode the same 4D amino acid
                if len(fdAminoAcidSet) == 1:
                    # print codonPos,fdAminoAcidList,fdAminoAcidSet
                    if orthogroupID not in orthogroupPosDict:
                        orthogroupPosDict[orthogroupID] = []
                    orthogroupPosDict[orthogroupID].append((codonPos,fdAminoAcidSet))
                    #if codonPos not in orthogroupPosDict[orthogroupID]:
                    #    orthogroupPosDict[orthogroupID][codonPos] = fdAminoAcidSet
            #else:
            #    print "incorrect number of species"
            #    sys.exit()
    return(orthogroupPosDict)
            

def createFDTvSequence(orthogroupPosDict,seqDict):
    # orthogroupPosDict[orthogroupID].append((codonPos,fdAminoAcidSet))
    speciesSpecificSeqDict = {}
    for orthogroupID in seqDict:
        for speciesID,sequence in seqDict[orthogroupID]:
            if speciesID not in speciesSpecificSeqDict:
                speciesSpecificSeqDict[speciesID] = []
            # ortholog sequence
            #sequence = seqDict[orthogroupID][speciesID]
            # initialize empty string to add 4D codons to
            if orthogroupID in orthogroupPosDict:
                catCodon = ''
                for codonPos,aminoAcid in orthogroupPosDict[orthogroupID]:
                    codon = sequence[codonPos:codonPos+3]
                    #print type(codon)
                    thirdCodonPosition = codon[2]
                    # print type(thirdCodonPosition)
                    catCodon += thirdCodonPosition
                    # catCodon += codon
                # add string of 4D codons to list
                #print speciesID,orthogroupID,catCodon
                #if orthogroupID not in speciesSpecificSeqDict[speciesID]:
                #     speciesSpecificSeqDict[speciesID] = []
                speciesSpecificSeqDict[speciesID].append((orthogroupID,catCodon))
    return(speciesSpecificSeqDict)


def createFDTvFastas(speciesSpecificSeqDict):
    speciesCoordDict = {}
    # speciesSpecificSeqDict[speciesID].append((orthogroupID,catCodon))
    for speciesID in speciesSpecificSeqDict:
        # this is going to be the full string of 4D codons from all orthogroups
        concatenatedSequence = ''
        concatenatedSequenceList = []
        outName = speciesID  + ".fasta"
        seqStart = 0
        seqStop = 0
        if speciesID not in speciesCoordDict:
            speciesCoordDict[speciesID] = []
        for orthogroupID,concatenated4DTvCodon in speciesSpecificSeqDict[speciesID]:
            #print speciesID,orthogroupID
            # concatenated4DTvCodon = speciesSpecificSeqDict[speciesID][orthogroupID]
            # print speciesID,orthogroupID,concatenated4DTvCodon
            # add 4D codons to full string
            concatenatedSequence += concatenated4DTvCodon
            seqLen = len(concatenated4DTvCodon)
            seqStart = seqStop + 1
            seqStop += (seqLen + 1)
            speciesCoordDict[speciesID].append((orthogroupID,seqStart,seqStop,seqLen))
        # creating a new SeqRecord object
        # https://biopython.org/wiki/SeqRecord
        newRecord = SeqRecord(
            Seq(concatenatedSequence),
            id=speciesID,
            name='',
            description='',
        )
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


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <fasta file list>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

fileList = sys.argv[1]

# A, G, L, P, R, S, T, V
codonTable = {'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 
              'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 
              'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 
              'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 
              'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
              'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 
              'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
              'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V'}


readFileList(fileList,codonTable)



