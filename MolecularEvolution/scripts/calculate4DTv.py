import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from Bio import AlignIO


###############
# SUBROUTINES #
###############


def calculateDistances(fastaFile,codonTable):
    purines = ['A','G']
    pyrimidines = ['C','T']
    baseName = os.path.basename(fastaFile)
    getFilePrefix = re.search('(.+)_aligned',baseName)
    filePrefix = getFilePrefix.group(1)
    sequences = SeqIO.parse(fastaFile,'fasta')
    pair = []
    for record in sequences:
        record.seq = record.seq.upper()
        pair.append(record)
    CDS1 = pair[0]
    CDS2 = pair[1]

    # kimura
    totalLen = 0
    transitions = 0
    transversions = 0
    for i in range(len(CDS1)):
        #print len(CDS1.seq)
        if CDS1[i] != '-' and CDS2[i] != '-':
            if CDS1[i] != CDS2[i]:
                if CDS1[i] in purines and CDS2[i] in pyrimidines:
                    transversions += 1
                elif  CDS1[i] in pyrimidines and CDS2[i] in purines:
                    transversions += 1
                else:
                    transitions += 1
            # count all ungapped positions
            totalLen += 1
    # print fastaFile,transitions,totalLen
    p = float(transitions)/float(totalLen)
    q = float(transversions)/float(totalLen)
    d = -(0.5)*np.log((1-2*p-q)*np.sqrt(1-2*q)) 
    kimura_out = open(filePrefix + '.kimura.out','w')
    kimura_out.write("fileID\ttransitions\ttransversions\ttotalSeqLen\tp\tq\tdistance\n")
    kimura_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (filePrefix,transitions,transversions,totalLen,p,q,d))

    # rate of transversions among four-fold degenerate synonymous sites 
    # FDSiteCount is the number of codons with a 4D synonymous site
    FDTv_out = open(filePrefix + '.4DTv.out','w')
    FDSiteCount = 0
    FDSiteThreshold = 25
    maximumFDTv = 0.5
    # https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-5-31
    # "We disregarded alignments with fewer than 25 conserved 4D sites as we cannot determine reliable 4 DTV distances for such proteins, and these are either too incomplete or evolving too fast at the peptide level for the purpose of our analysis."
    # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2593578/
    # "4DTV values are calculated for gene pairs having at least 10 fourfold degenerate sites."
    # "Fourfold degenerate sites are codons of amino acid residues G, A, T, P, V, and R, S, L"
    for i in range(0, len(CDS1.seq), 3):
        codonCDS1 = CDS1.seq[i:i+3]
        codonCDS2 = CDS2.seq[i:i+3]
        if '-' not in codonCDS1 and '-' not in codonCDS2:
            if codonCDS1 in codonTable and codonCDS2 in codonTable:
                if codonTable[codonCDS1] == codonTable[codonCDS2]:
                    FDSiteCount += 1
                    #print codonCDS1,codonCDS2
                    #print codonCDS1,codonTable[codonCDS1],codonCDS2,codonTable[codonCDS2]
    # Accelerated gene evolution and subfunctionalization in the pseudotetraploid frog Xenopus laevis // https://bmcbiol.biomedcentral.com/articles/10.1186/1741-7007-5-31
    
    '''
    "We identified the conserved four-fold degenerate amino acids within the alignments, extracted the corresponding codons in the underlying DNA sequence and calculated the 4 DTv distances (D 4DTv ) between each aligning pair as the fraction of four-fold degenerate (4D) third codon positions in which transversions are observed to have occurred."

    "4DTv ranges from 0 for recently duplicated peptides, to approx. 0.5 for paralogs that are so old that third codon nucleotides have essentially been randomized."
    '''

    if FDSiteCount >= FDSiteThreshold:
        # calculate 4DTv
        Tv = 0
        Ts = 0
        identicalSites = 0
        codonLength = 0
        for i in range(0, len(CDS1.seq), 3):
            codonCDS1 = CDS1.seq[i:i+3]
            codonCDS2 = CDS2.seq[i:i+3]
            if '-' not in codonCDS1 and '-' not in codonCDS2:
                # number of codons in sequence
                codonLength += 1
                if codonCDS1 in codonTable and codonCDS2 in codonTable:
                    if codonTable[codonCDS1] == codonTable[codonCDS2]:
                        fourD1 = codonCDS1[2]
                        fourD2 = codonCDS2[2]
                        #print fourD1,fourD2
                        if fourD1 != fourD2:
                            if fourD1 in purines and fourD2 in pyrimidines:
                                Tv += 1
                            elif fourD1 in pyrimidines and fourD2 in purines:
                                Tv += 1
                            else:
                                Ts += 1
                        else:
                            identicalSites += 1
                        #else:
                        #    print fourD1,fourD2
                    #else:
                    #    print codonCDS1,codonCDS2
        FDTv = float(Tv)/FDSiteCount
        FDTs = float(Ts)/FDSiteCount
        # "In addition, we calculated the fraction of 4D sites that had experienced any substitution, transition or transversion, D 4D . This distance measure gives better resolution for recent paralogs."
        D4D = float(Tv + Ts)/FDSiteCount
        # https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0218085#pone.0218085.s002
        # "De novo European eel transcriptome provides insights into the evolutionary history of duplicated genes in teleost lineages"
        # maximum FDTv of 0.5 is before correction
        if FDTv < maximumFDTv:
            # FDTv corrected
            # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2593578/ 
            FDTvCorrected = -0.5*np.log(1-2*FDTv)
            FDTv_out.write("fileID\tidenticalSites\ttransitions\ttransversions\t4DSiteCount\tcodonLength\tFDTs\tFDTv\tD4D\tFDTvCorrected\n")
            FDTv_out.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (filePrefix,identicalSites,Ts,Tv,FDSiteCount,codonLength,FDTs,FDTv,D4D,FDTvCorrected))
	else:
            FDTv_out.write("FDTv above threshold:\t%s\n" % (filePrefix))
    else:
	# if FDSiteCount >= FDSiteThreshold:
	FDTv_out.write("FDSiteCount below FDSiteThreshold:\t%s\n" % (filePrefix))
	

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <fasta file>\n"
if len(sys.argv) != 2:
    print usage
    sys.exit()


fastaFile = sys.argv[1]


# http://www.hgmd.cf.ac.uk/docs/cd_amino.html
# A, G, L, P, R, S, T, V
codonTable = {'GCT':'A', 'GCC':'A', 'GCA':'A', 'GCG':'A', 
              'GGT':'G', 'GGC':'G', 'GGA':'G', 'GGG':'G', 
              'CTT':'L', 'CTC':'L', 'CTA':'L', 'CTG':'L', 
              'CCT':'P', 'CCC':'P', 'CCA':'P', 'CCG':'P', 
              'CGT':'R', 'CGC':'R', 'CGA':'R', 'CGG':'R', 
              'TCT':'S', 'TCC':'S', 'TCA':'S', 'TCG':'S', 
              'ACT':'T', 'ACC':'T', 'ACA':'T', 'ACG':'T', 
              'GTT':'V', 'GTC':'V', 'GTA':'V', 'GTG':'V'}


calculateDistances(fastaFile,codonTable)

