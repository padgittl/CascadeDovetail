import sys,re,os
from Bio import SeqIO
from Bio.Seq import Seq


###############
# SUBROUTINES #
###############


# loop through gff file again and get mRNA coords based on gene ID associated with cds features
def getGeneCoords(gffFile,mrna2proteinDict,mrna2proteinDict_unverified,geneIDs):
    coordDict = {}
    coordDict_unverified = {}
    OUT = open('cannabisGenes.gff','w')
    with open(gffFile,'r') as GFF:
        for line in GFF:
            if not line.startswith('#'):
                if not line.isspace():
                    line = line.strip()
                    scaffoldID,source,feature,featureStart,featureStop,score,strand,frame,attribute = line.strip().split("\t")
                    sp = "cs" + scaffoldID
                    # NC_044371.1     Gnomon  CDS     40250   40402   .       +       0       ID=cds-XP_030489337.1;Parent=rna-XM_030633477.1;Dbxref=GeneID:115705987,Genbank:XP_030489337.1;Name=XP_030489337.1;gbkey=CDS;gene=LOC115705987;product=probable serine/threonine-protein kinase At1g09600;protein_id=XP_030489337.1
                    if feature == 'CDS':
                        if 'NC_029855.1' in scaffoldID:
                            if 'YP_' in attribute:
                                getProteinID = re.search('protein_id=(.+)',attribute)
                                proteinID = getProteinID.group(1)
                                #print proteinID
                                if proteinID in mrna2proteinDict_unverified:
                                    if proteinID not in coordDict_unverified:
                                        coordDict_unverified[proteinID] = []
                                    coordDict_unverified[proteinID].append((sp,featureStart,featureStop))
                                else:
                                    print "proteinID not in mrna2proteinDict_unverified"
                                    sys.exit()
                        else:
                            getProteinID = re.search('ID=cds\-(.+);Parent=',attribute)
                            proteinID = getProteinID.group(1)
                            #print proteinID
                            if proteinID not in coordDict:
                                coordDict[proteinID] = []
                            coordDict[proteinID].append((sp,featureStart,featureStop))

    for proteinID in coordDict_unverified:
        coordDict_unverified[proteinID].sort(key=lambda x:x[1], reverse=False)
        geneStart = coordDict_unverified[proteinID][0][1]
        coordDict_unverified[proteinID].sort(key=lambda x:x[2], reverse=True)
        geneStop = coordDict_unverified[proteinID][0][2]
        sp = coordDict_unverified[proteinID][0][0]
        #print proteinID,geneStart,geneStop
        if proteinID in geneIDs:
            OUT.write("%s\t%s\t%s\t%s\n" % (sp,proteinID,geneStart,geneStop))

    for proteinID in coordDict:
        coordDict[proteinID].sort(key=lambda x:x[1], reverse=False)
        geneStart = coordDict[proteinID][0][1]
        coordDict[proteinID].sort(key=lambda x:x[2], reverse=True)
        geneStop = coordDict[proteinID][0][2]
        sp = coordDict[proteinID][0][0]
        if proteinID in geneIDs:
            OUT.write("%s\t%s\t%s\t%s\n" % (sp,proteinID,geneStart,geneStop))


# search for 'mcscanx manual'
def readFasta(fastaFile):
    geneIDs = {}
    for record in SeqIO.parse(fastaFile,"fasta"):
        if record.id not in geneIDs:
            geneIDs[record.id] = 1
            #print record.id
    return(geneIDs)

    # sp# gene starting_position ending_position
    # sp#, sp is the two-letter short name for the species; # is the chromosome number.
    # species is cs
    

def readGenBankFile(genbankFile):
    mrna2proteinDict = {}
    # YP_009243681.1 - example from genbank file - all protein seqs with 'YP_'
    mrna2proteinDict_unverified = {}
    listFromFile = []
    with open(genbankFile,'r') as G:
        for line in G:
            line = line.strip()
            # append all lines of file to a list 
            listFromFile.append(line)
    for i in range(len(listFromFile)):
        if 'VERSION' in listFromFile[i] and 'REVERSION' not in listFromFile[i]:
            # VERSION     XP_030477616.1
            info,proteinID = listFromFile[i].strip().split('VERSION') 
            proteinID = proteinID.strip()
            # DBSOURCE    REFSEQ: accession XM_030621756.1
            accession_info,mrnaID = listFromFile[i+2].split('accession')
            mrnaID = mrnaID.strip()
            mrna2proteinDict[mrnaID] = proteinID
            if mrnaID == 'NC_029855.1':
                if proteinID not in mrna2proteinDict_unverified:
                    mrna2proteinDict_unverified[proteinID] = 1
                
    #for protID in mrna2proteinDict:
    #    print protID,mrna2proteinDict[protID]
    # there are 35 genes that don't have an XM or XP-based ID - they have a locus tag 
    # these look like 'unverified' genes from other projects/publications
    for record in SeqIO.parse(genbankFile, "genbank"):
        #print record.id
        for f in record.features:
            # NC_029855.1
            # https://www.biostars.org/p/110284/
            if "locus_tag" in f.qualifiers:
                locus_tag_id = f.qualifiers["locus_tag"][0]
                #if locus_tag_id not in mrna2proteinDict_unverified:
                #    mrna2proteinDict_unverified[locus_tag_id] = 1
            if "coded_by" in f.qualifiers:
                #if 'NC_029855.1' in f.qualifiers["coded_by"][0]:
                    #print f.qualifiers["coded_by"]
                continue
    return(mrna2proteinDict,mrna2proteinDict_unverified)



########
# MAIN #
########


usage = "Usage: " + sys.argv[0] + " <fasta file> <gff file> <genbank file>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

fastaFile = sys.argv[1]
gffFile = sys.argv[2]
genbankFile = sys.argv[3]


mrna2proteinDict,mrna2proteinDict_unverified = readGenBankFile(genbankFile)
geneIDs = readFasta(fastaFile)
getGeneCoords(gffFile,mrna2proteinDict,mrna2proteinDict_unverified,geneIDs)

