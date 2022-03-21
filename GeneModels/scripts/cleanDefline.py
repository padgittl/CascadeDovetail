import os,re,sys
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord


def cleanFastaFile(fastaFile,seqType):
    recordDict = {}
    recordList = []
    for record in SeqIO.parse(fastaFile,"fasta"):
        # print record.id
        if record.id not in recordDict:
            record.name = record.name.replace(record.name,'')
            record.description = record.description.replace(record.description,'')
            recordDict[record.id] = record
            recordList.append(record)
    SeqIO.write(recordList,"clean_transdecoder." + seqType + ".fasta", "fasta")
        

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <transdecoder protein fasta> <transdecoder cds fasta> \n"
if len(sys.argv) != 3:
    print usage
    sys.exit()

transdecoderProteinFastaFile = sys.argv[1]
transdecoderCDSFastaFile = sys.argv[2]

cleanFastaFile(transdecoderProteinFastaFile,'pep')
cleanFastaFile(transdecoderCDSFastaFile,'cds')
