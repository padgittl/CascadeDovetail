import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq

###############
# SUBROUTINES #
###############

#fileName = os.path.basename(fullPath)
#filePrefix,fileExt = os.path.splitext(fileName)
#baseName,extra = os.path.splitext(filePrefix)

def readFileList(fileList,originalConfig):
    directoryList = open('make_dirs.sh','w')
    with open(fileList,'r') as FL:
        for line in FL:
            fileName = line.strip()
            modifyConfig(fileName,originalConfig,directoryList)
            


'''
      seqfile = examples/abglobin.nuc * sequence data file name
      outfile = yn           * main result file
      verbose = 0  * 1: detailed output (list sequences), 0: concise output

        icode = 0  * 0:universal code; 1:mammalian mt; 2-10:see below

    weighting = 0  * weighting pathways between codons (0/1)?
   commonf3x4 = 0  * use one set of codon freqs for all pairs (0/1)? 
*       ndata = 1


* Genetic codes: 0:universal, 1:mammalian mt., 2:yeast mt., 3:mold mt.,
* 4: invertebrate mt., 5: ciliate nuclear, 6: echinoderm mt., 
* 7: euplotid mt., 8: alternative yeast nu. 9: ascidian mt., 
* 10: blepharisma nu.
* These codes correspond to transl_table 1 to 11 of GENEBANK.
'''



def modifyConfig(pal2nalFile,originalConfig,directoryList):
    fullPath = pal2nalFile.strip()
    fileName = os.path.basename(fullPath)
    # ../trimalOutput/XP_030477675.1_vs_XP_030490489.1_aligned_CDS_noFS_NT.trimal.fasta
    # fileNamePrefix,fileNameSuffix = fileName.split('_aligned_CDS_noFS_NT.fasta')
    getFileNamePrefix = re.search('(.+)_aligned',fileName)
    fileNamePrefix = getFileNamePrefix.group(1)
    #print fileNamePrefix
    make_directory_command = "mkdir " + fileNamePrefix
    directoryList.write("%s\n" % (make_directory_command))
    #print make_directory_command
    OUT = open(fileNamePrefix + '.yn00.ctl','w')
    with open(originalConfig,'r') as F:
        fileContent = F.read()
        #print fileContent
        oldSeqFileInfo = 'seqfile = examples/abglobin.nuc * sequence data file name'
        newSeqFileInfo = 'seqfile = ' + fullPath
        fileContent = fileContent.replace(oldSeqFileInfo,newSeqFileInfo)
        oldOutFileInfo = 'outfile = yn           * main result file'
        newOutFileInfo = 'outfile = ' + fileNamePrefix + ".yn00.txt"
        fileContent = fileContent.replace(oldOutFileInfo,newOutFileInfo)
        OUT.write(fileContent)
        OUT.close()
    F.close()


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <original config file> <file list> "
if len(sys.argv) != 3:
    print usage
    sys.exit()

originalConfig = sys.argv[1]
fileList = sys.argv[2] 

readFileList(fileList,originalConfig)



