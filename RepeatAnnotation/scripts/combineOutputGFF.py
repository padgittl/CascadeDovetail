import sys, re, os

################
# SUBROUTINES ##
################

def createFilteredRepMaskGFF(mipsREdatRepMaskGFF,outName):
    repDict = {}
    OUT = open(outName,'w')
    with open(mipsREdatRepMaskGFF,'r') as RMF:
        for line in RMF:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                if "LTR" not in line:
                    OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (contigID,source,feature,start,end,score,strand,frame,attribute))
    OUT.close()

def catFiles(denovoLTRGFF,outName):
    catCommand = "cat " + denovoLTRGFF + " >> " + outName
    os.system(catCommand)

def modifyOutputFile(outName):
    idCount = 0
    modifiedFileName = "modified_" + outName
    modifiedOUT = open(modifiedFileName,'w')
    with open(outName,'r') as OF:
        for line in OF:
            if '#' not in line:
                contigID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                idCount += 1
                if 'LTR' in feature:
                    # Scaffold_1000   RepeatMasker    similarity      20463   20580   23.5    -       .       Target "Motif:TXX_174008|TXX_Gh_DX401492.1_6900|MobileElement|02|34274|Gossypium" 537 661
                    modified_denovoLTR_attribute = "ID=rep" + str(idCount) + ";Name=" + feature + ";Alias=" + feature
                    #feature = "repeat"
                    modifiedOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (contigID,source,feature,start,end,score,strand,frame,modified_denovoLTR_attribute))
                elif 'similarity' in feature:
                    #feature = "repeat"
                    parseAttributeCol = re.search('Target \"(.+)\"',attribute)
                    attributeCol = parseAttributeCol.group(1)
                    modifiedMipsRepMaskAttribute = "ID=rep"+ str(idCount) + ";Name=" + attributeCol + ";Alias=" + attributeCol
                    modifiedOUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (contigID,source,feature,start,end,score,strand,frame,modifiedMipsRepMaskAttribute))
                else:
                    print "other element in feature"
                    sys.exit()
    modifiedOUT.close()
    return(modifiedFileName)

def sortFile(modifiedFileName):
    sortCommand = "sort -k1V -k4n -k5n " + modifiedFileName + " > sorted_" + modifiedFileName
    os.system(sortCommand)

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <mipsREdat RepMask GFF file> <de novo LTR GFF file> <out file name>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

mipsREdatRepMaskGFF = sys.argv[1]
denovoLTRGFF = sys.argv[2]
outName = sys.argv[3]

createFilteredRepMaskGFF(mipsREdatRepMaskGFF,outName)
catFiles(denovoLTRGFF,outName)
modifiedFileName = modifyOutputFile(outName)
sortFile(modifiedFileName)
