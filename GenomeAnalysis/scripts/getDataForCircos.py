import sys, re, os, math
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines


sys.path.append('/svgwrite/svgwrite-1.2.1/')
import svgwrite
from svgwrite import cm, mm


###############
# SUBROUTINES #
###############


def readDefenseGeneFile(defenseGeneFile):
    defenseGeneDict = {}
    with open(defenseGeneFile,'r') as F:
        for line in F:
            # Scaffold_73     HUMLU_CAS0053498.t1.p1  15108726        15113171        Q7XA40  GO:0006952      sp|Q7XA40|RGA3_SOLBU Putative disease resistance protein RGA3 OS=Solanum bulbocastanum OX=147425 GN=RGA3 PE=2 SV=2      defense response
            scaffoldID,geneID,geneStart,geneStop,uniprotID,goTermID,uniprotDescription,goTermDescription = line.strip().split('\t')
            if scaffoldID not in defenseGeneDict:
                defenseGeneDict[scaffoldID] = []
            defenseGeneDict[scaffoldID].append((geneID,int(geneStart),int(geneStop)))
    return(defenseGeneDict)


def defenseGeneDensity(longestScaffoldList,defenseGeneDict,windowSize,megabaseMultiplier):
    defenseGeneDensityDict = {}
    defenseGeneDensityPlotDict = {}
    defenseGeneDensityList = []
    DEFENSE_GENE_OUT = open('defenseGeneDensities_window' + str(windowSize) + '.txt','w')
    for scaffoldID,scaffoldLength in longestScaffoldList:
        if scaffoldID not in defenseGeneDensityDict:
            defenseGeneDensityDict[scaffoldID] = {}
        if scaffoldID not in defenseGeneDensityPlotDict:
            defenseGeneDensityPlotDict[scaffoldID] = []
        for i in range(0,scaffoldLength,windowSize):
            geneCount = 0
            position = float(i) + (windowSize / 2)
            windowInterval = i + windowSize
            for geneID,geneStart,geneStop in defenseGeneDict[scaffoldID]:
                if i <= geneStart and geneStart < windowInterval:
                    geneCount += 1
                    megabaseGeneCount = geneCount
            defenseGeneDensityList.append((position,megabaseGeneCount))
            defenseGeneDensityPlotDict[scaffoldID].append((position,megabaseGeneCount))
            if (i,windowInterval) not in defenseGeneDensityDict[scaffoldID]:
                defenseGeneDensityDict[scaffoldID][(i,windowInterval)] = megabaseGeneCount
            DEFENSE_GENE_OUT.write("%s\t%s\t%s\t%s\n" % (scaffoldID,i,windowInterval,megabaseGeneCount))
    return(defenseGeneDensityList,defenseGeneDensityDict,defenseGeneDensityPlotDict)



def readLinkageMapFile(linkageMapFile,largestScaffoldCountDict):
    markerDict = {}
    # Trait   Marker  Chr     Pos     df      GeneticMap      F       p       add_effect      add_F   add_p   dom_effect      dom_F   dom_p   errordf MarkerR2        Genetic Var     Residual Var    -2LnLikelihood
    # sex     S1_468115197    1       100000  2       0       1.1379600000000001      0.32196999999999998     0.1198  1.8134300000000001      0.1792  0.12751999999999999     2.2037499999999999      0.13882 277     8.1799999999999998E-3   0.39034999999999997     0.000004        239.86838
    with open(linkageMapFile,'r') as L:
        for line in L:
            if 'Marker' not in line:
                trait,linkageGroup_markerPos,chromID,chromPos,df,geneticMapPos,F,p,addEffect,addF,addP,domEffect,domF,domP,errorDF,markerR2,geneticVar,residualVar,minus2LogLikelihood = line.strip().split('\t')
                # print linkageGroup_markerPos,chromID,chromPos,geneticMapPos,p
                # largestScaffoldCountDict[largestScaffoldCount] = scaffoldID
                extra,scaffoldPos = linkageGroup_markerPos.split('_')
                chromID = int(chromID)
                if chromID in largestScaffoldCountDict:
                    scaffoldID = largestScaffoldCountDict[chromID]
                else:
                    print "chromID not in dict"
                    sys.exit()
                if scaffoldID not in markerDict:
                    markerDict[scaffoldID] = []
                markerDict[scaffoldID].append((int(scaffoldPos),float(p)))
    return(markerDict)


def computeSNPDensity(markerDict,scaffoldLengthDict,windowSize,megabaseMultiplier,pThreshold):
    snpDensityDict = {}
    snpDensityPlotDict = {}
    snpDensityList = []
    SNP_DENSITY_OUT = open('snpDensities_window' + str(windowSize) + '.txt','w')
    SIG_SNP_DENSITY_OUT = open('sigSNPDensities_window' + str(windowSize) + '.txt','w')
    SNP_SCATTER_OUT = open('snpScatterData_window' + str(windowSize) + '.txt','w')
    # significant p-values only
    SIG_SNP_SCATTER_OUT = open('sigSNPScatterData_window' + str(windowSize) + '.txt','w')
    for scaffoldID in markerDict:
        scaffoldLength = scaffoldLengthDict[scaffoldID]
        if scaffoldID not in snpDensityDict:
            snpDensityDict[scaffoldID] = {}
        if scaffoldID not in snpDensityPlotDict:
            snpDensityPlotDict[scaffoldID] = []
        for i in range(0,scaffoldLength,windowSize):
            snpCount = 0
            megabaseSNPCount = 0
            sigSNPCount = 0
            sigMegabaseSNPCount = 0
            position = float(i) + (windowSize / 2)
            windowInterval = i + windowSize
            for markerPos,p in markerDict[scaffoldID]:
                # http://www.circos.ca/documentation/tutorials/2d_tracks/scatter_plots/
                # chr start end value options
                # hs1 2001 2001 0.010
                # minusLogP = -math.log10(p)
                # print p
                SNP_SCATTER_OUT.write("%s\t%s\t%s\t%s\n" % (scaffoldID,markerPos,markerPos,p))
                if i <= markerPos and markerPos < windowInterval:
                    snpCount += 1
                    megabaseSNPCount = snpCount
                    if p < pThreshold:
                        SIG_SNP_SCATTER_OUT.write("%s\t%s\t%s\t%s\n" % (scaffoldID,markerPos,markerPos,p))
                        sigSNPCount += 1
                        sigMegabaseSNPCount = sigSNPCount
            snpDensityList.append((position,megabaseSNPCount))
            snpDensityPlotDict[scaffoldID].append((position,megabaseSNPCount))
            if (i,windowInterval) not in snpDensityDict[scaffoldID]:
                snpDensityDict[scaffoldID][(i,windowInterval)] = megabaseSNPCount
            SNP_DENSITY_OUT.write("%s\t%s\t%s\t%s\n" % (scaffoldID,i,windowInterval,megabaseSNPCount))
            SIG_SNP_DENSITY_OUT.write("%s\t%s\t%s\t%s\n" % (scaffoldID,i,windowInterval,sigMegabaseSNPCount))
    return(snpDensityList,snpDensityDict,snpDensityPlotDict)


def readScaffoldLengthsFile(scaffoldLengthsFile):
    scaffoldNumber = 10
    scaffoldLengthDict = {}
    scaffoldLenList = []
    with open(scaffoldLengthsFile,'r') as LEN:
        for line in LEN:
            scaffoldID,scaffoldLen = line.strip().split()
            scaffoldLenList.append((scaffoldID,int(scaffoldLen)))
    scaffoldLenList.sort(key=lambda x:x[1], reverse=True)
    longestScaffoldList = scaffoldLenList[0:scaffoldNumber]
    largestScaffoldCount = 0
    largestScaffoldCountDict = {}
    for scaffoldID,scaffoldLen in longestScaffoldList:
        largestScaffoldCount += 1
        if scaffoldID not in scaffoldLengthDict:
            scaffoldLengthDict[scaffoldID] = scaffoldLen
        if largestScaffoldCount not in largestScaffoldCountDict:
            largestScaffoldCountDict[largestScaffoldCount] = scaffoldID
    return(scaffoldLengthDict,longestScaffoldList,largestScaffoldCountDict)



def readGeneMapFile(geneMapFile,recordDict):
    geneData = {}
    geneCoordDict = {}
    with open(geneMapFile,'r') as F:
        for line in F:
            scaffoldID,oldGeneID,newGeneID,mRNAStart,mRNAStop,cdsStart,cdsStop = line.strip().split('\t')
            if newGeneID in recordDict:
                geneCoordDict[newGeneID] = (int(mRNAStart),int(mRNAStop))
                if scaffoldID not in geneData:
                    geneData[scaffoldID] = []
                geneData[scaffoldID].append((newGeneID,int(mRNAStart),int(mRNAStop)))
    return(geneData,geneCoordDict)


def readFasta(fastaFile):
    recordDict = {}
    for record in SeqIO.parse(fastaFile,"fasta"):
        record.name = record.name.replace(record.name,'')
        record.description = record.description.replace(record.description,'')
        if record.id not in recordDict:
            recordDict[record.id] = record
    return(recordDict)


def computeGeneDensity(longestScaffoldList,geneData,windowSize,megabaseMultiplier):
    geneDensityDict = {}
    geneDensityPlotDict = {}
    geneDensityList = []
    GENE_OUT = open('geneDensities_window' + str(windowSize) + '.txt','w')
    for scaffoldID,scaffoldLength in longestScaffoldList:
        if scaffoldID not in geneDensityDict:
            geneDensityDict[scaffoldID] = {}
        if scaffoldID not in geneDensityPlotDict:
            geneDensityPlotDict[scaffoldID] = []
        for i in range(0,scaffoldLength,windowSize):
            geneCount = 0
            position = float(i) + (windowSize / 2)
            windowInterval = i + windowSize
            #print i,position,windowInterval
            # geneData[scaffoldID].append((oldGeneID,newGeneID,int(mRNAStart),int(mRNAStop)))
            for geneID,geneStart,geneStop in geneData[scaffoldID]:
                if i <= geneStart and geneStart < windowInterval:
                    geneCount += 1
                    # megabaseGeneCount = megabaseMultiplier * geneCount
                    # 05/24/2021
                    # using larger windows, so didn't want to use multiplier to get in terms of 1Mb
                    megabaseGeneCount = geneCount
            geneDensityList.append((position,megabaseGeneCount))
            geneDensityPlotDict[scaffoldID].append((position,megabaseGeneCount))
            if (i,windowInterval) not in geneDensityDict[scaffoldID]:
                geneDensityDict[scaffoldID][(i,windowInterval)] = megabaseGeneCount
            # print position,megabaseGeneCount
            GENE_OUT.write("%s\t%s\t%s\t%s\n" % (scaffoldID,i,windowInterval,megabaseGeneCount))
    return(geneDensityList,geneDensityDict,geneDensityPlotDict)


def readLTRGFF(ltrGFF):
    ltrData = {}
    # 000000F RepeatMasker    LTR/Gypsy       2768    2924    31.9    +       418     003288F:275585..286585_INT   
    with open(ltrGFF,'r') as LTR:
        for line in LTR:
            if not line.startswith('#'):
                contigID,source,feature,start,end,score,strand,frame,attribute = line.strip().split("\t")
                # print contigID,source,feature,start,end,score,strand,frame,attribute
                if contigID not in ltrData:
                    ltrData[contigID] = []
                ltrData[contigID].append((feature,int(start),int(end)))
    return(ltrData)


def computeLTRDensity(longestScaffoldList,ltrData,windowSize,megabaseMultiplier):
    ltrDensityDict = {}
    ltrDensityPlotDict = {}
    ltrDensityList = []
    LTR_OUT = open('ltrDensities_window' + str(windowSize) + '.txt','w')
    for scaffoldID,scaffoldLength in longestScaffoldList:
        if scaffoldID not in ltrDensityDict:
            ltrDensityDict[scaffoldID] = {}
        if scaffoldID not in ltrDensityPlotDict:
            ltrDensityPlotDict[scaffoldID] = []
        for i in range(0,scaffoldLength,windowSize):
            ltrCount = 0
            position = float(i)
            windowInterval = i + windowSize
            for ltrType,ltrStart,ltrStop in ltrData[scaffoldID]:
                if i <= ltrStart and ltrStart < windowInterval:
                    ltrCount += 1
                    # megabaseMultiplier = 1000000.0 / windowSize
                    # megabaseLTRCount = megabaseMultiplier * ltrCount
                    # 05/24/2021 
                    # using larger windows, so didn't want to use multiplier to get in terms of 1Mb
                    megabaseLTRCount = ltrCount
            ltrDensityList.append((position,megabaseLTRCount))
            ltrDensityPlotDict[scaffoldID].append((position,megabaseLTRCount))
            if (i,windowInterval) not in ltrDensityDict[scaffoldID]:
                ltrDensityDict[scaffoldID][(i,windowInterval)] = megabaseLTRCount
            LTR_OUT.write("%s\t%s\t%s\t%s\n" % (scaffoldID,i,windowInterval,megabaseLTRCount))
            #print position,megabaseLTRCount
    return(ltrDensityList,ltrDensityDict,ltrDensityPlotDict)


def stackedHist(geneDensityDict,ltrDensityDict,snpDensityDict,windowSize):
    DATA_OUT = open('stackedHistData_window' + str(windowSize) + '.txt','w')
    # hs7 65000000 69999999 0.388090,0.070074,0.547485,0.842239,0.525658
    # ltrDensityDict[scaffoldID].append((position,megabaseLTRCount))
    # geneDensityDict[scaffoldID].append((position,megabaseGeneCount))
    # snpDensityDict[scaffoldID].append((position,megabaseSNPCount))
    # ltrDensityDict[scaffoldID][(i,windowInterval)] = megabaseLTRCount
    dataDict = {}
    for scaffoldID in geneDensityDict:
        if scaffoldID not in dataDict:
            dataDict[scaffoldID] = []
        for (windowStart,windowStop) in geneDensityDict[scaffoldID]:
            geneCount = geneDensityDict[scaffoldID][(windowStart,windowStop)]
            ltrCount = ltrDensityDict[scaffoldID][(windowStart,windowStop)]
            snpCount = snpDensityDict[scaffoldID][(windowStart,windowStop)]
            # print scaffoldID,windowStart,windowStop,geneCount,ltrCount,snpCount
            dataDict[scaffoldID].append((windowStart,windowStop,geneCount,ltrCount,snpCount))
    for scaffoldID in dataDict:
        dataDict[scaffoldID].sort(key=lambda x:x[0], reverse=False)
        for windowStart,windowStop,geneCount,ltrCount,snpCount in dataDict[scaffoldID]:
            DATA_OUT.write("%s\t%s\t%s\t%s,%s,%s\n" % (scaffoldID,windowStart,windowStop,snpCount,geneCount,ltrCount))


def plotGeneAndLTRDensity(longestScaffoldList,geneDensityPlotDict,ltrDensityPlotDict,windowSize,megabaseMultiplier):
    longestScaffold = longestScaffoldList[0][1]
    #fig, ax1 = plt.subplots(figsize=(7,0.75)) 
    for scaffoldID,scaffoldLength in longestScaffoldList:
        fractionalScaffoldLengthRelative2LongestScaffold = float(scaffoldLength)/longestScaffold
        figLength = 7 * fractionalScaffoldLengthRelative2LongestScaffold
        #print scaffoldLength,longestScaffold,scaffoldLengthRelative2LongestScaffold,figLength
        fig, ax1 = plt.subplots(figsize=(figLength,2))
        geneDensityArray = np.array(geneDensityPlotDict[scaffoldID])
        ltrDensityArray = np.array(ltrDensityPlotDict[scaffoldID])
        for i in range(len(ltrDensityArray)-1):
            #for position,ltrDensity in ltrDensityArray:
            position,ltrDensity = ltrDensityArray[i]
            ax1.bar(position,ltrDensity,windowSize, color='#d6604d', align='edge')
            #ax1.set_xlim(0,linkageGroupTotalLen)
            ax1.spines['left'].set_color('#d6604d')
            ax1.spines['right'].set_color('black')
            ax1.tick_params(axis='y', colors='#d6604d')
            #plt.ylim(0,900)
            plt.margins(0)
        ax2 = ax1.twinx()
        for j in range(len(geneDensityArray)-2):
            position,geneDensity = geneDensityArray[j]
            nextPosition,nextGeneDensity = geneDensityArray[j+1]
            # ax2.plot((position,nextPosition),(geneDensity,nextGeneDensity), color='black', marker='o', linewidth=0.25, markersize=0.5) 
            ax2.plot((position,nextPosition),(geneDensity,nextGeneDensity), color='black', marker='o', linewidth=1.25, markersize=1.75)
            # linewidth=0.25, markersize=0.5
            #ax2.set_xlim(0,linkageGroupTotalLen)
            #print linkageGroupTotalLen
            ax2.spines['left'].set_color('#d6604d')
            ax2.spines['right'].set_color('black')
            ax2.tick_params(axis='y', colors='black')
            plt.margins(0)
            #plt.ylim(0,50)
            # https://stackoverflow.com/questions/10171618/changing-plot-scale-by-a-factor-in-matplotlib
            ticks = ax2.get_xticks()/1000000
            #print ticks
            ax2.set_xticklabels(ticks.astype(int))

        legendGenes = mlines.Line2D([], [], color='black', marker='o',markersize=6,label='Gene Density')
        legendLTR = mlines.Line2D([], [], color='#d6604d', marker='s',linestyle="None",markersize=6, label='LTR Density')
        plt.legend(handles=[legendGenes,legendLTR], frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5,-0.1))
        plt.savefig(scaffoldID + "_geneRepeatDensityWindow" + str(windowSize) + ".svg", bbox_inches = 'tight', pad_inches=0)
        plt.savefig(scaffoldID + "_geneRepeatDensityWindow" + str(windowSize) + ".pdf", bbox_inches = 'tight', pad_inches=0)
        plt.close()
            

def plotSNPs(snpDensityPlotDict,longestScaffoldList,windowSize):
    longestScaffold = longestScaffoldList[0][1]
    for scaffoldID,scaffoldLength in longestScaffoldList:
        fractionalScaffoldLengthRelative2LongestScaffold = float(scaffoldLength)/longestScaffold
        figLength = 7 * fractionalScaffoldLengthRelative2LongestScaffold
        fig, ax1 = plt.subplots(figsize=(figLength,2))
        snpDensityArray = np.array(snpDensityPlotDict[scaffoldID])
        for i in range(len(snpDensityArray)-1):
            position,snpDensity = snpDensityArray[i]
            ax1.bar(position,snpDensity,windowSize, color='#d6604d', align='edge')
            ax1.spines['left'].set_color('#d6604d')
            ax1.spines['right'].set_color('black')
            ax1.tick_params(axis='y', colors='#d6604d')
        plt.margins(0)
        plt.tight_layout()
        plt.savefig(scaffoldID + "_snpScatterWindow" + str(windowSize) + ".svg")
        plt.savefig(scaffoldID + "_snpScatterWindow" + str(windowSize) + ".pdf")


############
# MAIN  ####
############


usage = "Usage: " + sys.argv[0] + " <scaffold lengths file> <gene fasta file> <gene map file> <ltr gff> <linkage map file> <defense gene file> <window size: e.g. 1000000>\n"
if len(sys.argv) != 8:
    print usage
    sys.exit()

scaffoldLengthsFile = sys.argv[1]
fastaFile = sys.argv[2]
geneMapFile = sys.argv[3]
ltrGFF = sys.argv[4]
linkageMapFile = sys.argv[5]
defenseGeneFile = sys.argv[6]
windowSize = sys.argv[7]

### thresholds 
windowSize = int(windowSize)
megabaseMultiplier = 1000000.0 / windowSize
pThreshold = 0.05

scaffoldLengthDict,longestScaffoldList,largestScaffoldCountDict = readScaffoldLengthsFile(scaffoldLengthsFile)

# defense genes
defenseGeneDict = readDefenseGeneFile(defenseGeneFile)
defenseGeneDensityList,defenseGeneDensityDict,defenseGeneDensityPlotDict = defenseGeneDensity(longestScaffoldList,defenseGeneDict,windowSize,megabaseMultiplier)

# all genes
recordDict = readFasta(fastaFile)
geneData,geneCoordDict = readGeneMapFile(geneMapFile,recordDict)
geneDensityList,geneDensityDict,geneDensityPlotDict = computeGeneDensity(longestScaffoldList,geneData,windowSize,megabaseMultiplier)

# ltrs
ltrData = readLTRGFF(ltrGFF)
ltrDensityList,ltrDensityDict,ltrDensityPlotDict = computeLTRDensity(longestScaffoldList,ltrData,windowSize,megabaseMultiplier)

# snps
markerDict = readLinkageMapFile(linkageMapFile,largestScaffoldCountDict)
snpDensityList,snpDensityDict,snpDensityPlotDict = computeSNPDensity(markerDict,scaffoldLengthDict,windowSize,megabaseMultiplier,pThreshold)

# get data for circos stacked hist
stackedHist(geneDensityDict,ltrDensityDict,snpDensityDict,windowSize)

# plot hists to get y-axis scale
plotGeneAndLTRDensity(longestScaffoldList,geneDensityPlotDict,ltrDensityPlotDict,windowSize,megabaseMultiplier)
plotSNPs(snpDensityPlotDict,longestScaffoldList,windowSize)
