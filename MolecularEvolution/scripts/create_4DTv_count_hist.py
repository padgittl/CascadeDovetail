import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


###############
# SUBROUTINES #
###############

def readCannabisOnlyFileList(cannabisFileList):
    dataList = []
    can_4dtv_out = open('cannabis_4dtv.txt','w')
    with open(cannabisFileList,'r') as CFL:
        for line in CFL:
            fullFileName = line.strip()
            fileName = os.path.basename(fullFileName)
            # print fileName
            filePrefix,extraStuff = fileName.split('.4DTv.out')
            # print filePrefix
            FDTvCorrected = readFile(fullFileName)
            if isinstance(FDTvCorrected, float):
                dataList.append(FDTvCorrected)
                can_4dtv_out.write("%s\n" % (FDTvCorrected))
            #else:
            #    print filePrefix
    return(dataList)


def readHopOnlyFileList(hopFileList):
    dataList = []
    hop_4dtv_out = open('hop_4dtv.txt','w')
    with open(hopFileList,'r') as HFL:
        for line in HFL:
            fullFileName = line.strip()
            fileName = os.path.basename(fullFileName)
            filePrefix,extraStuff = fileName.split('.4DTv.out')
            geneID1,geneID2 = filePrefix.split('_vs_')
            # print geneID1,geneID2
            FDTvCorrected = readFile(fullFileName)
            if isinstance(FDTvCorrected, float):
                dataList.append(FDTvCorrected)
                hop_4dtv_out.write("%s\n" % (FDTvCorrected))
            #else:
            #    print filePrefix
    return(dataList)


def readHopAndCannabisFileList(hop_vs_cannabisFileList):
    dataList = []
    hop_vs_can_4dtv_out = open('hop_vs_can_4dtv.txt','w')
    with open(hop_vs_cannabisFileList,'r') as HCFL:
        for line in HCFL:
            fullFileName = line.strip()
            fileName = os.path.basename(fullFileName)
            filePrefix,extraStuff = fileName.split('.4DTv.out')
            geneID1,geneID2 = filePrefix.split('_vs_')
            # print geneID1,geneID2
            FDTvCorrected = readFile(fullFileName)
            # print geneID1
            if isinstance(FDTvCorrected, float):
                dataList.append(FDTvCorrected)
                hop_vs_can_4dtv_out.write("%s\n" % (FDTvCorrected))
                #print FDTvCorrected
                #else:
                #    print FDTvCorrected
            #else:
            #    print filePrefix
    return(dataList)


def readFile(inputFile):
    with open(inputFile,'r') as I:
        for line in I:
            if not line.startswith('fileID'):
                #print line
                if not 'FDSiteCount below FDSiteThreshold' in line:
                    if not 'FDTv above threshold' in line:
                        fileID,identicalSites,transitions,transversions,FDSiteCount,codonLength,FDTs,FDTv,D4D,FDTvCorrected = line.strip().split('\t')
                        if FDTvCorrected == '-0.0':
                            FDTvCorrected = 0.0
                        FDTvCorrected = round(float(FDTvCorrected), 5)
                        return(FDTvCorrected)
                    else:
                        return(False)
                else:
                    return(False)


def createHist(cannabisDataList,hopDataList,hop_vs_cannabisDataList):
    binWidth = 0.05
    
    fig = plt.figure(figsize=(5,4))
    combinedLists = cannabisDataList + hopDataList + hop_vs_cannabisDataList
    bins = np.arange(min(combinedLists),max(combinedLists),binWidth)
    counts1, bins1, bars1 = plt.hist(cannabisDataList, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp, ' + str(len(cannabisDataList)) + ' gene pairs', alpha=0.75, linewidth=2)
    counts2, bins2, bars2 = plt.hist(hopDataList, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop, ' + str(len(hopDataList)) + ' gene pairs', alpha=0.75, linewidth=2)
    counts3, bins3, bars3 = plt.hist(hop_vs_cannabisDataList, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp, ' + str(len(hop_vs_cannabisDataList)) + ' gene pairs', alpha=0.75, linewidth=2)

    #plt.xlim((0, 6))
    #plt.xscale('log')
    #plt.yscale('log')
    # https://stackoverflow.com/questions/42045767/how-can-i-change-the-x-axis-in-matplotlib-so-there-is-no-white-space
    # plt.margins(x=0)
    plt.xlabel('4DTv',size=16)
    plt.ylabel('Count',size=16)
    #plt.legend(loc='upper right',frameon=False,fontsize=14)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('FDTvCountHist_binSize' + str(binWidth) + '_noFS.png', dpi=600)
    plt.savefig('FDTvCountHist_binSize' + str(binWidth) + '_noFS.pdf')
    plt.savefig('FDTvCountHist_binSize' + str(binWidth) + '_noFS.svg')
    plt.close()

    # density hist
    fig = plt.figure(figsize=(5,4))
    #binWidth = 0.05
    #binWidth = 0.025
    binWidth = 0.01
    combinedLists = cannabisDataList + hopDataList + hop_vs_cannabisDataList
    bins = np.arange(min(combinedLists),max(combinedLists),binWidth)
    counts1, bins1, bars1 = plt.hist(cannabisDataList, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp, ' + str(len(cannabisDataList)) + ' gene pairs', alpha=0.75, density=True, linewidth=1)
    counts2, bins2, bars2 = plt.hist(hopDataList, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop, ' + str(len(hopDataList)) + ' gene pairs', alpha=0.75, density=True, linewidth=1)
    counts3, bins3, bars3 = plt.hist(hop_vs_cannabisDataList, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp, ' + str(len(hop_vs_cannabisDataList)) + ' gene pairs', alpha=0.75, density=True, linewidth=1)
    
    plt.xlabel('4DTv',size=16)
    plt.ylabel('Density',size=16)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('FDTvDensityHist_binSize' + str(binWidth) + '_noFS.png', dpi=600)
    plt.savefig('FDTvDensityHist_binSize' + str(binWidth) + '_noFS.pdf')
    plt.savefig('FDTvDensityHist_binSize' + str(binWidth) + '_noFS.svg')
    plt.close()


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <cannabis file list> <hop file list> <hop vs cannabis file list>\n"
if len(sys.argv) != 4:
    print usage
    sys.exit()

cannabisFileList = sys.argv[1] 
hopFileList = sys.argv[2]
hop_vs_cannabisFileList = sys.argv[3]


cannabisDataList = readCannabisOnlyFileList(cannabisFileList)
hopDataList = readHopOnlyFileList(hopFileList)
hop_vs_cannabisDataList = readHopAndCannabisFileList(hop_vs_cannabisFileList)

'''
cannabisDataList = readFileList(cannabisFileList)
hopDataList = readFileList(hopFileList)
hop_vs_cannabisDataList = readFileList(hop_vs_cannabisFileList)
'''
createHist(cannabisDataList,hopDataList,hop_vs_cannabisDataList)



