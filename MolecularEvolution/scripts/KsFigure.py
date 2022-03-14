import sys,re,os
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt



###############
# SUBROUTINES #
###############


def readFile(fileName):
    valueList = []
    with open(fileName,'r') as F:
        for line in F:
            value = line.strip()
            valueList.append(float(value))
    return(valueList)


def zipLists(dataList,densityList):
    zippedList = []
    sortedDataList = []
    sortedDensityList = []
    for i in range(len(dataList)):
        zippedList.append((dataList[i],densityList[i]))
    zippedList.sort(key=lambda x: x[0], reverse=False)
    for datum,density in zippedList:
        #print j,k
        sortedDataList.append(datum)
        sortedDensityList.append(density)
    return(sortedDataList,sortedDensityList)
            

def createFigure(sortedCannabisDataList,sortedCannabisDensityList,sortedHopDataList,sortedHopDensityList,sorted_hop_vs_cannabisDataList,sorted_hop_vs_cannabisDensityList,cannabisMeansList,hopMeansList,hop_vs_cannabisMeansList,numberCannabisPairs,numberHopPairs,numberHop_vs_cannabisPairs,ksLowerThreshold,ksUpperThreshold,outPrefix):
    lambda_value = 6.1*10**-9    
    # numberCannabisPairs,numberHopPairs,numberHop_vs_cannabisPairs

    plt.plot(sortedCannabisDataList,sortedCannabisDensityList, label='Hemp vs Hemp, ' + str(numberCannabisPairs) + ' gene pairs', color='#0471A6', linewidth=3)
    plt.plot(sortedHopDataList,sortedHopDensityList, label='Hop vs Hop, ' + str(numberHopPairs) + ' gene pairs', color='#A41623', linewidth=3)
    plt.plot(sorted_hop_vs_cannabisDataList,sorted_hop_vs_cannabisDensityList, label='Hop vs Hemp, ' + str(numberHop_vs_cannabisPairs) + ' gene pairs', color='#F85E00', linewidth=3)
    
    ### hemp #731963 // hop #A41623 // hop vs hemp #F85E00

    # binWidth = 0.15
    ###binWidth = 0.04
    binWidth = 0.04
    # binWidth = 0.025
    plt.rcParams['patch.edgecolor'] = 'None'
    combinedLists = sortedCannabisDataList + sortedHopDataList + sorted_hop_vs_cannabisDataList
    # print sortedCannabisDataList
    bins = np.arange(min(combinedLists),max(combinedLists),binWidth)
    counts1, bins1, bars1 = plt.hist(sortedCannabisDataList, bins = bins, histtype = 'step', color ='#0471A6', alpha=0.5, normed=True, linewidth=3)
    counts2, bins2, bars2 = plt.hist(sortedHopDataList, bins = bins, histtype = 'step', color ='#A41623', alpha=0.5, normed=True, linewidth=3)
    counts3, bins3, bars3 = plt.hist(sorted_hop_vs_cannabisDataList, bins = bins, histtype = 'step', color ='#F85E00', alpha=0.5, normed=True, linewidth=3)

    allDensityLists = sortedCannabisDensityList + sortedHopDensityList + sorted_hop_vs_cannabisDensityList
    log_base = 10
    for i in cannabisMeansList:
        originalKs = log_base ** i
        time = (float(originalKs)/(2*lambda_value))/1000000
        time = round(time, 3)
        originalKs = round(originalKs, 3)
        # plt.annotate(text, (i,max(allDensityLists)))
        plt.annotate(str(originalKs), (i,0))
        plt.annotate(str(time), (i,2))
        plt.axvline(i, linestyle='--', c='#0471A6', alpha=0.5)
        print "hemp: " + str(originalKs) + ', time: ' + str(time)

    for j in hopMeansList:
        originalKs = log_base ** j
        time = (float(originalKs)/(2*lambda_value))/1000000
        time = round(time, 3)
        originalKs = round(originalKs, 3)
        plt.annotate(str(originalKs), (j,0))
        plt.annotate(str(time), (j,2))
        plt.axvline(j, linestyle='--', c='#A41623', alpha=0.5)
        print "hop: " + str(originalKs) + ', time: ' + str(time)

    for x in hop_vs_cannabisMeansList:
        originalKs = log_base ** x
        time = (float(originalKs)/(2*lambda_value))/1000000
        time = round(time, 3)
        originalKs = round(originalKs, 3)
        plt.annotate(str(originalKs), (x,0))
        plt.annotate(str(time), (x,2))
        plt.axvline(x, linestyle='--', c='#F85E00', alpha=0.5)
        print "hop vs hemp: " + str(originalKs) + ', time: ' + str(time)

    # plt.tight_layout()
    plt.legend(frameon=False, loc="upper left", fontsize=14)
    plt.xlabel('Log-transformed Ks values', size=16)
    plt.ylabel('Density', size=16)
    # ksLowerThreshold,ksUpperThreshold
    plt.savefig(outPrefix + '.png', dpi=600)
    plt.savefig(outPrefix + '.pdf')
    plt.savefig(outPrefix + '.svg')
    plt.close()


########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <cannabis data> <cannabis densities> <cannabis means> <hop data> <hop densities> <hop means> <hop vs cannabis data> <hop vs cannabis densities> <hop vs cannabis means> <ks lower threshold> <ks upper threshold> <out-file prefix>\n"
if len(sys.argv) != 13:
    print usage
    sys.exit()


cannabisData = sys.argv[1]
cannabisDensities = sys.argv[2]
cannabisMeans = sys.argv[3]
hopData = sys.argv[4]
hopDensities = sys.argv[5]
hopMeans = sys.argv[6]
hop_vs_cannabisData = sys.argv[7]
hop_vs_cannabisDensities = sys.argv[8]
hop_vs_cannabisMeans = sys.argv[9]
ksLowerThreshold = sys.argv[10]
ksUpperThreshold = sys.argv[11]
outPrefix = sys.argv[12]
ksLowerThreshold = float(ksLowerThreshold)
ksUpperThreshold = float(ksUpperThreshold)

cannabisDataList = readFile(cannabisData)
cannabisDensityList = readFile(cannabisDensities)
cannabisMeansList = readFile(cannabisMeans)

hopDataList = readFile(hopData)
hopDensityList = readFile(hopDensities)
hopMeansList = readFile(hopMeans)

hop_vs_cannabisDataList = readFile(hop_vs_cannabisData)
hop_vs_cannabisDensityList = readFile(hop_vs_cannabisDensities)
hop_vs_cannabisMeansList = readFile(hop_vs_cannabisMeans)

numberCannabisPairs = len(cannabisDataList)
numberHopPairs = len(hopDataList)
numberHop_vs_cannabisPairs = len(hop_vs_cannabisDataList)

sortedCannabisDataList,sortedCannabisDensityList = zipLists(cannabisDataList,cannabisDensityList)
sortedHopDataList,sortedHopDensityList = zipLists(hopDataList,hopDensityList)
sorted_hop_vs_cannabisDataList,sorted_hop_vs_cannabisDensityList = zipLists(hop_vs_cannabisDataList,hop_vs_cannabisDensityList)

createFigure(sortedCannabisDataList,sortedCannabisDensityList,sortedHopDataList,sortedHopDensityList,sorted_hop_vs_cannabisDataList,sorted_hop_vs_cannabisDensityList,cannabisMeansList,hopMeansList,hop_vs_cannabisMeansList,numberCannabisPairs,numberHopPairs,numberHop_vs_cannabisPairs,ksLowerThreshold,ksUpperThreshold,outPrefix)
