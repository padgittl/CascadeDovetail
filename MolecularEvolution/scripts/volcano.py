#!/bin/python
import sys, os, re, math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.path as mpath
import matplotlib.patches as mpatches
import matplotlib.lines as mlines
import matplotlib.cm as cm
import matplotlib.colors as colors


###############
# SUBROUTINES #
###############


def retrieveNonZeroQValues(goTermFile):
    pValueList = []
    qValueList = []
    with open(goTermFile,'r') as F:
        for line in F:
            if 'goTermID' not in line:
                # goTermID        pValue  qValue  obsK    expK    foldChange      goTermDescription
                goTermID,pValue,qValue,obsK,expK,foldChange,goTermDescription = line.strip().split('\t')
                pValue = float(pValue)
                qValue = float(qValue)
                if qValue > 0:
                    qValueList.append(qValue)
                if pValue > 0:
                    pValueList.append(pValue)
    return(pValueList,qValueList)


def readGOTermFile(goTermFile,pValue_addToZero,qValue_addToZero):
    obsDataList = []
    expDataList = []
    fullDataList = []
    fdrThreshold = 0.05
    xValues = []
    yValues = []
    qValues = []
    logQValues = []
    annotationList = []
    obsKThreshold = 6
    negLogPValueThreshold = 3
    with open(goTermFile,'r') as F:
        for line in F:
            if 'goTermID' not in line:
                # goTermID        pValue  qValue  obsK    expK    foldChange      goTermDescription
                goTermID,pValue,qValue,obsK,expK,foldChange,goTermDescription = line.strip().split('\t')
                pValue = float(pValue)
                qValue = float(qValue)
                obsK = int(obsK)
                expK = float(expK)
                qValues.append(qValue)
                foldChange = float(foldChange)
                if pValue == 0.0 or pValue == 0:
                    pValue = pValue + pValue_addToZero
                # if qValue == 0.0 or qValue == 0:
                #     qValue = qValue + qValue_addToZero
                    # https://discuss.analyticsvidhya.com/t/methods-to-deal-with-zero-values-while-performing-log-transformation-of-variable/2431/3
                    # -math.log(1e-323, 10)
                    # math.log(x, base)
                logFoldChange = math.log(foldChange, 2)
                #logQValue = math.log(qValue, 10)
                #logQValues.append(logQValue)
                #negLogQValue = -math.log(qValue, 10)
                negLogPValue = -math.log(pValue, 10)
                xValues.append(logFoldChange)
                yValues.append(negLogPValue)
                obsDataList.append(float(obsK))
                expDataList.append(float(expK))
                fullDataList.append((goTermID,float(pValue),float(qValue),negLogPValue,float(obsK),float(expK),foldChange,logFoldChange,goTermDescription))
                if negLogPValue > negLogPValueThreshold:
                    if obsK > expK:
                        if obsK >= obsKThreshold:
                            annotationList.append((goTermID,float(pValue),float(qValue),negLogPValue,float(obsK),float(expK),foldChange,logFoldChange,goTermDescription))
                    if obsK < expK:
                        annotationList.append((goTermID,float(pValue),float(qValue),negLogPValue,float(obsK),float(expK),foldChange,logFoldChange,goTermDescription))
                        #print goTermID,pValue,qValue,obsK,expK,foldChange,goTermDescription
    #fullDataList = zip(fullDataList,xValues,yValues)
    #fullDataList.sort(key = lambda x:x[0][2])
    fullDataList.sort(key = lambda x:x[2], reverse=False)
    annotationList.sort(key = lambda x:x[2], reverse=False)
    return(obsDataList,expDataList,fullDataList,xValues,yValues,qValues,logQValues,annotationList)


def volcano(xValues,yValues,qValues,logQValues,fullDataList,annotationList,goTermCategory):
    OUT = open(goTermCategory + '_volcano_data.txt','w')
    OUT.write("goTermID\tpValue\tqValue\tnegLogPValue\tobsK\texpK\tfoldChange\tlogFoldChange\tgoTermDescription\n")
    plt.rcParams["figure.figsize"] = (10,8)

    #color_norm_list = colors.Normalize(min(qValues), max(qValues))
    if goTermCategory == 'bp':
        label = 'Biological Processes'
        #color = '#7D2E68'
        color = 'black'
        marker = 'o'
        #color_map = cm.ScalarMappable(norm=color_norm_list, cmap='Reds')
        color_map = plt.cm.get_cmap('Reds')
        color_map = color_map.reversed()
    if goTermCategory == 'cc':
        label = 'Cellular Components'
        #color = '#005E7C'
        color = 'black'
        marker = 'o'
        color_map = plt.cm.get_cmap('Blues')
        color_map = color_map.reversed()
        #marker = 's'
        #color_map = 'Blues'
        #color_map = color_map.reversed()
        #color_palette = plt.get_cmap('Blues')
        #color_map = cm.ScalarMappable(norm=color_norm_list, cmap='Blues')
    if goTermCategory == 'mf':
        label = 'Molecular Functions'
        #color = '#BA5A31'
        color = 'black'
        #marker = 'D'
        marker = 'o'
        color_map = plt.cm.get_cmap('Purples')
        color_map = color_map.reversed()
        #color_map = 'Purples'
        #color_map = color_map.reversed()
        #color_map = cm.ScalarMappable(norm=color_norm_list, cmap='Oranges')
        # plt.scatter(xValues,yValues, 'bo', ms=5, c=qValues, cmap=color_map, alpha=0.5, label=label, markeredgewidth = 0.0, marker=marker)
    plt.scatter(xValues,yValues, c=qValues, cmap=color_map, label=label, norm=colors.SymLogNorm(1e-6, linscale=1.0))
    plt.xlabel("log$_{2}$(Fold-change)", size=16)
    plt.ylabel("-log$_{10}$(P-value)", size=16)

    # these lists are already sorted 
    # sortedFullDataList = fullDataList
    # annotationList
    N = 15
    topObservationList = [x for index, x in enumerate(annotationList) if index < N]
    for goTermID,pValue,qValue,negLogPValue,obsK,expK,foldChange,logFoldChange,goTermDescription in topObservationList:
        #print goTermID,pValue,qValue,negLogPValue,obsK,expK,foldChange,logFoldChange,goTermDescription
        fullGOInfo = goTermID + ', ' + goTermDescription
        plt.annotate(fullGOInfo,(logFoldChange,negLogPValue), color=color)
        OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (goTermID,pValue,qValue,negLogPValue,obsK,expK,foldChange,logFoldChange,goTermDescription))
    #plt.xlim(-4,4)
    #plt.ylim(0,7)
    plt.legend()
    plt.colorbar()
    plt.tight_layout()
    plt.savefig(goTermCategory + '_volcano.pdf')
    plt.savefig(goTermCategory + '_volcano.png')
    plt.savefig(goTermCategory + '_volcano.svg')
    

########
# MAIN #
########

usage = "Usage: python " + sys.argv[0] + " <go term file> <go term category: bp, cc, or mf>\n"
if len(sys.argv) != 3:
    print usage
    print "\n"
    sys.exit()

goTermFile = sys.argv[1]
goTermCategory = sys.argv[2]

# biologicalProcesses = bp
# cellularComponents = cc
# molecularFunction = mf

pValueList,qValueList = retrieveNonZeroQValues(goTermFile)

# non-zero p- and q-values
smallestNonZeroPValue =min(pValueList)
pValue_addToZero = float(smallestNonZeroPValue)/ 2

smallestNonZeroQValue = min(qValueList)
qValue_addToZero = float(smallestNonZeroQValue) / 2

obsDataList,expDataList,fullDataList,xValues,yValues,qValues,logQValues,annotationList = readGOTermFile(goTermFile,pValue_addToZero,qValue_addToZero)
volcano(xValues,yValues,qValues,logQValues,fullDataList,annotationList,goTermCategory)
