import sys, re, os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import math


###############
# SUBROUTINES #
###############


def readLengthsfile(lengthsFile):
    scaffoldLengthDict = {}
    with open(lengthsFile,'r') as F:
        # Scaffold_1531 476495644
        for line in F:
            scaffoldID,scaffoldLen = line.strip().split(' ')
            scaffoldLengthDict[scaffoldID] = int(scaffoldLen)
            #print scaffoldID,scaffoldLen
    return(scaffoldLengthDict)


# Scaffold_1531   MAKER0000001.t1 HUMLU_CAS0000001.t1.p1  13449   23487   13449   23487
def readGeneMapFile(geneMapFile,scaffoldLengthDict):
    filteredGeneMapDict = {}
    with open(geneMapFile,'r') as F:
        for line in F:
            # OUT.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (scaffoldID,transcriptID,newGeneID,geneStart,geneStop,cdsStart,cdsStop))
            scaffoldID,oldGeneID,newGeneID,geneStart,geneStop,cdsStart,cdsStop = line.strip().split('\t')
            if scaffoldID in scaffoldLengthDict:
                scaffoldLen = scaffoldLengthDict[scaffoldID]
                filteredGeneMapDict[newGeneID] = scaffoldID
    return(filteredGeneMapDict)


def readCannabisFileList(cannabisFileList):
    cannabisDataDict = {}
    with open(cannabisFileList,'r') as FL:
        for line in FL:
            fullFileName = line.strip()
            fileName = os.path.basename(fullFileName)
            # HUMLU_CAS0000006.t1_vs_XP_030499033.1.yn00.txt
            filePrefix,extraStuff = fileName.split('.yn00.txt')
            dataList = parse_yn00_file(fullFileName)
            publication = 'Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43'
            for i in range(len(dataList)):
                if publication in dataList[i]:
                    #print dataList[i]
                    # seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE
                    #print dataList[i+2]
                    #print dataList[i+3]
                    if filePrefix not in cannabisDataDict:
                        cannabisDataDict[filePrefix] = dataList[i+3]
            #else:
            #    print filePrefix
    #print len(cannabisDataDict)
    return(cannabisDataDict,len(cannabisDataDict))


def readHopFileList(hopFileList,filteredGeneMapDict):
    hopDataDict = {}
    with open(hopFileList,'r') as FL:
        for line in FL:
            fullFileName = line.strip()
            fileName = os.path.basename(fullFileName)
            # HUMLU_CAS0000006.t1_vs_XP_030499033.1.yn00.txt
            filePrefix,extraStuff = fileName.split('.yn00.txt')
            geneID1,geneID2 = filePrefix.split('_vs_')
            if geneID1 in filteredGeneMapDict and geneID2 in filteredGeneMapDict:
                dataList = parse_yn00_file(fullFileName)
                publication = 'Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43'
                for i in range(len(dataList)):
                    if publication in dataList[i]:
                        #print dataList[i]
                        # seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE
                        #print dataList[i+2]
                        #print dataList[i+3]
                        if filePrefix not in hopDataDict:
                            hopDataDict[filePrefix] = dataList[i+3]
            #else:
            #    print filePrefix
    #print len(hopDataDict)
    return(hopDataDict,len(hopDataDict))


def readHop_vs_CannabisFileList(hop_vs_cannabisFileList,filteredGeneMapDict):
    hop_vs_cannabisDataDict = {}
    with open(hop_vs_cannabisFileList,'r') as FL:
        for line in FL:
            fullFileName = line.strip()
            fileName = os.path.basename(fullFileName)
            # HUMLU_CAS0000006.t1_vs_XP_030499033.1.yn00.txt
            filePrefix,extraStuff = fileName.split('.yn00.txt')
            geneID1,geneID2 = filePrefix.split('_vs_')
            dataList = parse_yn00_file(fullFileName)
            publication = 'Yang Z, Nielsen R (2000) Estimating synonymous and nonsynonymous substitution rates under realistic evolutionary models. Mol. Biol. Evol. 17:32-43'
            if 'HUMLU_CAS' in geneID1:
                if geneID1 in filteredGeneMapDict:
                    for i in range(len(dataList)):
                        if publication in dataList[i]:
                            #print dataList[i]
                            # seq. seq.     S       N        t   kappa   omega     dN +- SE    dS +- SE
                            #print dataList[i+2]
                            #print dataList[i+3]
                            if filePrefix not in hop_vs_cannabisDataDict:
                                hop_vs_cannabisDataDict[filePrefix] = dataList[i+3]
            if 'HUMLU_CAS' in geneID2:
                if geneID2 in filteredGeneMapDict:
                    for i in range(len(dataList)):
                        if publication in dataList[i]:
                            if filePrefix not in hop_vs_cannabisDataDict:
                                hop_vs_cannabisDataDict[filePrefix] = dataList[i+3]
            #else:
            #    print filePrefix
    #print len(hop_vs_cannabisDataDict)
    return(hop_vs_cannabisDataDict,len(hop_vs_cannabisDataDict))
            

def parse_yn00_file(yn00File):
    dataList = []
    with open(yn00File,'r') as Y:
        for line in Y:
            line = line.strip()
            if not len(line.strip()) == 0:
                dataList.append(line)
    return(dataList)


def retrieve_dn_ds_values(dataDict,dsThresholdLower,dsThresholdUpper,description):
    lambda_value = 6.1*10**-9
    # lambda_value = 2.1*10**-9
    #timeKsLowerThreshold = 0.05
    #timeKsLowerThreshold = 0.0
    timeKsLowerThreshold = 0.01
    timeKsUpperThreshold1 = 1.0
    timeKsUpperThreshold2 = 2.0
    #timeKsUpperThreshold = 3.0
    pairCountForTimeCalculation1 = 0
    pairCountForTimeCalculation2 = 0

    dn_list = []
    ds_list = []
    log_ds_list = []
    dn_ds_list = []
    time_list1 = []
    time_list2 = []
    # '_dsLower' + str(dsThresholdLower) + "_dsUpper " + str(dsThresholdUpper) + ' 
    DNDS = open(description + '_kaks_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.txt','w')
    DS = open(description + '_ks_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.txt','w')
    LOG_DS = open(description + '_log_ks_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.txt','w')
    UNDEFINED_DN_DS = open(description + '_undefined_kaks_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.txt','w')

    # more-stringent thresholds for mixture model
    # timeKsLowerThreshold = 0.01
    # timeKsUpperThreshold = 2.0
    DNDS2 = open(description + '_kaks_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold2) + '.txt','w')
    DS2 = open(description + '_ks_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold2) + '.txt','w')
    LOG_DS2 = open(description + '_log_ks_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold2) + '.txt','w')
    TIME2 = open(description + '_time_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold2) + '.txt','w')
    
    DNDS1 = open(description + '_kaks_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold1) + '.txt','w')
    DS1 = open(description + '_ks_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold1) + '.txt','w')
    LOG_DS1 = open(description + '_log_ks_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold1) + '.txt','w')
    TIME1 = open(description + '_time_dsLower' + str(timeKsLowerThreshold) + "_dsUpper" + str(timeKsUpperThreshold1) + '.txt','w')

    for genePair in dataDict:
        dnds_info_list = []
        dnds_info = dataDict[genePair]
        dnds_info = dnds_info.split(' ')
        for i in dnds_info:
            if i != '':
                dnds_info_list.append(i)
        # ['2', '1', '438.3', '1409.7', '0.2533', '2.7644', '0.1714', '0.0393', '+-', '0.0054', '0.2296', '+-', '0.0274']
        dn = dnds_info_list[7]
        ds = dnds_info_list[10]
        dn = float(dn)
        ds = float(ds)
        if ds > 0:
            dn_ds = round(float(dn) / ds, 5)
            # dsThresholdLower,dsThresholdUpper
            if ds > dsThresholdLower and ds < dsThresholdUpper:
                dn_list.append(dn)
                ds_list.append(ds)
                log_ds = math.log10(ds)
                log_ds_list.append(log_ds)
                dn_ds_list.append(dn_ds)
                # Draft genome sequence of the mulberry tree Morus notabilis
                # Lynch, M. & Conery, J. S. The evolutionary fate and consequences of duplicate genes.
                # divergence date: T=Ks/2*lambda with lambda=6.1*10-9
                if ds >= timeKsLowerThreshold:
                    if ds <= timeKsUpperThreshold2:
                        pairCountForTimeCalculation2 += 1
                        time = (float(ds)/(2*lambda_value))/1000000
                        time_list2.append(time)
                        DS2.write("%s\t%s\n" % (genePair,ds))
                        LOG_DS2.write("%s\t%s\n" % (genePair,log_ds))
                        DNDS2.write("%s\t%s\n" % (genePair,dn_ds))
                        TIME2.write("%s\t%s\n" % (genePair,time))
                    if ds <= timeKsUpperThreshold1:
                        pairCountForTimeCalculation1 += 1
                        time = (float(ds)/(2*lambda_value))/1000000
                        time_list1.append(time)
                        DS1.write("%s\t%s\n" % (genePair,ds))
                        LOG_DS1.write("%s\t%s\n" % (genePair,log_ds))
                        DNDS1.write("%s\t%s\n" % (genePair,dn_ds))
                        TIME1.write("%s\t%s\n" % (genePair,time))
                #print ds
                DS.write("%s\t%s\n" % (genePair,ds))
                LOG_DS.write("%s\t%s\n" % (genePair,log_ds))
                DNDS.write("%s\t%s\n" % (genePair,dn_ds))
                #else:
                #    print genePair,ds
        else:
            UNDEFINED_DN_DS.write("%s\t%s\t%s\n" % (genePair,dn,ds))
            #    print genePair,dn,ds
    return(dn_list,ds_list,log_ds_list,dn_ds_list,time_list1,time_list2,pairCountForTimeCalculation1,pairCountForTimeCalculation2)
        

def createLogHist(can_log_ds_list,hop_log_ds_list,hop_vs_can_log_ds_list,dsBinSize,dsThresholdLower,dsThresholdUpper,numberCannabisGenePairs,numberHopGenePairs,numberHop_vs_CannabisGenePairs):
    binWidth = dsBinSize

    all_log_ds = can_log_ds_list + hop_log_ds_list + hop_vs_can_log_ds_list
    bins = np.arange(min(all_log_ds),max(all_log_ds),binWidth)
    counts1, bins1, bars1 = plt.hist(can_log_ds_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp')
    counts2, bins2, bars2 = plt.hist(hop_log_ds_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop')
    counts3, bins3, bars3 = plt.hist(hop_vs_can_log_ds_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp')
    plt.xlabel('Log-transformed Ks',size=16)
    plt.ylabel('Count',size=16)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    # _dsLower' + str(dsThresholdLower) + "_dsUpper " + str(dsThresholdUpper) + '
    plt.savefig('logKsCountHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('logKsCountHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('logKsCountHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()

    # density
    all_log_ds = can_log_ds_list + hop_log_ds_list + hop_vs_can_log_ds_list
    bins = np.arange(min(all_log_ds),max(all_log_ds),binWidth)
    counts1, bins1, bars1 = plt.hist(can_log_ds_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp', density=True, linewidth=2)
    counts2, bins2, bars2 = plt.hist(hop_log_ds_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop', density=True, linewidth=2)
    counts3, bins3, bars3 = plt.hist(hop_vs_can_log_ds_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp', density=True, linewidth=2)
    plt.xlabel('Log-transformed Ks',size=16)
    plt.ylabel('Density',size=16)
    plt.legend(frameon=False,fontsize=14, loc="upper left")
    plt.tight_layout()
    plt.savefig('logKsDensityHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('logKsDensityHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('logKsDensityHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()


def createHist(can_ds_list,can_dn_ds_list,hop_ds_list,hop_dn_ds_list,hop_vs_can_ds_list,hop_vs_can_dn_ds_list,can_time_list,hop_time_list,hop_vs_can_time_list,dsBinSize,dndsBinSize,dsThresholdLower,dsThresholdUpper,numberCannabisGenePairs,numberHopGenePairs,numberHop_vs_CannabisGenePairs,can_pairCountForTimeCalculation,hop_pairCountForTimeCalculation,hop_vs_can_pairCountForTimeCalculation):
    binWidth = dndsBinSize

    # dnds
    complete_dnds = can_dn_ds_list + hop_dn_ds_list + hop_vs_can_dn_ds_list
    bins = np.arange(min(complete_dnds),max(complete_dnds),binWidth)
    counts1, bins1, bars1 = plt.hist(can_dn_ds_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp')
    counts2, bins2, bars2 = plt.hist(hop_dn_ds_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop')
    counts3, bins3, bars3 = plt.hist(hop_vs_can_dn_ds_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp')
    # https://stackoverflow.com/questions/42045767/how-can-i-change-the-x-axis-in-matplotlib-so-there-is-no-white-space
    #plt.margins(x=0)
    plt.xlabel('Ka/Ks',size=16)
    plt.ylabel('Count',size=16)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('kaksCountHist_binSize' + str(dndsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('kaksCountHist_binSize' + str(dndsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('kaksCountHist_binSize' + str(dndsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()

    # density
    # numberCannabisGenePairs,numberHopGenePairs,numberHop_vs_CannabisGenePairs
    complete_dnds = can_dn_ds_list + hop_dn_ds_list + hop_vs_can_dn_ds_list
    bins = np.arange(min(complete_dnds),max(complete_dnds),binWidth)
    counts1, bins1, bars1 = plt.hist(can_dn_ds_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp', density=True)
    counts2, bins2, bars2 = plt.hist(hop_dn_ds_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop', density=True)
    counts3, bins3, bars3 = plt.hist(hop_vs_can_dn_ds_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp', density=True)
    plt.xlabel('Ka/Ks',size=16)
    plt.ylabel('Density',size=16)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('kaksDensityHist_binSize' + str(dndsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('kaksDensityHist_binSize' + str(dndsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('kaksDensityHist_binSize' + str(dndsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()


    # plot dS (same thing as Ks)
    dsBinWidth = dsBinSize
    complete_ds = can_ds_list + hop_ds_list + hop_vs_can_ds_list
    bins = np.arange(min(complete_ds),max(complete_ds),dsBinWidth)
    dsCounts1, dsBins1, dsBars1 = plt.hist(can_ds_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp', linewidth=2)
    dsCounts2, dsBins2, dsBars2 = plt.hist(hop_ds_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop', linewidth=2)
    dsCounts3, dsBins3, dsBars3 = plt.hist(hop_vs_can_ds_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp', linewidth=2)

    #plt.margins(x=0)
    plt.xlabel('Ks',size=16)
    plt.ylabel('Count',size=16)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('KsCountHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('KsCountHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('KsCountHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()

    # density
    # numberCannabisGenePairs,numberHopGenePairs,numberHop_vs_CannabisGenePairs
    dsBinWidth = dsBinSize
    complete_ds = can_ds_list + hop_ds_list + hop_vs_can_ds_list
    bins = np.arange(min(complete_ds),max(complete_ds),dsBinWidth)
    dsCounts1, dsBins1, dsBars1 = plt.hist(can_ds_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp, ' + str(numberCannabisGenePairs) + " gene pairs", density=True, linewidth=2)
    dsCounts2, dsBins2, dsBars2 = plt.hist(hop_ds_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop, ' + str(numberHopGenePairs) + " gene pairs", density=True, linewidth=2)
    dsCounts3, dsBins3, dsBars3 = plt.hist(hop_vs_can_ds_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp, ' + str(numberHop_vs_CannabisGenePairs) + " gene pairs", density=True, linewidth=2)
    plt.xlabel('Ks',size=16)
    plt.ylabel('Density',size=16)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('KsDensityHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('KsDensityHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('KsDensityHist_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()


    # plot dS (same thing as Ks) - x-axis log-scale, non-transformed Ks values
    # supplementary info for "A reference genome for pea provides insight into legume genome evolution"
    # https://www.nature.com/articles/s41588-019-0480-1
    dsBinWidth = dsBinSize
    complete_ds = can_ds_list + hop_ds_list + hop_vs_can_ds_list
    bins = np.arange(min(complete_ds),max(complete_ds),dsBinWidth)
    dsCounts1, dsBins1, dsBars1 = plt.hist(can_ds_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp')
    dsCounts2, dsBins2, dsBars2 = plt.hist(hop_ds_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop')
    dsCounts3, dsBins3, dsBars3 = plt.hist(hop_vs_can_ds_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp')
    plt.xscale('log')
    plt.xlabel('Log-scale x-axis, non-log Ks values',size=14)
    plt.ylabel('Count',size=14)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('KsCountHist_xaxisLog_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('KsCountHist_xaxisLog_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('KsCountHist_xaxisLog_binSize' + str(dsBinSize) + '_noFS_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()

    # count
    timeBinWidth = 1
    all_time_lists = can_time_list + hop_time_list + hop_vs_can_time_list
    bins = np.arange(min(all_time_lists),max(all_time_lists),timeBinWidth)
    timeCounts1, timeBins1, timeBars1 = plt.hist(can_time_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp', linewidth=1)
    timeCounts2, timeBins2, timeBars2 = plt.hist(hop_time_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop', linewidth=1)
    timeCounts3, timeBins3, timeBars3 = plt.hist(hop_vs_can_time_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp', linewidth=1)
    hopMaxCount = 0
    for i in range(len(timeCounts2)):
        if timeCounts2[i] > hopMaxCount:
            hopMaxCount = timeCounts2[i]
            hopMaxIndex = i
    hopAvgMaxBin = float(timeBins2[hopMaxIndex] + timeBins2[hopMaxIndex+1]) / 2
    plt.axvline(hopAvgMaxBin, linestyle='--', c='black')
    plt.text(hopAvgMaxBin+20, 100, str(round(hopAvgMaxBin, 1)))
    
    maxCount = 0
    for i in range(len(timeCounts3)):
        if timeCounts3[i] > maxCount:
            maxCount = timeCounts3[i]
            maxIndex = i
    avgMaxBin = float(timeBins3[maxIndex] + timeBins3[maxIndex+1]) / 2
    plt.axvline(avgMaxBin, linestyle='--', c='black')
    plt.text(avgMaxBin+20, 50, str(round(avgMaxBin, 1)))

    plt.xlabel('Time (mya)',size=14)
    plt.ylabel('Count',size=14)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    # '_dsLower' + str(dsThresholdLower) + "_dsUpper " + str(dsThresholdUpper) + '
    plt.savefig('KsBasedTimeCountHist_bin' + str(timeBinWidth) + '_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('KsBasedTimeCountHist_bin' + str(timeBinWidth) + '_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('KsBasedTimeCountHist_bin' + str(timeBinWidth) + '_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()

    # can_pairCountForTimeCalculation,hop_pairCountForTimeCalculation,hop_vs_can_pairCountForTimeCalculation
    # density
    # # numberCannabisGenePairs,numberHopGenePairs,numberHop_vs_CannabisGenePairs 
    timeBinWidth = 4
    all_time_lists = can_time_list + hop_time_list + hop_vs_can_time_list
    bins = np.arange(min(all_time_lists),max(all_time_lists),timeBinWidth)
    timeCounts1, timeBins1, timeBars1 = plt.hist(can_time_list, bins = bins, histtype = 'step', color ='#0471A6', label = 'Hemp vs Hemp, ' + str(can_pairCountForTimeCalculation) + " gene pairs", density=True, linewidth=2, alpha=0.8)
    timeCounts2, timeBins2, timeBars2 = plt.hist(hop_time_list, bins = bins, histtype = 'step', color ='#A41623', label = 'Hop vs Hop, ' + str(hop_pairCountForTimeCalculation) + " gene pairs", density=True, linewidth=2, alpha=0.8)
    timeCounts3, timeBins3, timeBars3 = plt.hist(hop_vs_can_time_list, bins = bins, histtype = 'step', color ='#F85E00', label = 'Hop vs Hemp, ' + str(hop_vs_can_pairCountForTimeCalculation) + " gene pairs", density=True, linewidth=2, alpha=0.8)

    hopMaxCount = 0
    for i in range(len(timeCounts2)):
        if timeCounts2[i] > hopMaxCount:
            hopMaxCount = timeCounts2[i]
            hopMaxIndex = i
    hopAvgMaxBin = float(timeBins2[hopMaxIndex] + timeBins2[hopMaxIndex+1]) / 2
    plt.axvline(hopAvgMaxBin, linestyle='--', c='black')
    # plt.text(hopAvgMaxBin+20, 0.1, str(round(hopAvgMaxBin, 1)))
    plt.text(hopAvgMaxBin+40, 0.04, str(round(hopAvgMaxBin, 1)))
    
    maxCount = 0
    for i in range(len(timeCounts3)):
        if timeCounts3[i] > maxCount:
            maxCount = timeCounts3[i]
            maxIndex = i
    avgMaxBin = float(timeBins3[maxIndex] + timeBins3[maxIndex+1]) / 2
    plt.axvline(avgMaxBin, linestyle='--', c='black')
    # plt.text(avgMaxBin+30, 0.1, str(round(avgMaxBin, 1)))
    plt.text(avgMaxBin+40, 0.04, str(round(avgMaxBin, 1))) 

    plt.xlabel('Time (mya)',size=14)
    plt.ylabel('Density',size=14)
    plt.legend(frameon=False,fontsize=14)
    plt.tight_layout()
    plt.savefig('KsBasedTimeDensityHist_bin' + str(timeBinWidth) + '_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.png', dpi=600)
    plt.savefig('KsBasedTimeDensityHist_bin' + str(timeBinWidth) + '_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.pdf')
    plt.savefig('KsBasedTimeDensityHist_bin' + str(timeBinWidth) + '_dsLower' + str(dsThresholdLower) + "_dsUpper" + str(dsThresholdUpper) + '.svg')
    plt.close()


############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <cannabis file list> <hop file list> <hop vs cannabis file list> <lengths file> <gene map file> <ds bin size> <dnds bin size> <ds lower threshold> <ds upper threshold>\n"
if len(sys.argv) != 10:
    print usage
    sys.exit()

cannabisFileList = sys.argv[1]
hopFileList = sys.argv[2]
hop_vs_cannabisFileList = sys.argv[3]
lengthsFile = sys.argv[4]
geneMapFile = sys.argv[5]
dsBinSize = sys.argv[6]
dndsBinSize = sys.argv[7]
dsThresholdLower = sys.argv[8]
dsThresholdUpper = sys.argv[9]

# https://www.protocols.io/view/introduction-to-calculating-dn-ds-ratios-with-code-qhwdt7e?step=4
# https://onlinelibrary.wiley.com/doi/full/10.1111/tpj.13801 "All Ks values less than or equal to 0.1 excluded" -- further, only WGD signals between 0.235 and 2 were investigated
# ds upper bound of 2 --> https://www.nature.com/articles/nature16548 // The genome of the seagrass Zostera marina reveals angiosperm adaptation to the sea
#dsThresholdLower = 0.05
#dsThresholdLower = 0.0
#dsThresholdUpper = 5
dsThresholdLower = float(dsThresholdLower)
dsThresholdUpper = float(dsThresholdUpper)
dsBinSize = float(dsBinSize)
dndsBinSize  = float(dndsBinSize)

scaffoldLengthDict = readLengthsfile(lengthsFile)
filteredGeneMapDict = readGeneMapFile(geneMapFile,scaffoldLengthDict)

cannabisDataDict,numberCannabisGenePairs = readCannabisFileList(cannabisFileList)
hopDataDict,numberHopGenePairs = readHopFileList(hopFileList,filteredGeneMapDict)
hop_vs_cannabisDataDict,numberHop_vs_CannabisGenePairs = readHop_vs_CannabisFileList(hop_vs_cannabisFileList,filteredGeneMapDict)

can_dn_list,can_ds_list,can_log_ds_list,can_dn_ds_list,can_time_list1,can_time_list2,can_pairCountForTimeCalculation1,can_pairCountForTimeCalculation2 = retrieve_dn_ds_values(cannabisDataDict,dsThresholdLower,dsThresholdUpper,'cannabis')
hop_dn_list,hop_ds_list,hop_log_ds_list,hop_dn_ds_list,hop_time_list1,hop_time_list2,hop_pairCountForTimeCalculation1,hop_pairCountForTimeCalculation2 = retrieve_dn_ds_values(hopDataDict,dsThresholdLower,dsThresholdUpper,'hop')
hop_vs_can_dn_list,hop_vs_can_ds_list,hop_vs_can_log_ds_list,hop_vs_can_dn_ds_list,hop_vs_can_time_list1,hop_vs_can_time_list2,hop_vs_can_pairCountForTimeCalculation1,hop_vs_can_pairCountForTimeCalculation2 = retrieve_dn_ds_values(hop_vs_cannabisDataDict,dsThresholdLower,dsThresholdUpper,'hop_vs_cannabis')

createHist(can_ds_list,can_dn_ds_list,hop_ds_list,hop_dn_ds_list,hop_vs_can_ds_list,hop_vs_can_dn_ds_list,can_time_list2,hop_time_list2,hop_vs_can_time_list2,dsBinSize,dndsBinSize,dsThresholdLower,dsThresholdUpper,numberCannabisGenePairs,numberHopGenePairs,numberHop_vs_CannabisGenePairs,can_pairCountForTimeCalculation2,hop_pairCountForTimeCalculation2,hop_vs_can_pairCountForTimeCalculation2)

createLogHist(can_log_ds_list,hop_log_ds_list,hop_vs_can_log_ds_list,dsBinSize,dsThresholdLower,dsThresholdUpper,numberCannabisGenePairs,numberHopGenePairs,numberHop_vs_CannabisGenePairs)
