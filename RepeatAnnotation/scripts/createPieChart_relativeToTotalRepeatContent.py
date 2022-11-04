import sys, re, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.mlab as mlab
import matplotlib.pyplot as plt
plt.switch_backend('agg')
import numpy as np

###############
# SUBROUTINES #
###############

#type    count   percentage
#DNA     7757346 0.292   
#LTR     2574127002      97.04   
#LINE    1457927 0.055   
#SINE    426     0.0     
#Satellite       0       0.0     
#Simple  62951315        2.373   
#Mobile Element  5369783 0.202   
#rRNA    247508  0.009   
#Other   743493  0.028   
#Total repeat content    2652654800      99.999 

def readValuesFile(valuesFile):
    with open(valuesFile,'r') as V:
        for line in V:
            if 'DNA' in line:
                # DNA     7757346 0.292
                repType,dnaCount,dnaPercentage = line.strip().split('\t')
                dnaPercentage = float(dnaPercentage)
            if 'LTR' in line:
                repType,ltrCount,ltrPercentage = line.strip().split('\t')
                ltrPercentage = float(ltrPercentage)
            if 'LINE' in line:
                repType,lineCount,linePercentage = line.strip().split('\t')
                linePercentage = float(linePercentage)
            if 'SINE' in line:
                repType,sineCount,sinePercentage = line.strip().split('\t')
                sinePercentage = float(sinePercentage)
            if 'Satellite' in line:
                repType,satelliteCount,satellitePercentage = line.strip().split('\t')
                satellitePercentage = float(satellitePercentage)
            if 'Simple' in line:
                repType,simpleCount,simplePercentage = line.strip().split('\t')
                simplePercentage = float(simplePercentage)
            if 'Mobile' in line:
                repType,mobileElementCount,mobileElementPercentage = line.strip().split('\t')
                mobileElementPercentage = float(mobileElementPercentage)
            if 'rRNA' in line:
                repType,rRNACount,rRNAPercentage = line.strip().split('\t')
                rRNAPercentage = float(rRNAPercentage)
            if 'Other' in line:
                repType,otherCount,otherPercentage = line.strip().split('\t')
                otherPercentage = float(otherPercentage)
            if 'Total Repeat Content' in line:
                repType,totalRepCount,totalRepPercentage = line.strip().split('\t')
                totalRepPercentage = float(totalRepPercentage)
    return(dnaPercentage,ltrPercentage,linePercentage,sinePercentage,satellitePercentage,simplePercentage,mobileElementPercentage,rRNAPercentage,otherPercentage,totalRepPercentage)
                

########
# MAIN #
########

usage = "Usage: " + sys.argv[0] + " <hop repeat file> \n"
if len(sys.argv) != 2:
    print usage
    sys.exit()

hopFile = sys.argv[1]

#fullFileName = arabidopsisGFF.strip()
#reducedFileName = os.path.basename(fullFileName)
#baseName,fileExt = os.path.splitext(reducedFileName)

dnaPercentageH,ltrPercentageH,linePercentageH,sinePercentageH,satellitePercentageH,simplePercentageH,mobileElementPercentageH,rRNAPercentageH,otherPercentageH,totalRepPercentageH = readValuesFile(hopFile)

pieRadius = 0.8

hopLabels = ['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile Element']
hopSizes = [ltrPercentageH,dnaPercentageH,linePercentageH,sinePercentageH,satellitePercentageH,simplePercentageH,mobileElementPercentageH]


print("Hop DNA\t%s\t" % (dnaPercentageH))
print("Hop LTR\t%s\t" % (ltrPercentageH))
print("Hop LINE\t%s\t" % (linePercentageH))
print("Hop SINE\t%s\t" % (sinePercentageH))
print("Hop Satellite\t%s\t" % (satellitePercentageH))
print("Hop Simple\t%s\t" % (simplePercentageH))
print("Hop Mobile Element\t%s\t" % (mobileElementPercentageH))
print("Hop rRNA\t%s\t" % (rRNAPercentageH))
print("Hop Other\t%s\t" % (otherPercentageH))


#colors = ['#313695','#4575b4','#74add1','#abd9e9','#ffffbf','#fee090','#fdae61','#f46d43','#d73027']
colors = ['#4575b4','#74add1','#abd9e9','#ffffbf','#fee090','#fdae61','#f46d43']
# colors = ['#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43']
#['LTR','DNA','LINE','SINE','Satellite','Simple','Mobile element']

#colors = ['#4575b4','#74add1','#abd9e9','#e0f3f8','#fee090','#fdae61','#f46d43']
 
#plt.rcParams["figure.frameon"] = False
plt.rcParams['axes.titlesize'] = 20
# plt.rcParams['font.size'] = 16
plt.rcParams['font.size'] = 16


#figWidth = 100
#figHeight = 30
figWidth = 20
figHeight = 10
fig, ax1 = plt.subplots(figsize=(10,5))
fig_size = plt.rcParams["figure.figsize"]
fig_size[0] = figWidth
fig_size[1] = figHeight
plt.rcParams["figure.figsize"] = fig_size

ax1 = plt.subplot(131)
ax1.set_aspect('equal')
        
hopExplode = 0
hopExplodeList = []
filteredHopSizes = []
filteredHopColors = []
hopStuff = []

for percValue,color in zip(hopSizes,colors):
    percValue = round(percValue,1)
    hopStuff.append((percValue,color))
hopStuff.sort(key=lambda x:x[0], reverse=True)
for percValue,color in hopStuff:
    if 0 < percValue and percValue < 50:
        #hopExplode += 0.05
        hopExplode += 0.075
        hopExplodeList.append(hopExplode)
        filteredHopSizes.append(percValue)
        filteredHopColors.append(color)
    else:
        hopExplode = 0.0
        hopExplodeList.append(hopExplode)
        filteredHopSizes.append(percValue)
        filteredHopColors.append(color)

# https://stackoverflow.com/questions/34035427/conditional-removal-of-labels-in-matplotlib-pie-chart
def filterAutopct(perc):
    return ('%1.1f%%' % perc) if perc > 0 else ''

hopWedges, hopTexts = ax1.pie(filteredHopSizes, colors=filteredHopColors, shadow=False, radius = pieRadius,startangle=-45, explode=hopExplodeList)

kw = dict(arrowprops=dict(arrowstyle="-"))

for i, p in enumerate(hopWedges):
    ang = (p.theta2 - p.theta1)/2. + p.theta1
    y = np.sin(np.deg2rad(ang))
    x = np.cos(np.deg2rad(ang))
    horizontalalignment = {-1: "right", 1: "left"}[int(np.sign(x))]
    connectionstyle = "angle,angleA=0,angleB={}".format(ang)
    kw["arrowprops"].update({"connectionstyle": connectionstyle})
    ax1.annotate(str(filteredHopSizes[i]) + "%", xy=(x, y), xytext=(0.9*np.sign(x), 1.2*y),
                 horizontalalignment=horizontalalignment,**kw)

plt.title('Hop Cascade Dovetail Assembly', fontsize=12)

plt.tight_layout()

#plt.title('Percentage of repeat types', fontsize=12)

plt.savefig("dovetailPercentRepeatPieChart.png", dpi=600)
plt.savefig("dovetailPercentRepeatPieChart.pdf")
plt.savefig("dovetailPercentRepeatPieChart.svg")



