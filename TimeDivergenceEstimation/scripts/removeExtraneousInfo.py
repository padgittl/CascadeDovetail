import sys, re, os


###############
# SUBROUTINES #
###############

'''
[Printing tree 1]

[TREE DESCRIPTION of tree tree_1]

tree tree_1 = ((((((cannabis_sativa:28.000000,humulus_lupulus:28.000000)N5:10.495251,(parasponia_andersonii:38.493576,trema_orientale:38.493576):0.001676):57.511430,morus_notabilis:96.006681)N3:0.000646,ziziphus_jujuba:96.007328)N2:0.016833,prunus_persica:96.024161)N1:27.168193,vitis_vinifera:123.192353)N0;
[Printing tree 2]

'''

def cleanFile(inputFile):
    OUT = open('r8s_bootstrap_cleaned_output.txt','w')
    with open(inputFile,'r') as F:
        for line in F:
            if 'TREE DESCRIPTION' not in line:
                if 'tree_' in line:
                    line = line.strip()
                    extraStuff,treeInfo = line.strip().split("=")
                    tree = treeInfo.strip()
                    OUT.write("%s\n" % (tree))
            

############
# MAIN  ####
############

usage = "Usage: " + sys.argv[0] + " <input file>"
if len(sys.argv) != 2:
    print usage
    sys.exit()

inputFile = sys.argv[1]

cleanFile(inputFile)
