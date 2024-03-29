# Pipeline for estimating time divergence using MCMCTree

## 1. Run OrthoFinder

<code>orthofinder -S diamond -M msa -f outputDirectory -n resultsName -t 8 -a 1</code>


## 2. Extract single-copy ortholog CDS sequences for all species

Orthogroups.tsv and Orthogroups_SingleCopyOrthologues.txt are generated from OrthoFinder

<code>python scripts/retrieveGenesFromOrthogroups.py Orthogroups.tsv Orthogroups_SingleCopyOrthologues.txt cannabis_sativa.fasta humulus_lupulus.fasta morus_notabilis.fasta parasponia_andersonii.fasta prunus_persica.fasta trema_orientale.fasta vitis_vinifera.fasta ziziphus_jujuba.fasta</code>

## 3. Generate multiple sequence alignment of single-copy orthologs with MACSE

--> orthogroup (OG)

<code>java -jar macse_v2.03.jar -prog alignSequences -seq OG.fasta -out_NT OG_NT.fasta -out_AA OG_AA.fasta</code>

<code>java -jar macse_v2.03.jar -prog exportAlignment -align OG_NT.fasta -codonForFinalStop --- -codonForInternalStop +++ -codonForInternalFS +++ -charForRemainingFS + -out_NT OG_NT.expAlign.fasta -out_AA OG_AA.expAlign.fasta</code>

## 4. Rename the geneID defline in each OG-specific fasta file to list the species name instead

<code>ls -1 ../macse/exportAlignments/*_NT.expAlign.fasta > fastaFileList.txt</code>

<code>python scripts/renameMSA.py fastaFileList.txt</code>

## 5. Concatenate the third position of four-fold degenerate codons

<code>ls -1 ../renamedMSA/*fasta > msaFileList.txt</code>

<code>python scripts/concatenateFDTvPositions.py msaFileList.txt</code>

<code>cat *fasta > allSpecies.thirdCodonPos.fa</code>

## 5. (alternate approaches using r8s and treePL) Concatenate full codons and create 'partition.txt' file for RAxML

<code>ls -1 ../renamedMSA/*fasta > msaFileList.txt</code>

<code>python scripts/concatenateFullCodons.py msaFileList.txt</code>

<code>cat *fasta > allSpecies.thirdCodonPos.fa</code>

### see 'partitions_example.txt' for an example of the format for the partitions file
