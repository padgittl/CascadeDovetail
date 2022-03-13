# Pipeline for estimating time divergence using MCMCTree

# Extract single-copy ortholog CDS sequences for all species
### Orthogroups.tsv and Orthogroups_SingleCopyOrthologues.txt are generated from OrthoFinder
<code>python scripts/retrieveGenesFromOrthogroups.py Orthogroups.tsv Orthogroups_SingleCopyOrthologues.txt cannabis_sativa.fasta humulus_lupulus.fasta morus_notabilis.fasta parasponia_andersonii.fasta prunus_persica.fasta trema_orientale.fasta vitis_vinifera.fasta ziziphus_jujuba.fasta</code>

# Generate multiple sequence alignment of single-copy orthologs with MACSE
MACSE alignSequences and exportAlignments
### orthogroup (OG)
<code>java -jar macse_v2.03.jar -prog alignSequences -seq OG.fasta -out_NT OG_NT.fasta -out_AA OG_AA.fasta</code>
<code>java -jar macse_v2.03.jar -prog exportAlignment -align OG_NT.fasta -codonForFinalStop --- -codonForInternalStop +++ -codonForInternalFS +++ -charForRemainingFS + -out_NT OG_NT.expAlign.fasta -out_AA OG_AA.expAlign.fasta</code>

# Rename the defline with the geneID for each OG-specific fasta file to contain the species name
<code>ls -1 ../macse/exportAlignments/*_NT.expAlign.fasta > fastaFileList.txt</code>
<code>python scripts/renameMSA.py fastaFileList.txt</code>

# Concatenate the third position of four-fold degenerate codons
ls -1 ../renamedMSA/*fasta > msaFileList.txt
python scripts/concatenateFDTvPositions.py msaFileList.txt
cat *fasta > allSpecies.thirdCodonPos.fa
