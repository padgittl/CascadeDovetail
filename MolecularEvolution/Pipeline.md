# Pipeline for assessing molecular evolution in genes from syntenic blocks

## blastall for MCScanX

<code>formatdb -i hop.pep.fasta -p T -o T</code>

<code>formatdb -i can.pep.fasta -p T -o T</code>

<code>blastall -p blastp -i hop.pep.fasta -d hopBlastpDB/hop.pep.fasta -e 1e-3 -b 5 -v 5 -m 8 -o hop_vs_hop.blast</code>

<code>blastall -p blastp -i can.pep.fasta -d hopBlastpDB/hop.pep.fasta -e 1e-3 -b 5 -v 5 -m 8 -o can_vs_hop.blast</code>

<code>blastall -p blastp -i hop.pep.fasta -d canBlastpDB/can.pep.fasta -e 1e-3 -b 5 -v 5 -m 8 -o hop_vs_can.blast</code>


## Prepare files for MCScanX
### To create hopGenes.gff
<code>python scripts/createMCScanXInputFiles.py combinedGeneModels.txt hop.pep.fasta</code>
### combinedGeneModels.txt is a mapping file that contains information about the gene models
### format: scaffoldID\toriginalGeneID\tnewGeneID\tgeneStart\tgeneStop\tcdsStart\tcdsStop\n

### To create hop.chr
<code>less hopGenes.gff | awk '{print $1 "\thop"}' | sort | uniq > hop.chr</code>

### To create cannabisGenes.gff
<code>python scripts/createMCScanXCannabisInputFiles.py can.pep.fasta GCF_900626175.2_cs10_genomic.gff GCF_900626175.2_cs10_protein.gpff</code>

### To create cannabis.chr
<code>less cannabisGenes.gff | awk '{print $1 "\tcan"}' | sort | uniq > cannabis.chr</code>

### To create hop_vs_cannabis.gff
<code>cat hopGenes.gff cannabisGenes.gff > hop_vs_cannabis.gff</code>

### To create hop_vs_cannabis.chr
<code>cat hop.chr cannabis.chr > hop_vs_cannabis.chr</code>

### To create hop_vs_cannabis.blast
<code>cat hop_vs_hop.blast can_vs_can.blast hop_vs_can.blast > hop_vs_cannabis.blast</code>

### To run MCScanX 
<code>MCScanX hop_vs_cannabis</code>


## Generate sequence alignments for anchor gene pairs from syntenic blocks

### Extract CDS for hop vs hop anchor gene pairs 

<code>python scripts/createGenePairFastaHopParalogs.py hop_vs_cannabis.collinearity hop.cds.fasta</code>

### Extract CDS for hemp vs hemp anchor gene pairs

<code>python scripts/createGenePairFastaCannabisParalogs.py hop_vs_cannabis.collinearity can.cds.fasta</code>

### Extract CDS for hop vs hemp anchor gene pairs

<code>cat hop.cds.fasta can.cds.fasta > hop_vs_cannabis.cds.fasta</code>

<code>python scripts/createGenePairFastaHopCannabisParalogs.py hop_vs_cannabis.collinearity hop_vs_cannabis.cds.fasta</code>

### MACSE alignSequences - same command format for each set of gene pairs

<code>java -jar macse_v2.03.jar -prog alignSequences -seq gene1_vs_gene2.cds.fasta -out_NT gene1_vs_gene2_NT.fasta -out_AA gene1_vs_gene2_AA.fasta</code>

<code>ls -1 *_NT.fasta > alignedNTFileList.txt</code>

<code>python scripts/identifyFrameshifts.py alignedNTFileList.txt</code>


### MACSE exportAlignment - same command format for each set of gene pairs

<code>java -jar macse_v2.03.jar -prog exportAlignment -align gene1_vs_gene2_NT.fasta -codonForFinalStop --- -codonForInternalStop ,,, -codonForInternalFS +++ -charForRemainingFS + -out_NT gene1_vs_gene2_expAlign_NT.fasta -out_AA gene1_vs_gene2_expAlign_AA.fasta</code>

## Calculate 4DTv and Kimura distance
<code>python scripts/calculate4DTv.py gene1_vs_gene2_expAlign_NT.fasta</code>

### visualize 4DTv distribution

<code>ls -1 ../../cannabis_paralogs/kimura_and_4dtv/*4DTv.out > cannabis4DTvFileList.txt</code>

<code>ls -1 ../../hop_paralogs/kimura_and_4dtv/*4DTv.out > hop4DTvFileList.txt</code>

<code>ls -1 ../../hop_vs_cannabis_paralogs/kimura_and_4dtv/*4DTv.out > hop_vs_cannabis4DTvFileList.txt</code>

<code>python scripts/create_4DTv_count_hist.py cannabis4DTvFileList.txt hop4DTvFileList.txt hop_vs_cannabis4DTvFileList.txt</code>

## Calculate Ks (synonymous substitution rate)
KsPipeline.md
MixtureModel.md

