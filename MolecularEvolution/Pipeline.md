# Pipeline for assessing molecular evolution in genes from syntenic blocks

<details>
<summary>blastall for MCScanX</summary>

<pre>formatdb -i hop.pep.fasta -p T -o T</pre>

<pre>formatdb -i can.pep.fasta -p T -o T</pre>

<pre>blastall -p blastp -i hop.pep.fasta -d hopBlastpDB/hop.pep.fasta -e 1e-3 -b 5 -v 5 -m 8 -o hop_vs_hop.blast</pre>

<pre>blastall -p blastp -i can.pep.fasta -d hopBlastpDB/hop.pep.fasta -e 1e-3 -b 5 -v 5 -m 8 -o can_vs_hop.blast</pre>

<pre>blastall -p blastp -i hop.pep.fasta -d canBlastpDB/can.pep.fasta -e 1e-3 -b 5 -v 5 -m 8 -o hop_vs_can.blast</pre>
</details>


<details>
<summary>Prepare files for MCScanX</summary>

### To create hopGenes.gff
<pre>python scripts/createMCScanXInputFiles.py combinedGeneModels.txt hop.pep.fasta</pre>
### combinedGeneModels.txt is a mapping file that contains information about the gene models
### format: scaffoldID\toriginalGeneID\tnewGeneID\tgeneStart\tgeneStop\tcdsStart\tcdsStop\n

### To create hop.chr
<pre>less hopGenes.gff | awk '{print $1 "\thop"}' | sort | uniq > hop.chr</pre>

### To create cannabisGenes.gff
<pre>python scripts/createMCScanXCannabisInputFiles.py can.pep.fasta GCF_900626175.2_cs10_genomic.gff GCF_900626175.2_cs10_protein.gpff</pre>

### To create cannabis.chr
<pre>less cannabisGenes.gff | awk '{print $1 "\tcan"}' | sort | uniq > cannabis.chr</pre>

### To create hop_vs_cannabis.gff
<pre>cat hopGenes.gff cannabisGenes.gff > hop_vs_cannabis.gff</pre>

### To create hop_vs_cannabis.chr
<pre>cat hop.chr cannabis.chr > hop_vs_cannabis.chr</pre>

### To create hop_vs_cannabis.blast
<pre>cat hop_vs_hop.blast can_vs_can.blast hop_vs_can.blast > hop_vs_cannabis.blast</pre>
</details>


<details>
<summary>Run MCScanX</summary>

<pre>MCScanX hop_vs_cannabis</pre>
</details>


<details>
<summary>Generate sequence alignments for anchor gene pairs from syntenic blocks</summary>

### Extract CDS for hop vs hop anchor gene pairs 

<pre>python scripts/createGenePairFastaHopParalogs.py hop_vs_cannabis.collinearity hop.cds.fasta</pre>

### Extract CDS for hemp vs hemp anchor gene pairs

<pre>python scripts/createGenePairFastaCannabisParalogs.py hop_vs_cannabis.collinearity can.cds.fasta</pre>

### Extract CDS for hop vs hemp anchor gene pairs

<pre>cat hop.cds.fasta can.cds.fasta > hop_vs_cannabis.cds.fasta</pre>

<pre>python scripts/createGenePairFastaHopCannabisParalogs.py hop_vs_cannabis.collinearity hop_vs_cannabis.cds.fasta</pre>

### MACSE alignSequences - same command format for each set of gene pairs

<pre>java -jar macse_v2.03.jar -prog alignSequences -seq gene1_vs_gene2.cds.fasta -out_NT gene1_vs_gene2_NT.fasta -out_AA gene1_vs_gene2_AA.fasta</pre>

<pre>ls -1 *_NT.fasta > alignedNTFileList.txt</pre>

<pre>python scripts/identifyFrameshifts.py alignedNTFileList.txt</pre>

### MACSE exportAlignment - same command format for each set of gene pairs

<pre>java -jar macse_v2.03.jar -prog exportAlignment -align gene1_vs_gene2_NT.fasta -codonForFinalStop --- -codonForInternalStop ,,, -codonForInternalFS +++ -charForRemainingFS + -out_NT gene1_vs_gene2_expAlign_NT.fasta -out_AA gene1_vs_gene2_expAlign_AA.fasta</pre>
</details>


<details>
<summary>Calculate 4DTv and Kimura distance</summary>

<pre>python scripts/calculate4DTv.py gene1_vs_gene2_expAlign_NT.fasta</pre>

### visualize 4DTv distribution

<pre>ls -1 ../../cannabis_paralogs/kimura_and_4dtv/*4DTv.out > cannabis4DTvFileList.txt</pre>

<pre>ls -1 ../../hop_paralogs/kimura_and_4dtv/*4DTv.out > hop4DTvFileList.txt</pre>

<pre>ls -1 ../../hop_vs_cannabis_paralogs/kimura_and_4dtv/*4DTv.out > hop_vs_cannabis4DTvFileList.txt</pre>

<pre>python scripts/create_4DTv_count_hist.py cannabis4DTvFileList.txt hop4DTvFileList.txt hop_vs_cannabis4DTvFileList.txt</pre>
</details>


<details>
<summary>Calculate Ks (synonymous substitution rate)</summary>
<pre>KsPipeline.md</pre>
<pre>MixtureModel.md</pre>
</details>
