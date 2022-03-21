# Pipeline for developing Transdecoder gene models

<pre>Transdecoder-v5.5.0</pre>
<pre>hisat2</pre>
<pre>StringTie v1.3.3b</pre>

<details><summary>Assemble RNA-seq transcripts</summary>
### Example hisat2 command
<pre>hisat2 --rna-strandness FR --no-discordant --no-mixed --dta -x HISAT2_INDEXES/polishedMaskedAssembly.fastsa -1 cascadeRNASeq_R1.fastq -2 cascadeRNASeq_R2.fastq -S alignedRNASeq.sam</pre>

### Convert sam to bam file
<pre>samtools view alignedRNASeq.sam -bS -o alignedRNASeq.bam</pre>

### Sort bam file
<pre>samtools sort -m 4G alignedRNASeq.bam -o sorted_alignedRNASeq.bam</pre>

### Index sorted bam file
<pre>samtools index sorted_alignedRNASeq.bam sorted_alignedRNASeq.bam.bai</pre>

### Assemble transcripts with stringtie
<pre>stringtie sorted_alignedRNASeq.bam -j 2 -o assembledRNASeqTranscripts.gtf --fr -A assembledRNASeqTranscripts.tab</pre>

### Create stringtie_merged.gtf
<pre>stringtie --merge hopConeTranscriptAssemblies.txt -o stringtie_merged.gtf</pre>
<p>hopConeTranscriptAssemblies.txt contains list of transcript assemblies from leaf, meristem, stem, gland, and cone tissue from three developmental time points</p>
</details>

<details><summary>Create gene models with Transdecoder</summary>
### Create transcripts.gff3
<pre>TransDecoder-v5.5.0/util/gtf_to_alignment_gff3.pl stringtie_merged.gtf > transcripts.gff3</pre>

### Create transcripts.fasta
<pre>TransDecoder-v5.5.0/util/gtf_genome_to_cdna_fasta.pl stringtie_merged.gtf combinedRepeatsDovetail.denovoLTRsOnly.singleScaffold.fasta > transcripts.fasta</pre>

### Get longest open reading frame with TransDecoder.LongOrfs
<pre>TransDecoder.LongOrfs -S -t transcripts.fasta -O longestORF</pre>

### Run hmmscan to get similarity to Pfam domains
<pre>hmmscan --cpu 16 --domtblout hopTranscripts.domtblout Pfam-A.hmm ../longestORF/longest_orfs.pep</pre>

### Run blastp to get similarity to UniProt plant genes
<pre>blastp -query ../longestORF/longest_orfs.pep -db uniprotPlants.fasta -outfmt 6 -evalue 1e-3 -num_threads 10 > hop_vs_uniprot.blastp.txt</pre>

### Run Transdecoder.Predict
<pre>TransDecoder.Predict -t transcripts.fasta --retain_pfam_hits hmmResults/hopTranscripts.domtblout --retain_blastp_hits blastpResults/hop_vs_uniprot.blastp.txt --output_dir longestORF</pre>

### Clean protein and CDS fasta file deflines
<pre>python scripts/cleanDefline.py transcripts.fasta.transdecoder.pep transcripts.fasta.transdecoder.cds</pre>

### Create transcripts.fasta.transdecoder.genomeCentric.gff3
<pre>TransDecoder-v5.5.0/util/cdna_alignment_orf_to_genome_orf.pl transcripts.fasta.transdecoder.gff3 transcripts.gff3 transcripts.fasta > transcripts.fasta.transdecoder.genomeCentric.gff3</pre>
</details>



