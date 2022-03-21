## Paired-end quality trimming for DNA short-read sequencing from Cascade
cutadapt --minimum-length=20 -q 30 --output=cascade_R1_qualFilter30.fastq --paired-output=cascade_R2_qualFilter30.fastq --pair-filter=any cascade_R1.fastq cascade_R2.fastq

## Run polca
<pre>MaSuRCA-3.4.2/bin/polca.sh -a unpolishedAssembly.fasta -r 'cascade_R1_qualFilter30.fastq cascade_R2_qualFilter30.fastq'</pre>
