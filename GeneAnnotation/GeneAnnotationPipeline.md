# Pipeline for assigning putative gene functions in the Cascade Hop Dovetail Assembly 

# Identify similarity to Pfam protein domains
<details>
<summary>Pfam command</summary>

Pfam release 33.1 (accessed 08/25/2020)

HMMER 3.3

<code>hmmscan --cpu 64 --domtblout combinedHopCascadeDovetail.domtblout Pfam-A.hmm geneModels.pep.fasta > combinedHopCascadeDovetail.err</code>

<code>python assignPfam.py combinedHopCascadeDovetail.domtblout pfamRepeatDomains.txt geneModels.pep.fasta</code>

</details>

# Identify similarity to UniProt genes
<details>
<summary>UniProt commands</summary>

# hop vs uniprot transposable element genes

<code>blastp -query geneModels.pep.fasta -db uniprot_transposable_element_KW0814.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_uniprotTEs.blastp -num_threads 16</code>

<code>blastp -query uniprot_transposable_element_KW0814.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out uniprotTEs_vs_hop.blastp -num_threads 16</code>

# hop vs uniprot bacteria genes
<code>blastp -query geneModels.pep.fasta -db uniprot_bacteria.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_bacteria.blastp -num_threads 16</code>

<code>blastp -query uniprot_bacteria.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out bacteria_vs_hop.blastp -num_threads 16</code>

# hop vs uniprot virus genes
<code>blastp -query geneModels.pep.fasta -db uniprot_viruses.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_viruses.blastp -num_threads 16</code>

<code>blastp -query uniprot_viruses.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out viruses_vs_hop.blastp -num_threads 16</code>

# hop vs uniprot plant genes
<code>blastp -query geneModels.pep.fasta -db uniprot-reviewed_yes+taxonomy_3193.08242020.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_uniprotPlants.blastp -num_threads 16</code>

<code>blastp -query uniprot-reviewed_yes+taxonomy_3193.08242020.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out uniprotPlants_vs_hop.blastp -num_threads 16</code>
</details>
