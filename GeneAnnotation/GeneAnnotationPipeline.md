# Pipeline for assigning putative gene functions in the Cascade Hop Dovetail Assembly 

## Run hmmscan to identify similarity to Pfam protein domains
<details>
<summary>Pfam command</summary>

Pfam release 33.1 (accessed 08/25/2020)

HMMER 3.3

<code>hmmscan --cpu 64 --domtblout combinedHopCascadeDovetail.domtblout Pfam-A.hmm geneModels.pep.fasta > combinedHopCascadeDovetail.err</code>

</details>

## Assign similarity to Pfam protein domains
<details>
<summary>Command</summary>
<code>python assignPfam.py combinedHopCascadeDovetail.domtblout pfamRepeatDomains.txt geneModels.pep.fasta</code>
</details>


## Identify similarity to UniProt genes with blastp
<details>
<summary>UniProt commands</summary>

### hop vs uniprot transposable element genes

<code>blastp -query geneModels.pep.fasta -db uniprotTEs.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_uniprotTEs.blastp -num_threads 16</code>

<code>blastp -query uniprotTEs.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out uniprotTEs_vs_hop.blastp -num_threads 16</code>

### hop vs uniprot bacteria genes
<code>blastp -query geneModels.pep.fasta -db uniprotBacteria.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_bacteria.blastp -num_threads 16</code>

<code>blastp -query uniprotBacteria.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out bacteria_vs_hop.blastp -num_threads 16</code>

### hop vs uniprot virus genes
<code>blastp -query geneModels.pep.fasta -db uniprotViruses.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_viruses.blastp -num_threads 16</code>

<code>blastp -query uniprotViruses.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out viruses_vs_hop.blastp -num_threads 16</code>

### hop vs uniprot plant genes
<code>blastp -query geneModels.pep.fasta -db uniprotPlants.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out hop_vs_uniprotPlants.blastp -num_threads 16</code>

<code>blastp -query uniprotPlants.fasta -db blastDB/geneModels.pep.fasta -evalue 1e-3 -outfmt '6 std qcovs' -out uniprotPlants_vs_hop.blastp -num_threads 16</code>
</details>

## Assign similarity to UniProt plant genes
<details>
<summary>Command</summary>
<code>python collectTopUniProtHits.py hop_vs_uniprotPlants.blastp uniprotPlants_vs_hop.blastp uniprotPlants.fasta uniprot_transposable_element_KW0814.fasta TEUniProtGenes.txt 20</code>
</details>

## Assign similarity to other UniProt genes
<details>
<summary>Command</summary>
<code>python getOtherTopUniprotHit.py hop_vs_bacteria.blastp bacteria_vs_hop.blastp uniprotBacteria.fasta Bacteria 30 hop</code>

<code>python getOtherTopUniprotHit.py hop_vs_uniprotTEs.blastp uniprotTEs_vs_hop.blastp uniprotTEs.fasta UniprotTE 30 hop</code>

<code>python getOtherTopUniprotHit.py hop_vs_viruses.blastp viruses_vs_hop.blastp uniprotViruses.fasta Virus 30 hop</code>  

</details>
