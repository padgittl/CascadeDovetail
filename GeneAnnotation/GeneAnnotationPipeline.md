# Pipeline for assigning putative gene functions in the Cascade Hop Dovetail Assembly 

# Identify similarity to Pfam protein domains
<details>
<summary>Pfam command</summary>

Pfam release 33.1 (accessed 08/25/2020)

HMMER 3.3

<code>hmmscan --cpu 64 --domtblout combinedHopCascadeDovetail.domtblout Pfam-A.hmm geneModels.pep.fasta > combinedHopCascadeDovetail.err</code>

</details>

# UniProt
