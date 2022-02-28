# Pipeline for assigning putative gene functions in the Cascade Hop Dovetail Assembly 

# Identify similarity to Pfam protein domains
<details>
<summary>Pfam information</summary>
Pfam release 33.1 (accessed 08/25/2020)\n
HMMER 3.3\n
<code>hmmscan --cpu 64 --domtblout combinedHopCascadeDovetail.domtblout Pfam-A.hmm geneModels.pep.fasta > combinedHopCascadeDovetail.err</code>\n
</details>

# UniProt
