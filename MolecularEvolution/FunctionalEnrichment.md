## Perform hypergeometric test

<details>
<summary>Hop vs hop anchor gene pairs</summary>
<pre>python scripts/hypergeometric_hop_vs_hop.py hop_fullNonRepeatAssociatedGeneList.tsv hop_vs_cannabis.collinearity uniprot-reviewed_yes+taxonomy_3193.allGOTerms.08242020.tab uniprot-reviewed_yes+taxonomy_3193.biologicalProcesses.08242020.tab biologicalProcesses</pre>

<pre>python scripts/hypergeometric_hop_vs_hop.py hop_fullNonRepeatAssociatedGeneList.tsv hop_vs_cannabis.collinearity uniprot-reviewed_yes+taxonomy_3193.allGOTerms.08242020.tab uniprot-reviewed_yes+taxonomy_3193.cellularComponents.08242020.tab cellularComponent</pre>

<pre>python scripts/hypergeometric_hop_vs_hop.py hop_fullNonRepeatAssociatedGeneList.tsv hop_vs_cannabis.collinearity uniprot-reviewed_yes+taxonomy_3193.allGOTerms.08242020.tab uniprot-reviewed_yes+taxonomy_3193.molecularFunction.08242020.tab molecularFunction</pre>

<p>hop_fullNonRepeatAssociatedGeneList.tsv contains information about similarity to UniProt genes and Pfam protein domains</p>
<pre>HUMLU_CAS0073067.t1.p1  uniprotPlantNonRepeat   Q9FX89  sp|Q9FX89|FB50_ARATH Putative F-box protein At1g49610 OS=Arabidopsis thaliana OX=3702 GN=At1g49610 PE=4 SV=2
HUMLU_CAS0073067.t1.p1  pfamNonRepeat   PF00646.34      F-box domain</pre>

## Volcano plot
<pre>python scripts/volcano.py hypergeometric_biologicalProcesses_allFDRValues.txt bp</pre>
<pre>python scripts/volcano.py hypergeometric_cellularComponent_allFDRValues.txt cc</pre>
<pre>python scripts/volcano.py hypergeometric_molecularFunction_allFDRValues.txt mf</pre>
</details>

<details>
<summary>Hemp vs hemp anchor gene pairs</summary>
<pre>python hypergeometric_hop_vs_hemp.py hop_fullNonRepeatAssociatedGeneList.tsv hop_vs_cannabis.collinearity uniprot-reviewed_yes+taxonomy_3193.allGOTerms.08242020.tab uniprot-reviewed_yes+taxonomy_3193.biologicalProcesses.08242020.tab biologicalProcesses</pre>
<pre>python hypergeometric_hop_vs_hemp.py hop_fullNonRepeatAssociatedGeneList.tsv hop_vs_cannabis.collinearity uniprot-reviewed_yes+taxonomy_3193.allGOTerms.08242020.tab uniprot-reviewed_yes+taxonomy_3193.cellularComponents.08242020.tab cellularComponent</pre>
<pre>python hypergeometric_hop_vs_hemp.py hop_fullNonRepeatAssociatedGeneList.tsv hop_vs_cannabis.collinearity uniprot-reviewed_yes+taxonomy_3193.allGOTerms.08242020.tab uniprot-reviewed_yes+taxonomy_3193.molecularFunction.08242020.tab molecularFunction</pre>

## Volcano plot
<pre>python scripts/volcano.py hypergeometric_biologicalProcesses_allFDRValues.txt bp</pre>
<pre>python scripts/volcano.py hypergeometric_cellularComponent_allFDRValues.txt cc</pre>
<pre>python scripts/volcano.py hypergeometric_molecularFunction_allFDRValues.txt mf</pre>
</details>