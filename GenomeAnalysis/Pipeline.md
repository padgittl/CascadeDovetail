## Create data files for circos
<pre>python scripts/getDataForCircos.py scaffoldLengths.txt hop.pep.fasta combinedGeneModels.txt denovoLTRsDovetail.gff geneticMap.tsv goTerms_biologicalProcesses.tsv 5000000 &</pre>

<details>
<summary>File descriptions</summary>

<details>
<summary>scaffoldLengths.txt</summary>
<pre>Scaffold_1531 476495644
Scaffold_19 434152558
Scaffold_1533 423633482
Scaffold_76 370465130
Scaffold_24 345299309</pre>
</details>

<details>
<summary>hop.pep.fasta</summary>
<pre>>HUMLU_CAS0090952.t1.p1
MDLTSPRYFHAPSSTFSSATGEPPLEANASFYGKTKNNPFAETFPDPLCKLNLKETSEFV
KSFPMPHGGTESNRVFRESSTQRRTEVGVNSVVTQRRFEAPPTPGRPVFSFSAGNLSRKG
FPSKWDDAEKWLISSSCHESPAHTIKPSESVRIAKPSDYNFKQQMEVFADKSRVTEEKVS
KKQFSSFHCSVSLDNHNSVRAFDGVSCSTDNVFLKDKFTNEIEPVLPNFRTSESTKEGFL
FKNSACEAMKDAGTEMVQHRDVGTEMTPLGSSTTSRCHTPFKISSPARHNTPANRSGPLG
LEHSNSISSTIDIAQLQECHLAKLQLGTHYDSVTSNWSSRQEEEEEISKSLRHFEIDNVN
CCQKNGPESRAVAWEEEEKTKCCLRYQREEAKIQAWVNLQSAKAEAQSRKLEVKIQKMRS
NLEEKLMKRMAVVHRKAEEWREAARQQHSDQIEKATVHAQKMVIRNNSHFSTATSCGCFP
CNNHFR*</pre>
</details>

<details>
<summary>combinedGeneModels.txt</summary>
<pre>scaffoldID	originalGeneID	newGeneID	geneStart	geneStop	cdsStart	cdsStop	
Scaffold_1531   MAKER0000005.t1 HUMLU_CAS0000005.t1.p1  43615   44430   43615   44430
Scaffold_1531   MSTRG.1239.1.p1 HUMLU_CAS0000006.t1.p1  50767   52290   51380   52057</pre>
</details>

<details>
<summary>denovoLTRsDovetail.gff</summary>
<pre>Scaffold_1531   RepeatMasker    LTR/Copia       34571   35520   8.6     +       6463    Scaffold_1531:202456750..202457742_LTR
Scaffold_1531   RepeatMasker    LTR/Copia       37287   40006   26.8    +       7858    Scaffold_1531:93765742..93770940_INT
Scaffold_1531   RepeatMasker    LTR/Copia       40089   41718   26.5    +       5781    Scaffold_1531:93765742..93770940_INT</pre>
</details>

<details>
<summary>geneticMap.tsv</summary>
<pre>Trait   Marker  Chr     Pos     df      GeneticMap      F       p       add_effect      add_F   add_p
   dom_effect      dom_F   dom_p   errordf MarkerR2        Genetic Var     Residual Var    -2LnLikelihood
sex     S1_468115197    1       100000  2       0       1.1379600000000001      0.32196999999999998
     0.1198  1.8134300000000001      0.1792  0.12751999999999999     2.2037499999999999      0.13882 277     8.1799999999999998E-3   0.39034999999999997     0.000004        239.86838</pre>
</details>

<details>
<summary>goTerms_biologicalProcesses.tsv</summary>
<pre>scaffoldID	geneID	geneStart	geneStop	uniprotID	goTerm	uniprotGeneName	goDescription	
Scaffold_73     HUMLU_CAS0053498.t1.p1  15108726        15113171        Q7XA40  GO:0006952      sp|Q7XA40|RGA3_SOLBU Putative disease resistance protein RGA3 OS=Solanum bulbocastanum OX=147425 GN=RGA3 PE=2 SV=2      defense response</pre>
</details>

</details>

## Create synteny file for circos
<pre>python scripts/prepareCircosSyntenyFile.py scaffoldLengths.txt combinedGeneModels.txt hop_vs_cannabis.match_size9.collinearity hopGenes.gff synteny.match_size9.txt</pre>

## Run circos
<pre>circos -conf circos.conf</pre>