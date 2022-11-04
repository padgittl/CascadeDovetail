# Pipeline for annotating repeats in the Cascade Hop Dovetail Assembly 


### *de novo* identification of long terminal retrotransposons (LTRs)

<details>
<summary>commands</summary>

gt suffixerator (GenomeTools) 1.6.1  
<code>gt suffixerator -db scaffoldID.fasta -indexname scaffoldID -tis -suf -lcp -des -ssp -dna</code>

gt ltrharvest (GenomeTools) 1.6.1  
<code>gt ltrharvest -index ../index/scaffoldID -out scaffoldID.ltrharvest.out -outinner scaffoldID.ltrharvest.outinner -gff3 scaffoldID.ltrharvest.gff3 -seqids yes > scaffoldID.ltrharvestScreen.out</code>

LTR_FINDER_parallel v1.1  
<code>LTR_FINDER_parallel -seq ../index/scaffoldID.fasta -threads 4 -harvest_out</code>

LTR_retriever v2.7 (wherever the index files are located will be the location where the LTR_retriever output files are generated!!)    
<code>LTR_retriever -genome ../index/scaffoldID.fasta -inharvest ../ltrharvest/scaffoldID.ltrharvestScreen.out -infinder ../LTR_finder_parallel/scaffoldID.fasta.finder.combine.scn -threads 4 > scaffoldID.ltr_retriever.out</code>

Combine LTR_retriever output files (if pipeline is performed on each scaffold separately)  
<code>cat *.out.gff > denovoLTRsDovetail.gff</code>
</details>



### Identification of non-LTR repeat sequences  

<details>
<summary>commands</summary>
RepeatMasker version 4.1.0  
Repeat library: [mipsREdat_9.3p_Eudicot_TEs.fasta](https://www.mmnt.net/db/0/0/ftp.mips.embnet.org/plants/REdat)  
<code>RepeatMasker -lib mipsREdat_9.3p_Eudicot_TEs.fasta -qq -pa 4 -cutoff 225 -norna -a -gff -dir outputDir/ scaffoldID.fasta</code>  

Combine RepeatMasker output files (if pipeline is performed on each scaffold separately)  
<code>cat outputDir/*.gff > mipsRERepMaskDovetail.gff</code>

### Combine RepeatMasker and *de novo* results

Create GFF file (after this command, we refer to the combined GFF file as combinedRepeats.gff, which is here denoted as outputFileName.gff, to signify that it is the output file of this script)  
<code>python scripts/combineOutputGFF.py mipsRERepMaskDovetail.gff denovoLTRsDovetail.gff outputFileName.gff</code>

Create masked fasta file (after this command, we refer to the combined fasta file as combinedRepeats.fasta, which is here denoted as outputFileName.fasta, to signify that it is the output file of this script)  
<code>bedtools maskfasta -fi assembly.fasta -bed combinedRepeats.gff -fo outputFileName.fasta</code>
</details>


### Analyze repeat results (this step was performed on individual scaffolds)

<details>
<summary>commands</summary>

Split GFF  
<code>python scripts/splitGFF.py combinedRepeats.gff</code>

Split fasta  
<code>python scripts/splitFasta.py combinedRepeats.fasta</code>

Calculate repeat percentages per scaffold  
<code>for f in Scaffold_*.gff; do echo python scripts/getRepeatPercentages.py $f polishedScaffoldLengths.txt '>' \`basename $f gff\`txt; done > getPerScaffoldRepeatContent.sh</code>

<details>
<summary>polishedScaffoldLengths.txt format</summary>
<pre>Scaffold_1531 476495644
Scaffold_19 434152558
Scaffold_1533 423633482
Scaffold_76 370465130
Scaffold_24 345299309
Scaffold_172 327882944
Scaffold_77 316519611
Scaffold_73 303741476
Scaffold_49 290858211
Scaffold_191 185200997</pre>
</details>

Create file list  
<code>ls -1 *txt > repeatCountFileList.txt</code>

Calculate repeat percentages for whole assembly relative to total repeat content  
<code>python scripts/getRepeatPercentageFromSingleGFFs4PieChart.py repeatCountFileList.txt > repeatPercentagesRelative2TotalRepeatContent.txt</code>

Calculate repeat percentages for whole assembly relative to assembly size (assembly size here is 3713677344 bp)    
<code>python scripts/getRepeatPercentageFromSingleGFFs4StackedBarChart.py repeatCountFileList.txt 3713677344 > repeatPercentagesRelative2AssemblySize.txt</code>
</details>


### Pie chart visualizations

<details>
<summary>commands</summary>

Repeat percentages relative to total repeat content  
<code>python scripts/createPieChart_relativeToTotalRepeatContent.py repeatPercentagesRelative2TotalRepeatContent.txt</code>

Repeat percentages relatives to assembly size  
<code>python scripts/createPieChart_relativetoAssemblySize.py repeatPercentagesRelative2AssemblySize.txt</code>

</details>