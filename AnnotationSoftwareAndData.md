# Software and Data for Annotating the Cascade Hop Dovetail Assembly

## Repeat Annotation
### *de novo* identification of long terminal retrotransposons (LTRs)
gt suffixerator (GenomeTools) 1.6.1  
gt ltrharvest (GenomeTools) 1.6.1  
LTR_FINDER_parallel v1.1  
LTR_retriever v2.7  

### Identification of non-LTR repeat sequences
RepeatMasker version 4.1.0  
Repeat library: mipsREdat_9.3p_Eudicot_TEs.fasta  
<code>RepeatMasker -lib mipsREdat_9.3p_Eudicot_TEs.fasta -qq -pa 4 -cutoff 225 -norna -a -gff -dir outputDir/ Scaffold.fasta</code>  

## Alignments for Gene Prediction
megablast 2.2.26  
blastx 2.10.0+  
exonerate version 2.3.0  

### EST data:  
*Humulus lupulus* ESTs from NCBI (25692 sequences, accessed 11/12/2018)  
*Humulus lupulus* ESTs from TrichOME (22959 sequences, accessed 03/28/2018)  

### Protein data:  
*Cannabis sativa* protein sequences from RefSeq (33639 sequences, accessed accessed 12/02/2020)  
*Prunus persica* protein sequences from Plaza (26843 sequences, accessed 10/13/2020)  
*Ziziphus jujuba* protein sequences from Plaza (28799 sequences, accessed 10/13/2020)  
UniProt Embryophyta protein sequences (38747 sequences, accessed 08/24/2020)  

## Gene Prediction  
BUSCO v4.1.1  
SNAP version 2006-07-28  
Augustus version 3.3.2  
MAKER 2.31.10  

### BUSCO is run in 'long' mode to train Augustus  
<pre>SGE_Batch -c "busco --in maskedGenomeAssembly.fasta --out outputDir --mode genome --config /path/busco_v4_config.ini --long" -r busco_v3_sge -q specified_queue -P 16</pre>  

### The masked assembly is directly provided to MAKER, and gene prediction proceeds in three rounds  
<pre>SGE_Batch -c "maker -RM_off" -r maker_round1_sge -q specified_queue</pre>  


### Web resources for guiding selection of MAKER parameters  
<details>
<summary>Resources</summary>

### MAKER forum on google
<https://groups.google.com/g/maker-devel?pli=1>

### Control files explained:  
<http://weatherby.genetics.utah.edu/MAKER/wiki/index.php/The_MAKER_control_files_explained>

### Training SNAP twice  
<https://biohpc.cornell.edu/doc/annotation_2019_exercises1.html>  
<https://robertslab.github.io/sams-notebook/2018/11/27/Annotation-Olurida_v081-MAKER-on-Mox.html>  
<https://robertslab.github.io/sams-notebook/2019/01/14/Annotation-Olurida_v081-MAKER-BUSCO-Augustus-Training.html>  
<https://gist.github.com/darencard/bb1001ac1532dd4225b030cf0cd61ce2>  

### "Generally there is little further improvement after 2 rounds of bootstrap training with the same evidence, and you run the risk of overtraining"      
[Genome Annotation and Curation Using MAKER and MAKER-P](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4286374)  

### Reference for training Augustus via BUSCO  
[Improving Illumina assemblies with Hi‐C and long reads: An example with the North African dromedary](https://onlinelibrary.wiley.com/doi/full/10.1111/1755-0998.13020)  

### Setting correct_est_fusion = 0 in the first round  
<https://github.com/wuying1984/MAKER2_PM_genome_annotation>  

### Setting correct_est_fusion = 1  
[The tea plant reference genome and improved gene annotation using long-read and paired-end sequencing data](https://www.nature.com/articles/s41597-019-0127-1)    
<https://github.com/ptranvan/Genome-annotation>  
<https://groups.google.com/g/maker-devel/c/J_ZLTFQ3xN4>  
<https://groups.google.com/g/maker-devel/c/tN-bxyC8IhQ>  

### Setting always_complete=1
<https://groups.google.com/g/maker-devel/c/tN-bxyC8IhQ>  
<https://github.com/wuying1984/MAKER2_PM_genome_annotation>  

</details>


## First round of MAKER
<details>
<summary>maker_opts.ctl</summary>

<pre>
#-----Genome (these are always required)
genome=maskedGenomeAssembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=allScaffolds_vs_ncbiESTs.modified.exonerate,allScaffolds_vs_trichomeESTs.modified.exonerate,hopCascadeDovetailMaskedStringtieTranscriptsHenning.gff3,hopCascadeDovetailMaskedStringtieTranscriptsMatthews.gff3 # aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein= #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=allScaffolds_vs_prunusPersica.modified.exonerate,allScaffolds_vs_uniprotEmbryophyta.modified.exonerate,allScaffolds_vs_ziziphusJujuba.modified.exonerate,allScaffolds_vs_refSeqCSativaProtein.modified.exonerate #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm= #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=1 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=1 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=0
# correct_est_fusion=1 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
</pre>
</details>

## First round of SNAP training
<details>
<summary>Commands</summary>

### create ZFF file  
<pre>SGE_Batch -c "maker2zff -d maskedGenomeAssembly.maker.output/maskedGenomeAssembly_master_datastore_index.log" -r maker2zff_sge -q specified_queue</pre>  

### fathom validate
<pre>SGE_Batch -c "fathom genome.ann genome.dna -validate > snap_validate_output.txt" -r fathomeValidate_sge -q specified_queue</pre>  

### grep for errors  
<pre>cat snap_validate_output.txt | grep "error" > fathomValidateErrors.txt</pre>  

### fathom categorize  
"break up the sequences into fragments with one gene per sequence"  
<https://vcru.wisc.edu/simonlab/bioinformatics/programs/snap/00README.txt>  
Why 1000? Seems to be standard practice; represents 1000 bp flanking gene  
<https://www.biostars.org/p/217144/>  
<https://reslp.github.io/blog/My-MAKER-Pipeline/>  
<pre>SGE_Batch -c "fathom genome.ann genome.dna -categorize 1000" -r fathomCategorize_sge -q specified_queue</pre>  

### fathom export
<pre>SGE_Batch -c "fathom uni.ann uni.dna -export 1000 -plus" -r fathomExport_sge -q specified_queue</pre>  

### forge  
<pre>SGE_Batch -c "forge export.ann export.dna" -r forge_sge -q specified_queue</pre>  

### create hmm file for MAKER  
<pre>/local/cluster/snap/hmm-assembler.pl maskedGenomeAssembly.fasta . > round1.hmm</pre>  
</details>


## Second round of MAKER
<details>
<summary>maker_opts.ctl</summary>

<pre>
#-----Genome (these are always required)
genome=maskedGenomeAssembly.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=allScaffolds_vs_ncbiESTs.modified.exonerate,allScaffolds_vs_trichomeESTs.modified.exonerate,hopCascadeDovetailMaskedStringtieTranscriptsHenning.gff3,hopCascadeDovetailMaskedStringtieTranscriptsMatthews.gff3 # aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein= #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=allScaffolds_vs_prunusPersica.modified.exonerate,allScaffolds_vs_uniprotEmbryophyta.modified.exonerate,allScaffolds_vs_ziziphusJujuba.modified.exonerate,allScaffolds_vs_refSeqCSativaProtein.modified.exonerate #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=round1.hmm  #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species= #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
# correct_est_fusion=0
correct_est_fusion=1 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary files
</pre>
</details>  


## Second round of SNAP training  
<details>  
<summary>Commands</summary>  

### create ZFF file
<pre>SGE_Batch -c "maker2zff -d maskedGenomeAssembly.maker.output/maskedGenomeAssembly_master_datastore_index.log" -r maker2zff_sge -q specified_queue</pre>  

### fathom validate
<pre>SGE_Batch -c "fathom genome.ann genome.dna -validate > snap_validate_output.txt" -r fathomeValidate_sge -q specified_queue</pre>  

### grep for errors
<pre>cat snap_validate_output.txt | grep "error" > fathomValidateErrors.txt</pre>  

### fathom categorize
<pre>SGE_Batch -c "fathom genome.ann genome.dna -categorize 1000" -r fathomCategorize_sge -q specified_queue</pre>  

### fathom export  
<pre>SGE_Batch -c "fathom uni.ann uni.dna -export 1000 -plus" -r fathomExport_sge -q specified_queue</pre>  

### forge  
<pre>SGE_Batch -c "forge export.ann export.dna" -r forge_sge -q specified_queue</pre>  

### create hmm file for MAKER  
<pre>/local/cluster/snap/hmm-assembler.pl maskedGenomeAssembly.fasta . > round2.hmm</pre>  
</details/>


## Third round of MAKER
<details>
<summary>maker_opts.ctl</summary>

<pre>
#-----Genome (these are always required)
genome=maskedScaffold.fasta #genome sequence (fasta file or fasta embeded in GFF3 file)
organism_type=eukaryotic #eukaryotic or prokaryotic. Default is eukaryotic

#-----Re-annotation Using MAKER Derived GFF3
maker_gff= #MAKER derived GFF3 file
est_pass=0 #use ESTs in maker_gff: 1 = yes, 0 = no
altest_pass=0 #use alternate organism ESTs in maker_gff: 1 = yes, 0 = no
protein_pass=0 #use protein alignments in maker_gff: 1 = yes, 0 = no
rm_pass=0 #use repeats in maker_gff: 1 = yes, 0 = no
model_pass=0 #use gene models in maker_gff: 1 = yes, 0 = no
pred_pass=0 #use ab-initio predictions in maker_gff: 1 = yes, 0 = no
other_pass=0 #passthrough anyything else in maker_gff: 1 = yes, 0 = no

#-----EST Evidence (for best results provide a file for at least one)
est= #set of ESTs or assembled mRNA-seq in fasta format
altest= #EST/cDNA sequence file in fasta format from an alternate organism
est_gff=allScaffolds_vs_ncbiESTs.modified.exonerate,allScaffolds_vs_trichomeESTs.modified.exonerate,hopCascadeDovetailMaskedStringtieTranscriptsHenning.gff3,hopCascadeDovetailMaskedStringtieTranscriptsMatthews.gff3 # aligned ESTs or mRNA-seq from an external GFF3 file
altest_gff= #aligned ESTs from a closly relate species in GFF3 format

#-----Protein Homology Evidence (for best results provide a file for at least one)
protein= #protein sequence file in fasta format (i.e. from mutiple oransisms)
protein_gff=allScaffolds_vs_prunusPersica.modified.exonerate,allScaffolds_vs_uniprotEmbryophyta.modified.exonerate,allScaffolds_vs_ziziphusJujuba.modified.exonerate,allScaffolds_vs_refSeqCSativaProtein.modified.exonerate #aligned protein homology evidence from an external GFF3 file

#-----Repeat Masking (leave values blank to skip repeat masking)
model_org= #select a model organism for RepBase masking in RepeatMasker
rmlib= #provide an organism specific repeat library in fasta format for RepeatMasker
repeat_protein= #provide a fasta file of transposable element proteins for RepeatRunner
rm_gff= #pre-identified repeat elements from an external GFF3 file
prok_rm=0 #forces MAKER to repeatmask prokaryotes (no reason to change this), 1 = yes, 0 = no
softmask=1 #use soft-masking rather than hard-masking in BLAST (i.e. seg and dust filtering)

#-----Gene Prediction
snaphmm=round2.hmm  #SNAP HMM file
gmhmm= #GeneMark HMM file
augustus_species=maskedHopCascadeDovetail_BUSCO #Augustus gene prediction species model
fgenesh_par_file= #FGENESH parameter file
pred_gff= #ab-initio predictions from an external GFF3 file
model_gff= #annotated gene models from an external GFF3 file (annotation pass-through)
est2genome=0 #infer gene predictions directly from ESTs, 1 = yes, 0 = no
protein2genome=0 #infer predictions from protein homology, 1 = yes, 0 = no
trna=0 #find tRNAs with tRNAscan, 1 = yes, 0 = no
snoscan_rrna= #rRNA file to have Snoscan find snoRNAs
unmask=0 #also run ab-initio prediction programs on unmasked sequence, 1 = yes, 0 = no

#-----Other Annotation Feature Types (features MAKER doesn't recognize)
other_gff= #extra features to pass-through to final MAKER generated GFF3 file

#-----External Application Behavior Options
alt_peptide=C #amino acid used to replace non-standard amino acids in BLAST databases
cpus=1 #max number of cpus to use in BLAST and RepeatMasker (not for MPI, leave 1 when using MPI)

#-----MAKER Behavior Options
max_dna_len=100000 #length for dividing up contigs into chunks (increases/decreases memory usage)
min_contig=10000 #skip genome contigs below this length (under 10kb are often useless)

pred_flank=200 #flank for extending evidence clusters sent to gene predictors
pred_stats=0 #report AED and QI statistics for all predictions as well as models
AED_threshold=1 #Maximum Annotation Edit Distance allowed (bound by 0 and 1)
min_protein=0
alt_splice=0 #Take extra steps to try and find alternative splicing, 1 = yes, 0 = no
always_complete=1 #extra steps to force start and stop codons, 1 = yes, 0 = no
map_forward=0 #map names and attributes forward from old GFF3 genes, 1 = yes, 0 = no
keep_preds=0 #Concordance threshold to add unsupported gene prediction (bound by 0 and 1)

split_hit=10000 #length for the splitting of hits (expected max intron size for evidence alignments)
single_exon=0 #consider single exon EST evidence when generating annotations, 1 = yes, 0 = no
single_length=250 #min length required for single exon ESTs if 'single_exon is enabled'
correct_est_fusion=1 #limits use of ESTs in annotation to avoid fusion genes

tries=2 #number of times to try a contig if there is a failure for some reason
clean_try=0 #remove all data from previous run before retrying, 1 = yes, 0 = no
clean_up=0 #removes theVoid directory with individual analysis files, 1 = yes, 0 = no
TMP= #specify a directory other than the system default temporary directory for temporary file
</pre>
</details>

