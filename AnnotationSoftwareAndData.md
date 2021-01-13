# Software and Data for Annotating the Cascade Hop Dovetail Assembly

## Repeat Annotation
### *de novo* identification of long terminal retrotransposons (LTRs)
<p>gt suffixerator (GenomeTools) 1.6.1</p>
gt ltrharvest (GenomeTools) 1.6.1 
LTR_FINDER_parallel v1.1
LTR_retriever v2.7

### identification of non-LTR repeat sequences
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
SGE_Batch -c "busco --in maskedGenomeAssembly.fasta --out outputDir --mode genome --config /path/busco_v4_config.ini --long" -r busco_v3_sge -q specified_queue -P 16

### The masked assembly is directly provided to MAKER
SGE_Batch -c "maker -RM_off" -r maker_round1_sge -q specified_queue


