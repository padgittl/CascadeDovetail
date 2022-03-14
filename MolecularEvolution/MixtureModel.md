## Generate files for mixture modeling

<pre>python scripts/create_dnds_hist.py cannabis_yn00OutFileList.txt hop_yn00OutFileList.txt hop_vs_cannabis_yn00OutFileList.txt tenLargestScaffoldLengths.txt combinedGeneModels.txt 0.0 5 0.03 0.0 5.0</pre>

### cannabis_yn00OutFileList.txt, hop_yn00OutFileList.txt, and hop_vs_cannabis_yn00OutFileList.txt are generated with this command: 
<pre>ls -1 /path/to/dnds/*_vs_*/*_vs_*yn00.txt > yn00OutFileList.txt</pre>

<pre>### tenLargestScaffoldLengths.txt contains this information:
Scaffold_1531 476495644
Scaffold_19 434152558
Scaffold_1533 423633482
Scaffold_76 370465130
Scaffold_24 345299309
Scaffold_172 327882944
Scaffold_77 316519611
Scaffold_73 303741476
Scaffold_49 290858211
Scaffold_191 185200997</pre>


## Perform mixture modeling with mclust in R
<pre>scripts/mclust_mixture_model.R</pre>

## Create Ks Figure
<pre>python scripts/KsFigure.py can_log_ks_mclust_data.txt can_log_ks_mclust_densities.txt can_log_ks_mclust_means.txt hop_log_ks_mclust_data.txt hop_log_ks_mclust_densities.txt hop_log_ks_mclust_means.txt hop_vs_can_log_ks_mclust_data.txt hop_vs_can_log_ks_mclust_densities.txt hop_vs_can_log_ks_mclust_means.txt 0.01 2.0 original_data_log_mclust_hist</pre>