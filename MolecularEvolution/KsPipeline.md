### Create fasta file list
<pre>ls -1 ../macse/exportAlignments/*_NT.fasta > fastaFileList.txt</pre>

### Copy fresh yn00 control file
<pre>cp ~/path/to/paml4.9j/yn00.ctl .</pre>

### Create "make_dirs.sh" file and a gene pair-specific control file for yn00
<pre>python scripts/prep_yn00_config.py yn00.ctl fastaFileList.txt</pre>
<pre>chmod +x make_dirs.sh</pre>
<pre>./make_dirs.sh</pre>

### Create file list of gene pair-specific control files from above command
<pre>ls -1 *_vs_*.yn00.ctl > genePairSpecificCtlFiles.txt</pre>

### Run yn00 and move output files to pair-specific directory
<pre>python scripts/run_yn00.py genePairSpecificCtlFiles.txt</pre>

### Create yn00 file list
<pre>ls -1 /path/to/dnds/*_vs_*/*_vs_*yn00.txt > yn00OutFileList.txt</pre>
