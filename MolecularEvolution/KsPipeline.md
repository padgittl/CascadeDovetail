# create fasta file list
<pre>ls -1 ../macse/exportAlignments/*_NT.fasta > fastaFileList.txt</pre>

# copy fresh control file
<pre>cp ~/path/to/paml4.9j/yn00.ctl .</pre>

# create "make_dirs.sh" file and create pair-specific control file for yn00
python scripts/prep_yn00_config.py yn00.ctl fastaFileList.txt
chmod +x make_dirs.sh
./make_dirs.sh

# create file list of gene pair-specific control files from above command
ls -1 *_vs_*.yn00.ctl > genePairSpecificCtlFiles.txt
genePairSpecificCtlFiles.txt.sh

# run yn00 and move output files to pair-specific dir
python scripts/run_yn00.py genePairSpecificCtlFiles.txt

# create yn00 file list
yn00OutFileList.txt.sh
ls -1 /path/to/dnds/*_vs_*/*_vs_*yn00.txt > yn00OutFileList.txt

less yn00OutFileList.txt | wc -l
3209

ls -1l *_vs_*/*_vs_*yn00.txt > fileSizes.txt
less fileSizes.txt  | sort -k5g | less
