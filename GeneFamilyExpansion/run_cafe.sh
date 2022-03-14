#! cafe 
# version 
# date 
load -i orthogroupCountFile.cafe.filtered.tsv -t 16 -l logfile.txt -p 0.05
tree (vitisVinifera:124.7894,(prunusPersica:119.1633,(ziziphusJujuba:103.5249,(morusNotabilis:82.6374,((cannabisSativa:22.6438,humulusLupulus:22.6438):36.6391,(tremaOrientale:7.6060,parasponiaAndersonii:7.6060):51.6769):23.3545):20.8875):15.6384):5.6261)
lambda -s
report resultfile
