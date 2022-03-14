# get tree figure
# rapid
python ../../../../../scripts/cafetutorial_draw_tree.v2.py -i reports/summary_run1_node.txt -t '(vitisVinifera:124.7894,(prunusPersica:119.1633,(ziziphusJujuba:103.5249,(morusNotabilis:82.6374,((cannabisSativa:22.6438,humulusLupulus:22.6438):36.6391,(tremaOrientale:7.6060,parasponiaAndersonii:7.6060):51.6769):23.3545):20.8875):15.6384):5.6261)' -d '(vitisVinifera<0>,(prunusPersica<2>,(ziziphusJujuba<4>,(morusNotabilis<6>,((cannabisSativa<8>,humulusLupulus<10>)<9>,(tremaOrientale<12>,parasponiaAndersonii<14>)<13>)<11>)<7>)<5>)<3>)<1>' -o reports/summary_run1_tree_rapid -y Rapid

# expansions
python ../../../../../scripts/cafetutorial_draw_tree.v2.py -i reports/summary_run1_node.txt -t '(vitisVinifera:124.7894,(prunusPersica:119.1633,(ziziphusJujuba:103.5249,(morusNotabilis:82.6374,((cannabisSativa:22.6438,humulusLupulus:22.6438):36.6391,(tremaOrientale:7.6060,parasponiaAndersonii:7.6060):51.6769):23.3545):20.8875):15.6384):5.6261)' -d '(vitisVinifera<0>,(prunusPersica<2>,(ziziphusJujuba<4>,(morusNotabilis<6>,((cannabisSativa<8>,humulusLupulus<10>)<9>,(tremaOrientale<12>,parasponiaAndersonii<14>)<13>)<11>)<7>)<5>)<3>)<1>' -o reports/summary_run1_tree_expansions -y Expansions

# contractions
python ../../../../../scripts/cafetutorial_draw_tree.v2.py -i reports/summary_run1_node.txt -t '(vitisVinifera:124.7894,(prunusPersica:119.1633,(ziziphusJujuba:103.5249,(morusNotabilis:82.6374,((cannabisSativa:22.6438,humulusLupulus:22.6438):36.6391,(tremaOrientale:7.6060,parasponiaAndersonii:7.6060):51.6769):23.3545):20.8875):15.6384):5.6261)' -d '(vitisVinifera<0>,(prunusPersica<2>,(ziziphusJujuba<4>,(morusNotabilis<6>,((cannabisSativa<8>,humulusLupulus<10>)<9>,(tremaOrientale<12>,parasponiaAndersonii<14>)<13>)<11>)<7>)<5>)<3>)<1>' -o reports/summary_run1_tree_contractions -y Contractions

