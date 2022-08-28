# priming run - optimization parameters go into the config file for step 2
$treePL best_tree_step1_priming_config.txt

# comment out priming command 
$treePL best_tree_step2_optimization_config.txt

# get smoothing value for the bootstrap run and comment out cross-validation command
$treePL bootstrap_config.txt

