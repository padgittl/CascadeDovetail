$ r8s

r8s>execute best_tree.nex
r8s>execute bootstrap_trees.nex

# remove characters and empty space from r8s output for input into the TreeAnnotator program -- this is specifically for the output from the bootstrap run
python ../scripts/removeExtraneousInfo.py r8s_bootstrap_output.txt
