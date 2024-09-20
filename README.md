# Testing for the Ortho/Para Labeling Tool

## Notes and assumptions
The namespaces are generated such that the number of species is the square-root of the total number of OTU's rounded up. 

All trees are in proper newick format when run.

For the sake of comparing directly to ete3 for each tree, no multifurcating trees were tested to evaluate runtime


## Conda environments
To build the exact conda environment run

`conda create --name [environment name] --file spec-file.txt`


## To re-run the original 

Change directories with `cd Automated\ Test/`

### Generate the balanced and pectinate trees run

Run `python ./generate_bal_pect_trees.py`

The for-loop can be modified to do a different range of OTU numbers, note that this may require generating new namespaces with the `generate_namespace.py` script. 

Note: In order to have balanced trees, this must be in powers of 2. However, random trees can be any number of OTUs.

### Generate 'random' trees with dendropy

Run `python ./generate_random_trees.py`

### Run the experiments and save runtimes of both algorithms to a csv

Run `python ./auto_tree_test.py`

The parameters can be changed at the top to specify output and input can be changed when the getRuntimes function is called
(see comments in the auto_tree_test file)