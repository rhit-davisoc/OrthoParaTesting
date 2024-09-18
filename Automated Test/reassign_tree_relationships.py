import dendropy
from ete3 import PhyloTree
import re
from statistics import mean, median

# Read in one tree from list and record it
treefname = "./test_output/high_ete_trials_5000.txt"

treefile = open(treefname,"r")
treestr = treefile.read()
treefile.close()

tree_arr = re.sub(r"Time:\s+[0-9]*\.[0-9]+\s+[A-Za-z0-9]+:\s",'',treestr).split('\n\n')

treenwk = tree_arr[0]

new_tree_file = open("./test_input/high_mix_5000.txt","w")
new_tree_file.write(treenwk)
new_tree_file.close()

# Convert tree to be all paralogous (all leaves are of species "AAA")
tree = dendropy.Tree.get(data=treenwk, schema="newick")

for i, leaf in enumerate(tree.leaf_nodes(),0):
    leaf.taxon.label = "AAA-" + str(i)

nwk = tree.as_string("newick")

new_tree_file = open("./test_input/high_para_5000.txt","w")
new_tree_file.write(nwk)
new_tree_file.close()

# Convert tree to be all orthologous (all species are different)
alph = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"]

for i, leaf in enumerate(tree.leaf_nodes(),0):
    leaf.taxon.label = alph[int(int(i/26)/26)%26] + alph[int(i/26)%26] + alph[i%26] + "-0"

nwk = tree.as_string("newick")

new_tree_file = open("./test_input/high_ortho_5000.txt","w")
new_tree_file.write(nwk)
new_tree_file.close()