import dendropy
import itertools
import argparse
import re

p = "./input/ncbi.tre"

tree = dendropy.Tree.get(path=p, schema="newick")

num_children = len(tree.leaf_nodes())

print("There are " + str(num_children) + " children")

tree.resolve_polytomies()

nwk = tree.as_string(schema="newick",suppress_internal_taxon_labels=True,
        suppress_internal_node_labels=True,suppress_rooting = False)

print("seed node children: " + str(tree.seed_node.num_child_nodes()))

# f = open("./tree.tre","w")
# f.write(nwk)

