import dendropy
import itertools
import argparse
import re
from statistics import mean, median

treefname = "./input/low_ete_trials_5000_2.txt"

treefile = open(treefname,"r")
treestr = treefile.read()
treefile.close()

tree_arr = re.sub(r"Time:\s+[0-9]*\.[0-9]+\s+[A-Za-z0-9]+:\s",'',treestr).split('\n\n')

for i in range(0,len(tree_arr)-1):
        tree = tree_arr[i]
        tree = dendropy.Tree.get(data=tree, schema="newick")

        for eg in tree.postorder_edge_iter():
                eg.length = 1

        dists = tree.calc_node_root_distances(return_leaf_distances_only=True)

        max_dist = max(dists)
        min_dist = min(dists)
        mean_dist = mean(dists)
        median_dist = median(dists)

        print("Min: %f\nMax: %f\n Median: %f\nMean: %f\n" % (min_dist, max_dist, median_dist, mean_dist))

