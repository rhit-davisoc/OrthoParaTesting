import dendropy
import itertools
import argparse
from ete3 import PhyloTree
import re
from statistics import mean, median

# treefname = "./test_output/low_1000.txt"

# treefile = open(treefname,"r")
# treestr = treefile.read()
# treefile.close()

# tree_arr = re.sub(r"Time:\s+[0-9]*\.[0-9]+\s+[A-Za-z0-9]+:\s",'',treestr).split('\n\n')

# for i in range(0,len(tree_arr)-1):
# treenwk = tree_arr[i]
# tree = dendropy.Tree.get(path='../input/low_1000.tre', schema="newick")
tree = dendropy.Tree.get(data='((A,B),(C,D));', schema="newick")

def getTotalChildrenCall(tree):
    total = 0

    to_visit = []
    current = tree.seed_node

    while current:
        # Gets childs and appends them to the To_visit list
        children = current.child_nodes()
        to_visit += children
        if len(children)>2:
            raise TypeError("nodes are expected to have two childs.")
        elif len(children)==0:
            total += 1
        else:
            total += len(children[0].leaf_nodes())
            total += len(children[1].leaf_nodes())
        try:
            current = to_visit.pop(0)
        except IndexError:
            current = None

    print("Total: %d\n" % total)

def calculateTotalDescendants(node, depth):
    children = node.child_nodes()

    if(len(children) < 2):
        return depth
    
    return depth + calculateTotalDescendants(children[0],depth + 1) + calculateTotalDescendants(children[1],depth + 1)

totalDesc = calculateTotalDescendants(tree.seed_node,0)

# getTotalChildrenCall(tree)

# print("Total Descendants: %d" % totalDesc)

si = dendropy.calculate.treemeasure.sackin_index(tree, normalize=False)

print(si)