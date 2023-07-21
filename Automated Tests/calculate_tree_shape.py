import dendropy
import itertools
import argparse
from ete3 import PhyloTree
import re
from statistics import mean, median

treefname = "./test_output/low_ete_trials_5000.txt"

treefile = open(treefname,"r")
treestr = treefile.read()
treefile.close()

tree_arr = re.sub(r"Time:\s+[0-9]*\.[0-9]+\s+[A-Za-z0-9]+:\s",'',treestr).split('\n\n')

for i in range(0,len(tree_arr)-1):
    treenwk = tree_arr[i]
    tree = dendropy.Tree.get(data=treenwk, schema="newick")

    total_A = 0
    total_B = 0
    visits = 0
    total_frac = 0
    
    to_visit = []
    current = tree.seed_node
    all_events = []
    while current:
        # Gets childs and appends them to the To_visit list
        childs = current.child_nodes()
        to_visit += childs
        visits += 1
        if len(childs)>2:
            raise TypeError("nodes are expected to have two childs.")
        elif len(childs)==0:
            pass # leaf
        else:
            total_A += len(childs[0].leaf_nodes())
            total_B += len(childs[1].leaf_nodes())
        try:
            current = to_visit.pop(0)
        except IndexError:
            current = None

    # t = PhyloTree(treenwk)

    # para = 0
    # ortho = 0

    # # Alternatively, you can scan the whole tree topology
    # events = t.get_descendant_evol_events()
    # # Print its orthology and paralogy relationships
    # print("Events detected from the root of the tree")
    # for ev in events:
    #     if ev.etype == "S":
    #         # print('   ORTHOLOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
    #         ortho += 1
    #     elif ev.etype == "D":
    #         # print('   PARALOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
    #         para += 1

    # print("Ortho: %d" % ortho)
    # print("Para: %d" % para)

    print("Total A: %d\nTotal B: %d\nTotal: %d\n" % (total_A,total_B,total_A + total_B))