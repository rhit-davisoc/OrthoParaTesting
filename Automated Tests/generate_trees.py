import dendropy
from dendropy.simulate import treesim
import sys
sys.path.append('../')
import phylo_label_class
import time
from ete3 import PhyloTree
import itertools
from statistics import mean, median
import pandas as pd
import functools
import gc

# Settings for namespace
namesp_fname = "./namespaces/namespace_10000.txt"
separator = "-"

# Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
namesp_file = open(namesp_fname,"r")
namesp_str = namesp_file.read()
namesp_file.close()

namespace = namesp_str.strip('][').strip("'").split(', ')
max_otu_num = len(namespace)
taxa = dendropy.TaxonNamespace(namespace)

# Generate a random tree with the given number of leaves (OTUs) and return the Newick format
@functools.lru_cache(maxsize = None)
def build_random_tree(otu_num):
    if(otu_num > max_otu_num):
        print("OTU number exceeds maximum limit of" + str(max_otu_num))

    if(verbose):
        print("Building random tree...")

    tree = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=otu_num)

    # for i, leaf in enumerate(tree.leaf_nodes(),0):
    #     leaf.taxon.label = "AAA-" + str(i)

    nwk = tree.as_string(schema="newick",suppress_edge_lengths=True).strip("[&R] ").replace("'","")

    return nwk

def nwk_from_file(fname):
    f = open(fname,"r")
    nwk = f.read()
    f.close()
    return nwk

# Build relationship dictionary using our tool and time the process
def build_reltree_dict(nwk):
    RelTree = phylo_label_class.RelTree(nwk,separator)

    if(verbose):
        print("Starting RelTree labeling...")
    start = time.time()
    dict = RelTree.get_relationship_dict()
    end = time.time()

    elapsed_time = end - start

    if(verbose):
        print('RelTree dictionary building time:', elapsed_time, 'seconds\n')

    return (elapsed_time,dict)

@functools.lru_cache(maxsize = None)
def build_reltree_compact(nwk):
    RelTree = phylo_label_class.RelTree(nwk,separator)
    if(verbose):
        print("Starting RelTree labeling...")

    elapsed_time = 0

    start = time.time()
    RelTree.label_tree_events_compact()
    end = time.time()
    elapsed_time = end - start

    if(verbose):
        print('RelTree dictionary building time:', elapsed_time, 'seconds\n')

    return (elapsed_time,None)

# Compute evolution history using ETE tool and time the process
def get_phylotree_events(nwk):
    t = PhyloTree(nwk)

    if(verbose):
        print("Starting ETE PhyloTree labeling...")

    start = time.time()
    events = t.get_descendant_evol_events()
    end = time.time()
    elapsed_time = end - start

    if(verbose):
        print('ETE PhyloTree event labeling time:', elapsed_time, 'seconds\n\n')

    return (elapsed_time,events)

# Compare results of our tool and the ETE tool. Report back mismatches/total and newick tree if results do not match.
def compare_results(nwk, RelTree_dict, ETE_events):
    mismatches = 0
    total = 0

    for ev in ETE_events:
        if ev.etype == "D":
            for tax_1, tax_2 in itertools.product(ev.in_seqs, ev.out_seqs):
                relationship = RelTree_dict[tax_1][tax_2]
                if not "paralogous" in relationship:
                    print("Mismatch in labeling between " + tax_1 + " and " + tax_2 + "! +  Reltree: " + relationship + " | ETE: paralogous")
                    mismatches += 1
                total += 1
        else:
            for tax_1, tax_2 in itertools.product(ev.in_seqs, ev.out_seqs):
                relationship = RelTree_dict[tax_1][tax_2]
                if not "orthologous" in relationship:
                    print("Mismatch in labeling between " + tax_1 + " and " + tax_2 + "! Reltree: " + relationship + " | ETE: orthologous")
                    mismatches += 1
                total += 1

    if(mismatches > 0):
        print("Mismatched/Total: " + str(mismatches) + "/" + str(total) + "\n")
        print(nwk)

def getHeight(root):
    children = root.child_nodes()

    if len(children) == 0:
        return 1
    else:
        left_height = getHeight(children[0])
        right_height = getHeight(children[1])

        return max(left_height + 1, right_height + 1)

def getDiameter(root, diameter):
    children = root.child_nodes()

    if len(children) == 0:
        return 1, diameter
    else:
        left_height, diameter = getDiameter(children[0], diameter)
        right_height, diameter = getDiameter(children[1], diameter)

        max_diameter = left_height + right_height + 1

        diameter = max(diameter, max_diameter)

        return max(left_height, right_height) + 1, diameter

def calculate_metrics(nwk,ctis,sis,small_outgs):
    tree = dendropy.Tree.get(data=nwk, schema="newick")

    cti = dendropy.calculate.treemeasure.colless_tree_imbalance(tree)
    ctis.append(cti)

    si = dendropy.calculate.treemeasure.sackin_index(tree,normalize=False)
    sis.append(si)

    r_children = tree.seed_node.child_nodes()
    small_outg = min(len(r_children[0].leaf_nodes()),len(r_children[1].leaf_nodes()))
    small_outgs.append(small_outg)


# Settings for trial runs
# otu_nums = [10,20,30,40,50,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000]
# otu_nums = [100,200,300,400,500,600,700,800,900,1000,1100,1200,1300,1400,1500,1600,1700,1800,1900,2000]
# otu_nums = [100,200,300]
# otu_nums = [100,200,300,400,500,600,700,800,900,1000]
otu_nums = [1100]
runs = 1
verbose = True

# print("Running " + str(runs) + " time(s) with " + str(otu_num) + " OTU tree(s)...\n")

rel_times = []
ete_times = []

time_returns = []

otus = []

b1s = []
ctis = []
sis = []
diams = []
heights = []
small_outgs = []

def time_check(nwk):
    start = time.time()
    tree = dendropy.Tree.get(data=nwk, schema="newick")
    end = time.time()
    elapsed = end - start
    return elapsed

for otu_num in otu_nums:
    print("Running trials for " + str(otu_num) + " OTU trees.")
    for i in range(0,runs):
        print("Run " + str(i + 1) + " out of " + str(runs) + ".\n")

        nwk = build_random_tree(otu_num)
