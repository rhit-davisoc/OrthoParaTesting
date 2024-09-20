import dendropy
from dendropy.simulate import treesim
import sys
sys.path.append('../')
from statistics import mean, median
import functools
import random
import argparse
import math

random.seed(1004)

# Settings for namespace
separator = "-"

def get_namespace(namesp_fname):
    # Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
    namesp_file = open(namesp_fname,"r")
    namesp_str = namesp_file.read()
    namesp_file.close()

    namespace = namesp_str.strip('][').strip("'").split(', ')
    random.shuffle(namespace)
    return namespace

# Generate a random tree with the given number of leaves (OTUs) and return the Newick format
@functools.lru_cache(maxsize = None)
def build_random_tree(namesp_fname,i,otu_num):

    # The seed for a tree with 1024 OTU's with ID 1 is 10240001, seed for the same size with ID 2 is 10240002.
    # Makes it easier to reproduce results
    # All ID's will be unique (given we stick with tests where ID's are under 1000, but this can be extended by adding my 0's)

    tree_seed = otu_num*10000 + i

    namespace = get_namespace(namesp_fname)
    max_otu_num = len(namespace)
    taxa = dendropy.TaxonNamespace(namespace)

    if(otu_num > max_otu_num):
        print("OTU number exceeds maximum limit of" + str(max_otu_num))

    if(verbose):
        print("Building random tree..." + str(tree_seed))

    r = random.Random(tree_seed)
    # r = random.seed()
    
    tree = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=otu_num,rng=r)

    nwk = tree.as_string(schema="newick",suppress_edge_lengths=True).strip("[&R] ").replace("'","")

    return nwk

def calculate_metrics(nwk):
    tree = dendropy.Tree.get(data=nwk, schema="newick")

    cti = dendropy.calculate.treemeasure.colless_tree_imbalance(tree)

    return cti

run_start = 0
run_end = 1000
verbose = True

tree_info_csv = "./random_trees_info/random_trees_info.csv"

for n in range(3,13): # Generate total number of OTU's to be powers of two, starting from 2^3 going to 2^12
    num_otu = pow(2,n)
    num_spec = math.ceil(math.sqrt(num_otu)) # Number of species is the ceiling
    filename = "./namespaces/" + str(num_otu) + "total_" + str(num_spec) + "spec.txt"
    for i in range(run_start,run_end):
        print("Creating random " + str(i) + " of size " + str(num_otu) + "\n")

        nwk = build_random_tree(filename,i,num_otu)

        cti = calculate_metrics(nwk)

        f = open(tree_info_csv, "a")
        f.write(str(num_otu) + "," + str(i) + "," + str(cti) + "\n")
        f.close()

        f = open("./random_trees/random" + str(num_otu) + "_ID" + str(i) + ".txt", "w")
        f.write(nwk)
        f.close()

# for otu_num in otu_nums:
#     print("Running trials for " + str(otu_num) + " OTU trees. Starting at run " + str(run_start))
#     for i in range(run_start,run_end):
#         runs = run_end - run_start
#         print("Run " + str(i + 1 - run_start) + " out of " + str(runs) + ".\n")

#         if(t==8):
#             tree_seed = otu_num*1000 + i

#             nwk = build_random_tree(tree_seed)

#             cti = calculate_metrics(nwk)

#             f = open(tree_info_csv, "a")
#             f.write(str(tree_seed) + "," + str(otu_num) + "," + str(cti) + "\n")
#             f.close()

#             f = open("./random_trees/random" + str(otu_num) + "_ID" + str(i) + ".txt", "w")
#             f.write(nwk)
#             f.close()
#         else:
#             tree_seed = otu_num*100 + i

#             nwk = build_random_tree(tree_seed)

#             cti = calculate_metrics(nwk)

#             f = open(tree_info_csv, "a")
#             f.write(str(tree_seed) + "," + str(otu_num) + "," + str(cti) + "\n")
#             f.close()

#             f = open("./test_trees/tree_" + str(tree_seed) + ".txt", "w")
#             f.write(nwk)
#             f.close()
        
