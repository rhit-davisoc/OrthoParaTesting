import dendropy
from dendropy.simulate import treesim
import sys
sys.path.append('../')
from statistics import mean, median
import functools
import random
import argparse

parser = argparse.ArgumentParser(
    description="Label orthologous/paralogous relationships of a given tree."
)
parser.add_argument(
    "--test",
    nargs="*",
    help="test number to run",
)
args = parser.parse_args()

# Settings for namespace
namesp_fname = "./namespaces/namespace_50000.txt"
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
def build_random_tree(tree_seed):
    if(otu_num > max_otu_num):
        print("OTU number exceeds maximum limit of" + str(max_otu_num))

    if(verbose):
        print("Building random tree...")

    r = random.Random(tree_seed)
    
    tree = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=otu_num,rng=r)

    nwk = tree.as_string(schema="newick",suppress_edge_lengths=True).strip("[&R] ").replace("'","")

    return nwk

def calculate_metrics(nwk):
    tree = dendropy.Tree.get(data=nwk, schema="newick")

    cti = dendropy.calculate.treemeasure.colless_tree_imbalance(tree)

    return cti

# Settings for trial runs
t = int(args.test[0])

print("t: " + str(t))

if (t == 1):
    otu_nums = [100,200,300,400,500,600,700,800,900,1000]
elif(t == 2):
    otu_nums = [1000,2000,3000]
elif(t == 3):
    otu_nums = [4000,5000]
elif(t == 4):
    otu_nums = [6000,7000]
elif(t == 5):
    otu_nums = [8000]
elif(t == 6):
    otu_nums = [9000]
elif(t == 7):
    otu_nums = [10000]
elif(t==8):
    otu_nums = [4096]

run_start = 0
run_end = 1000
verbose = True

for otu_num in otu_nums:
    print("Running trials for " + str(otu_num) + " OTU trees. Starting at run " + str(run_start))
    for i in range(run_start,run_end):
        runs = run_end - run_start
        print("Run " + str(i + 1 - run_start) + " out of " + str(runs) + ".\n")

        if(t==8):
            tree_seed = 4096*1000 + i

            nwk = build_random_tree(tree_seed)

            cti = calculate_metrics(nwk)

            f = open("./test_trees_info/tree_info_4096" + str(t) + ".csv", "a")
            f.write(str(tree_seed) + "," + str(otu_num) + "," + str(cti) + "\n")
            f.close()

            f = open("./test_trees/tree_" + str(tree_seed) + ".txt", "w")
            f.write(nwk)
            f.close()
        else:
            tree_seed = otu_num*100 + i

            nwk = build_random_tree(tree_seed)

            cti = calculate_metrics(nwk)

            f = open("./test_trees_info/tree_info_" + str(t) + ".csv", "a")
            f.write(str(tree_seed) + "," + str(otu_num) + "," + str(cti) + "\n")
            f.close()

            f = open("./test_trees/tree_" + str(tree_seed) + ".txt", "w")
            f.write(nwk)
            f.close()
        