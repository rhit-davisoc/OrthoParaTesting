import dendropy
from dendropy.simulate import treesim
import sys
sys.path.append('../')
import phylo_label_class
import time
from ete3 import PhyloTree
import itertools
from statistics import mean, median
import functools
import gc
import random


# Settings for experiments
otu_nums = [8,16,32,64,128,256,512,1024,2048,4096]
run_start = 0
run_end = 1000
verbose = True
output = "./results/runtime_comparison_8-8192.csv"

# Note: The input of the files can be changed when getRuntimes is called
# E.x. getRuntimes("random",otu_num,i,"./random_trees/random" + str(otu_num) + "_ID" + str(i) + ".txt",output)
# is how we get the runtimes of all the trees in the random_trees/ directory with naming convention 
# "random[otu_num]_ID[id].txt" (e.g "random128_ID574.txt" for the 578th random tree with 128 OTU's generated)


# Settings for namespace
separator = "-"

def get_namespace(namesp_fname):
    namesp_fname = "./namespaces/namespace_20000.txt"

    # Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
    namesp_file = open(namesp_fname,"r")
    namesp_str = namesp_file.read()
    namesp_file.close()

    namespace = namesp_str.strip('][').strip("'").split(', ')
    return namespace

# Generate a random tree with the given number of leaves (OTUs) and return the Newick format
@functools.lru_cache(maxsize = None)
def build_random_tree(otu_num,run,namespace_fname):

    namespace = get_namespace(namespace_fname)
    max_otu_num = len(namespace)
    taxa = dendropy.TaxonNamespace(namespace)

    if(otu_num > max_otu_num):
        print("OTU number exceeds maximum limit of" + str(max_otu_num))

    if(verbose):
        print("Building random tree...")

    s = otu_num + run
    r = random.Random(s)
    
    tree = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=otu_num,rng=r)

    nwk = tree.as_string(schema="newick",suppress_edge_lengths=True).strip("[&R] ").replace("'","")

    return nwk

def build_random_tree_species(otu_num,num_tax,s):
    namesp_fname = "./namespaces/" + str(otu_num) + "tot_" + str(num_tax) + "tax.txt"
    namesp_file = open(namesp_fname,"r")
    namesp_str = namesp_file.read()
    namesp_file.close()

    namespace = namesp_str.strip('][').strip("'").split(', ')
    max_otu_num = len(namespace)
    taxa = dendropy.TaxonNamespace(namespace)

    if(verbose):
        print("Building random tree...")

    r = random.Random(s)
    
    tree = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=otu_num,rng=r)

    nwk = tree.as_string(schema="newick",suppress_edge_lengths=True).strip("[&R] ").replace("'","")

    f = open("./species_num_tests/species" + str(num_tax) + "_id" + str(s) + ".txt","w")
    f.write(nwk)
    f.close()

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
    nwk=dendropy.Tree.get(path=nwk, schema="newick")
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
def get_phylotree_events(nwkfile):

    t = PhyloTree(nwkfile)

    if(verbose):
        print("Starting ETE PhyloTree labeling...")

    start = time.time()
    events = t.get_descendant_evol_events()
    end = time.time()
    elapsed_time = end - start

    # root = t.get_tree_root()
    # outg1, outg2=root.get_children()
    # smaller_out = outg1 if len(outg2) else outg2

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

def calculate_metrics(nwk):
    tree = dendropy.Tree.get(data=nwk, schema="newick")

    cti = dendropy.calculate.treemeasure.colless_tree_imbalance(tree)
    return cti

def time_check(nwk):
    start = time.time()
    tree = dendropy.Tree.get(data=nwk, schema="newick")
    end = time.time()
    elapsed = end - start
    return elapsed

def getRuntimes(type,otu,i,nwk_file,output):
    build_reltree_compact.cache_clear()
    rel_return = build_reltree_compact(nwk_file)
    build_reltree_compact.cache_clear()
    ete_return = get_phylotree_events(nwk_file)
    rel_time = rel_return[0]
    ete_time = ete_return[0]

    # cti = calculate_metrics(nwk)

    f = open(output,"a")
    f.write(type + "," + str(otu) + "," + str(i) + "," + str(rel_time) + "," + str(ete_time) + "\n")
    f.close()

def getRuntimesSpecies(species,otu,i,nwk,output):
    build_reltree_compact.cache_clear()
    rel_return = build_reltree_compact(nwk)
    build_reltree_compact.cache_clear()
    ete_return = get_phylotree_events(nwk)
    rel_time = rel_return[0]
    ete_time = ete_return[0]

    cti = calculate_metrics(nwk)

    f = open(output,"a")
    f.write(str(otu) + "," + str(species) + "," + str(i) + "," + str(cti) + "," + str(rel_time) + "," + str(ete_time) + "\n")
    f.close()

sys.setrecursionlimit(5000)

for otu_num in otu_nums:
    print("Running trials for " + str(otu_num) + " OTU trees. Starting at run " + str(run_start))
    for i in range(run_start,run_end):
        runs = run_end - run_start
        print("Run " + str(i + 1 - run_start) + " out of " + str(runs) + ".\n")

        getRuntimes("balanced",otu_num,i,"./balanced_and_pectinate_trees/balanced" + str(otu_num) + "_ID" + str(i) + ".txt",output)
        getRuntimes("pectinate_right",otu_num,i,"./balanced_and_pectinate_trees/pectinateR" + str(otu_num) + "_ID" + str(i) + ".txt",output)
        getRuntimes("random",otu_num,i,"./random_trees/random" + str(otu_num) + "_ID" + str(i) + ".txt",output)