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

# Settings for namespace
namesp_fname = "./namespaces/namespace_5000.txt"
separator = "-"

# Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
namesp_file = open(namesp_fname,"r")
namesp_str = namesp_file.read()
namesp_file.close()

namespace = namesp_str.strip('][').strip("'").split(', ')
max_otu_num = len(namespace)

# Generate a random tree with the given number of leaves (OTUs) and return the Newick format
def build_random_tree(otu_num):
    if(otu_num > max_otu_num):
        print("OTU number exceeds maximum limit of" + str(max_otu_num))

    if(verbose):
        print("Building random tree...")

    taxa = dendropy.TaxonNamespace(namespace)
    tree = treesim.birth_death_tree(0.4, 0.1, taxon_namespace=taxa,num_extant_tips=otu_num)
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


# Settings for trial runs
# otu_nums = [10,20,30,40,50,100,200,300,400,500,600,700,800,900,1000,2000,3000,4000,5000]
otu_nums = [3500,4500,5000]
runs = 100
verbose = True

# print("Running " + str(runs) + " time(s) with " + str(otu_num) + " OTU tree(s)...\n")

rel_times = []
ete_times = []

low_times = []
low_newick = []

high_times = []
high_newick = []

otus = []

for otu_num in otu_nums:
    print("Running trials for " + str(otu_num) + " OTU trees.")
    for i in range(0,runs):
        print("Run " + str(i + 1) + " out of " + str(runs) + ".\n")

        nwk = build_random_tree(otu_num)
        # nwk = nwk_from_file("./test_input/low_para_5000.txt")

        rel_return = build_reltree_dict(nwk)
        ete_return = get_phylotree_events(nwk)
        rel_time = rel_return[0]
        ete_time = ete_return[0]

        rel_times.append(rel_time)
        ete_times.append(ete_time)
        otus.append(otu_num)

        # df.loc[len(df.index)] = ['Amy', 89, 93] 

        # if (ete_return[0] < .04):
        #     low_times.append(ete_return[0])
        #     low_newick.append(nwk)

        # if (ete_return[0] > .3):
        #     high_times.append(ete_return[0])
        #     high_newick.append(nwk)

        # RelTree_dict = rel_return[1]
        # ETE_events = ete_return[1]

        # compare_results(nwk, RelTree_dict, ETE_events)

# Get time stats
# rel_min = min(rel_times)
# rel_max = max(rel_times)
# rel_mean = mean(rel_times)

# ete_min = min(ete_times)
# ete_max = max(ete_times)
# ete_mean = mean(ete_times)

# print("Our algorithm time statistics: ")
# print("Min: %f" % rel_min)
# print("Max: %f" % rel_max)
# print("Mean: %f" % rel_mean)

# print("\nETE algorithm time statistics: ")
# print("Min: %f" % ete_min)
# print("Max: %f" % ete_max)
# print("Mean: %f" % ete_mean)

times_dict = {"OTUs":otus,"Our Time":rel_times,"ETE Time":ete_times}

df = pd.DataFrame(times_dict)

df.to_pickle("./time_data/time_data_3.pkl") 

# f = open("./test_output/low_ete_trials_512.txt","w")

# for i in range(0,len(low_times)):
#     f.write("Time: %f\nNewick: %s\n" % (low_times[i], low_newick[i]))

# f.close()

# f = open("./test_output/high_ete_trials_512.txt","w")

# for i in range(0,len(high_times)):
#     f.write("Time: %f\nNewick: %s\n" % (high_times[i], high_newick[i]))

# f.close()