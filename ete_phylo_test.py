from ete3 import PhyloTree
import dendropy
import itertools
import argparse

# p = "./input/high_1000.tre"

# tree = dendropy.Tree.get(path=p, schema="newick")

# nwk = tree.as_string(schema="newick",suppress_internal_taxon_labels=True,
#         suppress_internal_node_labels=True,suppress_rooting = False)

nwk = '((A,B),(C,D));'

t = PhyloTree(nwk)

def get_species_name(label):
    i = 0
    return label.split('|')[i]
# t.set_species_naming_function(get_species_name)

para = 0
ortho = 0

# Alternatively, you can scan the whole tree topology
events = t.get_descendant_evol_events()
# Print its orthology and paralogy relationships
# print("Events detected from the root of the tree")
# for ev in events:
#     if ev.etype == "S":
#         # print('   ORTHOLOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
#         ortho += 1
#     elif ev.etype == "D":
#         # print('   PARALOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
#         para += 1
#     print(','.join(ev.inparalogs))

# print("Ortho: %d" % ortho)
# print("Para: %d" % para)