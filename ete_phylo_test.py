from ete3 import PhyloTree
import dendropy
import itertools
import argparse

p = "./tree.tre"

tree = dendropy.Tree.get(path=p, schema="newick")

nwk = tree.as_string(schema="newick",suppress_internal_taxon_labels=True,
        suppress_internal_node_labels=True,suppress_rooting = False)

t = PhyloTree(nwk)

def get_species_name(label):
    i = 0
    return label.split('|')[i]
# t.set_species_naming_function(get_species_name)

# Alternatively, you can scan the whole tree topology
events = t.get_descendant_evol_events()
# Print its orthology and paralogy relationships
print("Events detected from the root of the tree")
for ev in events:
    if ev.etype == "S":
        print('   ORTHOLOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
    elif ev.etype == "D":
        print('   PARALOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))