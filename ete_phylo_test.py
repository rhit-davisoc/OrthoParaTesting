from ete3 import PhyloTree
import dendropy
import itertools
import argparse

p = "./input/low_1000.tre"

tree = dendropy.Tree.get(path=p, schema="newick")

nwk = tree.as_string(schema="newick",suppress_internal_taxon_labels=True,
        suppress_internal_node_labels=True,suppress_rooting = False)

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
print("Events detected from the root of the tree")
for ev in events:
    if ev.etype == "S":
        # print('   ORTHOLOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
        ortho += 1
    elif ev.etype == "D":
        # print('   PARALOGY RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
        para += 1
    print(','.join(ev.inparalogs))

print("Ortho: %d" % ortho)
print("Para: %d" % para)

#LOW
#V83 1000 - O 705 - P 294 - .084
#S93 1000 - O 712 - P 287 - .071
#X30 1000 - O 719 - P 280 - .082
#Q231 1000 - O 708 - P 291 - .094


#HIGH
#M15 1000 - O 722 - P277 - .843