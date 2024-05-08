from ete3 import PhyloTree
import dendropy
import itertools
import argparse

# p = "./input/high_1000.tre"

# tree = dendropy.Tree.get(path=p, schema="newick")

# nwk = tree.as_string(schema="newick",suppress_internal_taxon_labels=True,
#         suppress_internal_node_labels=True,suppress_rooting = False)

# nwk = '((A,(B,E)),(C,D));'

# nwk = '(A,(B,(D,(E,(F,(G,(H,(I,J))))))));'

# nwk = '(((((((((A,B),C),D),E),F),G),H),I),J);'

f = open("./Automated Tests/balance_tests/pectinateR_1.txt")
nwk = f.read()
f.close()

t = PhyloTree(nwk)

def get_species_name(label):
    i = 0
    return label.split('|')[i]
# t.set_species_naming_function(get_species_name)

para = 0
ortho = 0

root = t.get_tree_root()

# Checks that is actually rooted
outgroups = root.get_children()
if len(outgroups) != 2:
    raise TypeError("Tree is not rooted")

# Cautch the smaller outgroup (will be stored as the tree outgroup)
o1 = set([n.name for n in outgroups[0].get_leaves()])
o2 = set([n.name for n in outgroups[1].get_leaves()])


if len(o2)<len(o1):
    smaller_outg = outgroups[1]
else:
    smaller_outg = outgroups[0]

print("Outg1: " + str(len(o1)) + " | Outg2: " + str(len(o2)))
print("Len of smaller_out: " + str(len(smaller_outg)))

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