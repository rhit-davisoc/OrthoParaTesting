from ete3 import PhyloTree
import sys
sys.path.append('../')
import phylo_label_class

# nwk = """
# ((Dme_001,Dme_002),(((Cfa_001,Mms_001),((Hsa_001,Hsa_003),Mmu_001)),
# (Ptr_002,(Hsa_002,Mmu_002))));
# """

nwk = "(HHH-1, (HHH-2,(HHH-3,BBB-3)));"

sep = "-"

RelTree = phylo_label_class.RelTree(nwk,sep)
RelTree.print_compact_relationship()

print("\n")

t = PhyloTree(nwk)
events = t.get_descendant_evol_events()
print("Events detected from the root of the tree")
for ev in events:
    if ev.etype == "S":
        print('   ORTHOLOGOUS RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
    elif ev.etype == "D":
        print('   PARALOGOUS RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))

matches = t.search_nodes(name="HHH-1")
human_seq = matches[0]
# Obtains its evolutionary events
events = human_seq.get_my_evol_events()
# Print its orthology and paralogy relationships
print("\nEvents detected that involve HHH-1:")
for ev in events:
    if ev.etype == "S":
        print('   ORTHOLOGOUS RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
        print('   IN-PARALOGS:', ','.join(ev.inparalogs))
    elif ev.etype == "D":
        print('   PARALOGOUS RELATIONSHIP:', ','.join(ev.in_seqs), "<====>", ','.join(ev.out_seqs))
        print('   IN-PARALOGS:', ','.join(ev.inparalogs))
        print('   OUT-PARALOGS:', ','.join(ev.outparalogs))