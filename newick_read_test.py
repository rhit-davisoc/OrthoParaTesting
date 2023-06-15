import dendropy
import itertools

# Tree to test with
s1 = "((man-a,man-b),((mouse-a,zebra-a),(mouse-c,deer-c)));"

tree = dendropy.Tree.get(
    data=s1,
    schema="newick")

rel_dict = {}

print(tree.as_string(schema="newick",) + "\n")
print(tree.as_ascii_plot())

def get_geneid(label):
    if "-" in label:
        return label.split('-')[1]
    else:
        return ""

def get_species(label):
    if "-" in label:
        return label.split('-')[0]
    else:
        return ""

# Label pairs of taxons paralogs or orthologs
# Orthologs -> Have a least common ancestor, is a result of speciation
# Paralogs -> Result of duplication (LCA is a duplication event)

def label_tree_events(tree, start=1.0):
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.species = [get_species(str(node.taxon))]
            node.gene = [get_geneid(str(node.taxon))]
            node.clade = [str(node.taxon)]
            rel_dict[str(node.taxon)] = {} 
        else:
            assign_event(node)
    return tree


def assign_event(node):
    if (node.num_child_nodes() == 2):
        children = node.child_nodes()
        event = ""

        #Speciation occured, label orthologs
        if not set(children[0].species) & set(children[1].species):
            event = "orthologs"
        else:
            event = "paralogs"

        for tax_1,tax_2 in itertools.product(children[0].clade,children[1].clade):
            print("%s and %s are %s" % (tax_1,tax_2,event))
            rel_dict[tax_1][tax_2] = event
            rel_dict[tax_2][tax_1] = event
        
        node.species = children[0].species + children[1].species
        node.clade = children[0].clade + children[1].clade


def print_2d_dict(dic,children):
    
    return


label_tree_events(tree)

# for node in tree.postorder_node_iter():
#     print("%s : %s" % (node.taxon, get_geneid(str(node.taxon))))
