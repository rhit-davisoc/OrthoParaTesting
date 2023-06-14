import dendropy

# Simple test tree in newick format (Begins with duplication of a and b gene, then two subsequent speciation events (man and mouse))
s1 = "((man-a,mouse-a),(man-b,mouse-b));"
tree = dendropy.Tree.get(
    data=s1,
    schema="newick")

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

# Label the events at each node (not including leaf nodes) to be either speciation or duplication
# Speciation -> Leads to different taxonomy (e.g. results in man and mouse split)
# Duplication -> Gene duplication (e.g. hemoglobin a and homoglobin b split)

def label_tree_events(tree, start=1.0):
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.species = get_species(str(node.taxon))
            node.gene = get_geneid(str(node.taxon))
            node.event = "leaf"
        else:
            assign_event(node)
    return tree


def assign_event(node):
    if (node.num_child_nodes() == 2):
        children = node.child_nodes()
        if (children[0].species == children[1].species):
            node.event = "duplication"
            node.species = ""
            node.gene = ""
        else:
            node.event = "speciation"
            node.species = ""
            node.gene = ""
    else:
           node.event = ""


label_tree_events(tree)

for node in tree.postorder_node_iter():
    print("%s : %s : %s" % (node.taxon, node.event, get_geneid(str(node.taxon))))
