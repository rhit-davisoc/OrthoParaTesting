import dendropy
import itertools

# Tree to test with
s1 = "((man-a,man-b),((mouse-a,zebra-a),(mouse-c,deer-c)));"

tree = dendropy.Tree.get(
    data=s1,
    schema="newick")

# print(tree.as_string(schema="newick",) + "\n")
# print(tree.as_ascii_plot())

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
# Orthologs -> Result of a speciation event
# Paralogs -> Result of a duplication event

def label_tree_events(tree, start=1.0):
    rel_dict = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.species = [get_species(str(node.taxon))]
            node.gene = [get_geneid(str(node.taxon))]
            node.clade = [str(node.taxon)]
            rel_dict[str(node.taxon)] = {} 
        else:
            assign_relationships(node, rel_dict)
    return rel_dict

# Assign paralogous or orthologous relationship(s) at a given node
def assign_relationships(node, rel_dict):
    if (node.num_child_nodes() == 2):
        children = node.child_nodes()
        event = ""

        #Speciation occured, label combination of species in each clade orthologous
        #(Ex. Child 1 has species {man-a, man-b} and child 2 has species {mouse-a,mouse-b}, 
        #so the combinations {man-a,mouse-a}, {man-a,mouse-b}, {man-b,mouse-a}, {man-b,mouse-b} would be orthologs)
        if not set(children[0].species) & set(children[1].species):
            event = "orthologous"
        else:
            #Duplication occured, label combinations of species at this node paralogous
            event = "paralogous"

        for tax_1,tax_2 in itertools.product(children[0].clade,children[1].clade):
            #Print out line by line as orthologous/paralogous relationships are found
            #print("%s and %s are %s" % (tax_1,tax_2,event))
            rel_dict[tax_1][tax_2] = event
            rel_dict[tax_2][tax_1] = event
        
        node.species = children[0].species + children[1].species
        node.clade = children[0].clade + children[1].clade


# Given a 2d dictionary of relationships, print out a the information in a matrix
def print_2d_dict(rel_dict):
    child_list = rel_dict.keys()

    print(f"{'':<10}", end = " ")
    for child in child_list:
        print("\t" + child, end = " ")

    print("\n")
    
    for child_row in child_list:
        print(f"{child_row:<15}", end = " ")
        for child_col in child_list:
            if(child_row != child_col):
                print(rel_dict[child_row][child_col] + "\t", end = " ")
            else:
                print("\t\t",end = " ")
        print("\n")
        

rel_dict = label_tree_events(tree)

print_2d_dict(rel_dict)