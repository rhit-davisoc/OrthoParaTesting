import dendropy
import itertools
import argparse

parser = argparse.ArgumentParser(description='Label orthologous/paralogous relationships of a given tree.')
parser.add_argument('--targets', nargs='*',
                    help='genes that will have a file output of their relationships')
parser.add_argument('--dir', nargs='?',
                    help='directory in which to output relationship information',default='')
args = parser.parse_args()


# Tree to test with
s1 = "(man-1,man-2,zebra-2);"
tree = dendropy.Tree.get(
    data=s1,
    schema="newick")

# print(tree.as_string(schema="newick",) + "\n")
print(tree.as_ascii_plot())

separator = input('Specify the separator: ')
print('')

def get_species(label):
    if ' ' in label:
        return label.split(' ')[0]
    
    elif separator not in label:
            print('Incorrect separator\n')
            exit()

    return label.split(separator)[0]

# Label pairs of taxons paralogs or orthologs
# Orthologs -> Result of a speciation event
# Paralogs -> Result of a duplication event
def label_tree_events(tree, start=1.0):
    rel_dict = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.species = [get_species(str(node.taxon))]
            node.sp_occured = False #Indicator for if speciation occured
            node.clade = [str(node.taxon.label)]
            rel_dict[node.taxon.label] = {}
        else:
            assign_relationships(node, rel_dict)
    return rel_dict

# Assign paralogous or orthologous relationship(s) at a given node
def assign_relationships(node, rel_dict):
    children = node.child_nodes()
    num_children = node.num_child_nodes()

    if (num_children == 2):
        event = ""

        #Speciation occured, label combination of species in each clade orthologous
        #(Ex. Child 1 has species {man-a, man-b} and child 2 has species {mouse-a,mouse-b}, 
        #so the combinations {man-a,mouse-a}, {man-a,mouse-b}, {man-b,mouse-a}, {man-b,mouse-b} would be orthologs)
        if not set(children[0].species) & set(children[1].species):
            event = "orthologous"
            node.sp_occured = True
        else:
            #Duplication occured, label combinations of species at this node paralogous
            if (children[0].sp_occured or children[1].sp_occured):
                event = "out-paralogous"
                node.sp_occured = True
            else:
                event = "in-paralogous"
                node.sp_occured = False

        for tax_1,tax_2 in itertools.product(children[0].clade,children[1].clade):
            #Print out line by line as orthologous/paralogous relationships are found
            #print("%s and %s are %s" % (tax_1,tax_2,event))
            rel_dict[tax_1][tax_2] = event
            rel_dict[tax_2][tax_1] = event
        
        node.species = children[0].species + children[1].species
        node.clade = children[0].clade + children[1].clade
    else:
        #There are more than 2 children (polytomy)
        speciation = False
        duplication = False
        event = ""

        for i in range(0, num_children):
            for k in range(i, num_children):
                if(i != k):
                    if speciation and duplication:
                        break
                    if not (set(children[i].species) & set(children[k].species)):
                        speciation = True
                    else:
                        duplication = True

        if speciation & duplication:
            event = "ambiguous"
        elif speciation:
            event = "orthologous"
        else:
            event = "paralogous"

        node.clade = []
        node.species = []

        for child in children:
            node.clade += child.clade
            node.species += child.species

        for tax_1,tax_2 in itertools.product(node.clade,node.clade):
            if(tax_1 != tax_2):
                if(event == "ambiguous" and (get_species(tax_1) == get_species(tax_2))):
                    rel_dict[tax_1][tax_2] = "paralogous"
                else:
                    rel_dict[tax_1][tax_2] = event


# Given a 2d dictionary of relationships, print out a the information in a matrix
def print_relationships_of_child(target_child,rel_dict):
    child_list = rel_dict.keys()
    
    print(target_child + ":")

    for child in child_list:
        if(target_child != child):
            print("\t" + child + ": " + rel_dict[target_child][child])

def print_all_relationships(child_list,rel_dict):
    for child in child_list:
        print_relationships_of_child(child,rel_dict)
        print("")


# Given a 2d dictionary of relationships, output all relationship for the target taxon to a text file
def write_relationships_of_child(target_child,rel_dict,dir):
    child_list = rel_dict.keys()

    f = open(dir + target_child + ".txt",'w')
    
    f.write('taxon\t\trelationship\textra information\n')

    for child in child_list:
        if(target_child != child):
            relationship = rel_dict[target_child][child]
            if(relationship == "in-paralogous" or relationship == 'out-paralogous'):
                f.write(child + '\t\t' + 'paralogous' + '\t\t' + relationship)
            else:
                f.write(child + "\t\t" + relationship)
            f.write('\n')

def write_all_relationships(child_list,rel_dict,dir):
    for child in child_list:
        write_relationships_of_child(child,rel_dict,dir)



rel_dict = label_tree_events(tree)

if(args.targets):
    print_all_relationships(args.targets,rel_dict)
    write_all_relationships(args.targets,rel_dict,args.dir)
else:
    print_all_relationships(rel_dict.keys(),rel_dict)