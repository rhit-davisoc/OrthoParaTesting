import dendropy
import itertools
import argparse

parser = argparse.ArgumentParser(
    description="Label orthologous/paralogous relationships of a given tree."
)
parser.add_argument(
    "--targets",
    nargs="*",
    help="genes that will have a file output of their relationships",
)
parser.add_argument(
    "--output",
    nargs="?",
    help="directory in which to output relationship information",
    default="",
)
parser.add_argument(
    "--input",
    nargs="?",
    help="newick tree file to identify relationships of",
    required=True,
)
parser.add_argument("--sep", nargs="?", help="separator of OTU and id")
parser.add_argument(
    "--display",
    nargs="?",
    help="indicator to display the tree in ascii format",
    default=False,
)
parser.add_argument(
    "--id_first",
    nargs="?",
    help="set true if the id is first in the separator",
    default=False,
)


args = parser.parse_args()

tree = dendropy.Tree.get(path=args.input, schema="newick")

if args.display:
    print(tree.as_ascii_plot())

if not args.sep:
    separator = input("Specify the separator: ")
    print("")
else:
    separator = args.sep


def get_species(label):
    i = 0
    if(args.id_first):
        i=1
    
    if " " in label:
        return label.split(" ")[i]

    elif separator not in label:
        print("Incorrect separator\n")
        exit()

    return label.split(separator)[i]


# Label pairs of taxons paralogs or orthologs
# Orthologs -> Result of a speciation event
# Paralogs -> Result of a duplication event
def label_tree_events(tree, start=1.0):
    rel_dict = {}
    for node in tree.postorder_node_iter():
        if node.is_leaf():
            node.species = [get_species(str(node.taxon))]
            node.taxon.sp_occurred = "NO"
            node.clade = [node.taxon]
            rel_dict[node.taxon.label] = {}
        else:
            assign_relationships(node, rel_dict)
    return rel_dict


# Assign paralogous or orthologous relationship(s) at a given node
def assign_relationships(node, rel_dict):
    children = node.child_nodes()
    num_children = node.num_child_nodes()

    if num_children == 2:
        node.species = children[0].species + children[1].species
        node.clade = children[0].clade + children[1].clade

        # Speciation occurred, label combination of species in each clade orthologous
        # (Ex. Child 1 has species {man-a, man-b} and child 2 has species {mouse-a,mouse-b},
        # so the combinations {man-a,mouse-a}, {man-a,mouse-b}, {man-b,mouse-a}, {man-b,mouse-b} would be orthologs)
        if not set(children[0].species) & set(children[1].species):
            event = "orthologous"

            for tax_1, tax_2 in itertools.product(children[0].clade, children[1].clade):
                tax_1.sp_occurred = "YES"
                tax_2.sp_occurred = "YES"
                rel_dict[tax_1.label][tax_2.label] = event
                rel_dict[tax_2.label][tax_1.label] = event

        else:
            # Duplication occurred, label combinations of species at this node paralogous
            for tax_1, tax_2 in itertools.product(children[0].clade, children[1].clade):
                if (tax_1.sp_occurred == "YES" or tax_2.sp_occurred == "YES") or not (
                    get_species(tax_1.label) == get_species(tax_2.label)
                ):
                    event = "out-paralogous"
                elif tax_1.sp_occurred == "UNKNOWN" or tax_2.sp_occurred == "UNKNOWN":
                    event = "paralogous"
                else:
                    event = "in-paralogous"

                rel_dict[tax_1.label][tax_2.label] = event
                rel_dict[tax_2.label][tax_1.label] = event
    else:
        # There are more than 2 children (polytomy)
        speciation = False
        duplication = False
        event = ""

        for i in range(0, num_children):
            for k in range(i, num_children):
                if i != k:
                    if speciation and duplication:
                        break
                    if not (set(children[i].species) & set(children[k].species)):
                        speciation = True
                    else:
                        duplication = True

        if speciation & duplication:
            poly_event = "ambiguous"
        elif speciation:
            poly_event = "orthologous"
        else:
            poly_event = "paralogous"

        node.clade = []
        node.species = []

        for child in children:
            node.clade += child.clade
            node.species += child.species

        for i in range(0, num_children):
            for k in range(i, num_children):
                if i != k:
                    if poly_event == "ambiguous":
                        if set(children[i].species) & set(children[k].species):
                            special_event = "paralogous"
                        else:
                            special_event = "ambiguous"

                        for tax_1, tax_2 in itertools.product(
                            children[i].clade, children[k].clade
                        ):
                            if tax_1.sp_occurred == "NO":
                                tax_1.sp_occurred = "UNKNOWN"
                            if tax_2.sp_occurred == "NO":
                                tax_2.sp_occurred = "UNKNOWN"

                            if special_event == "paralogous" and (
                                (
                                    tax_1.sp_occurred == "YES"
                                    or tax_2.sp_occurred == "YES"
                                )
                                or not (
                                    get_species(tax_1.label) == get_species(tax_2.label)
                                )
                            ):
                                event = "out-paralogous"
                            elif special_event == "paralogous":
                                event = "paralogous"
                            else:
                                event = "ambiguous"

                            rel_dict[tax_1.label][tax_2.label] = event
                            rel_dict[tax_2.label][tax_1.label] = event

                    elif poly_event == "orthologous":
                        for tax_1, tax_2 in itertools.product(
                            children[i].clade, children[k].clade
                        ):
                            event = "orthlogous"
                            tax_1.sp_occurred = "YES"
                            tax_2.sp_occurred = "YES"
                            rel_dict[tax_1.label][tax_2.label] = event
                            rel_dict[tax_2.label][tax_1.label] = event
                    else:
                        for tax_1, tax_2 in itertools.product(
                            children[i].clade, children[k].clade
                        ):
                            if (
                                tax_1.sp_occurred == "YES" or tax_2.sp_occurred == "YES"
                            ) or not (
                                get_species(tax_1.label) == get_species(tax_2.label)
                            ):
                                event = "out-paralogous"
                            elif (
                                tax_1.sp_occurred == "UNKNOWN"
                                or tax_2.sp_occurred == "UNKNOWN"
                            ):
                                event = "paralogous"
                            else:
                                event = "in-paralogous"

                            rel_dict[tax_1.label][tax_2.label] = event
                            rel_dict[tax_2.label][tax_1.label] = event


# Given a 2d dictionary of relationships, print out a the information in a matrix
def print_relationships_of_child(target_child, rel_dict):
    child_list = rel_dict.keys()

    print(target_child + ":")

    for child in child_list:
        if target_child != child:
            print("\t" + child + ": " + rel_dict[target_child][child])


def print_all_relationships(child_list, rel_dict):
    for child in child_list:
        print_relationships_of_child(child, rel_dict)
        print("")


# Given a 2d dictionary of relationships, output all relationship for the target taxon to a text file
def write_relationships_of_child(target_child, rel_dict, dir):
    child_list = rel_dict.keys()

    f = open(dir + target_child + ".txt", "w")

    f.write("taxon\t\trelationship\textra information\n")

    for child in child_list:
        if target_child != child:
            relationship = rel_dict[target_child][child]
            if relationship == "in-paralogous" or relationship == "out-paralogous":
                f.write(child + "\t\t" + "paralogous" + "\t\t" + relationship)
            else:
                f.write(child + "\t\t" + relationship)
            f.write("\n")


def write_all_relationships(child_list, rel_dict, dir):
    for child in child_list:
        write_relationships_of_child(child, rel_dict, dir)


rel_dict = label_tree_events(tree)

if args.targets:
    print_all_relationships(args.targets, rel_dict)
    write_all_relationships(args.targets, rel_dict, args.output)
else:
    print_all_relationships(rel_dict.keys(), rel_dict)
    write_all_relationships(rel_dict.keys(), rel_dict, args.output)