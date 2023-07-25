import dendropy
import itertools

class Event:
    def __init__(self):
        self.type = None
        self.leftSp = None
        self.leftNoSp = None
        self.rightSp = None
        self.rightNoSp = None

class RelTree:
    def __init__(self, nwk, sep, *, targets=None, id_first=False):
        self.tree = dendropy.Tree.get(data=nwk, schema="newick")
        self.separator = sep
        self.targets = targets
        self.id_first = id_first

    def display_tree(self):
        print(self.tree.as_ascii_plot())

    def get_species(self, label):
        i = 0
        if self.id_first:
            i = 1

        if " " in label:
            return label.split(" ")[i]

        elif self.separator not in label:
            print("Incorrect separator\n")
            exit()

        return label.split(self.separator)[i]

    # Label pairs of taxons paralogs or orthologs
    # Orthologs -> Result of a speciation event
    # Paralogs -> Result of a duplication event
    def label_tree_events(self, tree):
        rel_dict = {}
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                node.species = [self.get_species(str(node.taxon))]
                node.taxon.sp_occurred = "NO"
                node.clade = [node.taxon]
                rel_dict[node.taxon.label] = {}
            else:
                self.assign_relationships(node, rel_dict)
        return rel_dict

    # Speciation occurred, label combination of species in each clade orthologous
    # (Ex. Child 1 has species {man-a, man-b} and child 2 has species {mouse-a,mouse-b},
    # so the combinations {man-a,mouse-a}, {man-a,mouse-b}, {man-b,mouse-a}, {man-b,mouse-b} would be orthologs)
    def assign_orthologous(self, clade1, clade2, rel_dict):
        event = "orthologous"
        for tax_1, tax_2 in itertools.product(clade1, clade2):
            tax_1.sp_occurred = "YES"
            tax_2.sp_occurred = "YES"
            rel_dict[tax_1.label][tax_2.label] = event
            rel_dict[tax_2.label][tax_1.label] = event

    # Duplication occured, label combination of species in each clade paralogous
    def assign_paralogous(self, clade1, clade2, rel_dict):
        for tax_1, tax_2 in itertools.product(clade1, clade2):
            if tax_1.sp_occurred == "YES" or tax_2.sp_occurred == "YES":
                event = "out-paralogous"
            elif (
                tax_1.sp_occurred == "NO" and tax_2.sp_occurred == "NO"
            ):
                event = "in-paralogous"
            else:
                if(self.get_species(tax_1.label) == self.get_species(tax_2.label)):
                    event = "paralogous"
                else:
                    event = "out-paralogous"

            rel_dict[tax_1.label][tax_2.label] = event
            rel_dict[tax_2.label][tax_1.label] = event

    def assign_ambigious(self, child1, child2, rel_dict):
        if set(child1.species) & set(child2.species):
            special_event = "paralogous"
        else:
            special_event = "ambiguous"

        for tax_1, tax_2 in itertools.product(
            child1.clade, child2.clade
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
                    self.get_species(tax_1.label)
                    == self.get_species(tax_2.label)
                )
            ):
                event = "out-paralogous"
            elif special_event == "paralogous":
                event = "paralogous"
            else:
                event = "ambiguous"

            rel_dict[tax_1.label][tax_2.label] = event
            rel_dict[tax_2.label][tax_1.label] = event

    def get_poly_event(self,num_children, children):
        speciation = False
        duplication = False

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
            return "ambiguous"
        elif speciation:
            return "orthologous"
        else:
            return "paralogous"

    # Assign paralogous or orthologous relationship(s) at a given node and store it in a dictionary
    # Input:
    # node -> The node to assign relationships at. (Each children's clade will be assigned a relationship to the other clades).
    # rel_dict -> The dictionary to store the relationship information in.
    # Relationships can be orthologous, in-paralogous, out-paralogous, and ambigious.
    # Ambigious cases only happen after/during polytomies when speciation events are unknown between two OTU.
    def assign_relationships(self, node, rel_dict):
        children = node.child_nodes()
        num_children = node.num_child_nodes()

        if num_children == 2:
            node.species = children[0].species + children[1].species
            node.clade = children[0].clade + children[1].clade

            if not set(children[0].species) & set(children[1].species):
                self.assign_orthologous(children[0].clade, children[1].clade, rel_dict)
            else:
                self.assign_paralogous(children[0].clade, children[1].clade, rel_dict)
        else:
            # There are more than 2 children (polytomy)
            poly_event = self.get_poly_event(num_children, children)

            node.clade = []
            node.species = []

            for child in children:
                node.clade += child.clade
                node.species += child.species

            for i in range(0, num_children):
                for k in range(i, num_children):
                    if i != k:
                        if poly_event == "ambiguous":
                            self.assign_ambigious(children[i], children[k], rel_dict)
                        elif poly_event == "orthologous":
                            self.assign_orthologous(children[i].clade, children[k].clade, rel_dict)
                        else:
                            self.assign_paralogous(children[i].clade, children[k].clade, rel_dict)

    #Compact version of labeling a tree
    def label_tree_events_compact(self):
        tree = self.tree
        for node in tree.postorder_node_iter():
            if node.is_leaf():
                node.species = [self.get_species(str(node.taxon))]
                node.sp_occur = set()
                node.no_sp_occur = set([node.taxon.label])
            else:
                self.assign_relationships_compact(node)

    def assign_relationships_compact(self, node):
        children = node.child_nodes()
        num_children = node.num_child_nodes()

        if num_children == 2:
            node.species = children[0].species + children[1].species

            if not set(children[0].species) & set(children[1].species):
                node.event = "S" #-> Speciation
                node.sp_occur = children[0].sp_occur | children[1].sp_occur | children[0].no_sp_occur | children [1].no_sp_occur
                node.no_sp_occur = set()
            else:
                node.event = "D" #-> Duplication
                node.sp_occur = children[0].sp_occur | children[1].sp_occur
                node.no_sp_occur = children[0].no_sp_occur | children [1].no_sp_occur
        elif num_children == 1:
            print(num_children)

    def print_compact_relationship(self):
        self.label_tree_events_compact()

        tree = self.tree

        for node in tree.postorder_internal_node_iter():
            children = node.child_nodes()
            if node.event == "S":
                print('   ORTHOLOGY RELATIONSHIP:', ','.join((children[0].sp_occur | children[0].no_sp_occur)), "<====>", ','.join(children[1].sp_occur | children[1].no_sp_occur))
            else:
                if(children[0].no_sp_occur and children[1].no_sp_occur):
                    print('   IN-PARALOGOUS RELATIONSHIP:', ','.join(children[0].no_sp_occur), "<====>", ','.join(children[1].no_sp_occur))

                if(children[0].sp_occur and (children[0].no_sp_occur or children[1].no_sp_occur)):
                    print('   OUT-PARALOGOUS RELATIONSHIP:', ','.join(children[0].sp_occur), "<====>", ','.join(children[1].sp_occur | children[1].no_sp_occur))
                if(children[1].sp_occur):
                    print('   OUT-PARALOGOUS RELATIONSHIP:', ','.join(children[0].sp_occur | children[0].no_sp_occur), "<====>", ','.join(children[1].sp_occur))


    # Given a 2d dictionary of relationships, print out a the information in a matrix
    def print_relationships_of_child(self, target_child, rel_dict):
        child_list = rel_dict.keys()

        print(target_child + ":")

        for child in child_list:
            if target_child != child:
                print("\t" + child + ": " + rel_dict[target_child][child])

    def print_all_relationships(self, child_list, rel_dict):
        for child in child_list:
            self.print_relationships_of_child(child, rel_dict)
            print("")

    # Given a 2d dictionary of relationships, output all relationship for the target taxon to a text file
    def write_relationships_of_child(self, target_child, rel_dict, dir):
        child_list = rel_dict.keys()

        fname = dir + target_child.replace(self.separator, "") + ".txt"

        f = open(fname, "w")

        f.write("taxon\t\trelationship\textra information\n")

        for child in child_list:
            if target_child != child:
                relationship = rel_dict[target_child][child]
                if relationship == "in-paralogous" or relationship == "out-paralogous":
                    f.write(child + "\t\t" + "paralogous" + "\t\t" + relationship)
                else:
                    f.write(child + "\t\t" + relationship)
                f.write("\n")

    def write_all_relationships(self, child_list, rel_dict, dir):
        for child in child_list:
            self.write_relationships_of_child(child, rel_dict, dir)

    def get_relationship_dict(self):
        return self.label_tree_events(self.tree)
