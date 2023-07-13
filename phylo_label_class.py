import dendropy
import itertools

class RelTree:
    def __init__(self,nwk,sep,*,targets=None,id_first=False):
        self.tree = dendropy.Tree.get(data=nwk, schema="newick")
        self.separator = sep
        self.targets = targets
        self.id_first = id_first

    def get_species(self,label):
        i = 0
        if(self.id_first):
            i=1
        
        if " " in label:
            return label.split(" ")[i]

        elif self.separator not in label:
            print("Incorrect separator\n")
            exit()

        return label.split(self.separator)[i]


    # Label pairs of taxons paralogs or orthologs
    # Orthologs -> Result of a speciation event
    # Paralogs -> Result of a duplication event
    def label_tree_events(self,tree):
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


    # Assign paralogous or orthologous relationship(s) at a given node
    def assign_relationships(self, node, rel_dict):
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
                        self.get_species(tax_1.label) == self.get_species(tax_2.label)
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
                                        self.get_species(tax_1.label) == self.get_species(tax_2.label)
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
                                    self.get_species(tax_1.label) == self.get_species(tax_2.label)
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
    def print_relationships_of_child(self,target_child, rel_dict):
        child_list = rel_dict.keys()

        print(target_child + ":")

        for child in child_list:
            if target_child != child:
                print("\t" + child + ": " + rel_dict[target_child][child])


    def print_all_relationships(self,child_list, rel_dict):
        for child in child_list:
            self.print_relationships_of_child(child, rel_dict)
            print("")


    # Given a 2d dictionary of relationships, output all relationship for the target taxon to a text file
    def write_relationships_of_child(self,target_child, rel_dict, dir):
        child_list = rel_dict.keys()

        fname = dir + target_child.replace(self.separator,"") + ".txt"

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


    def write_all_relationships(self,child_list, rel_dict, dir):
        for child in child_list:
            self.write_relationships_of_child(child, rel_dict, dir)


    def get_relationship_dict(self):
        return self.label_tree_events(self.tree)