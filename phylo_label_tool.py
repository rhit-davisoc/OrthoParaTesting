""" Modules providing argument parsing (argparse) and our phylogenetic 
    tree labeling class (phylo_label_class)."""
import argparse
import phylo_label_class

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

with open(args.input,"r",encoding="utf-8") as nwk_file:
    nwk = nwk_file.read()

RelTree = phylo_label_class.RelTree(nwk, args.sep)

if args.display:
    RelTree.display_tree()

if not args.sep:
    separator = input("Specify the separator: ")
    print("")
else:
    separator = args.sep

rel_dict = RelTree.get_relationship_dict()

if args.targets:
    RelTree.print_all_relationships(args.targets, rel_dict)
    RelTree.write_all_relationships(args.targets, rel_dict, args.output)
else:
    RelTree.print_all_relationships(rel_dict.keys(), rel_dict)
    # write_all_relationships(rel_dict.keys(), rel_dict, args.output)
