import dendropy
import re

def getDiameter(root, diameter):
    children = root.child_nodes()

    if len(children) == 0:
        return 1, diameter
    else:
        left_height, diameter = getDiameter(children[0], diameter)
        right_height, diameter = getDiameter(children[1], diameter)

        max_diameter = left_height + right_height + 1

        diameter = max(diameter, max_diameter)

        return max(left_height, right_height) + 1, diameter

# Calculate metrics

def calculate_metrics(nwk,diams,b1s,nbars,ctis):
    tree = dendropy.Tree.get(data=nwk, schema="newick")

    x, diameter = getDiameter(tree.seed_node, 0)
    print("Diameter: %d" % diameter)
    diams.append(diameter)

    b1 = dendropy.calculate.treemeasure.B1(tree)
    print("B1: %f" % b1)
    b1s.append(b1)

    nbar = dendropy.calculate.treemeasure.N_bar(tree)
    print("N_Bar: %f" % nbar)
    nbars.append(nbar)

    cti = dendropy.calculate.treemeasure.colless_tree_imbalance(tree)
    print("CTI: %f" % cti)
    ctis.append(cti)


# Trees from time generated files
treefname = "./test_output/high_ete_trials_5000.txt"

treefile = open(treefname,"r")
treestr = treefile.read()
treefile.close()

tree_arr = re.sub(r"Time:\s+[0-9]*\.[0-9]+\s+[A-Za-z0-9]+:\s",'',treestr).split('\n\n')

diams = []
b1s = []
nbars = []
ctis = []

for i in range(0,len(tree_arr)-1):
    print("\nTree #%d:" %i)
    treenwk = tree_arr[i]

    calculate_metrics(treenwk,diams,b1s,nbars,ctis)

# nwk = '((A,B),(C,D));'

# calculate_metrics(nwk)