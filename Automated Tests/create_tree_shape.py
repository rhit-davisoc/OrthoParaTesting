import random

# Settings for namespace
namesp_fname = "./namespaces/namespace_5000.txt"

# Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
namesp_file = open(namesp_fname,"r")
namesp_str = namesp_file.read()
namesp_file.close()

namespace = namesp_str.strip('][').strip("'").split(', ')
random.shuffle(namespace)

# Size multiple of 2 and greater than or equal to 2.
tree_size = 512


# Build an unbalanced tree weighted to one side: e.g (((((a,b),c),d),e),f);
def build_unbalanced(size):
    nwk = "(" * int(tree_size-1)
    nwk += namespace[0] + "," + namespace[1] + ")"

    for i in range(2,tree_size):
        nwk += "," + namespace[i] + ")"

    nwk += ";"

    f = open("./test_input/unbalanced_tree.txt","w")
    f.write(nwk)
    f.close()


# Build balanced tree (((a,b),(c,d)),((e,f),(g,h)))
def build_balanced_helper(size,index):
    if size == 2:
        leaves = "(" + namespace[index] + "," + namespace[index+1] + ")"
        return leaves
    else:
        return "(" + build_balanced_helper(size/2,int(index)) + "," + build_balanced_helper(size/2,int(size/2) + index) + ")"
    
def build_balanced(size):
    nwk = build_balanced_helper(size,0) + ";"

    f = open("./test_input/balanced_tree.txt","w")
    f.write(nwk)
    f.close()

build_balanced(512)