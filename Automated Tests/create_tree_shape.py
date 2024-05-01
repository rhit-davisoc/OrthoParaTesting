import random

# Settings for namespace
namesp_fname = "./namespaces/namespace_10000.txt"

# Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
namesp_file = open(namesp_fname,"r")
namesp_str = namesp_file.read()
namesp_file.close()

namespace = namesp_str.strip('][').strip("'").split(', ')
random.shuffle(namespace)

# Size multiple of 2 and greater than or equal to 2.
tree_size = 4096


# Build a pectinate tree weighted to one side: e.g (((((a,b),c),d),e),f);
def build_unbalanced(size,seed,nmsp):
    nwk = "(" * int(size-1)
    nwk += nmsp[0] + "," + nmsp[1] + ")"

    for i in range(2,size):
        nwk += "," + nmsp[i] + ")"

    nwk += ";"

    f = open("./balance_tests/pectinate_" + str(seed) + ".txt","w")
    f.write(nwk)
    f.close()

def build_unbalanced_right(size,seed,nmsp):
    nwk = ""
    for i in range(0,size-2):
        nwk += "(" + nmsp[i] + ","

    nwk += "(" + nmsp[size-2] + "," + nmsp[size-1]

    nwk += ")" * int(size-1)

    nwk += ";"

    f = open("./balance_tests/pectinateR_" + str(seed) + ".txt","w")
    f.write(nwk)
    f.close()

# Build balanced tree e.g. (((a,b),(c,d)),((e,f),(g,h)))
def build_balanced_helper(nmsp,size,index):
    if size == 2:
        leaves = "(" + nmsp[index] + "," + nmsp[index+1] + ")"
        return leaves
    else:
        return "(" + build_balanced_helper(nmsp,size/2,int(index)) + "," + build_balanced_helper(nmsp,size/2,int(size/2) + index) + ")"
    
def build_balanced(size,seed,nmsp):
    nwk = build_balanced_helper(nmsp,size,0) + ";"

    f = open("./balance_tests/balanced_" + str(seed) +".txt","w")
    f.write(nwk)
    f.close()


def create_tests(seed):
    random.seed(seed)

    namespace = namesp_str.strip('][').strip("'").split(', ')
    random.shuffle(namespace)

    # build_balanced(tree_size,seed,namespace)
    # build_unbalanced(tree_size,seed,namespace)
    build_unbalanced_right(tree_size,seed,namespace)


for i in range(0,1000):
    print("Creating balanced/pectinate tree " + str(i) + "\n")
    create_tests(i)