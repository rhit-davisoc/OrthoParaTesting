import random

# Settings for namespace
namesp_fname = "./namespaces/namespace_10000.txt"

# Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
namesp_file = open(namesp_fname,"r")
namesp_str = namesp_file.read()
namesp_file.close()

namespace = namesp_str.strip('][').strip("'").split(', ')
random.shuffle(namespace)


# Build a pectinate tree weighted to one side: e.g (((((a,b),c),d),e),f);
def build_unbalanced(size,seed,nmsp):
    nwk = "(" * int(size-1)
    nwk += nmsp[0] + "," + nmsp[1] + ")"

    for i in range(2,size):
        nwk += "," + nmsp[i] + ")"

    nwk += ";"

    f = open("./balance_tests/pectinate" + str(size) + "_" + str(seed) + ".txt","w")
    f.write(nwk)
    f.close()

def build_unbalanced_right(size,seed,nmsp):
    nwk = ""
    for i in range(0,size-2):
        nwk += "(" + nmsp[i] + ","

    nwk += "(" + nmsp[size-2] + "," + nmsp[size-1]

    nwk += ")" * int(size-1)

    nwk += ";"

    f = open("./balance_tests/pectinateR" + str(size) + "_" + str(seed) + ".txt","w")
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

    f = open("./balance_tests/balanced" + str(size) + "_" + str(seed) +".txt","w")
    f.write(nwk)
    f.close()


def create_tests(seed,tree_size):
    random.seed(seed)

    namespace = namesp_str.strip('][').strip("'").split(', ')
    random.shuffle(namespace)

    build_balanced(tree_size,seed,namespace)
    build_unbalanced(tree_size,seed,namespace)
    build_unbalanced_right(tree_size,seed,namespace)

# for size in [2,4,8,16,32,64,128,256,512,1024,2048]:
#     for i in range(0,1000):
#         print("Creating balanced/pectinate tree " + str(i) + " of size " + str(size) + "\n")
#         create_tests(i,size)