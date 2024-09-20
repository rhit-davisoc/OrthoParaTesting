import random
import math

def get_namespace(namesp_fname):
    # Read in generated namespace file. Then create the namespace list and get total number of taxons (or OTUs)
    namesp_file = open(namesp_fname,"r")
    namesp_str = namesp_file.read()
    namesp_file.close()

    namespace = namesp_str.strip("][").replace("'",'').split(', ')

    random.shuffle(namespace)
    return namespace

# Build a pectinate tree weighted to one side: e.g (((((a,b),c),d),e),f);
def build_pectinate_left(size,seed,nmsp):
    nwk = "(" * int(size-1)
    nwk += nmsp[0] + "," + nmsp[1] + ")"

    for i in range(2,size):
        nwk += "," + nmsp[i] + ")"

    nwk += ";"

    f = open("./balanced_and_pectinate_trees/pectinate" + str(size) + "_ID" + str(seed) + ".txt","w")
    f.write(nwk)
    f.close()

def build_pectinate_right(size,seed,nmsp):
    nwk = ""
    for i in range(0,size-2):
        nwk += "(" + nmsp[i] + ","

    nwk += "(" + nmsp[size-2] + "," + nmsp[size-1]

    nwk += ")" * int(size-1)

    nwk += ";"

    f = open("./balanced_and_pectinate_trees/pectinateR" + str(size) + "_ID" + str(seed) + ".txt","w")
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

    f = open("./balanced_and_pectinate_trees/balanced" + str(size) + "_ID" + str(seed) +".txt","w")
    f.write(nwk)
    f.close()


def create_tests(namesp_fname,seed,tree_size):
    random.seed(seed)

    namespace = get_namespace(namesp_fname)

    build_balanced(tree_size,seed,namespace)
    build_pectinate_left(tree_size,seed,namespace)
    build_pectinate_right(tree_size,seed,namespace)

for n in range(3,13): # Generate total number of OTU's to be powers of two, starting from 2^3 going to 2^12
    num_otu = pow(2,n)
    num_spec = math.ceil(math.sqrt(num_otu)) # Number of species is the ceiling
    filename = "./namespaces/" + str(num_otu) + "total_" + str(num_spec) + "spec.txt"
    for i in range(0,1000):
        print("Creating balanced/pectinate tree " + str(i) + " of size " + str(num_otu) + "\n")
        create_tests(filename,i,num_otu)