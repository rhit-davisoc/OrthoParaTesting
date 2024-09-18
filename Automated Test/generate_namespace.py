import random
import math
from string import ascii_uppercase as alc

alph_array = ['A','B','C','D','E','F','G','H','I','J','K','L','M',
              'N','O','P','Q','R','S','T','U','V','W','X','Y','Z']

def generate_namespace(total):
    list = []
    count = 0
    for a in alc:
        if a == 'Z':
            repeat = total - count
        else:
            repeat = random.randint(1, math.ceil(total/13) - 1)

        for i in range(1,repeat + 1):
            list.append(a + a + a + "-" + str(i))
            count = count + 1
            if(count == total):
                break

        if(count == total):
            break
    
    f = open(filename,"w")
    out = str(list)

    f.write(out)
    f.close()

def num_to_name(num):
    name = ""

    while((int)(num/26) != 0):
        rem = num % 26
        name += alph_array[rem]
        num = (int)(num/26)

    rem = num % 26
    name += alph_array[rem]

    return name

def generate_namespace_num_species(total_tax,num_spec,filename):
    list = []
    count = 0

    for i in range(0,num_spec):
        if i == (num_spec-1):
            repeat = total_tax - count
        else:
            # repeat = random.randint(1, math.ceil(total_tax/(num_species/2)) - 1)
            repeat = math.ceil(total_tax/num_spec)

        for k in range(1,repeat + 1):
            name = num_to_name(i)
            list.append(name + "-" + str(k))
            count = count + 1
            if(count == total_tax):
                break

        if(count == total_tax):
            break
    
    f = open(filename,"w")
    out = str(list)

    f.write(out)
    f.close()

for i in range(3,14): # Generate total number of OTU's to be powers of two, starting from 2^3 going to 2^13
    num_otu = pow(2,i)
    num_spec = math.ceil(math.sqrt(num_otu)) # Number of species is the ceiling
    filename = "./namespaces/" + str(num_otu) + "total_" + str(num_spec) + "spec.txt"
    generate_namespace_num_species(num_otu,num_spec,filename)