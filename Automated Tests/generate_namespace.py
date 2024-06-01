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

def generate_namespace_num_species(total_tax,num_tax,filename):
    list = []
    count = 0

    for i in range(0,num_tax):
        if i == (num_tax-1):
            repeat = total_tax - count
        else:
            # repeat = random.randint(1, math.ceil(total_tax/(num_species/2)) - 1)
            repeat = math.ceil(total_tax/num_tax)

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

total_tax = 2000
num_tax_list = [1,2,5,10,20,30,40,50,100,200,300,400,500,1000,2000]

for num_tax in num_tax_list:
    filename = "./namespaces/2000tot_" + str(num_tax) + "tax.txt"
    generate_namespace_num_species(total_tax,num_tax,filename)