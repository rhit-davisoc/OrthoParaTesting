import random
import math
from string import ascii_uppercase as alc

filename = "./namespaces/namespace_50000.txt"
tax_num = 50000

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

generate_namespace(tax_num)