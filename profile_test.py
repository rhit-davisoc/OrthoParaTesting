import dendropy
from dendropy.simulate import treesim
import sys
sys.path.append('../')
import phylo_label_class
import time
from statistics import mean, median
import pandas as pd

nwk = "(c-1,(a-1,((b-1,(d-1,e-1)),c-2)));"
sep = "-"

RelTree = phylo_label_class.RelTree(nwk,sep)

start = time.time()
RelTree.label_tree_events_compact()
end = time.time()
elapsed_time = end - start

print("time: ", elapsed_time)