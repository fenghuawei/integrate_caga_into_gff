
import os
import pandas as pd
import sys

#f = open(sys.argv[1],'r')
cluster = sys.argv[2]
feature = sys.argv[3]
#f2 = open(sys.argv[3],'w')
#Tag_clusters_all.txt
data1 = pd.read_table(sys.argv[1],header=0)
sorted_df = data1.sort_values(by=[cluster, feature])

sorted_df.to_csv(sys.argv[4],sep='\t', index=False)

