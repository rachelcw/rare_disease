# Purpose: create a table with leafcutter results
# Input: leafcutter cluster_significance and effect sizes
# Output: leafcutter and spot results in bed format

import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("-c","--cluster", help="cluster_significance leafcutter")
parser.add_argument("-e","--effect", help="effect_sizes leafcutter")
parser.add_argument("-p","--pvalue", help="empercial_pvalue spot")
parser.add_argument("-md","--md", help="md spot - Mahalanobis distance")
parser.add_argument("-o","--output", help="output path")

args = parser.parse_args()

# # # LEAFCUTTER # # #

# # cluster_significance # #
cluster=pd.read_csv(args.cluster,sep="\t",usecols=["cluster","p.adjust","genes"])
cluster[["chr","cluster"]]=cluster['cluster'].str.split(':', expand=True)

# # effect sizes # #
ds=pd.read_csv(args.effect,sep="\t",usecols=["intron","deltapsi","healty","patient"])
ds[['chr', 'start', 'end','cluster']] = ds['intron'].str.split(':', expand=True)
ds.drop(['intron'], axis=1, inplace=True)

table=ds.merge(cluster, on=["chr","cluster"], how="left")
table[['cluster','strand']]=table['cluster'].str.rsplit('_', expand=True, n=1)
table["event"]='.'
table.rename(columns={"genes":"gene"}, inplace=True)

table=table[["chr","start","end","gene","event","strand","healthy","patient","deltapsi","p.adjust","cluster"]]

# # # SPOT # # #

# # empirical pvalue # #
pvalue=pd.read_csv(args.pvalue,sep="\t")

# # Mahalanobis distance # #
md=pd.read_csv(args.md,sep="\t")

data_merge=pd.merge(left=pvalue,right=md,how='outer',on=['CLUSTER_ID'],suffixes=('_spot.pv','_spot.md'))
data_merge.rename(columns={'CLUSTER_ID':'cluster'}, inplace=True)

spot_leafcutter=table.merge(data_merge, on=["cluster"], how="right")

spot_leafcutter.to_csv("/home/ls/rachelcw/projects/rare_disease/data/leafcutter_spot_results.bed",sep="\t", index=0)
