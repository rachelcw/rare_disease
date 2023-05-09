import sys
import pandas as pd
import pyensembl

args=sys.argv

# # cluster_significance # #
cluster=pd.read_csv(args[1],sep="\t",usecols=["cluster","p.adjust"])
cluster[["chr","cluster"]]=cluster['cluster'].str.split(':', expand=True)

# # effect sizes # #
ds=pd.read_csv(args[2],sep="\t",usecols=["intron","deltapsi","patient","healthy"])
ds[['chr', 'start', 'end','cluster']] = ds['intron'].str.split(':', expand=True)
ds.drop(['intron'], axis=1, inplace=True)

table=ds.merge(cluster, on=["chr","cluster"], how="left")
table[['cluster','strand']]=table['cluster'].str.rsplit('_', expand=True, n=1)

#   #   #   #   #   #   #   #
data = pyensembl.Genome(
    reference_name="GRCh38",
    annotation_name="GRCh38",
    gtf_path_or_url=args[3])
data.index()

table["gene"]=pd.NA

# how to relate to mitochondrial chromosome
for index, row in table.iterrows():
    gene_name = data.gene_names_at_locus(contig=row["chr"], position=int(row["start"]),end=int(row["end"]), strand=row["strand"])
    if gene_name != []:
        table.at[index,"gene"]=','.join(gene_name)

#   #   #   #   #   #   #   #

table[["chr","start","end","gene","strand","healthy","patient","deltapsi","p.adjust"]].to_csv("/home/ls/rachelcw/projects/rare_disease/data/leafcutter_results.bed",sep="\t", index=0)