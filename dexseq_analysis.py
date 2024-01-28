import pandas as pd
import glob

# read gtf file
# Define columns for the GTF file
columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

# Read GTF file into a DataFrame
df = pd.read_table("/data01/private/resources/GRCh38_hg38/hg38_gene.gtf", comment="#", header=None, names=columns, dtype={"seqname": str, "attribute": str})

# Extract gene information
df_genes = df[df["feature"] == "gene"]

# Extract gene_id and gene_name from the "attribute" column
df_genes["gene_id"] = df_genes["attribute"].str.extract(r'gene_id "([^"]+)"')
df_genes["gene_name"] = df_genes["attribute"].str.extract(r'gene_name "([^"]+)"')

# Select relevant columns
df_genes = df_genes[["gene_id", "gene_name"]]
#convert to dictionary
dict_genes = df_genes.set_index('gene_id').to_dict()['gene_name']

# Read in the DEXSeq results
dexseq_result = pd.DataFrame()
list_of_dexseq_counts = glob.glob('/home/ls/rachelcw/projects/rare_disease/data/dexseq/*.dexseq_counts')
for file in list_of_dexseq_counts:
    sample_name = file.split('/')[-1].split('.')[0]
    print(sample_name)
    df = pd.read_table(file,sep='\t', header=None, names=['gene_id_exon',sample_name])
    if dexseq_result.empty:
        dexseq_result = df
    else:
        dexseq_result = pd.merge(dexseq_result, df, how='outer', on='gene_id_exon')

dexseq_result['gene_id'] = dexseq_result['gene_id_exon'].str.split(':', expand=True)[0]

# add gene name
for row in dexseq_result.itertuples():
    if '+' in row.gene_id:
        #more than one gene
        list_gene = row.gene_id.split('+')
        dexseq_result.at[row.Index, 'gene_name'] = [dict_genes[gene] for gene in list_gene]
    elif '_' in row.gene_id:
        #gene not found in dexseq
        dexseq_result.at[row.Index, 'gene_name'] = pd.NA
    else:
        dexseq_result.at[row.Index, 'gene_name'] = dict_genes[row.gene_id]
