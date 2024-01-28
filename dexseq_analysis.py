import pandas as pd
import glob


def argparse():
    import argparse
    parser = argparse.ArgumentParser(description='Dexseq analysis')
    parser.add_argument('-d', '--dir', help='directory path of dexseq counts files')
    parser.add_argument('-o', '--output', help='output dexseq result file')
    parser.add_argument('-g', '--gtf', help='gtf file')
    return parser.parse_args()
    
# read gtf file
def read_gtf_to_dict(gtf_file):
    # Define columns for the GTF file
    columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

    # Read GTF file into a DataFrame
    df = pd.read_table(gtf_file, comment="#", header=None, names=columns, dtype={"seqname": str, "attribute": str})
    #/data01/private/resources/GRCh38_hg38/hg38_gene.gtf

    # Extract gene information
    df_genes = df[df["feature"] == "gene"]

    # Extract gene_id and gene_name from the "attribute" column
    df_genes["gene_id"] = df_genes["attribute"].str.extract(r'gene_id "([^"]+)"')
    df_genes["gene_name"] = df_genes["attribute"].str.extract(r'gene_name "([^"]+)"')

    # Select relevant columns
    df_genes = df_genes[["gene_id", "gene_name"]]
    #convert to dictionary
    dict_genes = df_genes.set_index('gene_id').to_dict()['gene_name']
    return dict_genes


# Read in the DEXSeq results

def read_dexseq_result(dexseq_counts_dir):
    dexseq_result = pd.DataFrame()
    list_of_dexseq_counts = glob.glob(dexseq_counts_dir+'/*.dexseq_counts')
    for file in list_of_dexseq_counts:
        sample_name = file.split('/')[-1].split('.')[0]
        print(sample_name)
        df = pd.read_table(file,sep='\t', header=None, names=['gene_id_exon',sample_name])
        if dexseq_result.empty:
            dexseq_result = df
        else:
            dexseq_result = pd.merge(dexseq_result, df, how='outer', on='gene_id_exon')
    return dexseq_result


def analysis_dexseq_result(dict_genes, dexseq_result, output):
    dexseq_result['gene_id'] = dexseq_result['gene_id_exon'].str.split(':', expand=True)[0]

    # add gene name
    for row in dexseq_result.itertuples():
        if '+' in row.gene_id:
            #more than one gene
            list_gene = row.gene_id.split('+')
            dexseq_result.at[row.Index, 'gene_id'] = list_gene
            list_gene_name = [dict_genes[gene] for gene in list_gene]
            dexseq_result.at[row.Index, 'gene_name'] = ','.join(list_gene_name)
        elif '_' in row.gene_id:
            #gene not found in dexseq
            dexseq_result.at[row.Index, 'gene_name'] = pd.NA
        else:
            dexseq_result.at[row.Index, 'gene_name'] = dict_genes[row.gene_id]
    # remove [] from gene_id and gene_name when saved to tsv
    dexseq_result['gene_id'] = dexseq_result['gene_id'].astype(str).str.replace(r'\[|\]', '')
    dexseq_result['gene_name'] = dexseq_result['gene_name'].astype(str).str.replace(r'\[|\]', '')
    # save to tsv file
    dexseq_result.to_csv(output ,sep='\t', index=False)
    # dexseq_result.to_csv('/home/ls/rachelcw/projects/rare_disease/data/dexseq/dexseq_result.tsv',sep='\t', index=False)

if __name__ == '__main__':
    args=argparse()
    dict_genes = read_gtf_to_dict(args.gtf)
    dexseq_result = read_dexseq_result(args.dir)
    analysis_dexseq_result(dict_genes, dexseq_result, args.output)