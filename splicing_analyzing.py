import pandas as pd
import argparse
import glob

FDR = 0.05

# get family name of the analyzed family, and the path to the analyzed family's data
# list of gene names to analyze - txt file
# genes from all three splicing tools 

def get_args():
    parser = argparse.ArgumentParser(description='Analyze splicing data')
    parser.add_argument('--family', type=str, help='Family name')
    parser.add_argument('--path', type=str, help='Path to the analyzed family data')
    parser.add_argument('--genes', type=str, help='List of genes to analyze')
    # parser.add_argument('--as', type=str, help='Path to the analyzed family data')
    # parser.add_argument('--mmsplice', type=str, help='Path to the analyzed family data')
    # parser.add_argument('--spidex', type=str, help='Path to the analyzed family data')
    args = parser.parse_args()
    return args

# def filter_genes_by_list():

def get_leafcutter_results(lc_path, list_of_genes):
    # cluster file
    cluster_path=lc_path + "/leafcutter_ds_cluster_significance.txt"
    cluster_df = pd.read_csv(cluster_path, sep='\t')
    cluster_df[["chr","cluster"]]=cluster_df['cluster'].str.split(':', expand=True)
    cluster_df[['cluster','strand']]=cluster_df['cluster'].str.rsplit('_', expand=True, n=1)
    # effect size file
    effect_size_path=lc_path + "/leafcutter_ds_effect_sizes.txt"
    effect_size_df = pd.read_csv(effect_size_path, sep='\t', usecols=["intron","deltapsi","patient","healthy"])
    effect_size_df[['chr', 'start', 'end','cluster']] = effect_size_df['intron'].str.split(':', expand=True)
    effect_size_df[['cluster','strand']]=effect_size_df['cluster'].str.rsplit('_', expand=True, n=1)
    effect_size_df['junction']=[f'{row["chr"]}:{row["start"]}:{row["end"]}:{row["strand"]}' for index, row in effect_size_df.iterrows() ]
    effect_size_df.drop(['intron'], axis=1, inplace=True)
    lc_results=effect_size_df.merge(cluster_df, on=["chr","cluster","strand"], how="outer")
    # filter the table by the list of genes
    lc_results = lc_results[lc_results['genes'].isin(list_of_genes)]
    # filter the table by FDR < 0.05
    lc_sig_results = lc_results[lc_results['p.adjust'] < 0.05]
    return lc_sig_results

def get_rmats_files(rmats_path):
    mats_files = glob.glob(rmats_path + "/*.MATS.JC.txt")
    for file in mats_files:
        if "A3SS" in file:
            a3ss_df = pd.read_csv(file, sep='\t')
            a3ss_df['junction']=[f'{row["chr"]}:{row["flankingEE"]}:{row["longExonStart_0base"]+1}:{row["strand"]}' for index, row in a3ss_df.iterrows() ]
            a3ss_df['event']="A3SS"
        elif "A5SS" in file:
            a5ss_df = pd.read_csv(file, sep='\t')
            a5ss_df['junction']=[f'{row["chr"]}:{row["longExonEnd"]}:{row["flankingES"]+1}:{row["strand"]}' for index, row in a5ss_df.iterrows() ]
            a5ss_df['event']="A5SS"
        elif "MXE" in file:
            mxe_df = pd.read_csv(file, sep='\t')
            mxe_df['event']="MXE"
        elif "RI" in file:
            ri_df = pd.read_csv(file, sep='\t')
            ri_df['junction']=[f'{row["chr"]}:{row["upstreamEE"]}:{row["downstreamES"]}:{row["strand"]}' for index, row in ri_df.iterrows() ]
            ri_df['event']="RI"
        elif "SE" in file:
            se_df = pd.read_csv(file, sep='\t')
            se_df['junction']=[f'{row["chr"]}:{row["upstreamEE"]}:{row["downstreamES"]+1}:{row["strand"]}' for index, row in se_df.iterrows() ]
            se_df['event']="SE"
    return a3ss_df, a5ss_df, mxe_df, ri_df, se_df

def get_rmats_results(rmats_path,list_of_genes):
    rmats_sig_results = pd.DataFrame()
    a3ss_df, a5ss_df, mxe_df, ri_df, se_df = get_rmats_files(rmats_path)
    #filter by gene list
    rmats_df_list = [df[df['geneSymbol'].isin(list_of_genes)] for df in [a3ss_df, a5ss_df, mxe_df, ri_df, se_df]]
    for df in rmats_df_list:
        if df.empty:
            continue
        else:
            df = df[['GeneID', 'geneSymbol', 'FDR', 'junction', 'event']]
            df = df.rename(columns={'geneSymbol': 'gene'})
            df = df[df['FDR'] < FDR]
            rmats_sig_results = pd.concat([rmats_sig_results, df])
    return rmats_sig_results


def get_majiq_results(majiq_path):
    # majiq_file = f"/home/ls/rachelcw/projects/rare_disease/data/{family}/healthy_patient.deltapsi_filtered_results.tsv"
    majiq_results = pd.read_csv(majiq_path, sep='\t')
    majiq_results.rename(columns={'#Gene Name': 'gene_name'}, inplace=True)
    return majiq_results

def get_common_junctions(lc_results, rmats_results, majiq_results):
    # get the common junction from rmats and leafcutter by the junction column
    common_junctions = pd.merge(lc_results, rmats_results, on='junction', how='inner')
    #check if the junction is found in majiq "junctions" column
    common_junctions["majiq"]=0
    for index, row in common_junctions.iterrows():
        junc_coords = f"{row["junction"].split(":")[1]}-{row["junction"].split(":")[2]}"
        for junction in majiq_results.loc[row["gene"]==majiq_results["gene_name"]]["Junctions coords"]:
            junction_list=junction.split(";")
            if junc_coords in junction_list:
                common_junctions.loc[index, "majiq"]=1
        
