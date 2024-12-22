import pandas as pd
import sys

if __name__ == "__main__":
    """
    Usage: python pathway_map_project.py <project_name>
    writes gene_to_pathways.json and pathways.txt to the project folder
    """
    args = sys.argv[1:]
    project_name = args[0]
    
    binary_input = '-count' not in args
    binary_output = '-no-binary' not in args
    # choose the input file
    if binary_input:
        df = pd.read_csv('../project_data/'+project_name+'/joined.csv')
    else:
        df = pd.read_csv('../project_data/'+project_name+'/joined_counts.csv')
    #df = df.set_index('Tumor_Sample_Barcode')

    # gene_to_id: a dictionary where key is name and value is ensembl id
    gene_to_id = pd.read_csv('./'+project_name+'/gene_to_id.csv')
    gene_to_id = gene_to_id.set_index('Gene Name')
    gene_to_id.dropna(subset=['HGNC Gene ID'], inplace=True)
    gene_to_id = gene_to_id['HGNC Gene ID'].to_dict()
    # make a dict that maps hgnc ids to pathways
    pathway_file = open('./'+project_name+'/pantherGeneList.txt', 'r')
    hgnc_id_to_pathway = {}
    while True:
        line = pathway_file.readline()
        if not line:
            break
        elif len(line) < 2:
            continue
        else:
            # sample line: 
            # P02755	HUMAN|HGNC=6936|UniProtKB=Q96RQ3	Methylmalonyl pathway	3	10	64
            line = line.strip().split()
            pathway_id = line[0]
            pathway_genes_list = line[1].split(',')
            hgnc_list = [token.split('|')[1].split('=')[1] for token in pathway_genes_list]
            for hgnc_id in hgnc_list:
                if hgnc_id in hgnc_id_to_pathway:
                    hgnc_id_to_pathway[hgnc_id].append(pathway_id)
                else:
                    hgnc_id_to_pathway[hgnc_id] = [pathway_id]
    # get the gene to pathway map by making a dictionary
    gene_pathway_map = {gene: hgnc_id_to_pathway[str(int(gene_to_id[gene]))] for gene in gene_to_id if str(int(gene_to_id[gene])) in hgnc_id_to_pathway}
    # remove low frequency genes
    low_freq_genes = [col for col in df.columns[3:] if df[col].sum() < 10]
    df = df.drop(low_freq_genes, axis = 1, errors='ignore')
    # pathway and gene lists
    pathway_list = list({pathway for pathways in gene_pathway_map.values() for pathway in pathways})
    gene_list = list(df.columns[1:])
    # map each patient's genes to pathways
    def map_pathways(row):
        pathway_count = {pathway: 0 for pathway in pathway_list}
        for gene in gene_list:
            if gene in gene_pathway_map:
                pathways = gene_pathway_map[gene]
                for pathway in pathways:
                    if binary_output:
                        pathway_count[pathway] = min(1, pathway_count[pathway] + row[gene])
                    else:
                        pathway_count[pathway] += row[gene]
        return pd.concat([row.iloc[:3], pd.Series(pathway_count)], axis = 0)
    # map the pathways
    new_df = df.apply(map_pathways, axis=1)
    new_df = new_df[new_df['months'] > 0]
    # save the new dataframe
    new_df.to_csv('./'+project_name+'/pathway_counts.csv', index=False)
    