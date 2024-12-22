import pandas as pd
import sys

if __name__ == "__main__":
    """
    Usage: python pathway_map_sample.py <input_file> <output_file> <geneid_file> <genelist_file> [-no-binary]
    writes gene_to_pathways.json and pathways.txt to the project folder
    """
    args = sys.argv[1:]
    input_file = args[0]
    output_file = args[1]
    geneid_file = args[2]
    genelist_file = args[3]
    # project_name = args[0]
    # sample_barcode = args[1]

    binary_output = '-no-binary' in args

    # df = pd.read_csv('../result/'+sample_barcode+'.csv')
    df = pd.read_csv(input_file)
    #df = df.set_index('Tumor_Sample_Barcode')

    # gene_to_id: a dictionary where key is name and value is ensembl id
    # gene_to_id = pd.read_csv('./'+project_name+'/gene_to_id.csv')
    gene_to_id = pd.read_csv(geneid_file)
    gene_to_id = gene_to_id.set_index('Gene Name')
    gene_to_id.dropna(subset=['HGNC Gene ID'], inplace=True)
    gene_to_id = gene_to_id['HGNC Gene ID'].to_dict()
    # make a dict that maps hgnc ids to pathways
    # pathway_file = open('./'+project_name+'/pantherGeneList.txt', 'r')
    pathway_file = open(genelist_file, 'r')
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
    # pathway and gene lists
    pathway_list = list({pathway for pathways in gene_pathway_map.values() for pathway in pathways})
    gene_list = list(df.columns[1:-1])

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
        return pd.concat([row.iloc[:1], pd.Series(pathway_count)], axis = 0)
    # map the pathways
    new_df = df.apply(map_pathways, axis=1)
    # save the new dataframe
    # new_df.to_csv('../result/'+project_name+'/'+sample_barcode+'/pathway_counts.csv', index=False)
    new_df.to_csv(output_file, index=False)
    