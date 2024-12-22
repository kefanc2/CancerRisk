import pandas as pd
import sys
import mygene


if __name__ == "__main__":
    """
    Usage: python save_gene_id.py <project_name>
    writes gene_to_id.csv and gene_ids.txt to the project folder
    """
    # Initialize the MyGeneInfo service
    mg = mygene.MyGeneInfo()
    # List of gene names to map
    args = sys.argv[1:]
    project_name = args[0]
    
    binary_input = '-count' not in args
    # choose the input file
    if binary_input:
        df = pd.read_csv('../project_data/'+project_name+'/joined.csv')
    else:
        df = pd.read_csv('../project_data/'+project_name+'/joined_counts.csv')
    df = df.set_index('Tumor_Sample_Barcode')
    # List of gene names to map
    gene_names = list(df.columns[2:])

    # Query MyGeneInfo to get human gene identifiers
    gene_info = mg.querymany(gene_names, scopes='symbol', fields='entrezgene,ensembl.gene,HGNC', species='human')
    
    # Convert the result to a DataFrame for easier viewing
    df = pd.DataFrame(gene_info)

    # Select relevant columns and rename them for clarity
    df = df[['query', 'entrezgene', 'ensembl', 'HGNC']].rename(columns={
        'query': 'Gene Name',
        'entrezgene': 'NCBI Gene ID',
        'ensembl': 'Ensembl Gene ID',
        'HGNC': 'HGNC Gene ID'
    })
    df['Ensembl Gene ID'] = df['Ensembl Gene ID'].apply(lambda x: x['gene'] if isinstance(x, dict) else None)

    # Drop duplicates and handle missing values
    df = df.drop_duplicates(subset=['Gene Name']).reset_index(drop=True)
    df.dropna(subset=['Ensembl Gene ID'], inplace=True)

    with open('./'+project_name+'/gene_ids.txt', 'w') as f:
        f.write(','.join(df['Ensembl Gene ID'].tolist()))

    df.to_csv('./'+project_name+'/gene_to_id.csv', index=False)

    

    

    

        
    
    

    

    