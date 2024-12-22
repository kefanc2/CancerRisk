import sys
import pandas as pd

if __name__ == "__main__":
    args = sys.argv[1:]
    project_name = args[0]
    df = pd.read_csv('../project_data/'+project_name+'/joined.csv')
    df = df.set_index('Tumor_Sample_Barcode')
    genes = list(df.columns[3:])
    with open('../project_data/'+project_name+'/gene_list.txt', 'w') as f:
        f.write(' '.join(genes) + '\n')