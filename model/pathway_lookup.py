import pandas as pd
import os

cwd = os.getcwd()
project_name = 'TCGA-LUAD'
pathway_file = open(f'{cwd}/../pathway_map/{project_name}/pantherGeneList.txt', 'r')
id_to_gene= {}
df = pd.read_csv(f'{cwd}/../pathway_map/{project_name}/gene_to_id.csv')
for i in range(len(df)):
    id_to_gene[df.iloc[i,3]] = df.iloc[i,0]

pathway_dict = {}
pathway_to_ids = {}
while True:
    line = pathway_file.readline()
    if not line:
        break
    elif len(line) < 2:
        continue
    else:
        line = line.strip().split()
        gene_list = [gene.split('|')[1].split('=')[1] for gene in line[1].split(',')]
        pathway_to_ids[line[0]] = [id_to_gene[int(gene)] for gene in gene_list if int(gene) in id_to_gene]
        pathway_dict[line[0]] = ' '.join(line[2:-3])

while True:
    line = input('type pathway id:\n>>>')
    if not line:
        break
    else:
        try:
            if line[0] != 'P':
                line = 'P' + '0'*(5-len(line)) + line
            print()
            print(pathway_dict[line])
            print(' '.join(pathway_to_ids[line]))
            print()
        except:
            print('Pathway not found')
            print()


