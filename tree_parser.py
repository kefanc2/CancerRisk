
# Python program to read
# json file
 
import json
import numpy as np
import pandas as pd
from extract_tree import *
 
"""
input:
    1. json file: tree example:
    {"dataset_name": "run_name", "mut_assignments": {"1": {"ssms": ["s0"], "cnvs": []}, "2": {"ssms": ["s2"], "cnvs": []}}}
    2. ssm file: ssm_id, gene_name, ref, total
    3. tree matrix: binary matrix, M[i,j]=1 if j is an ancestor of i
"""

def read_matrix(path):
    mfile = open(path, "r")
    m_array = []
    for line in mfile:
        m_array.append(list(map(lambda x: int(x), line.strip().split(" "))))
    mfile.close()
    M = np.array(m_array)
    mfile.close()
    return M

def read_ssm(path):
    ssm_in = open(path, "r")
    # skip header
    ssm_in.readline()
    # maps identifier to gene name
    gene_list_full = []
    ssm_dict = {}
    for line in ssm_in:
        fields = line.strip().split("\t")
        gene_name = fields[1]
        ssm_id = fields[0]
        ssm_dict[ssm_id] = gene_name
        gene_list_full.append(gene_name)
    ssm_in.close()
    return ssm_dict, gene_list_full

def read_json(path):
    f = open(path)
    data = json.load(f)
    tree = data['mut_assignments']
    f.close()
    return tree

def id2name(ssm_dict, ids):
    return list(map(lambda x: ssm_dict[x], ids))

def print_genes(result_dict):
    l = []
    for k in result_dict:
        l += result_dict[k]
    print(set(l))

def tree_to_matrix(barcode, output_folder, genes, ssm_path):
    # Opening JSON file
    idx = write_matrix(barcode, output_folder)
    # json_path = 'trees/'+ run_name + '.mutass/' + str(idx) + '.json'
    json_path = f'{output_folder}/{barcode}.mutass/{idx}.json'
    

    # binary matrix, M[i,j]=1 if j is an ancestor of i
    M = read_matrix(f'{output_folder}/{barcode}.matrix.txt')
    ssm_dict, gene_list_full = read_ssm(ssm_path)
    tree = read_json(json_path)
    nodes = list(tree.keys())
    # new emerging ssm in each population
    result_dict = dict.fromkeys(tree.keys(), [])
    for i in range(len(M)):
        nodei = nodes[i]
        for j in range(len(M)):
            nodej = nodes[j]
            if M[i,j] == 1 or i == j:
                result_dict[nodei] = result_dict[nodei] + id2name(ssm_dict, tree[nodej]['ssms']) 

    # print(result_dict)
    mut_matrix = np.zeros((len(result_dict), len(genes)), dtype=int)
    for i, population in enumerate(result_dict):
        # i: index, population: name of key
        mut_list = result_dict[population]
        for j in range(len(genes)):
            if genes[j] in mut_list:
                mut_matrix[i, j] = 1

    #add name of population
    output_df = pd.DataFrame(mut_matrix, columns=genes)
    output_df['Population'] = result_dict.keys()
    output_df.set_index('Population', inplace=True)
    time_stamp = dict([(population, len(result_dict[population])) for population in result_dict])

    output_df['time_stamp'] = output_df.apply(lambda x: time_stamp[x.name], axis=1)
    output_df.to_csv(f'{output_folder}/{barcode}.csv', index=True)

if __name__ == "__main__":
    args = sys.argv[1:]
    project_name = args[0]
    barcode = args[1]
    output_folder = args[2]
    gene_path = args[3]#'project_data/' + project_name + '/gene_list.txt'
    genes = open(gene_path).read().split(' ')
    ssm_path = args[4]#'result/' + barcode + '.txt'
    # tree_matrix = args[5]#f'project_data/{project_name}/trees/{barcode}.matrix.txt'
    tree_to_matrix(barcode, output_folder, genes, ssm_path)

            