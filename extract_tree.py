import json
import numpy as np
import sys
from pathlib import Path

def read_json(path):
    f = open(path)
    data = json.load(f)
    f.close()
    return data

def write_matrix(barcode, output_folder):
    # json_path = 'trees/'+ run_name + '.summ.json'
    json_path = f'{output_folder}/{barcode}.summ.json'
    j_summs = read_json(json_path)
    trees = j_summs['trees']
    llh_list = [(trees[i]['llh'],i) for i in trees]
    llh_list.sort()
    i = 0
    # FINDING THE FIRST TREE THAT HAS A FILE
    # SORTED IN ASCENDING ORDER OF LLH
    while True:
        best_idx = llh_list[i][1]
        # file = Path("trees/"+run_name+".mutass/"+str(best_idx)+".json")
        file = Path(f'{output_folder}/{barcode}.mutass/{best_idx}.json')
        if file.is_file():
            break
        i += 1
    # finding the best tree
    best_tree = trees[best_idx]
    structure = best_tree['structure']
    # structure: dictionary: example {'0': [1], '1': [2, 4], '2': [3]}
    l = len(best_tree['populations'])
    result_matrix = np.zeros((l-1,l-1))
    for i in range(1, len(structure)):
        parent = int(list(structure.keys())[i])
        children = structure[str(parent)]
        for child in children:
            result_matrix[child-1][parent-1] = 1
            result_matrix[child-1] += result_matrix[parent-1]
    result_matrix = np.array(result_matrix, dtype=int)
    # fout = open('trees/' + run_name + '.matrix.txt', "w")
    fout = open(f'{output_folder}/{barcode}.matrix.txt', "w")
    for i in result_matrix:
        fout.write(' '.join(list(map(str, i))) + '\n')
    fout.close()
    # fout = open('trees/' + run_name + '.dot', "w")
    fout = open(f'{output_folder}/{barcode}.dot', "w")
    fout.write('digraph G {\n')
    for i in range(len(structure)):
        if str(i) in structure:
            for j in structure[str(i)]:
                fout.write(str(i) + ' -> ' + str(j) + " [label=" + str(best_tree['populations'][str(j)]['num_ssms']) + '];\n')
    fout.write('}')
    fout.close()
    return best_idx

if __name__ == "__main__":
    """
    usage: python tree_parser.py <barcodr> <output folder>
    """
    args = sys.argv[1:]
    barcode = args[0]
    output_folder = args[1]
    write_matrix(barcode, output_folder)