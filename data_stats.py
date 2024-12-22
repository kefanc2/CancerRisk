import os
import sys
import subprocess
from tqdm import tqdm
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np

def show_project(project_name):
    list_dir = list(map(lambda x: x[:-4],os.listdir('./data/' + project_name)))
    list_result = os.listdir('./result/' + project_name)
    result = []
    for sample_barcode in tqdm(list_dir):
        maf_path = 'data/'+project_name+'/'+sample_barcode+'.maf'
        file = Path(maf_path)
        if file.is_file():
            s = os.path.getsize(maf_path)
            result.append(s)
    min_value = min(result)
    max_value = max(result)
    num_bins = 100
    log_space = np.logspace(np.log10(min_value), np.log10(max_value), num=num_bins)
    plt.vlines(100000, 0, 30, colors='r', linestyles='dashed')
    plt.vlines(2000000, 0, 30, colors='r', linestyles='dashed')
    plt.xscale('log')
    plt.hist(result, bins=log_space)
    plt.title('Histogram of file sizes')
    plt.xlabel('File size')
    plt.ylabel('Frequency')
    plt.show()
    
            
                

if __name__ == "__main__":
    """
    Usage: python main.py project_name
    or:
    Usage: python main.py project_name (-remap) (-resume)
    """
    args = sys.argv[1:]
    project_name = args[0]
    show_project(project_name)
    