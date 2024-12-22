import os
import sys
import subprocess
from tqdm import tqdm
from pathlib import Path

def run_project(project_name, resume=False):
    list_dir = list(map(lambda x: x[:-4],os.listdir('./data/' + project_name)))
    list_result = os.listdir('./result/' + project_name)
    for sample_barcode in tqdm(list_dir):
        if resume and sample_barcode in list_result:
            continue
        maf_path = 'data/'+project_name+'/'+sample_barcode+'.maf'
        file = Path(maf_path)
        if file.is_file():
            s = os.path.getsize(maf_path)
            if s > 100000 and s < 2000000:
                subprocess.call(['bash', './infer_sample.sh', project_name, sample_barcode])

if __name__ == "__main__":
    """
    Usage: python main.py project_name
    or:
    Usage: python main.py project_name (-remap) (-resume)
    """
    args = sys.argv[1:]
    project_name = args[0]
    resume = '-resume' in args
    if project_name not in os.listdir('./data'):
        print('Project not found')
        sys.exit()
    if project_name not in os.listdir('./result'):
        os.mkdir('./result/' + project_name)
    if '-remap' not in args:
        run_project(project_name, resume=resume)
    elif '-remap' in args:
        list_dir = list(filter(lambda x: len(x) == 12,os.listdir('./result/' + project_name)))
        for sample_barcode in tqdm(list_dir):
            subprocess.call(['bash', './remap.sh', project_name, sample_barcode])
    