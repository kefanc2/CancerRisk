import os
import sys
from tqdm import tqdm

if __name__ == "__main__":
    project_name = sys.argv[1]
    fname = f'cohort_data/{project_name}.maf'

    with open(fname, 'r') as f:
        header = f.readline()
        barcode_index = header.split('\t').index('Tumor_Sample_Barcode')
        if project_name not in os.listdir(f'.'):
            os.mkdir(project_name)
        currently_writing = ''
        while True:
            line = f.readline()
            if not line:
                break
            barcode = line.split('\t')[barcode_index][:12]
            if barcode != currently_writing:
                if currently_writing:
                    g.close()
                currently_writing = barcode
                g = open(f'{project_name}/{barcode}.maf', 'w')
                g.write(header)
            g.write(line)
