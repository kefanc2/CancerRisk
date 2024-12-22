import pandas as pd
import numpy as np
import sys


def vital_status(status):
    if status == 'Alive':
        return 0
    return 1

def parse_month(s):
    try:
        s['months'] = max(float(s['months'])/30, 0)
    except:
        s['months'] = 0
    return s

def trim_barcode(s):
    s['Tumor_Sample_Barcode'] = s['Tumor_Sample_Barcode'][:12]
    return s

def join_df(mutations_path, clinical_path, out_path, verbose=False):
    df = pd.read_csv(mutations_path, index_col = 0) # mutations data
    df.index.name = 'Tumor_Sample_Barcode'
    try:
        clinical = pd.read_csv(clinical_path) # clinical data
    except:
        print('Encoding error, trying windows-1254')
        clinical = pd.read_csv(clinical_path,encoding='windows-1254')
    clinical['days_to_death'].replace('[Not Applicable]', np.nan, inplace=True)
    clinical['months'] = clinical['days_to_death'].fillna(clinical['days_to_last_followup'])
    clinical['status'] = clinical['vital_status'].apply(vital_status)
    clinical = clinical[['Tumor_Sample_Barcode', 'status', 'months']]
    clinical.set_index('Tumor_Sample_Barcode', inplace=True)
    clinical = clinical.apply(parse_month, axis=1)
    result = clinical.join(df)
    print(f'Joined data has shape {result.shape}')
    result.to_csv(out_path)


if __name__ == "__main__":
    """
    Usage: python join_clinical.py mutations_path clinical_path out_path
    Alternative: python join_clinical.py project_name
    """
    args = sys.argv[1:]
    if len(args) >= 3:
        join_df(args[0], args[1], args[2])
    else:
        mutations_path = './'+args[0]+'/'+'mutations.csv'
        clinical_path = './'+args[0]+'/'+'clinical.csv'
        out_path = './'+args[0]+'/'+'joined.csv'
        join_df(mutations_path, clinical_path, out_path)
        counts_path = './'+args[0]+'/'+'counts.csv'
        out_path = './'+args[0]+'/'+'joined_counts.csv'
        join_df(counts_path, clinical_path, out_path)
    