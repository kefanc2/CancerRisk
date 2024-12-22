import os
import sys
import pandas as pd

def status_to_bool(status):
    if status == 'Alive':
        return 0
    elif status == 'Dead':
        return 1
    else:
        return None

def get_status_clinical(df, sample_barcode):
    try:
        return status_to_bool(df.loc[sample_barcode]['vital_status'])
    except Exception as e:
        print(e)
        return None

def get_status(df, sample_barcode):
    try:
        df['Tumor_Sample_Barcode'] = [x[:12] for x in df['Tumor_Sample_Barcode']]
        return [df.loc[df['Tumor_Sample_Barcode'] == sample_barcode]['status'].values[0], df.loc[df['Tumor_Sample_Barcode'] == sample_barcode]['months'].values[0]]
    except Exception as e:
        return None

def get_censor(clinical, sample_barcode):
    try:
        clinical['Tumor_Sample_Barcode'] = [x[:12] for x in clinical['Tumor_Sample_Barcode']]
        dtd = clinical.loc[clinical['Tumor_Sample_Barcode'] == sample_barcode]['days_to_death'].values[0]
        if dtd == '[Not Applicable]':
            return 1
        else:
            return 0
    except Exception as e:
        return 1

if __name__ == "__main__":
    args = sys.argv[1:]
    project_name = args[0]
    sample_barcode = args[1]

    df = pd.read_csv('../project_data/'+project_name+'/clinical.csv')
    df.set_index('Tumor_Sample_Barcode_min', inplace=True)
    print(get_censor(df, sample_barcode))