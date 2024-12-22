import pandas as pd
import numpy as np
from sksurv.preprocessing import OneHotEncoder
import joblib
from get_status import get_status
import os
import random
import utils
from tqdm import tqdm

class Analyzer():
    def __init__(self, project_name, test=False):
        self.project_name = project_name
        self.df = pd.read_csv('../pathway_map/'+project_name+'/pathway_counts.csv')
        X = self.df.iloc[:,3:] # select columns with features (currently gene info)
        surv_info = self.df[['status','months']]  # columns with survival info
        self.y = np.array(list(zip(surv_info.status, surv_info.months)) , dtype =[('Status', '?'), ('Survival_in_months', '<f8')])
        self.Xt = OneHotEncoder().fit_transform(X) # one hot encoder w/ categorical variables
        
        model_path = f'rf/{project_name}_rf.pkl'
        self.model = joblib.load(model_path)
        
        self.cached_trees = []
        self.cached_statuses = []
        self.test = test
        
    def process(self):
        sample_list = os.listdir('../result/'+self.project_name)
        
        if self.test: # for testing purposes
            random.shuffle(sample_list)
            sample_list = sample_list[:10]
            
        for sample_barcode in tqdm(sample_list):
            if not utils.is_valid_barcode(sample_barcode):
                continue
            status, populations = self.process_sample(sample_barcode)
            yield status, populations
        
    def process_sample(self, sample_barcode):
        """
        Args:
            sample_barcode (str): sample barcode

        Returns:
            status: [status(0/1), months of survival]
            populations: tree data in a list of all paths to leaf nodes, each node is [num_mutations, score, clone id]
        """
        phylo_pathway_data = pd.read_csv('../result/'+self.project_name+'/'+sample_barcode+'/pathway_counts.csv')
        phylo_data = pd.read_csv('../result/'+self.project_name+'/'+sample_barcode+'/'+sample_barcode+'.csv')
        X_phylo = phylo_pathway_data.iloc[:,1:] # select columns with features (currently gene info)
        new_order = self.Xt.columns
        X_phylo = X_phylo[new_order]
        Xt_phylo = OneHotEncoder().fit_transform(X_phylo) # one hot encoder w/ categorical variables

        scores_phylo = self.model.predict(Xt_phylo)

        phylo_result = list(zip(phylo_data['time_stamp'], scores_phylo))

        m = open('../result/'+self.project_name+'/'+sample_barcode+'/'+sample_barcode+'.matrix.txt')
        M = []
        while line := m.readline():
            M.append(list(map(int,line.strip().split(' '))))
        M = np.array(M)
        xy = [[entry[0],entry[1]] for entry in phylo_result]
        leaves = np.where(M.sum(axis=0)==0)[0]
        num_leaves = len(leaves)
        populations = [[] for i in range(num_leaves)]
        
        for i in range(num_leaves):
            for j in range(M.shape[0]):
                if M[leaves[i],j] == 1 or leaves[i] == j:
                    populations[i].append([xy[j][0], xy[j][1], j])
        populations[i].sort()

        status = get_status(self.df, sample_barcode)
        self.cached_trees.append(populations)
        self.cached_statuses.append(status) #[status, months]
        return status, populations
    
    def analyze(self):
        return self.cached_statuses, self.cached_trees
