from analyzer import Analyzer
import pandas as pd
import numpy as np
from sksurv.preprocessing import OneHotEncoder
import joblib
from get_status import get_status
import os
import random
import utils
from tqdm import tqdm

def get_max_by_time(tree, num_timepoints):
    # tree: a list of paths, each path is a list of nodes [num_mutations, score, clone id]
    max_t = max([tree[i][-1][0] for i in range(len(tree))])
    out = [tree[0][0][1]] * num_timepoints
    max_clones = [set()] * num_timepoints
    for path in tree:
        for node in path[::-1]:
            t = min(node[0] * num_timepoints // max_t, num_timepoints-1)
            score = node[1]
            for i in range(int(t), num_timepoints):
                if score > out[i]:
                    out[i] = score
                    max_clones[i] = set([node[2]])
                elif score == out[i]:
                    max_clones[i].add(node[2])
                else:
                    break
    return out, max_clones

class EightPlotTimeline(Analyzer):
    def __init__(self, project_name, test, num_timepoints=100):
        super().__init__(project_name, test)
        self.num_timepoints = num_timepoints
        
    def process(self):
        sample_list = os.listdir('../result/'+self.project_name)
        
        if self.test: # for testing purposes
            random.shuffle(sample_list)
            sample_list = sample_list[:10]
            
        for sample_barcode in tqdm(sample_list):
            if not utils.is_valid_barcode(sample_barcode):
                continue
            status, populations, max_pathways, max_mutations, initial_score, max_score = self.process_sample(sample_barcode)
            yield status, populations, max_pathways, max_mutations, initial_score, max_score

    def process_sample(self, sample_barcode):
        """
        Args:
            sample_barcode (str): sample barcode

        Returns:
            status: [status(0/1), months of survival]
            populations: tree data in a list of all paths to leaf nodes, each node is [num_mutations, score, clone id]
        """
        status = get_status(self.df, sample_barcode)
        if not status:
            return None, None, None, None, None, None
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
        
        max_pathways = list(X_phylo.columns[X_phylo.loc[np.argmax(scores_phylo)] >= 1])
        max_mutations = list(phylo_data.columns[phylo_data.loc[np.argmax(scores_phylo)] >= 1])
        if 'Population' in max_mutations:
            max_mutations.remove('Population')
        if 'Population' in max_pathways:
            max_pathways.remove('Population')
        if 'time_stamp' in max_mutations:
            max_mutations.remove('time_stamp')
        
        for i in range(num_leaves):
            for j in range(M.shape[0]):
                if M[leaves[i],j] == 1 or leaves[i] == j:
                    populations[i].append([xy[j][0], xy[j][1], j])
        populations[i].sort()

        self.cached_trees.append(populations)
        self.cached_statuses.append(status) #[status, months]
        initial_score = populations[0][0][1]
        max_score = max([max([node[1] for node in population]) for population in populations])
        return status, populations, max_pathways, max_mutations, initial_score, max_score
    
    