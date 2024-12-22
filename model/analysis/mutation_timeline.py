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

def add_to_dict(d, key):
    if isinstance(key, list) or isinstance(key, set):
        for k in key:
            add_to_dict(d, k)
    else:
        if key in d:
            d[key] += 1
        else:
            d[key] = 1

def get_max_by_time(tree, num_timepoints):
    # tree: a list of paths, each path is a list of nodes [num_mutations, score, clone id]
    max_t = max([tree[i][-1][0] for i in range(len(tree))])
    out = [float('-inf')] * num_timepoints # initialize to negative infinity
    max_clones = [set()] * num_timepoints
    for path in tree:
        for node in path[::-1]:
            # t = min(node[0] * num_timepoints // max_t, num_timepoints-1)
            # print(node[0], num_timepoints, max_t, node[0] * num_timepoints // max_t)
            t = min(node[0] * num_timepoints // max_t, num_timepoints-1) # timing from the first node
            score = node[1]
            for i in range(int(t), num_timepoints):
                o = max_clones[0]
                if score > out[i]:
                    out[i] = score
                    max_clones[i] = set([node[2]])
                elif score == out[i]:
                    max_clones[i].add(node[2])
                else:
                    break
    return out, max_clones

class MutationTimeline(Analyzer):
    def __init__(self, project_name, test, num_timepoints=100):
        super().__init__(project_name, test)
        self.num_timepoints = num_timepoints
        self.mutation_timeline_counts_alive = [{} for i in range(self.num_timepoints)]
        self.mutation_timeline_counts_dead = [{} for i in range(self.num_timepoints)]
        self.pathway_timeline_counts_alive = [{} for i in range(self.num_timepoints)]
        self.pathway_timeline_counts_dead = [{} for i in range(self.num_timepoints)]
        
    def analyze(self):
        statuses = [status[0] for status in self.cached_statuses]
        return self.mutation_timeline_counts_alive, self.mutation_timeline_counts_dead, \
            self.pathway_timeline_counts_alive, self.pathway_timeline_counts_dead, len(statuses) - sum(statuses), sum(statuses)

    def process(self):
        sample_list = os.listdir('../result/'+self.project_name)
        
        if self.test: # for testing purposes
            random.shuffle(sample_list)
            sample_list = sample_list[:10]
            
        for sample_barcode in tqdm(sample_list):
            if not utils.is_valid_barcode(sample_barcode):
                continue
            status, populations, mutation_timeline, pathway_timeline = self.process_sample(sample_barcode)
            if status == None:
                continue
            for i in range(self.num_timepoints):
                if status[0] == 0:
                    add_to_dict(self.mutation_timeline_counts_alive[i], mutation_timeline[i])
                    add_to_dict(self.pathway_timeline_counts_alive[i], pathway_timeline[i])
                else:
                    add_to_dict(self.mutation_timeline_counts_dead[i], mutation_timeline[i])
                    add_to_dict(self.pathway_timeline_counts_dead[i], pathway_timeline[i])
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
        if status is not None:
            self.cached_trees.append(populations)
            self.cached_statuses.append(status) #[status, months]
        scores, clones = get_max_by_time(populations, self.num_timepoints)
        
        mutation_timeline = [set() for i in range(self.num_timepoints)]
        pathway_timeline = [set() for i in range(self.num_timepoints)]
        
        # make a dictionary from node to score
        node_to_score = {}
        node_to_time = {}
        for path in populations:
            for node in path:
                node_to_score[node[2]] = node[1]
                node_to_time[node[2]] = node[0]
        
        # make a timeline of mutations and pathways
        for i, clones_i in enumerate(clones):
            for clone in clones_i:
                pathways = list(X_phylo.columns[X_phylo.iloc[clone] >= 1])
                mutations = list(phylo_data.columns[phylo_data.iloc[clone] >= 1])
                if 'Population' in mutations:
                    mutations.remove('Population')
                if 'Population' in pathways:
                    pathways.remove('Population')
                if 'time_stamp' in mutations:
                    mutations.remove('time_stamp')
                if 'time_stamp' in pathways:
                    pathways.remove('time_stamp')
                mutation_timeline[i].update(mutations)
                pathway_timeline[i].update(pathways)
                
        return status, populations, mutation_timeline, pathway_timeline
            
            
