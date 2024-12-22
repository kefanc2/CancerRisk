import sys
import joblib
import pandas as pd
import numpy as np
from sksurv.preprocessing import OneHotEncoder
import matplotlib.pyplot as plt
import os
from get_status import get_status
import random
from tqdm import tqdm
from scipy.stats import ttest_ind

def get_max_by_time(scores, num_timepoints):
    max_t = max([scores[i][-1][0] for i in range(len(scores))])
    indices = [0] * len(scores)
    result = []
    for n in range(num_timepoints + 1):
        curr_t = int(max_t / num_timepoints * n)
        max_score = float('-inf')
        for i in range(len(scores)):
            while indices[i] < len(scores[i]) - 1 and scores[i][indices[i]][0] < curr_t:
                indices[i] += 1
            if indices[i] < len(scores[i]) and scores[i][indices[i]][1] > max_score:
                max_score = scores[i][indices[i]][1]
            elif indices[i] > len(scores[i]) - 1:
                print('sss')
        result.append(max_score)
    return result
    

def process_sample(sample_barcode, project_name, df, model, num_timepoints):
    phylo_pathway_data = pd.read_csv('../result/'+project_name+'/'+sample_barcode+'/pathway_counts.csv')
    phylo_data = pd.read_csv('../result/'+project_name+'/'+sample_barcode+'/'+sample_barcode+'.csv')

    X_phylo = phylo_pathway_data.iloc[:,1:] # select columns with features (currently gene info)
    new_order = Xt.columns
    X_phylo = X_phylo[new_order]
    Xt_phylo = OneHotEncoder().fit_transform(X_phylo) # one hot encoder w/ categorical variables

    scores_phylo = model.predict(Xt_phylo)
    # scores_phylo = coxnet_pred.predict(Xt_phylo)

    phylo_result = list(zip(phylo_data['time_stamp'], scores_phylo))

    m = open('../result/'+project_name+'/'+sample_barcode+'/'+sample_barcode+'.matrix.txt')
    M = []
    while line := m.readline():
        M.append(list(map(int,line.strip().split(' '))))
    M = np.array(M)
    xy = [[entry[0],entry[1]] for entry in phylo_result]
    num_population = M.shape[0]
    populations = [[] for i in range(num_population)]
    max_pathways = list(X_phylo.columns[X_phylo.loc[np.argmax(scores_phylo)] == 1])
    max_mutations = list(phylo_data.columns[phylo_data.loc[np.argmax(scores_phylo)] == 1])
    if 'Population' in max_mutations:
        max_mutations.remove('Population')
    if 'Population' in max_pathways:
        max_pathways.remove('Population')

    for i in range(num_population):
        for j in range(num_population):
            if M[i,j] == 1 or i == j:
                populations[i].append(xy[j])
    populations[i].sort()
    # populations = np.array(populations)
    # populations = populations[M.sum(axis=0) == 0]

    status = get_status(df, sample_barcode)
    if status:
        month = status[1]
        status = status[0]
        initial_score = populations[0][0][1]
        max_score = max([max([node[1] for node in population]) for population in populations])
        msbt = get_max_by_time(populations, num_timepoints)
    else:
        status = None
        month = None
        initial_score = None
        max_score = None
        msbt = None

    return [status, month, initial_score, max_score], [max_mutations, max_pathways], msbt


if __name__ == "__main__":
    plt.figure(figsize=(15,6))
    args = sys.argv[1:]
    test = '-t' in args
    project_name = args[0]
    df = pd.read_csv('../pathway_map/'+project_name+'/pathway_counts.csv')
    X = df.iloc[:,3:] # select columns with features (currently gene info)
    surv_info = df[['status','months']]  # columns with survival info
    y = np.array(list(zip(surv_info.status, surv_info.months)) , dtype =[('Status', '?'), ('Survival_in_months', '<f8')])
    Xt = OneHotEncoder().fit_transform(X) # one hot encoder w/ categorical variables
    
    model_path = 'rf/' + project_name + '_rf.pkl'#args[1]
    # model_path = 'rf/TCGA-BRCA_rf_0.6.pkl'#args[1]
    #model_path = project_name + '.pkl'#args[1]
    model = joblib.load(model_path)

    NUM_MEASURES = 2

    sample_list = os.listdir('../result/'+project_name)
    if test:
        random.shuffle(sample_list)
        sample_list = sample_list[:10]
    result = np.zeros((len(sample_list), NUM_MEASURES + 2))
    mutation_count_alive = {}
    mutation_count_dead = {}
    pathway_count_alive = {}
    pathway_count_dead = {}
    num_timepoints = 10
    max_scores_alive = np.array([])
    max_scores_dead = np.array([])
    for i,sample_barcode in tqdm(enumerate(sample_list)):
        if len(sample_barcode) != 12: # skip non-sample files
            continue
        measures, mutation_data, msbt = process_sample(sample_barcode, project_name, df, model, num_timepoints)
        result[i] = measures
        # count the number of pathways
        if measures[0] == 0:
            for mutation in mutation_data[0]:
                if mutation not in mutation_count_alive:
                    mutation_count_alive[mutation] = 0
                mutation_count_alive[mutation] += 1
            for pathway in mutation_data[1]:
                if pathway not in pathway_count_alive:
                    pathway_count_alive[pathway] = 0
                pathway_count_alive[pathway] += 1
            if max_scores_alive.shape[0] < 1:
                max_scores_alive = np.array([msbt])
            else:
                max_scores_alive = np.append(max_scores_alive, np.array([msbt]), axis=0)
        elif measures[0] == 1:
            for mutation in mutation_data[0]:
                if mutation not in mutation_count_dead:
                    mutation_count_dead[mutation] = 0
                mutation_count_dead[mutation] += 1
            for pathway in mutation_data[1]:
                if pathway not in pathway_count_dead:
                    pathway_count_dead[pathway] = 0
                pathway_count_dead[pathway] += 1
            if max_scores_dead.shape[0] < 1:
                max_scores_dead = np.array([msbt])
            else:
                max_scores_dead = np.append(max_scores_dead, np.array([msbt]), axis=0)
    mutation_count_alive = sorted(mutation_count_alive.items(), key=lambda x: x[1], reverse=True)
    mutation_count_dead = sorted(mutation_count_dead.items(), key=lambda x: x[1], reverse=True)
    pathway_count_alive = sorted(pathway_count_alive.items(), key=lambda x: x[1], reverse=True)
    pathway_count_dead = sorted(pathway_count_dead.items(), key=lambda x: x[1], reverse=True)
    non_nan_rows = ~np.isnan(result).any(axis=1)

    # Filter the matrix to remove rows with NaN values
    result = result[non_nan_rows]
    result0 = result[result[:, 0] == 0]
    result1 = result[result[:, 0] == 1]
    result = np.array(result)
    # Get the indices that would sort the array by max score
    sorted_indices = np.argsort(result[:, 3])

    # Use the indices to sort the rows of the matrix
    sorted_result = result[sorted_indices]
    print()
    for ind in sorted_indices[:10]:
        print(sample_list[ind], result[ind])
    print()
    for ind in sorted_indices[-10:]:
        print(sample_list[ind], result[ind])

    # plt.subplot(2, 5, 1)
    # plt.boxplot([result0[:,2], result1[:,2]], labels=['Alive', 'Dead'])
    # plt.title('Initial Score')
    # plt.subplot(2, 5, 2)
    # plt.boxplot([result0[:,3], result1[:,3]], labels=['Alive', 'Dead'])
    # plt.title('Max Score')
    # plt.subplot(2, 5, 3)
    # plt.scatter(result1[:,1], result1[:,3], label='Dead', color='red', s=0.5)
    # plt.scatter(result0[:,1], result0[:,3], label='Alive', color='blue', s=0.5)
    # plt.xlabel('Months')
    # plt.ylabel('Max Score')
    # plt.legend()
    # plt.subplot(2, 5, 4)
    # plt.scatter(result1[:,2], result1[:,3], label='Dead', color='red', s=0.5)
    # plt.scatter(result0[:,2], result0[:,3], label='Alive', color='blue', s=0.5)
    # plt.ylabel('Max Score')
    # plt.xlabel('Initial Score')
    # plt.legend()
    # plt.subplot(2, 5, 5)
    # plt.bar([x[0] for x in mutation_count_alive[:10]], np.array([x[1] for x in mutation_count_alive[:10]])/len(result[result[:,0] == 0]), color='blue')
    # plt.xticks(rotation=90)
    # plt.title('Top 10 Mutations Alive')
    # plt.subplot(2, 5, 6)
    # plt.bar([x[0] for x in mutation_count_dead[:10]], np.array([x[1] for x in mutation_count_dead[:10]])/len(result[result[:,0] == 1]), color='red')
    # plt.xticks(rotation=90)
    # plt.title('Top 10 Mutations Dead')
    # plt.subplot(2, 5, 7)
    # plt.bar([x[0] for x in pathway_count_alive[:10]], np.array([x[1] for x in pathway_count_alive[:10]])/len(result[result[:,0] == 0]), color='blue')
    # plt.xticks(rotation=90)
    # plt.title('Top 10 Pathways Alive')
    # plt.subplot(2, 5, 8)
    # plt.bar([x[0] for x in pathway_count_dead[:10]], np.array([x[1] for x in pathway_count_dead[:10]])/len(result[result[:,0] == 1]), color='red')
    # plt.xticks(rotation=90)
    # plt.title('Top 10 Pathways Dead')
    plt.subplot(1,2,1)
    avg_scores_alive = np.mean(max_scores_alive, axis=0)
    avg_scores_dead = np.mean(max_scores_dead, axis=0)
    stdev_scores_alive = np.std(max_scores_alive, axis=0)
    stdev_scores_dead = np.std(max_scores_dead, axis=0)
    plt.plot(avg_scores_alive, label='Surviving', color='blue', marker='o')
    plt.plot(avg_scores_dead, label='Deceased', color='red', marker='o')
    plt.fill_between(range(num_timepoints + 1), (avg_scores_alive - stdev_scores_alive), (avg_scores_alive + stdev_scores_alive), color='blue', alpha=0.2)
    plt.fill_between(range(num_timepoints + 1), (avg_scores_dead - stdev_scores_dead), (avg_scores_dead + stdev_scores_dead), color='red', alpha=0.2)
    plt.xlabel('Relative Time')
    plt.ylabel('Max Score')
    plt.legend()
    plt.title('Max Score by Time')
    plt.subplot(1,2,2)
    pvals = [ttest_ind(max_scores_alive[:,i], max_scores_dead[:,i])[1] for i in range(num_timepoints + 1)]
    plt.plot(pvals, marker='o')
    plt.xlabel('Relative Time')
    plt.ylabel('P-Value')
    plt.title('T-Test P-Values')
    plt.show()



    
    