import numpy as np
import analysis.eightplot_analyzer
import matplotlib.pyplot as plt
from scipy.stats import chi2_contingency

def add_to_dict(d, key):
    if key in d:
        d[key] += 1
    else:
        d[key] = 1

def get_scores_timeline(project_name, test=False, num_timepoints=10):

    a = analysis.eightplot_analyzer.EightPlotTimeline(project_name, test, num_timepoints)

    mutation_count_alive = {}
    mutation_count_dead = {}
    pathway_count_alive = {}
    pathway_count_dead = {}
    months_alive = []
    months_dead = []

    for vals in a.process():
        status = vals[0]
        if status == None:
            continue
        max_pathways = vals[2]
        max_mutations = vals[3]
        if (status[0] == 0):
            for mutation in max_mutations:
                add_to_dict(mutation_count_alive, mutation)
            for pathway in max_pathways:
                add_to_dict(pathway_count_alive, pathway)
            months_alive.append(status[1])
        else:
            for mutation in max_mutations:
                add_to_dict(mutation_count_dead, mutation)
            for pathway in max_pathways:
                add_to_dict(pathway_count_dead, pathway)
            months_dead.append(status[1])
                
    mutation_count_alive = sorted(mutation_count_alive.items(), key=lambda x: x[1], reverse=True)
    mutation_count_dead = sorted(mutation_count_dead.items(), key=lambda x: x[1], reverse=True)
    pathway_count_alive = sorted(pathway_count_alive.items(), key=lambda x: x[1], reverse=True)
    pathway_count_dead = sorted(pathway_count_dead.items(), key=lambda x: x[1], reverse=True)
                
    return mutation_count_alive, mutation_count_dead, pathway_count_alive, pathway_count_dead, months_alive, months_dead

def get_diff_pathway_or_mutation(names_list, dict_alive, dict_dead, alive_count, dead_count):
    result = {}
    for i in names_list:
        # surviving and present
        try:
            sp = dict_alive[i]
        except:
            sp = 0
        # surviving and absent
        sa = alive_count - sp
        # deceased and present
        try:
            dp = dict_dead[i]
        except:
            dp = 0
        # deceased and absent
        da = dead_count - dp
        
        if sp == 0 and dp == 0:
            continue
        # contingency table
        table = np.array([[sp, sa], [dp, da]])
        # chi-squared test
        stat, p, dof, expected = chi2_contingency(table)
        if p < 0.05:
            result[i] = p
    return result

if __name__ == "__main__":
    project_name = 'TCGA-HNSC'
    test = False
    result = get_scores_timeline(project_name, test=test, num_timepoints=10)
    mutation_count_alive, mutation_count_dead, pathway_count_alive, pathway_count_dead, months_alive, months_dead = result
    pathway_dict_alive = {}
    pathway_dict_dead = {}
    mutation_dict_alive = {}
    mutation_dict_dead = {}
    alive_count = len(months_alive)
    dead_count = len(months_dead)
    
    for i in range(len(pathway_count_alive)):
        pathway_dict_alive[pathway_count_alive[i][0]] = pathway_count_alive[i][1]
    for i in range(len(pathway_count_dead)):
        pathway_dict_dead[pathway_count_dead[i][0]] = pathway_count_dead[i][1]
    for i in range(len(mutation_count_alive)):
        mutation_dict_alive[mutation_count_alive[i][0]] = mutation_count_alive[i][1]
    for i in range(len(mutation_count_dead)):
        mutation_dict_dead[mutation_count_dead[i][0]] = mutation_count_dead[i][1]
        
    pathway_set = set(pathway_dict_alive.keys()).union(set(pathway_dict_dead.keys()))
    mutation_set = set(mutation_dict_alive.keys()).union(set(mutation_dict_dead.keys()))
    
    diff_genes = get_diff_pathway_or_mutation(mutation_set, mutation_dict_alive, mutation_dict_dead, alive_count, dead_count)
    diff_pathways = get_diff_pathway_or_mutation(pathway_set, pathway_dict_alive, pathway_dict_dead, alive_count, dead_count)
    print(diff_genes)
    print(diff_pathways)
    
    
        
            
        
    
        
    
