import analyzer
import numpy as np
import analysis.eightplot_analyzer
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def add_to_dict(d, key):
    if key in d:
        d[key] += 1
    else:
        d[key] = 1

def get_scores_timeline(project_name, test=False, num_timepoints=10):

    a = analysis.eightplot_analyzer.EightPlotTimeline(project_name, test, num_timepoints)

    init_scores_alive = []
    init_scores_dead = []
    max_scores_alive = []
    max_scores_dead = []
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
        init_score = vals[4]
        max_score = vals[5]
        if (status[0] == 0):
            init_scores_alive.append(init_score)
            max_scores_alive.append(max_score)
            for mutation in max_mutations:
                add_to_dict(mutation_count_alive, mutation)
            for pathway in max_pathways:
                add_to_dict(pathway_count_alive, pathway)
            months_alive.append(status[1])
        else:
            init_scores_dead.append(init_score)
            max_scores_dead.append(max_score)
            for mutation in max_mutations:
                add_to_dict(mutation_count_dead, mutation)
            for pathway in max_pathways:
                add_to_dict(pathway_count_dead, pathway)
            months_dead.append(status[1])
                
    mutation_count_alive = sorted(mutation_count_alive.items(), key=lambda x: x[1], reverse=True)
    mutation_count_dead = sorted(mutation_count_dead.items(), key=lambda x: x[1], reverse=True)
    pathway_count_alive = sorted(pathway_count_alive.items(), key=lambda x: x[1], reverse=True)
    pathway_count_dead = sorted(pathway_count_dead.items(), key=lambda x: x[1], reverse=True)
                
    return init_scores_alive, init_scores_dead, max_scores_alive, max_scores_dead, mutation_count_alive, mutation_count_dead, pathway_count_alive, pathway_count_dead, months_alive, months_dead

if __name__ == "__main__":
    project_name = 'TCGA-COAD'
    test = False
    result = get_scores_timeline(project_name, test=test, num_timepoints=10)
    init_scores_alive, init_scores_dead, max_scores_alive, max_scores_dead, mutation_count_alive, mutation_count_dead, pathway_count_alive, pathway_count_dead, months_alive, months_dead = result
    
    fig, axs = plt.subplots(2, 4, figsize=(16, 8))
    fig.tight_layout(pad=3)  # Space between plots
    
    
    axs[0,0].boxplot([init_scores_alive, init_scores_dead], tick_labels=['Surviving', 'Deceased'])
    axs[0,0].set_title('(a) Risk Score of Root Nodes')
    axs[0,0].set_ylabel('Risk Score')
    axs[0,1].boxplot([max_scores_alive, max_scores_dead], tick_labels=['Surviving', 'Deceased'])
    axs[0,1].set_title('(b) Risk Score of the Most Dangerous Clone')
    axs[0,1].set_ylabel('Risk Score')
    axs[0,2].scatter(months_alive, max_scores_alive, label='Surviving', color='blue', s=0.5)
    axs[0,2].scatter(months_dead, max_scores_dead, label='Deceased', color='red', s=0.5)
    axs[0,2].set_title('(c) Survival Time vs. Max Score')
    axs[0,2].set_xlabel('Months')
    axs[0,2].set_ylabel('Risk Score')
    axs[0,3].scatter(init_scores_alive, max_scores_alive, label='Surviving', color='blue', s=0.5)
    axs[0,3].scatter(init_scores_dead, max_scores_dead, label='Deceased', color='red', s=0.5)
    axs[0,3].set_title('(d) Score of Root vs. Max')
    axs[0,3].set_xlabel('Root Score')
    axs[0,3].set_ylabel('Max Score')
    plt.subplot(2, 4, 5)
    
    ef_up = max([x[1]/len(init_scores_alive) for x in mutation_count_alive[:10]] + [x[1]/len(init_scores_dead) for x in mutation_count_dead[:10]]) * 1.1
    
    axs[1,0].bar([x[0] for x in mutation_count_alive[:10]], np.array([x[1] for x in mutation_count_alive[:10]])/len(init_scores_alive), color='blue')
    plt.xticks(rotation=90)
    plt.ylim(0, ef_up)
    axs[1,0].set_title('(e) Top 10 Mutations in Surviving Patients')
    plt.ylabel('Frequency')
    
    plt.subplot(2, 4, 6)
    axs[1,1].bar([x[0] for x in mutation_count_dead[:10]], np.array([x[1] for x in mutation_count_dead[:10]])/len(init_scores_dead), color='red')
    plt.xticks(rotation=90)
    plt.ylim(0, ef_up)
    axs[1,1].set_title('(f) Top 10 Mutations in Deceased Patients')
    plt.ylabel('Frequency')
    
    gh_up = max([x[1]/len(init_scores_alive) for x in pathway_count_alive[:10]] + [x[1]/len(init_scores_dead) for x in pathway_count_dead[:10]]) * 1.1
    
    plt.subplot(2, 4, 7)
    axs[1,2].bar([x[0] for x in pathway_count_alive[:10]], np.array([x[1] for x in pathway_count_alive[:10]])/len(init_scores_alive), color='blue')
    plt.xticks(rotation=90)
    plt.ylim(0, gh_up)
    axs[1,2].set_title('(g) Top 10 Pathways in Surviving Patients')
    plt.ylabel('Frequency')
    
    plt.subplot(2, 4, 8)
    axs[1,3].bar([x[0] for x in pathway_count_dead[:10]], np.array([x[1] for x in pathway_count_dead[:10]])/len(init_scores_dead), color='red')
    plt.xticks(rotation=90)
    plt.ylim(0, gh_up)
    axs[1,3].set_title('(h) Top 10 Pathways in Deceased Patients')
    plt.ylabel('Frequency')
    
    
    # plt.savefig(f'../report/img/{project_name[5:].lower()}.png')
    plt.show()