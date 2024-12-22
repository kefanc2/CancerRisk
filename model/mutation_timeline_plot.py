import analyzer
import numpy as np
import analysis.mutation_timeline
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def add_to_dict(d, key):
    if key in d:
        d[key] += 1
    else:
        d[key] = 1

def get_mutation_timeline(project_name, test=False, num_timepoints=10):

    a = analysis.mutation_timeline.MutationTimeline(project_name, test, num_timepoints)

    for vals in a.process():
        pass
    
    return a.analyze()
    
    
if __name__ == "__main__":
    project_name = 'TCGA-BRCA'
    test = False
    num_timepoints = 100
    mutation_timeline_alive, mutation_timeline_dead, pathway_timeline_alive, \
        pathway_timeline_dead, alive_count, dead_count = get_mutation_timeline(project_name, test=test, num_timepoints=num_timepoints)
    k = 5
    
    mutation_timeline_topk_alive = sorted(mutation_timeline_alive[-1].items(), key=lambda x: x[1], reverse=True)[:k]#[sorted(mutation_timeline_alive[i].items(), key=lambda x: x[1], reverse=True)[:k] for i in range(len(mutation_timeline_alive))]
    mutation_timeline_topk_dead = sorted(mutation_timeline_alive[-1].items(), key=lambda x: x[1], reverse=True)[:k]#[sorted(mutation_timeline_alive[i].items(), key=lambda x: x[1], reverse=True)[:k] for i in range(len(mutation_timeline_dead))]
    pathway_timeline_topk_alive = sorted(mutation_timeline_alive[-1].items(), key=lambda x: x[1], reverse=True)[:k]#[sorted(mutation_timeline_alive[i].items(), key=lambda x: x[1], reverse=True)[:k] for i in range(len(pathway_timeline_alive))]
    pathway_timeline_topk_dead = sorted(pathway_timeline_dead[-1].items(), key=lambda x: x[1], reverse=True)[:k]#[sorted(pathway_timeline_dead[i].items(), key=lambda x: x[1], reverse=True)[:k] for i in range(len(pathway_timeline_dead))]
    # gene_name_timeline_alive_alltime = set([pair[0] for x in mutation_timeline_topk_alive for pair in x])
    # gene_name_timeline_dead_alltime = set([pair[0] for x in mutation_timeline_topk_dead for pair in x])
    # pathway_name_timeline_alive_alltime = set([pair[0] for x in pathway_timeline_topk_alive for pair in x])
    # pathway_name_timeline_dead_alltime = set([pair[0] for x in pathway_timeline_topk_dead for pair in x])
    gene_name_timeline_alive_alltime = set([x[0] for x in mutation_timeline_topk_alive])
    gene_name_timeline_dead_alltime = set([x[0] for x in mutation_timeline_topk_dead])
    pathway_name_timeline_alive_alltime = set([x[0] for x in pathway_timeline_topk_alive])
    pathway_name_timeline_dead_alltime = set([x[0] for x in pathway_timeline_topk_dead])
    
    genes_to_plot = gene_name_timeline_alive_alltime.union(gene_name_timeline_dead_alltime)
    fig, axs = plt.subplots(1, 2, figsize=(16, 7))
    fig.tight_layout(pad=3)  # Space between plots
    for i in genes_to_plot:
        timeline = [mutation_timeline_alive[j][i]/alive_count if i in mutation_timeline_alive[j] else 0 for j in range(len(mutation_timeline_alive))]
        axs[0].plot(timeline, label=i)
    for i in genes_to_plot:
        timeline = [mutation_timeline_dead[j][i]/dead_count if i in mutation_timeline_dead[j] else 0 for j in range(len(mutation_timeline_dead))]
        axs[1].plot(timeline, label=i)
    axs[0].legend()
    axs[1].legend()
    axs[0].set_title('(a) Mutation Timeline in Surviving Patients')
    axs[1].set_title('(b) Mutation Timeline in Deceased Patients')
    axs[0].set_xlabel('Relative Time')
    axs[0].set_ylabel('Frequency')
    axs[1].set_xlabel('Relative Time')
    axs[1].set_ylabel('Frequency')
    upper = max([max(axs[0].get_ylim()), max(axs[1].get_ylim())]) * 1.1
    axs[0].set_ylim(0, upper)
    axs[1].set_ylim(0, upper)
    plt.show()
        
    