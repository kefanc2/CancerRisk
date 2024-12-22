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


def plot_gene(gene, ax, title=''):
    timeline_alive = [mutation_timeline_alive[j][gene]/alive_count if gene in mutation_timeline_alive[j] else 0 for j in range(len(mutation_timeline_alive))]
    ax.plot(timeline_alive, label='Surviving')
    timeline_dead = [mutation_timeline_dead[j][gene]/dead_count if gene in mutation_timeline_dead[j] else 0 for j in range(len(mutation_timeline_dead))]
    ax.plot(timeline_dead, label='Deceased')
    ax.legend()
    if title == '':
        ax.set_title('Mutation Timeline of ' + gene)
    else:
        ax.set_title(title)
    ax.set_xlabel('Relative Time')
    ax.set_ylabel('Frequency')
    
    
if __name__ == "__main__":
    project_name = 'TCGA-LUAD'
    test = False
    num_timepoints = 100
    mutation_timeline_alive, mutation_timeline_dead, pathway_timeline_alive, \
        pathway_timeline_dead, alive_count, dead_count = get_mutation_timeline(project_name, test=test, num_timepoints=num_timepoints)
    k = 5
    
    # plot_gene('TP53', axs[0])
    while True:
        gene = input('Enter gene name to plot: ')
        if gene == 'exit' or gene == 'quit' or gene == 'q' or not gene:
            break
        title = input('Enter title (optional): ')
        fig, axs = plt.subplots(1, 1, figsize=(8, 6))
        fig.tight_layout(pad=3)  # Space between plots
        plot_gene(gene, axs, title)
        plt.show()
    # plot_gene('EGFR', axs, 'EGFR in LUAD')
    # plt.show()
        
    