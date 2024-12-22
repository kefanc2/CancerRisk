import analyzer
import numpy as np
import analysis.score_timeline
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind

def get_scores_timeline(project_name, test=False, num_timepoints=10):

    a = analysis.score_timeline.ScoreTimeline(project_name, test=test, num_timepoints=num_timepoints)

    max_scores_alive = np.array([])
    max_scores_dead = np.array([])

    for status, tree in a.process():
        if status == None:
            continue
        scores, clones = a.analyze()
        if status[0] == 0:
            if max_scores_alive.shape[0] < 1:
                max_scores_alive = np.array([scores])
            else:
                max_scores_alive = np.append(max_scores_alive, np.array([scores]), axis=0)
        elif status[0] == 1:
            if max_scores_dead.shape[0] < 1:
                max_scores_dead = np.array([scores])
            else:
                max_scores_dead = np.append(max_scores_dead, np.array([scores]), axis=0)
                
    return max_scores_alive, max_scores_dead

if __name__ == "__main__":
    num_timepoints = 10
    test = False
    
    cohorts = ['TCGA-COAD']
    n = len(cohorts)
    
    fig, axs = plt.subplots(n, 2, figsize=(12, 5*n))
    fig.tight_layout(pad=4)  # Space between plots
    
    for i in range(n):
        max_scores_alive, max_scores_dead = get_scores_timeline(cohorts[i], test = test, num_timepoints=num_timepoints)
        avg_scores_alive = np.mean(max_scores_alive, axis=0)
        avg_scores_dead = np.mean(max_scores_dead, axis=0)
        stdev_scores_alive = np.std(max_scores_alive, axis=0)
        stdev_scores_dead = np.std(max_scores_dead, axis=0)
        plt.subplot(n,2,i*2+1)
        plt.plot(avg_scores_alive, label='Surviving', color='blue', marker='o')
        plt.plot(avg_scores_dead, label='Deceased', color='red', marker='o')
        plt.fill_between(range(num_timepoints), (avg_scores_alive - stdev_scores_alive), (avg_scores_alive + stdev_scores_alive), color='blue', alpha=0.2)
        plt.fill_between(range(num_timepoints), (avg_scores_dead - stdev_scores_dead), (avg_scores_dead + stdev_scores_dead), color='red', alpha=0.2)
        plt.xlabel('Relative Time')
        plt.ylabel('Risk Score of Most Dangerous Clone')
        plt.legend()
        plt.title(f'({chr(96+i*2+1)}) Max Risk Score of {cohorts[i]} by Time')
        plt.subplot(n,2,i*2+2)
        pvals = [ttest_ind(max_scores_alive[:,i], max_scores_dead[:,i])[1] for i in range(num_timepoints)]
        plt.plot(pvals, marker='o', color='blue')
        plt.xlabel('Relative Time')
        plt.ylabel('P-Value')
        plt.title(f'({chr(96+i*2+2)}) T-Test P-Values of {cohorts[i]}')
    
    
    plt.show()
    
    