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
    
    max_scores_alive_luad, max_scores_dead_luad = get_scores_timeline('TCGA-LUAD', test = test, num_timepoints=num_timepoints)
    avg_scores_alive_luad = np.mean(max_scores_alive_luad, axis=0)
    avg_scores_dead_luad = np.mean(max_scores_dead_luad, axis=0)
    stdev_scores_alive_luad = np.std(max_scores_alive_luad, axis=0)
    stdev_scores_dead_luad = np.std(max_scores_dead_luad, axis=0)
    
    max_scores_alive_coad, max_scores_dead_coad = get_scores_timeline('TCGA-COAD', test = test, num_timepoints=num_timepoints)
    avg_scores_alive_coad = np.mean(max_scores_alive_coad, axis=0) 
    avg_scores_dead_coad = np.mean(max_scores_dead_coad, axis=0)
    stdev_scores_alive_coad = np.std(max_scores_alive_coad, axis=0)
    stdev_scores_dead_coad = np.std(max_scores_dead_coad, axis=0)
    
    max_scores_alive_hnsc, max_scores_dead_hnsc = get_scores_timeline('TCGA-HNSC', test = test, num_timepoints=num_timepoints)
    avg_scores_alive_hnsc = np.mean(max_scores_alive_hnsc, axis=0) 
    avg_scores_dead_hnsc = np.mean(max_scores_dead_hnsc, axis=0)
    stdev_scores_alive_hnsc = np.std(max_scores_alive_hnsc, axis=0)
    stdev_scores_dead_hnsc = np.std(max_scores_dead_hnsc, axis=0)
    
    
    fig, axs = plt.subplots(3, 2, figsize=(10, 12))
    fig.tight_layout(pad=4)  # Space between plots
    
    
    plt.subplot(3,2,1)
    plt.plot(avg_scores_alive_luad, label='Surviving', color='blue', marker='o')
    plt.plot(avg_scores_dead_luad, label='Deceased', color='red', marker='o')
    plt.fill_between(range(num_timepoints), (avg_scores_alive_luad - stdev_scores_alive_luad), (avg_scores_alive_luad + stdev_scores_alive_luad), color='blue', alpha=0.2)
    plt.fill_between(range(num_timepoints), (avg_scores_dead_luad - stdev_scores_dead_luad), (avg_scores_dead_luad + stdev_scores_dead_luad), color='red', alpha=0.2)
    plt.xlabel('Relative Time')
    plt.ylabel('Risk Score of Most Dangerous Clone')
    plt.legend()
    plt.title('(a) Max Risk Score of TCGA-LUAD by Time')
    plt.subplot(3,2,2)
    pvals = [ttest_ind(max_scores_alive_luad[:,i], max_scores_dead_luad[:,i])[1] for i in range(num_timepoints)]
    plt.plot(pvals, marker='o', color='blue')
    plt.xlabel('Relative Time')
    plt.ylabel('P-Value')
    plt.title('(b) T-Test P-Values of TCGA-LUAD')
    
    plt.subplot(3,2,3)
    plt.plot(avg_scores_alive_coad, label='Surviving', color='blue', marker='o')
    plt.plot(avg_scores_dead_coad, label='Deceased', color='red', marker='o')
    plt.fill_between(range(num_timepoints), (avg_scores_alive_coad - stdev_scores_alive_coad), (avg_scores_alive_coad + stdev_scores_alive_coad), color='blue', alpha=0.2)
    plt.fill_between(range(num_timepoints), (avg_scores_dead_coad - stdev_scores_dead_coad), (avg_scores_dead_coad + stdev_scores_dead_coad), color='red', alpha=0.2)
    plt.xlabel('Relative Time')
    plt.ylabel('Risk Score of Most Dangerous Clone')
    plt.legend()
    plt.title('(c) Max Risk Score by TCGA-COAD Time')
    plt.subplot(3,2,4)
    pvals_coad = [ttest_ind(max_scores_alive_coad[:,i], max_scores_dead_coad[:,i])[1] for i in range(num_timepoints)]
    plt.plot(pvals_coad, marker='o', color='blue')
    plt.xlabel('Relative Time')
    plt.ylabel('P-Value')
    plt.title('(d) T-Test P-Values of TCGA-COAD')
    
    plt.subplot(3,2,5)
    plt.plot(avg_scores_alive_hnsc, label='Surviving', color='blue', marker='o')
    plt.plot(avg_scores_dead_hnsc, label='Deceased', color='red', marker='o')
    plt.fill_between(range(num_timepoints), (avg_scores_alive_hnsc - stdev_scores_alive_hnsc), (avg_scores_alive_hnsc + stdev_scores_alive_hnsc), color='blue', alpha=0.2)
    plt.fill_between(range(num_timepoints), (avg_scores_dead_hnsc - stdev_scores_dead_hnsc), (avg_scores_dead_hnsc + stdev_scores_dead_hnsc), color='red', alpha=0.2)
    plt.xlabel('Relative Time')
    plt.ylabel('Risk Score of Most Dangerous Clone')
    plt.legend()
    plt.title('(e) Max Risk Score by TCGA-HNSC Time')
    plt.subplot(3,2,6)
    pvals_hnsc = [ttest_ind(max_scores_alive_hnsc[:,i], max_scores_dead_hnsc[:,i])[1] for i in range(num_timepoints)]
    plt.plot(pvals_hnsc, marker='o', color='blue')
    plt.xlabel('Relative Time')
    plt.ylabel('P-Value')
    plt.title('(f) T-Test P-Values of TCGA-HNSC')
    
    plt.savefig('max_score_by_time.png')
    
    