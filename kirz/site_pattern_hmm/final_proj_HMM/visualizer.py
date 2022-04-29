import sys
from matplotlib import pyplot as plt
import numpy as np

# TODO: Include introgressed segment areas
def visualizer(naive_gamma, bw_gamma, true_intro_states):
    ng = naive_gamma[:, 1]
    bw = bw_gamma[:, 1]
    tis = true_intro_states
    baum_welch = np.ones((2, 2))

    plt.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(1, 3, figsize=(12, 9), sharex=False, sharey=False, dpi=100)
    axes[0].plot(np.linspace(0, 40, ng.shape[0]), ng, marker='.', color='tab:blue')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].set_title('Naive HMM')
    axes[1].plot(np.linspace(0, 40, bw.shape[0]), bw, marker='.', color='tab:green')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].set_title('Baum-Welch HMM')
    axes[2].plot(np.linspace(0, 40, tis.shape[0]), tis, marker='.', color='tab:orange')
    axes[2].spines['top'].set_visible(False)
    axes[2].spines['right'].set_visible(False)
    axes[2].set_title('True introgressed states')
    fig.supxlabel('500bp window, measured in thousands', size=14, weight='bold')
    fig.supylabel('Percent of introgression in window according to model', size=14, weight='bold')
    fig.suptitle('Identifying Introgressed Segments', weight='bold')
    fig.tight_layout()
    plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)