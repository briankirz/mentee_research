import sys
from matplotlib import pyplot as plt
import numpy as np

# TODO: Include introgressed segment areas
def compare_3(naive_gamma, bw_gamma, true_intro_states):
    ng = naive_gamma[:, 1]
    bw = bw_gamma[:, 1]
    tis = true_intro_states
    x_axis_start = 0
    x_axis_end = 40

    # Zoom in on introgressed states
    # zoom_start = 8500
    # zoom_stop = 8620
    # ng = naive_gamma[zoom_start:zoom_stop, 1]
    # bw = bw_gamma[zoom_start:zoom_stop, 1]
    # tis = true_intro_states[zoom_start:zoom_stop]
    # x_axis_start = zoom_start
    # x_axis_end = zoom_stop

    plt.rcParams.update({'font.size': 12})
    fig, axes = plt.subplots(1, 3, figsize=(12, 9), sharex=False, sharey=False, dpi=100)
    axes[0].plot(np.linspace(x_axis_start, x_axis_end, ng.shape[0]), ng, marker='.', color='tab:blue')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)
    axes[0].set_title('Naive HMM')
    axes[1].plot(np.linspace(x_axis_start, x_axis_end, bw.shape[0]), bw, marker='.', color='tab:green')
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)
    axes[1].set_title('Baum-Welch HMM')
    axes[2].plot(np.linspace(x_axis_start, x_axis_end, tis.shape[0]), tis, marker='.', color='tab:orange')
    axes[2].spines['top'].set_visible(False)
    axes[2].spines['right'].set_visible(False)
    axes[2].set_title('True introgressed states')
    # Delete "Measured in thousands" when zooming
    fig.supxlabel('500bp window measured in thousands', size=14, weight='bold')
    fig.supylabel('Percent of introgression in window according to model', size=14, weight='bold')
    fig.suptitle('Identifying Introgressed Segments', weight='bold')
    fig.tight_layout()
    plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)

    # ZOOMING
    # plt.rcParams.update({'font.size': 12})
    # fig, axes = plt.subplots(1, 2, figsize=(12, 3), sharex=False, sharey=False, dpi=100)
    # axes[0].plot(np.linspace(x_axis_start, x_axis_end, bw.shape[0]), bw, marker='.', color='tab:green')
    # axes[0].spines['top'].set_visible(False)
    # axes[0].spines['right'].set_visible(False)
    # axes[0].set_title('Baum-Welch HMM')
    # axes[1].plot(np.linspace(x_axis_start, x_axis_end, ng.shape[0]), ng, marker='.', color='tab:blue')
    # axes[1].spines['top'].set_visible(False)
    # axes[1].spines['right'].set_visible(False)
    # axes[1].set_title('Naive HMM')
    # # Delete "Measured in thousands" when zooming
    # fig.supxlabel('500bp window', size=14, weight='bold')
    # fig.supylabel('Percent of introgression in window according to model', size=14, weight='bold')
    # fig.suptitle('Identifying Introgressed Segments', weight='bold')
    # fig.tight_layout()
    # plt.savefig('/Users/briankirz/Downloads/gamma_test.png', facecolor='white', dpi=300)