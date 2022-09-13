import sys
import numpy as np
import gzip
import matplotlib.pyplot as plt

# The purpose of this function is to take the avg_tpr and avg_fpr filepath out of the ./performance_results folder based on parameter inputs
# Then, generate a graph and put the results into the ./performance_graphs folder

def graph_performance(total_reps, g_threshold, bw_limit, w_threshold, pattern, dxy):
    
    # Set the ./performance_results filepaths
    avg_fpr_filepath = './performance_results/avg_fpr_gt{0}_BW{1}_wt{2}_{3}_{4}_prufer_{5}_reps.csv.gz'.format(str(g_threshold), str(bw_limit), str(w_threshold), pattern, dxy, total_reps)
    avg_tpr_filepath = './performance_results/avg_tpr_gt{0}_BW{1}_wt{2}_{3}_{4}_prufer_{5}_reps.csv.gz'.format(str(g_threshold), str(bw_limit), str(w_threshold), pattern, dxy, total_reps)
    
    #download the average true positive and false positive rates from the filepaths
    avg_fpr = np.genfromtxt(avg_fpr_filepath, delimiter='\t')
    avg_tpr = np.genfromtxt(avg_tpr_filepath, delimiter='\t')
    
    # ROC CURVE
    plt.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(9, 6), dpi=300)
    plt.title('Receiver Operating Characteristic: Baum-Welch Performance')
    plt.plot(avg_fpr, avg_tpr, marker='o', markersize='3', markerfacecolor='r', color='b')
    plt.legend(["{0}-Rep Average".format(total_reps)], loc="upper left")
    plt.plot([0, .5], [0, .5],'r--')
    plt.xlim([0, .2])
    plt.ylim([0, .75])
    plt.ylabel('Average True Positive Rate')
    plt.xlabel('Average False Positive Rate')

    plt.savefig('./performance_graphs/ROC_gt{0}_BW{1}_wt{2}_{3}_{4}_prufer_{5}_reps.png'.format(str(g_threshold), str(bw_limit), str(w_threshold), pattern, dxy, total_reps), facecolor='white', dpi=300)
    
    # TPR/FPR To BW Iteration Ratio
    plt.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(9, 6), dpi=300)
    plt.title('Average TPR/FPR Ratio as Baum-Welch Iterates')
    plt.plot(np.arange(bw_limit), list(avg_tpr/avg_fpr), marker='o', markersize='3', markerfacecolor='r', color='b')
    plt.legend(["{0}-Rep Average".format(total_reps)], loc="upper left")
    plt.plot([0, .5], [0, .5],'r--')
    plt.xlim([0, bw_limit])
    plt.ylim([0, 100])
    plt.ylabel('True Positive Rate / False Positive Rate')
    plt.xlabel('Baum-Welch Iteration Number')
    plt.show()
    plt.savefig('./performance_graphs/ratio_gt{0}_BW{1}_wt{2}_{3}_{4}_prufer_{5}_reps.png'.format(str(g_threshold), str(bw_limit), str(w_threshold), pattern, dxy, total_reps), facecolor='white', dpi=300)
    
    return None

# Read in sys args
total_reps = str(sys.argv[1])
g_threshold = float(sys.argv[2])
bw_limit = int(sys.argv[3])
w_threshold = float(sys.argv[4])
pattern = str(sys.argv[5])
dxy = str(sys.argv[6])

graph_performance(total_reps, g_threshold, bw_limit, w_threshold, pattern, dxy)