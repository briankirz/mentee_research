import sys
import numpy as np
import gzip
import matplotlib.pyplot as plt

def calc_rep_performance(rep_num, Rep_performances, g_threshold, w_threshold, pattern, dxy):
    
    if dxy:
        distance = 'dxy'
    else:
        distance = ''
        
    # download the array from the rep_results filepath
    rep = np.genfromtxt(rep_filepath, delimiter='\t')
    # max number of BW iterations, should be 100
    max_iter = len(rep[0])-4
    fpr_array = np.zeros(max_iter)
    tpr_array = np.zeros(max_iter)
    
    
    # set the OSCAR filepath 
    rep_filepath = './hmm_results/BW{0}_wthreshold{1]_{2}_{3}_gthreshold{4}_prufer_results_rep_id_{5}.csv.gz'.format(str(max_iter), str(w_threshold), pattern, distance, str(g_threshold), rep_num+1)


    
    
    # Starting in the fourth column (naive gamma) to the end
    for gamma in range(4, max_iter+4):
        # initilaize false positives and false negatives for this gamma
        fp = 0  # number of false positives
        fn = 0  # number of false negatives
        tp = 0  # number of true postivies
        tn = 0  # number of true negatives
        # for each window in the rep (should be 40k)
        for w in range(0, len(rep)):
            # % of true introgression at a window
            true_val = rep[w][2]
            # the percentage introgression guessed by the model at that gamma
            gamma_val = rep[w][gamma]
            if 0 < true_val <= 1.:
                # true positive
                if gamma_val >= threshold:
                    tp += 1
                # false negative
                else:
                    fn += 1
            # Underlying window is not 100% introgressed
            elif true_val == 0:
                # false positive
                if gamma_val >= threshold:
                    fp += 1
                # true negative
                else:
                    tn += 1
            else:
                print("Error in eval: window shows introgression percentage below zero or above 1")      
        # once every window is tallied (loop is over), calculate TPR/FPR
        fpr = fp / (fp + tn)
        tpr = tp / (tp + fn)
        fpr_array[gamma-4] = fpr
        tpr_array[gamma-4] = tpr        
    # once every gamma is calculated, the fpr/tpr matrices are complete
    # add them to dictionary and return it
    R_p = Rep_performances
    R_p[rep_num] = fpr_array, tpr_array
    return R_p

# Creates two tables of performances (all true/false positive rates)
def table_performances(Rep_performances):
    # extract dimensions
    num_reps = len(Rep_performances) # should be 1000
    max_BW = len(Rep_performances[list(Rep_performances.keys())[0]][0]) # should be 100
    # initialize result matrices
    all_fprs = np.empty((num_reps, max_BW))
    all_tprs = np.empty((num_reps, max_BW))
    
    for key in Rep_performances:
        all_fprs[int(key)-1] = Rep_performances[key][0]
        all_tprs[int(key)-1] = Rep_performances[key][1]
    return all_fprs, all_tprs

def avg_performance(total_reps=1000, g_threshold=.9, w_threshold=1., pattern, dxy):
    
    # guess_threshold = g_threshold
    # window_threshold = w_threshold
    # pattern = pattern
    # dxy = dxy
    
    
    # Initialize the performance dictionary
    Rep_performances = {}
    for rep_num in range(total_reps):
        
        # add the performance for rep #rep_num to the dictionary
        Rep_performances = calc_rep_performance(rep_num, Rep_performances, g_threshold, w_threshold, pattern, dxy)
        
        
    # turn the dictionary of results into two nparrays (fpr, tpr)
    perf_tables = table_performances(Rep_performances)
    fpr = perf_tables[0]
    tpr = perf_tables[1]
    
    #Set dimensions of final arrays
    num_reps = total_reps
    max_BW = len(Rep_performances[list(Rep_performances.keys())[0]][0])
    
    # Initialize empty 2D arrays.
    avg_fpr = np.empty((num_reps, max_BW)) # 1000 reps x 100 BW iters
    avg_tpr = np.empty((num_reps, max_BW)) # 1000 reps x 100 BW iters
    
    # create mask arrays of the same size
    fpr_mask = np.ones(fpr.shape, dtype=bool)
    tpr_mask = np.ones(tpr.shape, dtype=bool)
    # initialize all pre-convergence indices as false
    for row in range(num_reps):
        fpr_mask[row, :np.nonzero(fpr[row])[0][-1]+1] = False
        tpr_mask[row, :np.nonzero(fpr[row])[0][-1]+1] = False
    # use np.where to quickly convert to nans: where mask is false fill with an np.nan
    nan_fpr = np.where(fpr_mask, np.nan, fpr)
    nan_tpr = np.where(tpr_mask, np.nan, tpr)
    
    # Average performances for both arrays
    avg_fpr = np.nanmean(nan_fpr, axis=0)
    avg_tpr = np.nanmean(nan_tpr, axis=0)    
    
    
    
    
    # HAVE TO INCLUDE INPUTS HERE
    
    
    # TODO: Label output plots with modular input
    
    # ROC CURVE
    plt.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(9, 6), dpi=300)
    plt.title('Receiver Operating Characteristic')
    plt.plot(avg_fpr, avg_tpr, marker='o', markersize='3', markerfacecolor='r', color='b')
    plt.legend(["1000-Rep Average"], loc="upper left")
    plt.plot([0, .5], [0, .5],'r--')
    plt.xlim([0, .1])
    plt.ylim([0, .5])
    plt.ylabel('Average True Positive Rate')
    plt.xlabel('Average False Positive Rate')
    plt.savefig('./hmm_results/{0}_rep_ROC'.format(total_reps), facecolor='white', dpi=300)
    # local filepath
    # plt.savefig('./{0}_rep_ROC'.format(total_reps), facecolor='white', dpi=300)
    
    
    # TPR/FPR To BW Iteration Ratio
    plt.rcParams.update({'font.size': 12})
    fig = plt.figure(figsize=(9, 6), dpi=300)
    plt.title('Average TPR/FPR Ratio as Baum-Welch Iterates')
    plt.plot(np.arange(max_BW), list(avg_tpr/avg_fpr), marker='o', markersize='3', markerfacecolor='r', color='b')
    plt.legend(["1000-Rep Average"], loc="upper left")
    plt.plot([0, .5], [0, .5],'r--')
    plt.xlim([0, max_BW])
    plt.ylim([0, 200])
    plt.ylabel('True Positive Rate / False Positive Rate')
    plt.xlabel('Baum-Welch Iteration Number')
    plt.show()    
    plt.savefig('./hmm_results/{0}_rep_pr_ratio'.format(total_reps), facecolor='white', dpi=300)
    # local filepath
    #plt.savefig('./{0}_rep_pr_ratio'.format(total_reps), facecolor='white', dpi=300)
    
    return avg_fpr, avg_tpr
    

    
# Read in sys args
total_reps = int(sys.argv[1])
threshold = float(sys.argv[2])

avg_performance(total_reps, threshold)



