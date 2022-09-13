import sys
import numpy as np
import gzip

def calc_rep_performance(rep_num, Rep_performances, bw_limit, g_threshold, w_threshold, pattern, dxy):

    # max number of BW iterations, should be 100
    fpr_array = np.zeros(bw_limit)
    tpr_array = np.zeros(bw_limit)
    
    # set the results filepath 
    rep_filepath = './hmm_results/results_BW{0}_wt{1}_{2}_{3}_prufer_rep_id_{4}.csv.gz'.format(str(bw_limit), str(w_threshold), pattern, dxy, rep_num+1)

    # download the array from the rep_results filepath
    rep = np.genfromtxt(rep_filepath, delimiter='\t')
    
    
    # Starting in the fourth column (naive gamma) to the end
    for gamma in range(4, bw_limit+4):
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
                if gamma_val >= g_threshold:
                    tp += 1
                # false negative
                else:
                    fn += 1
            # Underlying window is not 100% introgressed
            elif true_val == 0:
                # false positive
                if gamma_val >= g_threshold:
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

def avg_performance(total_reps, bw_limit, g_threshold=.9, w_threshold=1., pattern='patterna', dxy='nodxy'):

    # Initialize the performance dictionary
    Rep_performances = {}
    for rep_num in range(total_reps):
        # add the performance for rep #rep_num to the dictionary
        Rep_performances = calc_rep_performance(rep_num, Rep_performances, bw_limit, g_threshold, w_threshold, pattern, dxy)
        
        
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
    
    np.savetxt('./performance_results/avg_fpr_gt{0}_BW{1}_wt{2}_{3}_{4}_prufer_{5}_reps.csv.gz'.format(str(g_threshold), str(bw_limit), str(w_threshold), pattern, dxy, str(total_reps)), avg_fpr, fmt='%1.3f', delimiter='\t', newline='\n')
    np.savetxt('./performance_results/avg_tpr_gt{0}_BW{1}_wt{2}_{3}_{4}_prufer_{5}_reps.csv.gz'.format(str(g_threshold), str(bw_limit), str(w_threshold), pattern, dxy, str(total_reps)), avg_tpr, fmt='%1.3f', delimiter='\t', newline='\n')
    
    return avg_fpr, avg_tpr
    

    
# Read in sys args
total_reps = int(sys.argv[1])
bw_limit = int(sys.argv[2])
g_threshold = float(sys.argv[3])
w_threshold = float(sys.argv[4])
pattern = str(sys.argv[5])
dxy = str(sys.argv[6])

avg_performance(total_reps, bw_limit, g_threshold, w_threshold, pattern, dxy)


BW100_wt1.0_patternc_prufer_results_rep_id_794.csv.gz
