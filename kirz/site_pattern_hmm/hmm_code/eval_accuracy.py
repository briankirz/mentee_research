import sys
import numpy as np


def eval_accuracy(true_windows, tested_HMM, normalized=True, threshold=.9):
    # matrix of windows containing introgressed sites (1. = 100% coverage)
    tiw = true_windows
    # takes the second column of gamma, which measures the likelihood of Neanderthal ancestry
    gamma = tested_HMM[:, 1]

    # performance values preset to -1 to easily show errors
    false_pos_r = -1
    miss_rate = -1
    sensitivity = -1
    specificity = -1


    # If we need to normalize, change values of gamma such that the highest-scoring site is set to 100%
    if normalized:
        # Find the value of the highest-scoring site in gamma
        max_prob = np.amax(gamma)
        # Multiply every probability in gamma (including the highest-scoring site) by 1/max_prob
        gamma = gamma * (1 / max_prob)

    # Compare gamma and tiw to see if they have the same number of elements (windows)
    if len(tiw) == len(gamma):

        fp = 0  # number of false positives
        fn = 0  # number of false negatives
        tp = 0  # number of true postivies
        tn = 0  # number of true negatives

        # loop through gamma and record the number of true/false positives/negatives
        for w in range(len(tiw)):
            # Underlying window is partially or completely introgressed
            if 0 < tiw[w] <= 1.:
                # true positive
                if gamma[w] >= threshold:
                    tp += 1 # tiw[w]
                # false negative
                else:
                    fn += 1 # tiw[w]
            # Underlying window is not 100% introgressed
            elif tiw[w] == 0:
                # false positive
                if gamma[w] >= threshold:
                    fp += 1 # tiw[w]
                # true negative
                else:
                    tn += 1 # tiw[w]
            # something went wrong and tiw shows a value below zero or above 1
            else:
                print("ERROR in eval: window shows introgression percentage below zero or above 1")

        # check at the end to see if tp + tn + fp + fn == number of windows
        if tp + tn + fp + fn != len(tiw):
            print("ERROR in eval: true positives and true negatives do not sum to the number of windows")
            print(tp + tn + fp + fn)

        # Calculate rates for performance
        print("#tp: " + str(tp) + "\t#tn: " + str(tn), "\t#fp: " + str(fp), "\t#fn: " + str(fn))
        false_pos_r = fp / (fp + tn)
        miss_rate = fn / (fn + tp)
        sensitivity = tp / (tp + fn)
        specificity = tn / (tn + fp)

    else:
        print("ERROR IN EVAL: C column of gamma and true introgressed window array have different number of windows")

    performance = np.array([false_pos_r, miss_rate, sensitivity, specificity])
    return performance