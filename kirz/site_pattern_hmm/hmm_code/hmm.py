import sys
import numpy as np
import extract_obs
import support_code as supp
import time
# Commented out when we don't need to use performance visualization
import visualizer as vis
import eval_accuracy as eval

# Set "log_zero" to negative infinity to help with log sums. It will be taken in as an argument in base cases
log_zero = np.NINF


# TODO: Should I reverse the given | syntax for "given" probabilities?
def calc_alpha(A, B, pi, Ob, N, T):
    '''
    Calculates the alpha matrix with the Forward Algorithm
    INPUT:
        - A, a 2x2 array of transition probabilities
        - B, a 2x2 array of emission probabilities
        - pi, a 2x1 array of of initial state probabilities
        - Ob, a string representing the observed sequence (binary-encoded numbers indexed to labels: 0:N::1:C)
        - N, the number of states
        - T, the length of the sequence
    ______________________________________________________
    OUTPUT: ALPHA MATRIX
        - Dimensions: (T+1) x N
        - Definition: The forward variable (a) gives the probability of observing a prefix of the emission sequence and
          being in some given state at the end of the prefix or a(t)(i) = P(o1, o2, ... ot, q(t) = Si)
          for all i in {1...N} and 1<=t<=T where T is the number of observations in the sequence
        - Explanation: The first/second column represents the probability that the model is in the Species/Introgressed
          state after the first t characters in the sequence. This information is found in the t^th row of the matrix.

        Observation example: NN...C
        ---------------------------------------------------------------------------
        | ALPHA     | State 0, or Species           | State 1, or Introgressed    |
        ---------------------------------------------------------------------------
        | t=0       | P(first state is S)           | P(first state is I)         |
        | t=1 (N)   | P(seen N | ends w state S)    | P(seen N | ends w state I)  |
        | t=2 (NN)  | P(seen NN | ends w state S)   | P(seen NN | ends w state I) |
        ...
        | T=len(NN...C) | P((NN...C) | ends up at S)| P((NN...C) | ends w state I)|
        ---------------------------------------------------------------------------

        - Both values in the last row comprise the total probability of the observed sequence being produced
          by the HMM.
    '''

    # Must be T+1 columns in the alpha matrix bc the first one is the state after 0 prefix characters
    # This should be the same as the initial distribution likelihood found in pi
    alpha = np.zeros((T + 1, N))

    # initialize the first row to be the initial distribution values (the pi matrix of the HMM)
    # represents the probabilities of being in some state (S/1st or I/2nd) before seeing any (t=0) observed emissions
    alpha[0, :] = pi

    # Filled in one row at a time, starting with the 2nd row (we've already filled in the 1st row in the last step)
    # T counts the character number in the sequence.
    for t in range(1, T + 1):

        # k stores the character of the previous observed emission
        k = Ob[t - 1]

        # Filling in one column at a time, starting with column 1 (Species state) then column 2 (Introgression state)
        for j in range(N):

            # This placeholder is set to negative infinity the first time each cell is encountered, resetting it.
            # We need it to store a zero value when logged by logaddexp when calculating the first logsum
            lprob = log_zero

            # The i loop occurs in a single cell.
            # In a single cell, we calculate the sum of probabilities (in log form) of the transitions from all possible
            # previous states in time t-1 (the previous row) into the new one.
            # In this case, N=2, meaning there were 2 possible previous states that could have led to the current one
            # This code asks: what is the probability that each possible scenario (previous state being S or I) led to
            # our current state j? It then combines those probabilities and closes both timelines.
            for i in range(N):
                # lp represents a sum of log probabilities:
                # (forward variable at time t-1 for state i)
                # + likelihood that last row's state i transitioned to this state j using the transition matrix A
                # + the likelihood that state i emitted this observed character k using the emission matrix B
                lp = alpha[t - 1][i] + A[i][j] + B[i][k]

                # during the first iteration, lprob is reset as equal to lp, as lprob starts set to NINF
                # the second time around, lp is recalculated and represents the probability that we got to this state
                # j from the Introgressed state. Now, calling logaddexp(lprob, lp) represents the sum of these:
                # (prob we're in state j if the last state was S + prob we're in state j if the last state was I)
                lprob = logaddexp(lprob, lp)

            # After the probabilities based on both of the cells in the previous row were treated and combined,
            # we take that final number and set it as the forward variable:
            # the likelihood we observe prefix (...t) of the observe sequence and end up in state j
            alpha[t][j] = lprob
    return alpha


def calc_beta(A, B, Ob, N, T):
    beta = np.zeros((T + 1, N))
    for t in range(T - 1, -1, -1):
        k = Ob[t]
        for j in range(N):
            lprob = log_zero
            for i in range(N):
                lp = beta[t + 1][i] + A[j][i] + B[j][k]
                #    print(np.exp(lprob))
                #    print(np.exp(lp))
                #    print(np.exp(lprob) + np.exp(lp))
                #    print('--------------------')
                lprob = logaddexp(lprob, lp)
            # print('--------------------')
            beta[t][j] = lprob
    return beta


def calc_xi(A, B, Ob, N, T, alpha, beta):
    xi = np.zeros((T, N, N))
    for t in range(T):
        k = Ob[t]
        lp_traverse = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
                # print("t = " + str(t))
                # print("k = " + str(k))
                # print("i = " + str(i))
                # print("j = " + str(j))
                # P(getting to this arc)
                # P(making this transition)
                # P(emitting this character)
                # P(going to the end)
                lp = (
                        alpha[t][i]
                        + A[i][j]
                        + B[i][k]
                        + beta[t + 1][j]
                )
                # print("alpha[" + str(t) + "][" + str(i) + "] = log(" + str(np.exp(alpha[t][i])) + ")")
                # print("A[" + str(i) + "][" + str(j) + "] = log(" + str(np.exp(A[i][j])) + ")")
                # print("B[" + str(i) + "][" + str(k) + "] = log(" + str(np.exp(B[i][k])) + ")")
                # print("beta[" + str(t+1) + "][" + str(j) + "] = log(" + str(np.exp(beta[t+1][j])) + ")")
                lp_traverse[i][j] = lp
                # print("lp_traverse[" + str(i) + "][" + str(j) + "] = log(" + str(np.exp(lp)) + ")")
                # print("---------------------")
        # Normalize the probability for this time step
        xi[t, :, :] = lp_traverse - supp.logsum(lp_traverse)
        # print("---------------------")
        # print("logsum of floor = " + str(np.exp(supp.logsum(lp_traverse))))
        # print("---------------------")
        # print("lp_traverse = " + str(np.exp(lp_traverse)))
        # print("---------------------")
        # print("final xi floor = " + str(np.exp(xi[t, :, :])))
        # print("---------------------")
    return xi


def calc_gamma(xi, N, T):
    gamma = np.zeros((T, N))
    for t in range(T):
        for i in range(N):
            gamma[t][i] = supp.logsum(xi[t, i, :])
    return gamma


# iteratively update A
# A (transition) [i][j] is the sum of all the
# transitions from i to j, normalized by the sum
# of the transitions out of i
def update_A(N, xi, gamma):
    A = np.zeros((N, N))
    # Sum of all the transitions out of state i
    trans_out = np.zeros(N)
    for i in range(N):
        trans_out[i] = supp.logsum(gamma[:, i])
    for i in range(N):
        for j in range(N):
            A[i][j] = supp.logsum(xi[:, i, j]) - trans_out[i]
    return A


# iteratively update B
# B (emission) [i][k] is the sum of all the
# transitions out of i when k is observed
# divided by the sum of transitions out of i
# CHANGED GAMMA INPUT TO XI
def update_B(Ob, N, M, T, xi):
    B = np.zeros((N, M))
    for i in range(N):
        ksum = np.zeros(M) + log_zero  # ksum[k] is the sum of all i with k
        for t in range(T):
            k = Ob[t]
            for j in range(N):
                ksum[k] = logaddexp(ksum[k], xi[t, i, j])
        ksum = ksum - supp.logsum(ksum)  # Normalize
        B[i, :] = ksum
    return B


# iteratively update pi
def update_pi(N, gamma):
    pi = np.zeros(N)
    for i in range(N):
        pi[i] = gamma[0][i]
    return pi


# Creates a Hidden Markov Model to detect Neanderthal Introgression in modern haplotypes
# Input:
#   loci, a list of variant positions measured in kb
#   ancestries, a list of haplotype lists from 4 ancestries (AFR1, AFR2, TEST, NEAN) with 1s (derived) / 0s (ancestral)
# Output:
#   o_results: name of text file
#   TODO: Explain
def hmm(i_loci, i_ancestries, i_true_states, o_results):
    loci = i_loci
    ancestries = i_ancestries
    true_states = i_true_states

    # PRUFER'S PARAMETERS
    # Ancestral switch rate
    # TEST s = .25
    s = 0.0005
    # Prior probability for archaic ancestry at any locus
    # TEST p = .25
    p = 0.01
    # Probability of archaic ancestry conditional on all SNPs in the window being of state "C"
    # TEST u = .9
    u = 0.99
    # Probability cutoff for HMM's "guess" at a true state (HMM must be >=threshold% sure hidden state introgressed)
    threshold = .9
    # State space: state 0 = 'S' for Species (sapiens), state 1 = 'I' for Introgressed (neandertalensis)
    N = 2
    # Observation space:  observation 0 = 'N' or not consistent, observation 1 = 'C' or consistent with introgression
    M = 2
    # Log-likelihood convergence threshold - used to tell when Baum-Welch has gone far enough
    convergence_threshold = 0.01
    # Intialize the start time.
    start = time.time()

    # KIRZ's PARAMETERS (not specified by Prufer)
    # Should results be normalized based on relative probability? Prufer leaves this unclear
    normalized = False
    # Primary Baum-Welch adjustment parameter, to make sure it doesn't go on too long
    optimization_limit = 20
    # Remember that since we count the Naive HMM as BW# = 0, this will result in 4 optimization rounds.
    # If you want to run 5 rounds, put 6

    # PREPROCESSING
    # We begin by extracting the sequence:
    extraction = extract_obs.extract_O(loci, ancestries, true_states)
    # extraction is a tuple made up of an Observation Sequence, which is a string of letters ("NNC..CN")...
    O = extraction[0]
    # O = "NNCCN"  # Dummy Observed Sequence for testing/explanation
    # T is equal to the length of the sequence
    T = len(O)
    # ... and Win_intro_percent, a Dictionary of 500bp-bins and their contents that I included to keep track of the true
    # introgression state windows, and how "covered" each is by introgressed segments. This is crucial for evaluation.
    # It has the structure (Window # -> Percentage of Introgression)
    Win_intro_percent = extraction[1]
    # creates an list in which to store the window numbers of loci that have true introgressed hidden states
    true_intro_windows = []
    for key in Win_intro_percent: true_intro_windows.append(Win_intro_percent[key])
    # transposes the true introgression site list and stores it in a numpy array for the purposes of visual display
    tiw = np.array([true_intro_windows]).T
    # easy way to keep track of window stops and starts
    Windows = extraction[2]

    # index letter observations for future use
    observations = ['N', 'C']
    # Ob is the same as the observation sequence, but with 'N'-> 0 and 'C'-> 1 for quick referencing.
    Ob = [observations.index(label) for label in O]
    # Stage 1: Checkpoint that marks time after windows binned / sequence generated
    stage1 = time.time()

    # SETTING UP THE HMM

    # Initialize A (the Transition Matrix), B (the Emission Matrix), and pi (the Starting Distribution Matrix)
    # All calculations are done in log-space to prevent point-underflows
    # One potential improvement(?) to this method is to do some random massages here to find new local maxima

    # TODO: Draw these arrays out visually to prevent mistakes from switching
    # Transition Array (2x2)
    A = np.array(((1 - s, s), (s, 1 - s)))
    lp_A = np.log(A)
    # Emission Probabilities (2x2)
    B = np.array(((u, 1 - u), (1 - u, u)))
    lp_B = np.log(B)
    # Initial State Distribution (2x1)
    # pi = np.array((p, 1 - p)) this is an example of a switching mistake
    pi = np.array((1 - p, p))
    lp_pi = np.log(pi)

    # Initialize log-likelihood trackers and print initial inference
    logP_old = log_zero
    alpha = calc_alpha(lp_A, lp_B, lp_pi, Ob, N, T)
    logP_new = supp.logsum(alpha[T, :])

    # NAIVE HMM matrices (no Baum-Welch)
    beta = calc_beta(lp_A, lp_B, Ob, N, T)
    xi = calc_xi(lp_A, lp_B, Ob, N, T, alpha, beta)
    gamma = calc_gamma(xi, N, T)

    # Stage 2: Checkpoint that marks time Naive HMM matrices have been generated
    stage2 = time.time()

    # BAUM-WELCH OPTIMIZATION

    # Initializing a dictionary of gammas: this will allow the comparison of estimated likelihoods over rounds of B/W
    # It has the structure (current_optimization or algorithm step number -> tuple (gamma matrix, performance))
    All_gammas = {}

    optimization_count = 0

    # Iterate until convergence is reached between results, performance decreases, or the hard cap is met
    while logP_new - logP_old > convergence_threshold and optimization_count < optimization_limit:

        # calculate variables / fill out matrices
        bw_alpha = calc_alpha(lp_A, lp_B, lp_pi, Ob, N, T)
        bw_beta = calc_beta(lp_A, lp_B, Ob, N, T)
        bw_xi = calc_xi(lp_A, lp_B, Ob, N, T, bw_alpha, bw_beta)
        bw_gamma = calc_gamma(bw_xi, N, T)

        # recording optimization count / performance progress
        if optimization_count >= 1:
            print("Optimization count " + str(optimization_count))
            print("Improvement of " + str(logP_new - logP_old) + " from last model")
            All_gammas[optimization_count] = bw_gamma
        # we set it to just run once
        elif optimization_count == 0:
            All_gammas[optimization_count] = gamma

        # once variables have been calculated and progress displayed, the counter ticks up
        optimization_count += 1

        # update lambda, the underlying assumptions of the HMM
        new_A = update_A(N, bw_xi, bw_gamma)
        new_B = update_B(Ob, N, M, T, bw_xi)
        new_pi = update_pi(N, bw_gamma)

        # recalculate the forward variable (alpha matrix) from the new lambda
        bw_alpha = calc_alpha(new_A, new_B, new_pi, Ob, N, T)

        # continue iterating only if performance improves, or
        # the likelihood of seeing this sequence given this new HMM increases
        logP_old = logP_new
        # compares last two probabilities of the alpha matrix (%chance of seeing the complete prefix)
        # to the old log-probability of seeing the complete prefix given the HMM parameters
        if supp.logsum(bw_alpha[T, :]) > logP_old:
            lp_A, lp_B, lp_pi = new_A, new_B, new_pi
            logP_new = supp.logsum(bw_alpha[T, :])
    # Stage 3: Checkpoint that marks time after B-W is complete
    stage3 = time.time()

    # EXPRESS THE RESULTS IN MATPLOTLIB

    # Makes sure All_gammas is runnign properly
    # print(len(All_gammas.keys()))
    # print(All_gammas[0])
    # print(All_gammas[0].shape)
    # print(All_gammas[0] == gamma)
    # print(All_gammas[optimization_limit - 1] == bw_gamma)

    # TODO: PERFORMANCES
    # Creates Y rows X 4 columns array to record performances, for each Yth step in Baum-Welch
    performances = np.empty(shape=(len(All_gammas.keys()), 4), dtype=float)
    # For each gamma array in all_gammas (there are a number equal to the runs of baum-welch)
    for key in All_gammas:
        # initialize each row of performances to its respective version of the gamma array
        performances[key] = eval.eval_accuracy(tiw, np.exp(All_gammas[key]), normalized, threshold)

    # print(performances)
    # print(performances.shape)
    # # print(performances.shape[0])
    # # print(performances.shape[1])
    # # print(performances[0].shape())

    # TODO: VISUALIZE DATA
    # vis.compare_3(np.exp(gamma), np.exp(bw_gamma), tiw)
    # vis.display_performance(performances)

    # TODO: MEASURE THE PERFORMANCE OF AN HMM'S GAMMA VS THE REAL THING
    # TODO: print(tiw)
    # print(eval.eval_accuracy(tiw, np.exp(bw_gamma), normalized, threshold))
    # print("False Positive Rate, False Negative Rate,
    # True Positive Rate (Sensitivity), True Negative Rate (Specificity)")

    # # Commented code here used to export the gamma matrix for the purposes of displaying it
    # TODO: not showing up for some reason
    # np.savetxt(
    #     '/Users/briankirz/Downloads/temp_gamma_matrix.csv.gz',
    #     np.exp(gamma),
    #     fmt='%1.3f',
    #     delimiter=',',
    # )

    # list = [tp, tn, fp, fn]
    # np.asarray(list)
    # make sure this is a numpy array
    # list = np.array[list]

    # TODO: Create Results Numpyarray
    # Results is a numpy array that will be filled and exported with all the results of a single rep id
    # It has columns containing the following information:
    # Window Start position | Window Stop position | True Introgression % | BW{X} gamma | BW{X+1} gamma | BW{X+2} gamma
    num_windows = len(Windows)


    # adding extra column at the end to show window labels
    results = np.zeros((num_windows, optimization_limit + 4))

    # recording results
    for key in Windows:
        # initializing starts
        results[key-1][0] = Windows[key][0]
        # initializing stops
        results[key-1][1] = Windows[key][1]
        # initializing true introgression percentages
        results[key-1][2] = Win_intro_percent[key]
        # indicating window labels (1 = C, 0 = N)
        results[key-1][3] = Ob[key-1]
    # iterating through all baum-welch gamma matrices
    for g in range(0, optimization_limit):
        # for each particular window position in gamma, what is the percentage change of introgression?
        for w in range(0, num_windows):
            results[w][g + 4] = np.exp(All_gammas[g][w][1])


    np.savetxt('/Users/briankirz/Downloads/testing_results.csv.gz',
               results,
               fmt='%1.3f',
               delimiter=',',
               newline='\n',
               )

    # TODO: Writing to output textfile

#     results_txt = o_results
#     with open(results_txt, 'w') as out:

#         # TODO: Complex regex for deriving the rep id number from the filepath loci
#         rep_id_number = 1
#         out.write('Rep ID #' + str(rep_id_number) + " results:\n\n")

#         out.write('There are {0} consistent sites in the observed sequence'.format(np.count_nonzero(O == 'C')) + '\n')
#         # Commented out because I think it's unnecessary
#         # out.write('The consistent sites observations occur in window(s)\n{0}'.format(np.where(O == 'C')) + '\n')

#         # TODO: Runtime analysis

#         # makes manipulating format easier
#         stage1_time = "{:.2f}".format(stage1 - start)
#         out.write('\nRuntime for generating observation sequence: {0} seconds'.format(stage1_time))

#         stage2 = stage2 - start
#         stage2_minutes = str("{:.0f}".format((stage2 - stage2 % 60) / 60)) + ' minutes '
#         stage2_seconds = str("{:.2f}".format(stage2 % 60)) + ' seconds'
#         stage2_time = stage2_minutes + stage2_seconds
#         out.write('\nRuntime for running Naive HMM: ' + stage2_time)

#         stage3 = stage3 - start
#         stage3_minutes = str("{:.0f}".format((stage3 - stage3 % 60) / 60)) + ' minutes '
#         stage3_seconds = str("{:.2f}".format(stage3 % 60)) + ' seconds'
#         stage3_time = stage3_minutes + stage3_seconds
#         out.write('\nRuntime for ' + str(optimization_count) + ' steps of Baum-Welch: ' + stage3_time)

#         # TODO: True Introgressed Windows
        # out.write('\nTrue Introgressed Windows: ' + )

        # TODO: Windows with >90%  chance of being introgressed according to HMMs %5

        # Naive HMM windows

        # Loop that prints each 5th HMM after 5

        # TODO: RESULTS

        # Naive HMM results

        # HMM % 5 results

        # Best HMM results

        # False Positive Rate:
        # True Positive Rate (sensitivity):
        # False Negative Rate (miss rate):
        # True Negative Rate (specificity):

    return np.exp(gamma)


# test command
# sh hmm.sh var_pos pol_geno_mat intro_pos results.txt
# sh hmm.sh ../cs282_sim_data/rep_id_1_var_pos.csv.gz ../cs282_sim_data/rep_id_1_geno_mat.csv.gz ../cs282_sim_data/rep_id_1_intro_pos.csv.gz rep_id_1_results.txt

hmm(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
