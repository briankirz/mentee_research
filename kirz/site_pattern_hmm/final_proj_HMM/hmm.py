import sys
import numpy as np
<<<<<<< HEAD
import extract_obs
import support_code as supp
import visualizer as vis
import eval_accuracy as eval
import warnings

warnings.filterwarnings('ignore', message='divide by zero encountered in log')
warnings.filterwarnings('ignore', message='invalid value encountered in double_scalars')
try:
    logaddexp = np.logaddexp
except AttributeError:
    def logaddexp(logx, logy):
        if logy - logx > 100:
            return logy
        elif logx - logy > 100:
            return logx
        minxy = min(logx, logy)
        return minxy + np.log(np.exp(logx - minxy) + np.exp(logy - minxy))

# TODO: NORMALIZATION NECESSARY
LOG0 = np.NINF


def old_calc_alpha(A, B, pi, Ob, N, T):
    alpha = np.zeros((T, N))
    alpha[0, :] = pi * B[:, Ob[0]]
    for t in range(1, T):
        for j in range(N):
            alpha[t, j] = alpha[t - 1].dot(A[:, j]) * B[j, Ob[t]]
    return alpha


# Calculates the alpha matrix with the forward algorithm
# Matrix is (T+1)xN, where the last column is the total probability of the output
def calc_alpha(A, B, pi, Ob, N, T):
    alpha = np.zeros((T + 1, N))
    # initialize the first row to be the initial values
    alpha[0, :] = pi
    for t in range(1, T + 1):
        # SHOULD THIS BE OB T??
        k = Ob[t - 1]
        for j in range(N):
            # The probability of the state is the sum of the
            # transitions from all the states from time t-1
            lprob = LOG0
            for i in range(N):
                lp = alpha[t - 1][i] + A[i][j] + B[i][k]
                lprob = logaddexp(lprob, lp)
            alpha[t][j] = lprob
    return alpha


def old_calc_beta(A, B, Ob, N, T):
=======
# does this work?
import extract_obs
import support_code as supp
import warnings
warnings.filterwarnings('ignore', message='divide by zero encountered in log')
warnings.filterwarnings('ignore', message='invalid value encountered in double_scalars')

# Calculates the alpha matrix
def calc_alpha(A, B, pi, Ob, N, T):
    alpha = np.zeros((T, N))

    # Make A more readable
    ocn_from_ocn = A[0][0]
    isl_from_ocn = A[0][1]
    ocn_from_isl = A[1][0]
    isl_from_isl = A[1][1]

    # Make states more readable
    ocn = 0
    isl = 1

    # base cases
    # First state was S0 ocean, resulting in first character of sequence:
    # a0(S0) = P(S0) * P([observed char] emitted from S0)
    # P(S0) given by pi matrix: first element is probability of starting in ocean
    # second element is probability of starting in island
    # First state was S1 island, resulting in first character of sequence:
    # a0(S1) = P(S1) * P([observed char] emitted from S1)

    # a0(Ocean) = P(started on ocean) * P(emitted first letter in obs from ocean)
    alpha[0][0] = pi[0] * B[0][Ob[0]]
    # a0(Island) = P(started on island) * P(emitted fist letter in ob from island)
    alpha[0][1] = pi[1] * B[1][Ob[0]]

    # each subsequent entry in alpha represents the likelihood that the HMM
    # started by producing the same first T letters in the sequence
    # T corresponds with row in alpha

    # PROB TOTAL PREFIX WHEN STATE IS OCEAN
    # a(t)(S0) =
    # a(t-1)(S0) * P(last state gave ocn from ocn) * P(this state emitted char)
    # +
    # a(t-1)(S1) * P(last state gave ocn from isl) * P(this state emitted char)

    # PROB TOTAL PREFIX WHEN STATE IS ISLAND
    # a(t)(S1) =
    # a(t-1)(S0) * P(last state gave isl from ocn) * P(this state emitted char)
    # +
    # a(t-1)(S1) * P(last state gave isl from isl) * P(this state emitted char)

    # ALPHA AT CERTAIN STATE IS THE SUM OF THESE TWO

    for row_char in range(1, T):
        latest_char = Ob[row_char]  # points to nucleotide number code of recent char

        # Probability latest char was result of perfect prefix when last state ocn
        a_ocean = alpha[row_char - 1][ocn] * ocn_from_ocn * B[ocn][latest_char] + alpha[row_char - 1][
            isl] * ocn_from_isl * B[ocn][latest_char]

        # Probability latest char was result of perfect prefix when last state isl
        a_island = alpha[row_char - 1][isl] * isl_from_isl * B[isl][latest_char] + alpha[row_char - 1][
            ocn] * isl_from_ocn * B[isl][latest_char]

        # STATE IS OCEAN
        alpha[row_char][0] = a_ocean
        # STATE IS ISLAND
        alpha[row_char][1] = a_island
    return alpha


# Calculates the beta matrix
def calc_beta(A, B, Ob, N, T):
>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa
    beta = np.zeros((T, N))
    # The beta array at row t column N
    # represents the probability that we're going to see the rest of the observed
    # sequence from Ob[t] onward to Ob[T] (its full length) given that we are
    # now on state N
<<<<<<< HEAD
    # Setting last line of Beta to 1
    beta[T - 1] = np.ones((N))
    # Move backwards starting at T-1
=======

    # Setting last line of Beta to 1
    beta[T - 1] = np.ones((N))
    # Move backwards starting at T-1

>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa
    for t in range(T - 2, -1, -1):
        for j in range(N):
            beta[t, j] = (beta[t + 1] * B[:, Ob[t + 1]]).dot(A[j, :])
    return beta


<<<<<<< HEAD
def calc_beta(A, B, Ob, N, T):
    beta = np.zeros((T + 1, N))
    # Setting last line of Beta to 1
    beta[T - 1] = np.ones((N))
    for t in range(T - 1, -1, -1):
        k = Ob[t]
        for i in range(N):
            # The probability of the state is the sum of the
            # transitions from all the states from time t+1.
            lprob = LOG0
            for j in range(N):
                lp = beta[t + 1][j] + A[i][j] + B[i][k]
                lprob = logaddexp(lprob, lp)
            beta[t][i] = lprob
    return beta


def old_calc_xi(A, B, Ob, N, T, alpha, beta):
    xi = np.zeros((T - 1, N, N))
    for t in range(T - 1):
        denominator = np.dot(np.dot(alpha[t, :].T, A) * B[:, Ob[t + 1]].T, beta[t + 1, :])
        # TODO: zeroes debug
        for i in range(N):
            # for j in range(N):
            #: here is standing in for j

            numerator = alpha[t, i] * A[i, :] * B[:, Ob[t + 1]].T * beta[t + 1, :].T
            # TODO: zeroes debug
=======
# Calculates the xi matrix
def calc_xi(A, B, Ob, N, T, alpha, beta):
    xi = np.zeros((T - 1, N, N))

    for t in range(T - 1):
        denominator = np.dot(np.dot(alpha[t, :].T, A) * B[:, Ob[t + 1]].T, beta[t + 1, :])
        for i in range(N):
            # for j in range(N):
            #: here is standing in for j
            numerator = alpha[t, i] * A[i, :] * B[:, Ob[t + 1]].T * beta[t + 1, :].T
>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa
            xi[t, i, :] = numerator / denominator
    return xi


<<<<<<< HEAD
def calc_xi(A, B, Ob, N, T, alpha, beta):
    xi = np.zeros((T, N, N))
    for t in range(T):
        k = Ob[t]
        lp_traverse = np.zeros((N, N))
        for i in range(N):
            for j in range(N):
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
                lp_traverse[i][j] = lp
        # Normalize the probability for this time step
        xi[t, :, :] = lp_traverse - supp.logsum(lp_traverse)
    return xi


def old_calc_gamma(xi, N, T, alpha, beta):
=======
# Calculates the gamma matrix
def calc_gamma(xi, N, T, alpha, beta):
>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa
    gamma = np.zeros((T, N))
    for t in range(T - 1):
        gamma[t][0] = xi[t][0][0] + xi[t][0][1]
        gamma[t][1] = xi[t][1][0] + xi[t][1][1]

    # input last instance of gamma
    denominator = alpha[T - 1][0] * beta[T - 1][0] + alpha[T - 1][1] * beta[T - 1][1]
    gamma[T - 1][0] = alpha[T - 1][0] * beta[T - 1][0] / denominator
    gamma[T - 1][1] = alpha[T - 1][1] * beta[T - 1][1] / denominator

    return gamma


<<<<<<< HEAD
def calc_gamma(xi, N, T, alpha, beta):
    gamma = np.zeros((T, N))
    # TODO: DISCUSS: Sum of all transitions out of state i at time t, Size is T
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
=======
# iteratively update A
def update_A(N, T, xi, gamma):
    A = np.zeros((N, N))

    # The new A[i][j] is equal to a quotient for all i,j {1,...,N}
    # Numerator - the sum of all xi[t][i][j] from t to T-1
    # Denominator - the sum of all gamma[t][i] from t to T-1

    for i in range(N):
        for j in range(N):
            numerator = 0
            denominator = 0
            for t in range(T - 1):
                numerator = numerator + xi[t][i][j]
                denominator = denominator + gamma[t][i]
            A[i][j] = numerator / denominator

>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa
    return A


# iteratively update B
<<<<<<< HEAD
# B (emission) [i][k] is the sum of all the
# transitions out of i when k is observed
# divided by the sum of transitions out of i
# CHANGED GAMMA INPUT TO XI
def update_B(Ob, N, M, T, xi):
    B = np.zeros((N, M))
    for i in range(N):
        ksum = np.zeros(M) + LOG0  # ksum[k] is the sum of all i with k
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
=======
def update_B(Ob, N, M, T, gamma):
    B = np.zeros((N, M))

    # The new B[i][j] is equal to a quotient for all j {1,...,N} and k {1,...,M}
    # Numerator - the sum of all gamma[t][j] from t to T-1 and ot=vk????
    # Denominator - the sum of all gamma[t][i] from t to T-1

    for j in range(N):
        denominator = np.sum(gamma, axis=0)[j]
        for k in range(M):
            numerator = 0
            for t in range(T - 1):
                if Ob[t] == k:
                    numerator += gamma[t][j]
            B[j][k] = numerator / denominator

    return B


# iteratively update pi
def update_pi(N, gamma):
    pi = np.zeros(N)

    for i in range(N):
        pi[i] = gamma[0][i]

>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa
    return pi


# Creates a Hidden Markov Model to detect Neanderthal Introgression in modern haplotypes
# Input:
#   loci, a list of variant positions measured in kb
#   ancestries, a list of haplotype lists from 4 ancestries (AFR1, AFR2, TEST, NEAN) with 1s (derived) / 0s (ancestral)
# Output:
#   TODO: Explain
def hmm(i_loci, i_ancestries):
    loci = i_loci
    ancestries = i_ancestries

    # ASSUMED PARAMETERS, taken from Prufer:
    # Ancestral switch rate
    s = 0.0005
    # Prior probability for archaic ancestry at any locus
    p = 0.01
    # Probability of archaic ancestry conditional on all SNPs in the window being of state "C"
    u = 0.99

<<<<<<< HEAD
    # Uncomment for Simplified Example
    # s = .25
    # p = .25
    # u = .9

    # Probability cutoff for HMM's "guess" at a true state (HMM must be >=threshold% sure hidden state introgressed)
    threshold = .9
    # Should results be normalized based on relative probability? Prufer leaves this unclear
    normalized = True
    # Primary Baum-Welch adjustment parameter, makes sure it doesn't go on too long
    optimization_limit = 10

    # TODO: Extract the sequence here and assign it to O and win_intro_percent
    # TODO: NAIVE IMPLEMENTATION HAS UNDERFLOW BUG - [:19,000] works, [:20_000] does not
    extraction = extract_obs.extract_O(1, 2)
    # Takes the observation sequence, the first element of the tuple that extract_O returns
    O = extraction[0]
    # Dummy observed sequence
    # O = "NNCCN"
    # Takes the Window dictionary where window number corresponds to the percent of introgression, the second elem
    Win_intro_percent = extraction[1]
    # Creates an array to store those states for data display
    true_intro_windows = []
    for key in Win_intro_percent: true_intro_windows.append(Win_intro_percent[key])
    # Transposes the true introgression site array and stores them in a numpy array for the purposes of display
    tiw = np.array([true_intro_windows]).T

    # index letter observations for future use
    observations = ['N', 'C']
    Ob = [observations.index(label) for label in O]
=======
    # TODO: Extract the sequence here and assign it to O and B_Dict
    # TODO: For now, this is where I'll hardcode the sequence
    O = "CCCCNCCNCNCCCCNNNNNNNNNCNCCCCNNNCCCCCCCC"
>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa

    # IMMUTABLE PARAMETERS
    N = 2  # state 0 = 'species', state 1 = 'introgressed'
    M = 2  # observation 0 = 'N' or not consistent, observation 1 = 'C' or consistent with introgression
    T = len(O)
<<<<<<< HEAD
    # set threshold for log-likelihood convergence - BAUM-WELCH ONLY
    convergence_threshold = 0.01
    All_gammas = {} # key is current_optimization and value is (gamma matrix, performance)

    # TODO: Initialize A, B, and pi
    # TODO: DO ALL THE CALCULATIONS IN LOG SPACE TO AVOID UNDERFLOWS
    # TODO: Randomize these slightly???
    # Transition Array (2x2)
    A = np.array(((1 - s, s), (s, 1 - s)))
    lp_A = np.log(A)
    # Emission Probabilities (2x2)
    B = np.array(((u, 1 - u), (1 - u, u)))
    lp_B = np.log(B)
    # Initial State Distribution (2x1)
    # pi = np.array((p, 1 - p))
    pi = np.array((1 - p, p))
    lp_pi = np.log(pi)
    print(pi)

    # TODO: Initialize log-likelihood trackers and print initial inference
    logP_old = LOG0
    alpha = calc_alpha(lp_A, lp_B, lp_pi, Ob, N, T)
    logP_new = supp.logsum(alpha[T, :])

    # NAIVE HMM
    beta = calc_beta(lp_A, lp_B, Ob, N, T)
    xi = calc_xi(lp_A, lp_B, Ob, N, T, alpha, beta)
    gamma = calc_gamma(xi, N, T, alpha, beta)

    # SUB-19k Observation Sequence ONLY
    # old_alpha = old_calc_alpha(A, B, pi, Ob, N, T)
    # old_beta = old_calc_beta(A, B, Ob, N, T)
    # old_xi = old_calc_xi(A, B, Ob, N, T, old_alpha, old_beta)
    # old_gamma = old_calc_gamma(old_xi, N, T, old_alpha, old_beta)

    # TODO: Iterate until convergence is reached or performance decreases
    optimization_count = 0
    while logP_new - logP_old > convergence_threshold and optimization_count < optimization_limit:

        # calculate variables
        bw_alpha = calc_alpha(lp_A, lp_B, lp_pi, Ob, N, T)
        bw_beta = calc_beta(lp_A, lp_B, Ob, N, T)
        bw_xi = calc_xi(lp_A, lp_B, Ob, N, T, bw_alpha, bw_beta)
        bw_gamma = calc_gamma(bw_xi, N, T, bw_alpha, bw_beta)

        if optimization_count >= 1:
            print("Optimization count " + str(optimization_count))
            print("Improvement of " + str(logP_new - logP_old) + " from last model")
            All_gammas[optimization_count] = bw_gamma
        elif optimization_count == 0:
            All_gammas[optimization_count] = gamma

        optimization_count += 1

        # update lambda
        new_A = update_A(N, T, bw_xi, bw_gamma)
        new_B = update_B(Ob, N, M, T, bw_xi)
        new_pi = update_pi(N, bw_gamma)

        # recalculate alpha from lambda'
        bw_alpha = calc_alpha(new_A, new_B, new_pi, Ob, N, T)

        # only continue iterating if performance improves
        logP_old = logP_new
        if supp.logsum(bw_alpha[T, :]) > logP_old:
            lp_A, lp_B, lp_pi = new_A, new_B, new_pi
            logP_new = supp.logsum(bw_alpha[T, :])

    # check to see if there was any improvement
    if optimization_count > 0:
        print('\nadjusted A\n', np.exp(lp_A))
        print('\nadjusted B\n', np.exp(lp_B))
        print('\nadjusted pi\n', np.exp(lp_pi))
        print('______________________________')
        print('\nnaive alpha\n', np.exp(alpha))
        print('\nBW alpha\n', np.exp(bw_alpha))
        print('\nnaive beta\n', np.exp(beta))
        print('\nBW beta\n', np.exp(bw_beta))
        print('\nnaive xi\n', np.exp(xi))
        print('\nBW xi\n', np.exp(bw_xi))
        print('______________________________')
        print('\nnaive gamma\n', np.exp(gamma))
        print('\nnaive gamma shape\n', np.exp(gamma).shape)
        print('\narray of where unlogged gamma has nonzero values\n', np.where(np.exp(gamma)[:, 0] > 0)[0])
        print('\narray of where unlogged gamma has introgression chances above 1%\n', np.where(np.exp(gamma)[:, 1] > .001)[0])
        print('\nBW gamma\n', np.exp(bw_gamma))
        print('\nBW gamma shape\n', np.exp(bw_gamma).shape)

    # # TESTING GAMMA VS INTROGRESSED SEQUENCES
    # iw1_start = 8534
    # iw1_stop = 8614
    # print('\nObserved sequence around introgressed segment 1, covering 500bp windows '
    #       + str(iw1_start) + " to " + str(iw1_stop) + '\n', O[iw1_start:iw1_stop])
    # print('\ngamma at introgressed segment 1:\n', np.exp(gamma[iw1_start:iw1_stop]))
    # iw2_start = 13696
    # iw2_stop = 13802
    # print('\nObserved sequence around introgressed segment 2, covering 500bp windows '
    #       + str(iw2_start) + " to " + str(iw2_stop) + '\n', O[iw2_start:iw2_stop])
    # print('\ngamma at introgressed segment 2:\n', np.exp(gamma[iw2_start:iw2_stop]))
    # iw3_start = 23152
    # iw3_stop = 23170
    # print('\nObserved sequence around introgressed segment 3, covering 500bp windows '
    #       + str(iw3_start) + " to " + str(iw3_stop) + '\n', O[iw3_start:iw3_stop])
    # print('\ngamma at introgressed segment 3:\n', np.exp(gamma[iw3_start:iw3_stop]))
    # iw4_start = 25678
    # iw4_stop = 25842
    # print('\nObserved sequence around introgressed segment 4, covering 500bp windows '
    #       + str(iw4_start) + " to " + str(iw4_stop) + '\n', O[iw4_start:iw4_stop])
    # print('\ngamma at introgressed segment 4:\n', np.exp(gamma[iw4_start:iw4_stop]))
    # iw5_start = 34097
    # iw5_stop = 34249
    # print('\nObserved sequence around introgressed segment 5, covering 500bp windows '
    #       + str(iw5_start) + " to " + str(iw5_stop) + '\n', O[iw5_start:iw5_stop])
    # print('\ngamma at introgressed segment 5:\n', np.exp(gamma[iw5_start:iw5_stop]))
    # iw6_start = 35230
    # iw6_stop = 35276
    # print('\nObserved sequence around introgressed segment 6, covering 500bp windows '
    #       + str(iw6_start) + " to " + str(iw6_stop) + '\n', O[iw6_start:iw6_stop])
    # print('\ngamma at introgressed segment 6:\n', np.exp(gamma[iw6_start:iw6_stop]))

    # Tried to see what was going on at 29k-30k
    # print('\nObserved sequence around incorrect estimation 29k-30k, covering 500bp windows ' +
    #       '\n', O[29_000:30_000])

    # TODO: EXPRESS THE RESULTS IN MATPLOTLIB

    # print(len(All_gammas.keys()))
    # print(All_gammas[0])
    # print(All_gammas[0].shape)
    # print(All_gammas[0] == gamma)
    # print(All_gammas[optimization_limit-1] == bw_gamma)

    # TODO: PERFORMANCES
    # performances = np.empty(shape=(len(All_gammas.keys()), 4), dtype=float)
    # for key in All_gammas:
    #     performances[key] = eval.eval_accuracy(tiw, np.exp(All_gammas[key]), normalized, threshold)

    # print(performances)
    # print(performances.shape)
    # # print(performances.shape[0])
    # # print(performances.shape[1])
    # # print(performances[0].shape())

    # TODO: VISUALIZE DATA
    vis.compare_3(np.exp(gamma), np.exp(bw_gamma), tiw)
    #vis.display_performance(performances)

    # TODO: MEASURE THE PERFORMANCE OF AN HMM'S GAMMA VS THE REAL THING
    #print(eval.eval_accuracy(tiw, np.exp(bw_gamma), normalized, threshold))
    #print("False Postive Rate, False Negative Rate, True Positive Rate (Sensitivity), True Negative Rate (Specificity)")


    # # Commented code here used to export the gamma matrix for the purposes of displaying it
    # np.savetxt(
    #     '/Users/briankirz/Downloads/temp_gamma_matrix.csv.gz',
    #     gamma,

    # list = [tp, tn, fp, fn]
    # np.asarray(list)
    # make sure this is a numpy array
    # list = np.array[list]

    #
    #     fmt='%1.3f',
    #     delimiter=',',
    # )
    return np.exp(gamma)
=======

    # set threshold for log-likelihood convergence
    convergence_threshold = 0.01

    # index letter observations for future use
    observations = ['N', 'C']
    Ob = [observations.index(label) for label in O]

    # TODO: Initialize A, B, and pi
    # Transition Array (2x2)
    A = np.array(((1 - s, s), (s, 1 - s)))
    # Emission Probabilities (2x2)
    B = np.array(((u, 1 - u), (1 - u, u)))
    # Initial State Distribution (2x1)
    pi = np.array((p, 1 - p))

    # TODO: Initialize log-likelihood trackers and print initial inference
    logP_old = np.NINF
    alpha = calc_alpha(A, B, pi, Ob, N, T)
    logP_new = supp.compute_logP(alpha)
    supp.print_results(A, B, pi, Ob, N, T, logP_new)

    # TODO: Temporary Naive Matrices (Comment out when Baum-Welch is implemented)
    # beta = calc_beta(A, B, Ob, N, T)
    # xi = calc_xi(A, B, Ob, N, T, alpha, beta)
    # gamma = calc_gamma(xi, N, T, alpha, beta)

    # TODO: Iterate until convergence is reached or performance decreases
    while logP_new - logP_old > convergence_threshold:

        # calculate variables
        beta = calc_beta(A, B, Ob, N, T)
        xi = calc_xi(A, B, Ob, N, T, alpha, beta)
        gamma = calc_gamma(xi, N, T, alpha, beta)

        # update lambda
        new_A = update_A(N, T, xi, gamma)
        new_B = update_B(Ob, N, M, T, gamma)
        new_pi = update_pi(N, gamma)

        # recalculate alpha from lambda'
        alpha = calc_alpha(new_A, new_B, new_pi, Ob, N, T)

        # only continue iterating if performance improves
        logP_old = logP_new
        if supp.compute_logP(alpha) > logP_old:
            A, B, pi = new_A, new_B, new_pi
            logP_new = supp.compute_logP(alpha)

            # print current inference
            supp.print_results(A, B, pi, Ob, N, T, logP_new)

    # print lambda*
    print('\nA\n', A)
    print('\nB\n', B)
    print('\npi\n', pi)
    print('\nbinned sequence\n', O)
    # print('\ngamma\n', gamma)

    # Test for printing
    # comparison = []
    # for t in range(T):
    #     comparison.append([Ob[t], gamma[t][1]])
    # for c in range(T):
    #     print(comparison[c])
>>>>>>> 5d466732456cdf91ab356cfe485ca8850a68defa


hmm(sys.argv[1], sys.argv[2])
