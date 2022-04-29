import sys
import numpy as np
# does this work?
import extract_obs
import support_code as supp
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


# TODO: IS THIS NECESSARY?
VERY_SMALL_NUMBER = 1e-300
LOG0 = np.log(VERY_SMALL_NUMBER)

# OLD IMPLEMENTATION
# Calculates the alpha matrix
# def calc_alpha(A, B, pi, Ob, N, T):
#     alpha = np.zeros((T, N))
#     alpha[0, :] = pi * B[:, Ob[0]]
#     for t in range(1, T):
#         for j in range(N):
#             # TODO: Solve float point underflow issue here
#             alpha[t, j] = alpha[t-1].dot(A[:, j]) * B[j, Ob[t]]
#     return alpha
# NEW IMPLEMENTATION WITH LOGS
# Calculates the alpha matrix
def calc_alpha(A, B, pi, Ob, N, T):
    alpha = np.zeros((T, N))
    # initialize the first row to be the initial values
    alpha[0, :] = pi * B[:, Ob[0]]
    for t in range(1, T):
        # SHOULD THIS BE OB T-1??
        k = Ob[t]
        for j in range(N):
            # The probability of the state is the sum of the
            # transitions from all the states from time t-1
            lprob = LOG0
            for i in range(N):
                lp = alpha[t-1][i] + A[i][j] + B[i][k]
                lprob = logaddexp(lprob, lp)
            alpha[t][j] = lprob
    return alpha

# OLD IMPLEMENTATION
# Calculates the beta matrix
# def calc_beta(A, B, Ob, N, T):
#     beta = np.zeros((T, N))
#     # The beta array at row t column N
#     # represents the probability that we're going to see the rest of the observed
#     # sequence from Ob[t] onward to Ob[T] (its full length) given that we are
#     # now on state N
#
#     # Setting last line of Beta to 1
#     beta[T - 1] = np.ones((N))
#     # Move backwards starting at T-1
#     for t in range(T - 2, -1, -1):
#         for j in range(N):
#             # TODO: Solve float point underflow issue
#             beta[t, j] = (beta[t + 1] * B[:, Ob[t + 1]]).dot(A[j, :])
#     return beta
# TODO: NEW IMPLEMENTATION WITH LOGS
# Calculates the beta matrix
def calc_beta(A, B, Ob, N, T):
    beta = np.zeros((T, N))
    # Setting last line of Beta to 1
    beta[T - 1] = np.ones((N))
    for t in range(T - 2, -1, -1):
        k = Ob[t+1]
        for i in range(N):
            # The probability of the state is the sum of the
            # transitions from all the states from time t+1.
            lprob = LOG0
            for j in range(N):
                lp = beta[t+1][j] + A[i][j] + B[i][k]
                lprob = logaddexp(lprob, lp)
            beta[t][i] = lprob
    return beta



# OLD IMPLEMENTATION
# Calculates the xi matrix
def calc_xi(A, B, Ob, N, T, alpha, beta):
    xi = np.zeros((T - 1, N, N))
    for t in range(T - 1):
        denominator = np.dot(np.dot(alpha[t, :].T, A) * B[:, Ob[t + 1]].T, beta[t + 1, :])
        # TODO: zeroes debug
        for i in range(N):
            # for j in range(N):
            #: here is standing in for j

            numerator = alpha[t, i] * A[i, :] * B[:, Ob[t + 1]].T * beta[t + 1, :].T
            # TODO: zeroes debug
            xi[t, i, :] = numerator / denominator
    return xi
# NEW IMPLEMENTATION
# Calculates the xi matrix
# def calc_xi(A, B, Ob, N, T, alpha, beta):
#     xi = np.zeros((T - 1, N, N))
#     for t in range(T - 1):
#         normalizer = LOG0
#         for i in range(N):
#             for j in range(N):
#
#
#
#     return xi



# Calculates the gamma matrix
def calc_gamma(xi, N, T, alpha, beta):
    gamma = np.zeros((T, N))
    for t in range(T - 1):
        gamma[t][0] = xi[t][0][0] + xi[t][0][1]
        gamma[t][1] = xi[t][1][0] + xi[t][1][1]

    # input last instance of gamma
    denominator = alpha[T - 1][0] * beta[T - 1][0] + alpha[T - 1][1] * beta[T - 1][1]
    gamma[T - 1][0] = alpha[T - 1][0] * beta[T - 1][0] / denominator
    gamma[T - 1][1] = alpha[T - 1][1] * beta[T - 1][1] / denominator

    return gamma


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

    return A


# iteratively update B
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

    return pi


# Creates a Hidden Markov Model to detect Neanderthal Introgression in modern haplotypes
# Input:
#   loci, a list of variant positions measured in kb
#   ancestries, a list of haplotype lists from 4 ancestries (AFR1, AFR2, TEST, NEAN) with 1s (derived) / 0s (ancestral)
# Output:
#   TODO: Explain
def hmm_underflow(i_loci, i_ancestries):
    loci = i_loci
    ancestries = i_ancestries

    # ASSUMED PARAMETERS, taken from Prufer:
    # Ancestral switch rate
    s = 0.0005
    # Prior probability for archaic ancestry at any locus
    p = 0.01
    # Probability of archaic ancestry conditional on all SNPs in the window being of state "C"
    u = 0.99

    # TODO: Extract the sequence here and assign it to O and B_Dict
    # TODO: For now, this is where I'll hardcode the sequence
    # DEBUG EXTRA [:19,000] works, [:20_000] does not
    # TODO: Fix input
    O = extract_obs.extract_O(1, 2)[:19_000]

    # IMMUTABLE PARAMETERS
    N = 2  # state 0 = 'species', state 1 = 'introgressed'
    M = 2  # observation 0 = 'N' or not consistent, observation 1 = 'C' or consistent with introgression
    T = len(O)
    # set threshold for log-likelihood convergence - BAUM-WELCH ONLY
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
    # logP_old = np.NINF
    alpha = calc_alpha(A, B, pi, Ob, N, T)
    # logP_new = supp.compute_logP(alpha)
    # supp.print_results(A, B, pi, Ob, N, T, logP_new)

    # TODO: Temporary Naive Matrices (Comment out when Baum-Welch is implemented)
    beta = calc_beta(A, B, Ob, N, T)
    xi = calc_xi(A, B, Ob, N, T, alpha, beta)
    gamma = calc_gamma(xi, N, T, alpha, beta)

    # TODO: Iterate until convergence is reached or performance decreases
    # while logP_new - logP_old > convergence_threshold:
    #
    #     # calculate variables
    #     beta = calc_beta(A, B, Ob, N, T)
    #     xi = calc_xi(A, B, Ob, N, T, alpha, beta)
    #     gamma = calc_gamma(xi, N, T, alpha, beta)
    #
    #     # update lambda
    #     new_A = update_A(N, T, xi, gamma)
    #     new_B = update_B(Ob, N, M, T, gamma)
    #     new_pi = update_pi(N, gamma)
    #
    #     # recalculate alpha from lambda'
    #     alpha = calc_alpha(new_A, new_B, new_pi, Ob, N, T)
    #
    #     # only continue iterating if performance improves
    #     logP_old = logP_new
    #     if supp.compute_logP(alpha) > logP_old:
    #         A, B, pi = new_A, new_B, new_pi
    #         logP_new = supp.compute_logP(alpha)
    #
    #         # print current inference
    #         supp.print_results(A, B, pi, Ob, N, T, logP_new)

    # print lambda*
    # print('\nA\n', A)
    # print('\nB\n', B)
    # print('\npi\n', pi)
    # print('\nbinned sequence\n', O)
    print('\nalpha\n', alpha)
    print('\nbeta\n', beta)
    print('\ngamma\n', gamma)
    print('gamma shape\n', gamma.shape)
    print('array of where gamma has nonzero values\n', np.where(gamma[:, 0] > 0)[0])

    # Test for printing
    # comparison = []
    # for t in range(T):
    #     comparison.append([Ob[t], gamma[t][1]])
    # for c in range(T):
    #     print(comparison[c])

    # Commented code here used to export the gamma matrix for the purposes of displaying it
    # np.savetxt(
    #     '/Users/briankirz/Downloads/temp_gamma_matrix.csv.gz',
    #     gamma,
    #     fmt='%1.3f',
    #     delimiter=',',
    # )


hmm_underflow(sys.argv[1], sys.argv[2])
