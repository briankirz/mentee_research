import sys
import numpy as np
import warnings

warnings.filterwarnings('ignore', message='divide by zero encountered in log')
warnings.filterwarnings('ignore', message='invalid value encountered in double_scalars')


# SUPPORT CODE
# Stochastic Parameters
# def initialize():
#     std = 0.05
#
#     A = np.array([[np.random.normal(0.7, std), np.random.normal(0.3, std)],
#                   [np.random.normal(0.4, std), np.random.normal(0.6, std)]])
#
#     B = np.array([[np.random.normal(0.3, std), np.random.normal(0.2, std),
#                    np.random.normal(0.2, std), np.random.normal(0.3, std)],
#                   [np.random.normal(0.2, std), np.random.normal(0.3, std),
#                    np.random.normal(0.3, std), np.random.normal(0.2, std)]])
#     pi = [np.random.normal(0.8, std), np.random.normal(0.2, std)]
#
#     A /= np.tile(np.sum(A, axis=1), (2, 1)).T
#     B /= np.tile(np.sum(B, axis=1), (4, 1)).T
#     pi /= np.sum(pi)
#     return A, B, pi

# Viterbi Algorithm (Most likely sequence)
def viterbi(A, B, pi, Ob, N, T):
    V = np.zeros((T, N))
    tb = np.zeros((T, N), dtype='int')

    for i in range(N):
        V[0, i] = np.log(pi[i]) + np.log(B[i, Ob[0]])
        tb[0, i] = -1

    for t in range(1, T):
        for j in range(N):
            V[t, j] = np.NINF
            for i in range(N):
                p = V[t - 1, i] + np.log(A[i, j]) + np.log(B[j, Ob[t]])
                if p > V[t, j]:
                    V[t, j] = p
                    tb[t, j] = i

    logP_max = np.NINF

    for j in range(N):
        if V[-1, j] > logP_max:
            logP_max = V[-1, j]
            Q = [j]
            i = j
    for t in range(T - 2, -1, -1):
        i = tb[t + 1, i]
        Q.insert(0, i)

    return Q, logP_max


def compute_logP(alpha):
    if alpha[-1][0] >= 0:
        return np.log(np.sum(alpha[-1]))
    else:
        return np.log(np.sum(np.exp(alpha[-1])))


def print_results(A, B, pi, Ob, N, T, logP_new):
    Q, logP_max = viterbi(A, B, pi, Ob, N, T)
    print('logP(O|lambda): {0:.2f}'.format(logP_new),
          ''.join([str(i) for i in Q]),
          'logP(O|Q*): {0:.2f}'.format(logP_max))


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
    beta = np.zeros((T, N))
    # The beta array at row t column N
    # represents the probability that we're going to see the rest of the observed
    # sequence from Ob[t] onward to Ob[T] (its full length) given that we are
    # now on state N

    # Setting last line of Beta to 1
    beta[T - 1] = np.ones((N))
    # Move backwards starting at T-1

    for t in range(T - 2, -1, -1):
        for j in range(N):
            beta[t, j] = (beta[t + 1] * B[:, Ob[t + 1]]).dot(A[j, :])
    return beta


# Calculates the xi matrix
def calc_xi(A, B, Ob, N, T, alpha, beta):
    xi = np.zeros((T - 1, N, N))

    for t in range(T - 1):
        denominator = np.dot(np.dot(alpha[t, :].T, A) * B[:, Ob[t + 1]].T, beta[t + 1, :])
        for i in range(N):
            # for j in range(N):
            #: here is standing in for j
            numerator = alpha[t, i] * A[i, :] * B[:, Ob[t + 1]].T * beta[t + 1, :].T
            xi[t, i, :] = numerator / denominator
    return xi


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


# Extracts the observed sequence (binned)
# Input:
#   loci, a list of variant positions measured in kb
#   ancestries, a list of haplotype lists from 4 ancestries (AFR1, AFR2, TEST, NEAN) with 1s (derived) / 0s (ancestral)
#   s, the ancestry switch rate
# Output: the observed sequence and a dictionary containing information about it
def extract_O(i_loci, i_ancestries, i_s):
    O = 'NCNCN'
    Bin_dict = {}
    return O, Bin_dict


# Creates a Hidden Markov Model to detect Neanderthal Introgression in modern haplotypes
# Input:
#   loci, a list of variant positions measured in kb
#   ancestries, a list of haplotype lists from 4 ancestries (AFR1, AFR2, TEST, NEAN) with 1s (derived) / 0s (ancestral)
# Output:
#   ??
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

    # TODO: Extract the sequence here and assign it to O and B_Dict
    # TODO: For now, this is where I'll hardcode the sequence
    O = "CCCCNCCNCNCCCCNNNNNNNNNCNCCCCNNNCCCCCCCC"

    # IMMUTABLE PARAMETERS
    N = 2  # state 0 = 'species', state 1 = 'introgressed'
    M = 2  # observation 0 = 'N' or not consistent, observation 1 = 'C' or consistent with introgression
    T = len(O)

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
    logP_new = compute_logP(alpha)
    print_results(A, B, pi, Ob, N, T, logP_new)

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
        if compute_logP(alpha) > logP_old:
            A, B, pi = new_A, new_B, new_pi
            logP_new = compute_logP(alpha)

            # print current inference
            print_results(A, B, pi, Ob, N, T, logP_new)

    # print lambda*
    print('\nA\n', A)
    print('\nB\n', B)
    print('\npi\n', pi)
    print('\nbinned sequence\n', O)
    # print('\ngamma\n', gamma)

    # Test for printing


hmm(sys.argv[1], sys.argv[2])
