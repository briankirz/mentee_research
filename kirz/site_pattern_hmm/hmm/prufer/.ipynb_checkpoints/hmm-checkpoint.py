import sys
import re
import gzip
import numpy as np
import extract_obs

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

# Set "log_zero" to negative infinity to help with log sums. It will be taken in as an argument in base cases
log_zero = np.NINF

def logsum(array):
    # If the array is of multiple dimensions, it is 'flattened' along one dimension
    if len(array.shape) > 1:
        vec = np.reshape(array, (np.product(array.shape),))
    else:
        vec = array
    # the recurrence relation has to include a base case
    # before the is initialized is initialized the base case is negative infinity,
    # which has an underlying probability of zero from a logaddexp perspective
    sum = np.NINF
    for num in vec:
        sum = logaddexp(sum, num)
    return sum

def calc_alpha(A, B, pi, Ob, N, T):
    # Must be T+1 columns in the alpha matrix bc the top one is the state after 0 prefix characters
    # This is the same as the initial distribution likelihood found in pi
    alpha = np.zeros((T + 1, N))
    # initialize the first row to be the initial distribution values
    # represents the probabilities of being in some state (S/1st or I/2nd) before seeing any (t=0) observed emissions
    alpha[0, :] = pi
    
    # Compute each row, starting with 2nd row. 1st row filled in last step.
    # t counts the character number in the sequence.
    for t in range(1, T + 1):
        
        # k stores the character of the previous observed emission
        k = Ob[t - 1]
        # Compute each column, starting with 1st (Species state) then 2nd (Introgression state)
        for j in range(N):
            
            # Placeholder is set to negative infinity the first time each cell is encountered, resetting it.
            # It stores a probability interpreted by logaddexp as zero when calculating the first logsum
            lprob = np.NINF
            # The i loop occurs in a single cell.
            # Inside the cell, calculate the sum of probabilities (in log form) of the transitions from all possible
            # previous states in time t-1 (the previous row) into the new state j.
            # In this case, N=2, meaning there were 2 possible previous states that could have led to the current one
            # This code answers: "What is the probability that each possible scenario (previous state being S or I) led to
            # our current state j?" When the loop is finished, the value of the cell is set to the combination of those probabilities.
            for i in range(N):
                
                # lp represents a sum of log probabilities:
                # (forward variable at time t-1 for state i)
                # + likelihood that last row's state i transitioned to this state j using the transition matrix A
                # + the likelihood that state i emitted this observed character k using the emission matrix B
                lp = alpha[t - 1][i] + A[i][j] + B[i][k]
                # during the first iteration, lprob is reset as equal to lp, as lprob starts set to NINF
                # the second time around, lp is recalculated and represents the probability that the current state j
                # was reached from the Introgressed state. Now, calling logaddexp(lprob, lp) represents the sum of these:
                # (prob we're in state j if the last state was S + prob we're in state j if the last state was I)
                lprob = logaddexp(lprob, lp)
                
            # After the probabilities based on both of the cells in the previous row were treated and combined,
            # the final number is set as the forward variable:
            # the likelihood we observe prefix (...t) of the observed sequence and end up in state j
            alpha[t][j] = lprob
    return alpha


def calc_beta(A, B, Ob, N, T):
    # Must be T+1 columns in the beta matrix because the bottom one is the state before a 0-character suffix
    # This is given as 100% in the base case, so we still initialize the matrix to zeroes.
    # This is because the underlying proability assumed by the logaddexp occurrence is 1 (log(1) = 0).
    beta = np.zeros((T + 1, N))
    
    # Compute each row, starting with the 2nd from the bottom. The bottom row was filled out during initialization.
    # t counts the position of the state relative to the sequence
    for t in range(T - 1, -1, -1):
        
        # k stores the character just after (emitted by) the state being investigated
        k = Ob[t]
        # Compute each column, starting with 1st (Species state) then 2nd (Introgression state)
        for j in range(N):
            
            # Placeholder is set to negative infinity the first time each cell is encountered, resetting it.
            # It stores a probability interpreted by logaddexp as zero when calculating the first logsum
            lprob = np.NINF
            # The i loop occurs in a single cell, the variable iterating over the states in the previously-calculated row
            # Inside the cell, calculate the sum of probabilities (in log form) of the transitions from all possible
            # previous states in time t+1 (the previous/lower row) to the current row t.
            # This code answers: "What is the probability that each state was arrived at through the emission of
            # the most recent suffix character k from our current state in column j and a subsequent transition
            # from state j to i?" When the loop is finished, the value of the cell is set to the combination of those probabilities.
            for i in range(N):
                
                # lp represents a sum of log probabilities:
                # (backward variable at time t+1 for state i)
                # + likelihood that the lower row's state i transitioned to this state j using the transition matrix A
                # + the likelihood that state i emitted this observed character k using the emission matrix B
                lp = beta[t + 1][i] + A[j][i] + B[j][k]
                # during the first iteration, lprob is reset as equal to lp, as lprob starts set to NINF
                # the second time around, lp is recalculated and represents the probability that the current state j
                # transitioned the Introgressed state. Now, calling logaddexp(lprob, lp) represents the sum of these:
                # (prob we're in state j if the next state is S + prob we're in state j if the next state is I)
                lprob = logaddexp(lprob, lp)
                
            # After the proababilities based on both of the cells in the lower row were treated and combined,
            # the final number is set as the backward variable:
            # the likelihood we observe suffix(t...) of the observed sequence as a result of state j
            beta[t][j] = lprob
    return beta


def calc_xi(A, B, Ob, N, T, alpha, beta):
    # Must be T columns in the xi matrix because there are T-1 transitions between observed characters,
    # plus one state change from the state that emitted the last character to final state.
    xi = np.zeros((T, N, N))
    
    # Compute each 2x2 row or "floor" of the matrix from top to bottom. t=0 represents the first transition between observations
    for t in range(T):
        k = Ob[t]
        lp_traverse = np.zeros((N, N))
        
        # These loops will circle each "floor" and calculate each cell at the [i, j]th coordiante of that floor based
        # on the corresponding alpha and beta matrix positions and the transition and emission matrices
        for i in range(N):
            for j in range(N):
                
                # lp, or the probability of this transition, is equal to the sum of
                # P(getting to this state)
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
        # Each "room" on floor t has been calculated. Now that we have the values of all four cells, we can calculate
        # the total probability of all cases on the top floor as the sum of logarithm probabilities within it.
        # When the "floor" loop is over, this next step "subtracts the logs" (divides the probabilities) of each cell
        # by the total probability of floor T.
        # Normalize the probability for this time step (divide by P(O|lambda))
        xi[t, :, :] = lp_traverse - logsum(lp_traverse)
    return xi


def calc_gamma(xi, N, T):
    # Must be T columns in the gamma matrix because there are T observed loci
    gamma = np.zeros((T, N))
    
    # Compute each row, starting with the first and going down. Each corresponds to a locus
    for t in range(T):
        
        # Compute each column, starting with the Species state (i=0) and then the Introgressed state (i=1)
        for i in range(N):
            
            # Sum up the probabilities for state i at this position t by combining all relevant instances
            # where the hidden state could be i at time t
            gamma[t][i] = logsum(xi[t, i, :])
    return gamma


def update_A(N, xi, gamma):
    # Initialize a blank new transition matrix
    A = np.zeros((N, N))
    # Initialize the sum of all transitions out of i
    trans_out = np.zeros(N)
    
    # for every state i in gamma (Species or Neanderthal)
    for i in range(N):
        # how many transitions out of state i were there
        # summing probabilities because a confidence of 1 counts as 1 transition, 50% confidence counts as half, etc.
        trans_out[i] = logsum(gamma[:, i])
    
    # for every starting state i in xi
    for i in range(N):
        # for every receiving state j in xi
        for j in range(N):
            # A (transition) [i][j] is the sum of all the transitions from i to j
            # This normalized by the previously-calculated sum of the total number of inferred transitions from state i
            A[i][j] = logsum(xi[:, i, j]) - trans_out[i]
            
    return A


def update_B(Ob, N, M, T, xi):
    # Initialize a blank new emission matrix
    B = np.zeros((N, M))
    # For every state i
    for i in range(N):
        # Initialize the matrix of all emissions from state i
        # ksum[k] is the sum of all i with k
        ksum = np.zeros(M) + np.NINF
        # for every observed locus t in the sequence
        for t in range(T):
            # set k to the observation at the current locus
            k = Ob[t]
            # for every state j
            for j in range(N):
                # find the sum of all emissions of k from state i when transitioning to each state j and add them
                ksum[k] = logaddexp(ksum[k], xi[t, i, j])
        # Normalize the sum of all emissions of k from that state i by the sum of all emissions at that position
        ksum = ksum - logsum(ksum)
        # Set the new emission matrix to the normalized probability of every type of emission k from state i
        B[i, :] = ksum
    return B


# iteratively update pi
def update_pi(N, gamma):
    # Initialize a blank new initial distribution matrix
    pi = np.zeros(N)
    # for every state i
    for i in range(N):
        # The adjusted chances that a observed sequence will start on state i
        # are set to the probability that the first locus was i in the last iteration of the model
        pi[i] = gamma[0][i]
    return pi


# Creates a Hidden Markov Model to detect Neanderthal Introgression in modern haplotypes
def hmm(i_loci, i_ancestries, i_true_states, rep_id, opt_limit=100, w_threshold = 1., pattern = "patterna", dxy = "nodxy"):
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
    # State space: state 0 = 'S' for Species (sapiens), state 1 = 'I' for Introgressed (neandertalensis)
    N = 2
    # Observation space:  observation 0 = 'N' or not consistent, observation 1 = 'C' or consistent with introgression
    M = 2
    # Log-likelihood convergence threshold - used to tell when Baum-Welch has gone far enough
    convergence_threshold = 0.01
    # NOT NECESSARY TO INCLUDE
    # Probability cutoff for HMM's "guess" at a true state (HMM must be >=threshold% sure hidden state introgressed)
    # threshold = .9

    # KIRZ's PARAMETERS (not specified by Prufer)
    # Primary Baum-Welch adjustment parameter, to make sure it doesn't go on too long
    optimization_limit = int(opt_limit)
    # Remember that since we count the Naive HMM as BW# = 0, this will result in 4 optimization rounds.
    # If you want to run 5 rounds, put 6
    
    # minimum percentage of consistent sites necessary to qualify a window for overall consistency
    window_threshold = float(w_threshold)
    # setting site patterns (patterna = 0011, patternb = 1100, patternc = 0011 or 1100
    if dxy == "dxy":
        booldxy = True
    elif dxy == "nodxy":
        booldxy = False
    else:
        booldxy = False
        print("ERROR: invalid dxy (use 'dxy' or 'nodxy')")

    # PREPROCESSING
    # We begin by extracting the sequence:
    extraction = extract_obs.extract_O(loci, ancestries, true_states, w_threshold, pattern, booldxy)
    
    # extraction is a tuple made up of an Observation Sequence, which is a string of letters ("NNC..CN")...
    O = extraction[0]
    # O = "NNCCN"  # Dummy Observed Sequence for testing/explanation
    # T is equal to the length of the sequence
    T = len(O)
    # ... and Win_intro_percent, a Dictionary of 500bp-bins and their contents included to keep track of the true
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

    # SETTING UP THE HMM

    # Initialize A (the Transition Matrix), B (the Emission Matrix), and pi (the Starting Distribution Matrix)
    # All calculations are done in log-space to prevent point-underflows
    
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
    logP_old = np.NINF
    alpha = calc_alpha(lp_A, lp_B, lp_pi, Ob, N, T)
    logP_new = logsum(alpha[T, :])

    # NAIVE HMM matrices (no Baum-Welch)
    beta = calc_beta(lp_A, lp_B, Ob, N, T)
    xi = calc_xi(lp_A, lp_B, Ob, N, T, alpha, beta)
    gamma = calc_gamma(xi, N, T)

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
        if logsum(bw_alpha[T, :]) > logP_old:
            lp_A, lp_B, lp_pi = new_A, new_B, new_pi
            logP_new = logsum(bw_alpha[T, :])
            

    # CREATING RESULTS
    # Results is a numpy array that will be filled and exported with all the results of a single rep id
    num_windows = len(Windows)
    # adding extra column to show observation labels
    results = np.zeros((num_windows, optimization_limit + 4))
    # The observations are found in column index 3
    observation_col_index = 3

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
    for g in range(0, len(All_gammas)):
        # for each particular window position in gamma, what is the percentage change of introgression?
        for w in range(0, num_windows):
            results[w][g + 4] = np.exp(All_gammas[g][w][1])


# OUTPUTTING RESULTS
    np.savetxt('./hmm_results/results_BW{0}_wt{1}_{2}_{3}_prufer_rep_id_{4}.csv.gz'.format(str(optimization_limit), str(window_threshold), pattern, dxy, rep_id),
               results,
               fmt='%1.3f',
               delimiter='\t',
               newline='\n',
               )
    
    # IF YOU WANT TO USE THIS UPDATE THE FILEPATH        
    # Optional: Change the textfile to reflect 'N' or 'C'
    # instead of 0 or 1 in the Observed label column for better searching
    # with gzip.open('./hmm_results/prufer_results_rep_id_{0}.csv.gz'.format(str(sys.argv[1])), 'rt') as f_in:
    #     with open('./hmm_results/prufer_results_rep_id_{0}_CN.csv.gz'.format(str(sys.argv[1])), 'w') as f_out:
    #         lines = f_in.readlines()
    #         for line in lines:
    #             split_line = re.split('\t', line)
    #             # initialize new line with first element of original 
    #             newline = split_line[0]
    #             # if the current column index is the same as the observation
    #             # we exclude the first element because it doesn't fit the \t recurrence
    #             for column in range (1, len(split_line)):
    #                 if column == observation_col_index:
    #                     # replace the element with C
    #                     if split_line[observation_col_index] == '1.000':
    #                         newline += ('\tC')
    #                     # replace the element with N
    #                     elif split_line[observation_col_index] == '0.000':
    #                         newline += ('\tN')
    #                     else:
    #                         print("ERROR: value in observation column not 1 or 0")
    #                 # copy all other columns normally
    #                 else:
    #                     newline += ('\t'+split_line[column])
    #             # copy complete new line in new file
    #             f_out.write(newline)

                
                
                
# TREATING INPUT

# Read in sys args.
rep = str(sys.argv[1]) # ${REP}
bw_limit = int(sys.argv[2]) # 100
w_threshold = float(sys.argv[3]) # 1.
patterns = str(sys.argv[4]) # patterna, patternb, or patternc
dxy = str(sys.argv[5]) # dxy or nodxy

# Load the simulated data.
var_pos = './sim_data/rep_id_{0}_var_pos.csv.gz'.format(rep)
geno_mat = './sim_data/rep_id_{0}_geno_mat.csv.gz'.format(rep)
intro_pos = './sim_data/rep_id_{0}_intro_pos.csv.gz'.format(rep)



hmm(var_pos, geno_mat, intro_pos, rep, bw_limit, w_threshold, patterns, dxy)    

# hmm(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
