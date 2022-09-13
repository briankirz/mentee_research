import sys
import gzip
import numpy as np


# Define a function to split a genotype matrix into non-overlapping windows.
def genotype_matrix_windows(
        variant_positions,
        polarized_genotype_matrix,
        window_size=500,
        sequence_length=20000000,
):
    # Intialize a dictionary with the start and stop position for each window.
    windows = {}
    index = 1
    for window_start in range(0, int(sequence_length), int(window_size)):
        windows[index] = [window_start, (window_start + window_size)]
        index += 1
    # Locate what window each variant is in.
    # windows dictionary is now: window # (1-40,000) -> [start (0), stop (500)] with optional variant_index ]
    index = 0
    pos = variant_positions[index]
    for key in windows:
        start, stop = windows[key]
        while start <= pos < stop:
            windows[key].append(index)
            index += 1
            if index < len(variant_positions):
                pos = variant_positions[index]
            else:
                break
    return windows

def calc_window_intro_percent(Binned_windows, true_introgression_positions):
    Windows = Binned_windows
    true_intro_pos = true_introgression_positions
    
    # Initializing dictionary of Window Introgression Percentages
    Win_intro_percent = {}
    # Extract the columns into numpy arrays and round.
    # Sorting makes iterating easier. Not changing any start positions. intro_starts is 'official' starting position
    intro_starts = np.sort(np.round(true_intro_pos[:, 0]))
    intro_stops = np.sort(np.round(true_intro_pos[:, 1]))
    intro_sizes = intro_stops - intro_starts
    start_mods = intro_starts % 500
    stop_mods = intro_stops % 500
    start_windows = ((intro_starts - start_mods) / 500) + 1
    stop_windows = ((intro_stops - stop_mods) / 500) + 1
    
    # Initialize all windows to 0% first
    for key in Windows:
        Win_intro_percent[key] = 0.
    
    # For each segment
    for t in range(intro_sizes.shape[0]):
        start_win = int(start_windows[t])
        stop_win = int(stop_windows[t])
        
        # if segment within a single window
        if start_win == stop_win:
            Win_intro_percent[start_win] += intro_sizes[t] / 500
        # if segment not in a single window
        else:
            # all fully introgressed segments
            for f in range(start_win + 1, stop_win):
                Win_intro_percent[f] = 1.
            
            # clean up starting window
            Win_intro_percent[start_win] += (Windows[start_win][1] - intro_starts[t]) / 500
            # clean up stopping window
            if stop_win <= len(Windows):
                Win_intro_percent[stop_win] += (intro_stops[t] - Windows[stop_win][0]) / 500
            
    return Win_intro_percent

# Extracts the observed sequence (binned)
def extract_O(variable_positions, polarized_genotype_matrix, true_introgression_positions, w_threshold, pattern, dxy):
    

    # load the variant positions
    var_pos = np.loadtxt(variable_positions, delimiter=',')
    # Load the genotype matrix.
    pol_geno_mat = np.loadtxt(polarized_genotype_matrix, dtype=int, delimiter=',')
    # Load the introgressed region dataframe.
    true_intro_pos = np.loadtxt(true_introgression_positions, delimiter=',')
    # set the window threshold, or the proportion of consistent sites necessary to label C
    window_threshold = float(w_threshold)
    # Define what C, a pattern consistent with introgression, would look like.

    # Indexed from 1 - 400
    # Windows is of the format key -> value
    # Window # (1-400) -> [Start position, stop position, (optional var_pos positions)]
    Windows = genotype_matrix_windows(var_pos, pol_geno_mat, window_size=500, sequence_length=20_000_000)
    Wip = calc_window_intro_percent(Windows, true_intro_pos)

    # EXTRACTING OBSERVED SEQUENCE
    # Intialize observed sequence.
    obs_seq = []

    # Iterate through all the windows by key.
    for key in Windows:
        # Extract the values for the window key.
        window_vals = Windows[key]
        
        
        # Typically Windows[key] starts with [start, stop, ...].
        # If there are 1 or more variants then the length is greater than 2
        if len(window_vals) > 2:

            # Extract variable positions in that window. [2:] excludes start pos and end pos
            variants = np.asarray(window_vals[2:], dtype=np.int32)
            # Subset the genotype matrix for that window.
            window_geno_mat = pol_geno_mat[variants, :]
            
            # DXY CALCULATION
            # if we're including a lower Dxy distance as a potential source of window consistency
            # all we're doing here is determining a boolean - whether the window should be considered consistent
            c_by_dxy = False
            if dxy:
                num_pops = len(window_geno_mat[0])
                num_vars = len(window_geno_mat)
                # This is the index of the column population we're testing/referencing: In this case, EUR/NEAN
                # AFR1 | AFR2 | EUR | NEAN
                test_pop_index = 2
                ref_pop_index = 3
                # create a matrix of neanderthal alleles at each variant site (the last column)
                nean = window_geno_mat[:, ref_pop_index]
                
                # initialize the resulting matrix for all non-reference populations observed
                dxy_distances = np.zeros(num_pops-1)
                # the last column is the archaic pop, so we leave it zero
                for col in range(num_pops-1):
                    # create 1-D matrix. elem is allele at each respective variant site for that pop
                    # [[1 1 0 0]
                    #  [0 0 1 1]
                    #  [1 1 1 0]
                    #  [1 1 1 0]]
                    # For example, for col 0, or the first population, which is AFR1, would be the first column
                    # pop = [1 0 1 1]
                    pop = window_geno_mat[:, col]
                    dxy_distances[col] = np.sum((np.multiply(pop, 1-nean) + np.multiply(nean, 1-pop)) / num_vars)
                # dxy_distances now should look something like this:
                # [AFR1dxy, AFR2dxy, EURdxy]
                # If the minimum value is the test_pop_index (test EUR is closest), then c_by_dxy is True
                if np.argmin(dxy_distances) == test_pop_index:
                    c_by_dxy = True
                    
            # TRADITIONAL WINDOW_THRESHOLD
            # Keeping tally of consistent sites so we determine if the window is above threshold
            # if c_by_dxy is true, we don't need to caclulate this and can save time (we know it'll be closer)?
            # if not c_by_dxy:
            c_sites_tally = 0
            total_sites = len(window_vals)-2

            c_pattern_a = np.array([0, 0, 1, 1])
            c_pattern_b = np.array([1, 1, 0, 0])

            # Checking all of the sites in a single window
            for site in window_geno_mat:
                if pattern == "patterna": #0011
                    # If the C matrix is equal to the windowed matrix declare it consistent.
                    if np.array_equal(c_pattern_a, site):
                        c_sites_tally += 1
                elif pattern == "patternb": #1100
                    if np.array_equal(c_pattern_b, site):
                        c_sites_tally += 1
                elif pattern == "patternc": #0011 or 1100
                    if np.array_equal(c_pattern_a, site) or np.array_equal(c_pattern_b, site):
                        c_sites_tally += 1
                else:
                    print("ERROR: Invalid Pattern")
            c_site_proportion = c_sites_tally / total_sites
            c_by_threshold = (c_site_proportion >= window_threshold)
                    
            # DETERMINE WINDOW LABEL

            # if it's not above threshold, but dxy is active and EUR are closer than AFR still label C
            if c_by_threshold or c_by_dxy:
                # print('C')
                # print('C site proportion: ' + str(c_site_proportion*100) + '%')
                obs_seq.append('C')
            else:
                # print('N')
                # print('C site proportion: ' + str(c_site_proportion*100) + '%')
                obs_seq.append('N')
                    
            # print(window_geno_mat)
            # print('---------------')

        # If there are no variants in the window declare in non-consistent.
        else:
            # print('N')
            obs_seq.append('N')



    # Convert the observation sequence list to an array.
    obs_seq_array = np.asarray(obs_seq)

    # print('there are {0} many consistent observations'.format(np.count_nonzero(obs_seq_array == 'C')))
    # print('the consistent observations occur in window(s) {0}'.format(np.where(obs_seq_array == 'C')))
    # print('the run time for generating one observed sequence is {0} minutes'.format((end - start) / float(60)))

    return obs_seq_array, Wip, Windows


# extract_O(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
