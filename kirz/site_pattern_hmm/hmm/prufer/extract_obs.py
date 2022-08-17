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


# def calc_window_intro_percent(Binned_windows, true_introgression_positions):
#     # Creates a dictionary of windows that correspond to the % they are covered by introgressed segments
#     # INPUT:
#     # Dictionary variant_Windows = {} where 1 -> [0, 500), 2 -> [500, 100), ..., 40_000 -> [19_999_500, 20_000_000)
#     #                             The key is the window number from 1 to 40,000 and the value is an array.
#     #                             The first two elements of this array are the start and stop positions of the window
#     #                             The following elements of the array are the positions of the variants, if any
#     # nparray true_introgression_positions = an nparray with the start and stop locations of the introgressed segments
#     # Output:
#     # - Wip (window_introgression_percent), a Dictionary of 500bp-bins and their contents that I included to keep track of
#     # the true introgression state windows, and how "covered" each is by introgressed segments. This is crucial for evaluation.
#     # It has the structure (Window # -> Percentage of Introgression as a float between 0 and 1

#     Windows = Binned_windows
#     true_intro_pos = true_introgression_positions
#     # Initializing dictionary of Window Introgression Percentages
#     Win_intro_percent = {}

#     # Extract the columns into numpy arrays and round.
#     # Sorting makes iterating easier. Not changing any start positions. intro_starts is official starting positions
#     intro_starts = np.sort(np.round(true_intro_pos[:, 0]))
#     # print('Starts: {0}'.format(intro_starts))
#     intro_stops = np.sort(np.round(true_intro_pos[:, 1]))
#     # print('Stops: {0}'.format(intro_stops))
#     intro_sizes = np.sort(intro_stops - intro_starts)
#     # print('Sizes: {0}'.format(intro_sizes))

#     # The index of the true introgression segment in start/stop/sizes
#     intro_index = 0
#     for key in Windows:
#         # if intro_index is the same as the number of true introgressed segments, we can end and assign the rest 0
#         if intro_index == intro_sizes.shape[0]:
#             Win_intro_percent[key] = 0.
#         else:
#             # Tracking indices
#             # integer starting and ending positions of the true introgressed segments
#             curr_start = int(intro_starts[intro_index])
#             curr_stop = int(intro_stops[intro_index])
#             # integer offset of curr_start and curr_stop from most recent window
#             curr_start_mod = int(intro_starts[intro_index] % 500)
#             curr_stop_mod = int(intro_stops[intro_index] % 500)
#             # current window that contains the beginning or end of the current segment
#             curr_start_window = int(((curr_start - curr_start_mod) / 500) + 1)
#             curr_stop_window = int(((curr_stop - curr_stop_mod) / 500) + 1)            
#             # boolean that tracks whether the segment falls completely within a window
#             tiny_intro = curr_stop - curr_start < 500
            
            
#             if key == 2229:
#                 print(tiny_intro)
#                 print(key)
#                 print(curr_start)
#                 print(curr_stop)
#                 print(curr_start_mod)
#                 print(curr_stop_mod)
#                 print(curr_start_window)
#                 print(curr_stop_window)
#                 break
                
                
#             # skips windows that come before the current start window
#             if key < curr_start_window:
#                 Win_intro_percent[key] = 0.
#             elif key == curr_start_window:
#                 # If the introgressed segment is less than 500, we need to do a special case to find the percentage
#                 if tiny_intro:
#                     Win_intro_percent[key] = (curr_stop - curr_start) / 500
#                     # print("Tiny intro (<500bp) of length " + str(curr_stop - curr_start) +
#                     #       " located at " + str(curr_start) + " to " + str(curr_stop) +
#                     #       " in window " + str(curr_start_window) + " covers " + str(Win_intro_percent[key]) + "%")
#                     intro_index += 1
#                 else:  # normal case
#                     # calculates the % of the window that is covered by the segment from curr_start to the window's end
#                     Win_intro_percent[key] = (Windows[key][1] - curr_start) / 500
#                     # print("key " + str(key) + " is in current start window")
#                     # print("Checking proper window placement:\n" + "Introgression event " + str(
#                     #     intro_index + 1) + " is placed in " +
#                     #       "Window " + str(curr_start_window) + ":")
#                     # print("[" + str(Windows[curr_start_window][0]) +
#                     #       " ... start: " + str(curr_start) + " ... " +
#                     #       str(Windows[curr_stop_window][1]) + ")")
#             # If we get here, we're in the middle of the introgressed segment
#             elif curr_start_window < key < curr_stop_window:
#                 Win_intro_percent[key] = 1.
#             # If we get here, we're in the last window containing the segment. It should be partially introgressed.
#             elif key == curr_stop_window:
#                 Win_intro_percent[key] = (curr_stop - Windows[key][0]) / 500
#                 # print("key " + str(key) + " is in current stop window")
#                 # print("Checking proper window placement:\n" + "Introgression event " + str(
#                 #     intro_index + 1) + " is placed in " +
#                 #       "Window " + str(curr_stop_window) + ":")
#                 # print("[" + str(Windows[curr_stop_window][0]) +
#                 #       " ... stop: " + str(curr_stop) + " ... " +
#                 #       str(Windows[curr_stop_window][1]) + ")")

#                 # stop window is initiated, intro_index now goes to the next one
#                 intro_index += 1
#                 # check to make sure that we record the same number of windows as there are segments
#                 if intro_index > intro_sizes.shape[0]:
#                     print("ERROR: Recorded more windows than there are segments")
#                     break
#             else:  # Error check
#                 print("----------------------")
#                 print("ERROR: bug in key iteration for calculation of introgression percentages")
#                 print("intro index is " + str(intro_index))
#                 print("key " + str(key))
#                 print("curr_start " + str(curr_start))
#                 print("curr_stop " + str(curr_stop))
#                 print("curr_start_window " + str(curr_start_window))
#                 print("curr_stop_window " + str(curr_stop_window))
#                 print("----------------------")
#                 break

    # CHECKING WORK
    # tracking first window in 100% run
    # seg_start = 0
    # for key in Win_intro_percent:
    #     # looking for the beginnings of segments
    #     if Win_intro_percent[key] != 0:
    #         # not completely covered
    #         if Win_intro_percent[key] != 1:
    #             # initialize segment tracker
    #             # if we've already initialized a segment
    #             if seg_start != 0:
    #                 print("Windows {0} to {1} are 100% introgressed".format(seg_start, key - 1))
    #                 print("Window {0} is {1}% introgressed".format(key, Win_intro_percent[key] * 100))
    #                 print("-------------------------------")
    #                 seg_start = 0
    #             else:
    #                 print("Window {0} is {1}% introgressed".format(key, Win_intro_percent[key] * 100))
    #         # if we're in a segment and the percentage is 1
    #         else:
    #             if seg_start == 0:
    #                 seg_start = key
    #
    # return Win_intro_percent


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

    
        
    # TODO: dxy and window_threshold
    

    # Iterate through all the windows by key.
    for key in Windows:
        # Extract the values for the window key.
        window_vals = Windows[key]
        
        # TODO IF TIME: Make a little graph of the distribution of the number of variant sites per window
        
        # Typically Windows[key] starts with [start, stop, ...].
        # If there are 1 or more variants then the length is greater than 2
        if len(window_vals) > 2:
            # Extract variable positions in that window. [2:] excludes start pos and end pos
            variants = np.asarray(window_vals[2:], dtype=np.int32)
            # Subset the genotype matrix for that window.
            window_geno_mat = pol_geno_mat[variants, :]
            # Keeping tally of consistent sites so we determine if the window is above threshold
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
            
            # Determines the window label
            c_site_proportion = c_sites_tally / total_sites
            if c_site_proportion >= window_threshold:
                print('C')
                print('C site proportion: ' + str(c_site_proportion*100) + '%')
                obs_seq.append('C')
                
            else:
                print('N')
                print('C site proportion: ' + str(c_site_proportion*100) + '%')
                obs_seq.append('N')
                    
            print(window_geno_mat)
            print('---------------')

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


extract_O(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6])
