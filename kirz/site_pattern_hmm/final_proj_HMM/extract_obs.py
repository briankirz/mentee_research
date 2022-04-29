import sys
import gzip
import numpy as np
import pandas as pd
import time


# Define a function to split a genotype matrix into non-overlapping windows.
def genotype_matrix_windows(
        variant_positions,
        polarized_genotype_matrix,
        window_size=50_000,
        sequence_length=20_000_000,
):
    # Intialize a dicctionary with the start and stop position for each window.
    windows = {}
    index = 1
    for window_start in range(0, int(sequence_length), int(window_size)):
        windows[index] = [window_start, (window_start + window_size)]
        index += 1
    # Locate what window each variant is in.
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


# Dave was here for a sanity check.
# Extracts the observed sequence (binned)
# Variable positions corresponds to the first index in the genotype matrix
# Different from the intro positions, which are just the start and stop positions of all introgressed segements
# Input: var_pos - filepath to the array of all variable positions (filename rep_id_{REP}_var_pos.csv.gz)
#        pol_gen_mat - filepath to the polarized genotype matrix (filename rep_id_{REP}_polarized_geno_mat.csv.gz)
# Output: the observed sequence and a dictionary containing information about it
def extract_O(variable_positions, polarized_genotype_matrix):
    # var_pos = variable_positions
    # pgm = polarized_genotype_matrix

    # Load the mutated tree sequence.
    # rep_id_1_mts = tskit.load('./cs282_sim_data/rep_id_1_mut_tree_seq.ts')
    # Load the variable positions.
    rep_id_1_var_pos = np.loadtxt('../cs282_sim_data/rep_id_1_var_pos.csv.gz', delimiter=',')
    # Load the genotype matrix.
    rep_id_1_polarized_geno_mat = np.loadtxt('../cs282_sim_data/rep_id_1_geno_mat.csv.gz', dtype=int,
                                             delimiter=',')
    # Load the introgressed region dataframe.
    rep_id_1_intro_pos_df = pd.read_csv('../cs282_sim_data/rep_id_1_intro_pos.csv.gz', float_precision='round_trip')
    # Inspect the tree-sequence summary.
    # rep_id_1_mts

    # Typically 50k long
    # print(rep_id_1_var_pos)
    # print(rep_id_1_polarized_geno_mat)

    # Indexed from 1 - 400

    Windows = genotype_matrix_windows(rep_id_1_var_pos, rep_id_1_polarized_geno_mat, window_size=500)

    # print(Windows[1])
    # print(len(Windows[1])-2)
    # print('This should be smaller than 50,000 {0}'.format(rep_id_1_var_pos[1352]))
    # print('This should be larger than 50,000 {0}'.format(rep_id_1_var_pos[1353]))
    # print(Windows[1][2])
    # print(Windows[1][-1])
    # wind_1_idx = np.asarray(Windows[400][2:], dtype=np.int32)
    # wind_1_geno_mat = rep_id_1_polarized_geno_mat[wind_1_idx,:]
    # print(wind_1_idx)
    # print(wind_1_idx.shape)
    # print(wind_1_geno_mat)
    # print(wind_1_geno_mat.shape)

    # Intialize observed sequence.
    obs_seq = []
    # Define what C would look like.
    c_pattern = np.array([0, 0, 1, 1])
    # Intialize the start time.
    start = time.time()
    # Iterate through all the windows by key to discover variant sites
    for key in Windows:
        # Extract the values for the window key.
        window_vals = Windows[key]
        # Print the tracker for me.
        print('there are {0} variants in window {1}'.format(len(window_vals), key))
        # If there are variants in that window. Does this mean a window with a single 'C' in it gets left out?
        if len(window_vals) > 2:
            # Extract variable positions in that window. [2:] excludes start pos and end pos
            variants = np.asarray(window_vals[2:], dtype=np.int32)
            # Subset the genotype matrix for that window.
            window_geno_mat = rep_id_1_polarized_geno_mat[variants, :]
            print(window_geno_mat)
            # Define what C matrix would look like given an arbitrary number of variants.
            c_mat = np.tile(c_pattern, (window_geno_mat.shape[0], 1))
            # If the C matrix is equal to the windowed matrix declare it consistent.
            if np.array_equal(c_mat, window_geno_mat) == True:
                print('C')
                obs_seq.append('C')
            # Else declare the window non-cosistent.
            else:
                print('N')
                obs_seq.append('N')
        # If there are no variants in the window declare in non-consistent.
        else:
            print('N')
            obs_seq.append('N')


    # Intialize the end time.
    end = time.time()
    ## Convert the observation sequence list to an array.
    obs_seq_array = np.asarray(obs_seq)
    print('there are {0} many consistent observations'.format(np.count_nonzero(obs_seq_array == 'C')))
    print('the consistent observations occur in window(s) {0}'.format(np.where(obs_seq_array == 'C')))
    print('the run time for generating one observed sequence is {0} seconds'.format((end - start) / float(60)))
    # Convert the observed sequence list into one string.
    obs_seq_str = ''.join(obs_seq)
    # print(obs_seq_str)

    # Testing the transformation of intro_pos into windowed percentages for confirmation

    # print(rep_id_1_intro_pos_df)
    # Extract the panda columns into numpy arrays and round.
    # sorting makes iterating easier
    rep_id_1_intro_starts = np.sort(np.round(rep_id_1_intro_pos_df['start'].values))
    rep_id_1_intro_stops = np.sort(np.round(rep_id_1_intro_pos_df['stop'].values))
    rep_id_1_intro_sizes = np.sort(rep_id_1_intro_stops - rep_id_1_intro_starts)

    print("\nIntrogression event start/end positions, or TRUE STATES:")
    print(rep_id_1_intro_starts)
    print(rep_id_1_intro_stops)
    print("\nWindow sizes of each introgression position:")
    print(rep_id_1_intro_sizes)
    print("----------------------")

    # TODO: NEW! USING INTRO_POS, DETERMINE THE PERCENT OF REAL INTROGRESSION IN EACH WINDOW.
    #       Windows = {} where 1 -> [0, 500), 2-> [500, 1000), ..., 40_000 -> [19_999_500, 20_000_000)
    #       Windows_intro_percent = {} where 1 -> [0.], 2-> [.23454], 3-> [1.], 4 -> [.40234], ...
    Windows_intro_percent = {}

    # The index of the true introgression segment in start/stop/sizes
    intro_index = 0
    for key in Windows:
        # we're all done with the last segment and the rest can be called zero
        if intro_index == rep_id_1_intro_sizes.shape[0]:
            Windows_intro_percent[key] = 0.
        else:
            curr_start = int(rep_id_1_intro_starts[intro_index])
            curr_stop = int(rep_id_1_intro_stops[intro_index])
            curr_start_mod = int(rep_id_1_intro_starts[intro_index] % 500)
            curr_stop_mod = int(rep_id_1_intro_stops[intro_index] % 500)
            curr_start_window = int(((curr_start - curr_start_mod) / 500) + 1)
            curr_stop_window = int(((curr_stop - curr_stop_mod) / 500) + 1)
            tiny_intro = curr_stop - curr_start < 500

            if key < curr_start_window:
                Windows_intro_percent[key] = 0.
            elif key == curr_start_window:
                # print(str(curr_start), str(curr_stop), str(curr_start_window * 500))
                # If the introgressed segment is less than 500, (it shouldn't be)
                # Store stop - start in the dictionary and print the error
                if tiny_intro:
                    Windows_intro_percent[key] = curr_stop - curr_start
                    print("Tiny intro (<500bp) of length " + str(Windows_intro_percent[key]) +
                          " located at " + str(curr_start) + " to " + str(curr_start) +
                          " in window " + str(curr_start_window))
                else: # normal case
                    Windows_intro_percent[key] = (Windows[key][1] - curr_start) / 500
            # We're in the middle of the introgressed segment
            elif curr_start_window < key < curr_stop_window:
                Windows_intro_percent[key] = 1.
            elif key == curr_stop_window:
                Windows_intro_percent[key] = (curr_stop - Windows[key][0]) / 500
                # print("key " + str(key) + " is in current stop window")
                # print("Checking proper window placement:\n" + "Introgression event " + str(intro_index+1) + " is placed in " +
                #       "Window " + str(curr_stop_window) + ":" )
                # print("[" + str(Windows[curr_stop_window][0]) +
                #       " ... start: " + str(curr_stop) + " ... " +
                #       str(Windows[curr_stop_window][1]) + ")")

                # stop window is initiated, intro_index now goes to the next one
                intro_index += 1
                # check to make sure that we record the same number of windows as there are segments
                if intro_index > rep_id_1_intro_sizes.shape[0]:
                    print("ERROR: Recorded more windows than there are segments")
                    break
            else: # Error check
                print("----------------------")
                print("ERROR: bug in key iteration for calculation of introgression percentages")
                print("intro index is " + str(intro_index))
                print("key " + str(key))
                print("curr_start " + str(curr_start))
                print("curr_stop " + str(curr_stop))
                print("curr_start_window " + str(curr_start_window))
                print("curr_stop_window " + str(curr_stop_window))
                print("----------------------")
                break

    # CHECK MY WORK
    intro_index = 0
    inside_segment = False
    for key in Windows_intro_percent:
        # Are we at the first or last affected window of an introgressed segment?
        window_edge = 0 < Windows_intro_percent[key] < 1
        if window_edge:
            if not inside_segment: # Segment Start
                inside_segment = True
                print("Introgression segment " + str(intro_index + 1) + " begins at window " +
                      str(key) + " with " + str(Windows_intro_percent[key] * 100) + ' % coverage')
            else: # Segment End
                inside_segment = False
                print("Introgression segment " + str(intro_index + 1) + " ended at window " +
                      str(key) +
                      " with " + str(Windows_intro_percent[key] * 100) + ' % coverage\n')
                intro_index += 1

        # checking whether places outside segments are 0 and places inside segments are 1
        elif inside_segment and Windows_intro_percent[key] != 1:
            print("ERROR: Inside segment and intro % not 1")
        elif not inside_segment and  Windows_intro_percent[key] != 0:
            print("ERROR: OUtside segment and intro % not 0")







    # intro_index = 0
    # # print("Current Start: " + str(curr_start))
    # # print("Current Start Mod: " + str(curr_start_mod))
    # # print("Proper Window for Current Start: " + str(curr_start_window))
    # print("Checking proper window placement:\n" + "Introgression event " + str(intro_index+1) + " is placed in " +
    #       "Window " + str(curr_start_window) + ":" )
    # print("[" + str(Windows[curr_start_window][0]) +
    #       " ... start: " + str(curr_start) + " ... " +
    #       str(Windows[curr_start_window][1]) + ")")






    return obs_seq_str, Windows_intro_percent  # , Bin_dict


# extract_O(sys.argv[1], sys.argv[2])
