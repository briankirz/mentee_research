import sys
import gzip
import numpy as np
import pandas as pd


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


# Extracts the observed sequence (binned)
# Variable positions corresponds to the first index in the genotype matrix
# Different from the intro positions, which are just the start and stop positions of all introgressed segements
# Input: var_pos - filepath to the array of all variable positions (filename rep_id_{REP}_var_pos.csv.gz)
#        pol_gen_mat - filepath to the polarized genotype matrix (filename rep_id_{REP}_polarized_geno_mat.csv.gz)
# Output: the observed sequence and a dictionary containing information about it
def extract_O(variable_positions, polarized_genotype_matrix):

    # Load the mutated tree sequence.
    # rep_id_1_mts = tskit.load('./cs282_sim_data/rep_id_1_mut_tree_seq.ts')
    # Load the variable positions.
    rep_id_1_var_pos = np.loadtxt('../cs282_sim_data/rep_id_1_var_pos.csv.gz', delimiter=',')
    # Load the genotype matrix.
    rep_id_1_polarized_geno_mat = np.loadtxt('../cs282_sim_data/rep_id_1_polarized_geno_mat.csv.gz', dtype=int,
                                             delimiter=',')
    # Load the introgressed region dataframe.
    rep_id_1_intro_pos_df = pd.read_csv('../cs282_sim_data/rep_id_1_intro_pos.csv.gz', float_precision='round_trip')
    # Inspect the tree-sequence summary.
    # rep_id_1_mts

    # Typically 50k long
    # print(rep_id_1_var_pos)
    # print(rep_id_1_polarized_geno_mat)

    # Indexed from 1 - 400
    Windows = genotype_matrix_windows(rep_id_1_var_pos, rep_id_1_polarized_geno_mat)

    print(Windows[1])
    print(len(Windows[1])-2)
    print('This should be smaller than 50,000 {0}'.format(rep_id_1_var_pos[1352]))
    print('This should be larger than 50,000 {0}'.format(rep_id_1_var_pos[1353]))
    print(Windows[1][2])
    print(Windows[1][-1])
    wind_1_idx = np.asarray(Windows[1][2:], dtype=np.int32)
    print(wind_1_idx)
    print(wind_1_idx.shape)
    print(rep_id_1_polarized_geno_mat)
    print(rep_id_1_polarized_geno_mat.shape)



    # throw vp and pgm into genotype_matrix_windows

    O = 'NCNCN'
    Bin_dict = {}
    print("test")
    return O, Bin_dict


extract_O(sys.argv[1], sys.argv[2])
