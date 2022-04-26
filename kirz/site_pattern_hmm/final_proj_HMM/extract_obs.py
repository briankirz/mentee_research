import sys
import gzip

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
# Input: var_pos - an array of all variable positions (filename rep_id_{REP}_var_pos.csv.gz)
#        pol_gen_mat - a polarized genotype matrix (filename rep_id_{REP}_polarized_geno_mat.csv.gz)
# Output: the observed sequence and a dictionary containing information about it
def extract_O(variable_positions, polarized_genotype_matrix):
    with gzip.open(variable_positions, 'rb') as var_pos, open(polarized_genotype_matrix, 'w') as pol_gen_mat:
        vp = []
        for line in var_pos:
            vp.append(line)

        # pgm = []
        # for line in pol_gen_mat:
        #     pgm.append(line)

    print(vp)
    O = 'NCNCN'
    Bin_dict = {}
    print("test")
    return O, Bin_dict


extract_O(sys.argv[1], sys.argv[2])
