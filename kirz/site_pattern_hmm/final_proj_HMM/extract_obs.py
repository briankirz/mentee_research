import sys


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
# Input:
#   loci, a list of variant positions measured in kb
#   ancestries, a list of haplotype lists from 4 ancestries (AFR1, AFR2, TEST, NEAN) with 1s (derived) / 0s (ancestral)
#   ????? s, the ancestry switch rate
# Output: the observed sequence and a dictionary containing information about it
def extract_O(i_loci, i_ancestries):
    loci = i_loci
    ancestries = i_ancestries

    afr1 = i_ancestries[0]
    afr2 = i_ancestries[1]
    test = i_ancestries[2]
    nean = i_ancestries[3]

    O = 'NCNCN'
    Bin_dict = {}
    print("test")
    return O, Bin_dict


extract_O(sys.argv[1], sys.argv[2])
