# Import packages.
import msprime
import numpy as np
import sys

# Define IUA model of introgression.
def iua_human_model(f):
    # Intialize demographic model.
    iua_model = msprime.Demography()
    # We assume constant and equal effective population sizes for
    # all lineages.
    iua_model.add_population(name='AFR', initial_size=10_000)
    iua_model.add_population(name='EUR', initial_size=10_000)
    iua_model.add_population(name='NEA', initial_size=10_000)
    iua_model.add_population(name='AMH', initial_size=10_000)
    iua_model.add_population(name='HUM', initial_size=10_000)
    # Introgression from the Neanderthal to the Eurasian lineage
    # occuring 1,600 generations ago with a probability of f.
    iua_model.add_mass_migration(
        time=1_600, source='EUR', dest='NEA', proportion=f,
    )
    # The African and Eurasian lineages merge into the anatomically
    # modern human lineage 4,000 generations ago.
    iua_model.add_population_split(
        time=4_000, derived=['AFR', 'EUR'], ancestral='AMH',
    )
    # The anatomically modern human and Neanderthal lineages merge
    # into the ancestral human lineage 16,000 generations ago.
    iua_model.add_population_split(
        time=16_000, derived=['AMH', 'NEA'], ancestral='HUM',
    )
    return iua_model

# Define a function to extract introgressed regions.
def intro_ls_tracts(
    ts,
    donor_pop_id,
    donor_samp_id,
    recipient_samp_id,
):
    """
    ###########################################################################
    INPUT
        ts: Simulated tree sequence.
        donor_pop_id: ID of the donor population.
        donor_samp_id: ID of the donor sample.
        recipient_samp_id: ID of the recipient sample.
    ---------------------------------------------------------------------------
    OUTPUT: NumPy array of start and end positions of tracts with a history of 
            introgression plus lineage sorting.
    ###########################################################################
    """
    # Intialize and empty list to store the output and tract length stopping
    # point.
    intro_ls_tracts = []
    tract_left = None
    # For all trees in the tree sequence.
    for tree in ts.trees():
        # Record the population of the MRCA node of our two target populations.
        mrca_pop = ts.node(tree.mrca(recipient_samp_id, donor_samp_id)).population
        # Record the left interval.
        left = tree.interval[0]
        # If the two target populations find their MRCA in the donor population,
        # i.e., if there is introgression + lineage sorting, and if we are not
        # already in an introgressed tract.
        if mrca_pop == donor_pop_id and tract_left is None:
            # Record the left position as the start of the introgressed tract.
            tract_left = left
        # Else if the two target populations do not find their MRCA in the
        # donor population, and if the last tree was within and introgressed
        # tract.
        elif mrca_pop != donor_pop_id and tract_left is not None:
            # Record the introgressed tract from its intial start point to the
            # left position of the current tree.
            intro_ls_tracts.append((tract_left, left))
            tract_left = None
    # If the last tree in the tree sequence contains an introgressed tract.
    if tract_left is not None:
        # Then make the right value of the last introgressed tract the sequence
        # length.
        intro_ls_tracts.append((tract_left, ts.sequence_length))
    return np.array(intro_ls_tracts)

# For 1000 replicates.
for rep in list(range(1, 1001)):
    # Simulate a tree sequence.
    ts_f_03 = msprime.sim_ancestry(
        samples=[
            msprime.SampleSet(2, ploidy=1, population='AFR'), # Sample two haplotypes from AFR for our HMM.
            msprime.SampleSet(1, ploidy=1, population='EUR'),
            msprime.SampleSet(1, ploidy=1, population='NEA'),
        ],
        demography=iua_human_model(0.03), # Rate of introgression at 0.03.
        sequence_length=20_000_000, # Generate 20 Mb haplotypes per samples.
        recombination_rate=1e-8,
        record_migrations=True, # Need this to keep track of what segments are introgressed.
        discrete_genome=False,
        random_seed=rep,
    )
    # Overlay mutations.
    mts_f_03 = msprime.sim_mutations(
        tree_sequence=ts_f_03, rate=1.5e-8,
        model='jc69', random_seed=rep,
        discrete_genome=False,
    )
    # Extract the genotype matrix.
    genotype_matrix = mts_f_03.genotype_matrix()
    # Save the genotype matrix.
    np.savetxt(
        './sim_data/rep_id_{0}_geno_mat.csv.gz'.format(rep),
        genotype_matrix,
        fmt='%d',
        delimiter=',',
        )
    # Extract the variable positions.
    variable_positions = mts_f_03.tables.sites.position
    # Save the variable positions.
    np.savetxt(
        './sim_data/rep_id_{0}_var_pos.csv.gz'.format(rep),
        [variable_positions],
        fmt='%1.15f',
        delimiter=',',
        newline='\n',
    )
    # Identify introgressed tracts.
    ex_intro_tracts_array = intro_ls_tracts(
        ts=mts_f_03,
        donor_pop_id=2,
        donor_samp_id=3,
        recipient_samp_id=2,
    )
    # Save the introgressed segments.
    np.savetxt(
        './sim_data/rep_id_{0}_intro_pos.csv.gz'.format(rep),
        ex_intro_tracts_array,
        fmt='%1.15f',
        delimiter=',',
        newline='\n',
    )
    # Save the tree-sequence just for good measures.
    mts_f_03.dump('./sim_data/rep_id_{0}_mut_tree_seq.ts'.format(rep))
    # Progress Checker
    print('Finished simulating rep {0}!'.format(rep))