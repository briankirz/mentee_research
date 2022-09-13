import sys
import numpy as np
import gzip

# takes in the highest rep number and ranks all of the reps
# by percent introgressed and number of introgressed segs
# outputs each in a separate file
def rank_rep_introgression(total_reps):
    # assumes genome length of 20,000,000 base pairs
    
    # [[rep1, %], [rep2, %], ... [rep1000, %]]
    intro_percent = np.zeros((total_reps, 2))
    intro_num = np.zeros((total_reps, 2))
    
    # Create tables for numbers and percentages
    for rep in range(1, total_reps+1):
        # set the filepath
        rep_filepath = './sim_data/rep_id_{0}_intro_pos.csv.gz'.format(str(rep))
        # download the true introgressed position array from the rep filepath
        true_intro_pos = np.loadtxt(rep_filepath, delimiter=',')
        
        # assign the rep and number of segments
        intro_num[rep-1][0] = rep
        intro_num[rep-1][1] = len(true_intro_pos)
        
        # find the percentage of the rep that's introgressed
        intro_starts = np.sort(np.round(true_intro_pos[:, 0]))
        intro_stops = np.sort(np.round(true_intro_pos[:, 1]))
        # calculate introgression segment sizes
        intro_sizes = intro_stops - intro_starts
        # sum them all and divide by the size of the genome
        intro_percent[rep-1][0] = rep
        intro_percent[rep-1][1] = np.sum(intro_sizes) / 20_000_000
    
    
    # order both tables
    intro_num = intro_num[intro_num[:, 1].argsort()]
    intro_percent = intro_percent[intro_percent[:, 1].argsort()]
    
    # save number table
    np.savetxt('./ranked_rep_intro_num.csv.gz',
               intro_num,
               fmt='%1.3f',
               delimiter='\t',
               newline='\n',
               )
    np.savetxt('./ranked_rep_intro_percent.csv.gz',
               intro_percent,
               fmt='%1.5f',
               delimiter='\t',
               newline='\n',
               )
        
                   
    # save percentages table
        

    return None

total_reps = int(sys.argv[1])

rank_rep_introgression(total_reps)