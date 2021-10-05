import pandas as pd
import sys

# Create windows centered around the x-values
# Windows at the beginning and end of the chromosomes will be set as either [0, window width] or [end-window width, end]
# for the beginning and end, respectively
# These windows will overlap
def create_x_val_windows(in_bed, power, outfn):
    # Load x-values
    bed = pd.read_table(snakemake.input['bed'], header=None, names=['chr', 'start', 'end', 'length'])

    # Window width is 10^X, where X is the power being processed
    win = int(10 ** power) # Truncate window size as int for now - this handles non-integer window sizes

    # Handle windows that run outside the chromosome bounds
    bed['zeroes']    = 0                                        # 0 is the lowest the beginning can be
    bed['check_beg'] = bed['start'] - win // 2                  # Set current beginning value
    bed['win_beg']   = bed[['check_beg', 'zeroes']].max(axis=1) # Check if the beginning runs below 0, set accordingly

    bed['check_end'] = bed['win_beg'] + win                   # Set current end value
    bed['win_end'] = bed[['check_end', 'length']].min(axis=1) # Check if end runs past chromosome length, fix as needed
    bed['win_beg'] = bed['win_end'] - win                     # Fix beginning value based on the final end value

    # Sort output
    bed = bed.sort_values(['chr', 'win_beg', 'win_end'], ascending=(True, True, True))

    bed.to_csv(outfn, sep='\t', columns=['chr', 'win_beg', 'win_end', 'start', 'end'], header=False, index=False)

create_x_val_windows(snakemake.input['bed'], float(snakemake.wildcards['pow']), snakemake.output['bed'])

