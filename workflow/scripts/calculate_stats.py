from functools import reduce
import pandas as pd
import sys

# Extract name of sample from filename
def extract_sample(fname, samp_list):
    for s in samp_list:
        if s in fname:
            return s

# Calculate mean and standard deviation of the methylation values for each x-value
def calculate_stats(in_beds, out_mean, out_stdv, samp_list):
    # Collect all input files
    dfs = []
    for b in in_beds:
        sm = extract_sample(b, samp_list)
        df = pd.read_table(b, header=None, names=['chr', 'start', 'end', 'avg_'+sm], na_values='.')
        dfs.append(df)

    # Collapse to a single data frame
    out = reduce(lambda df1, df2: df1.merge(df2, on=['chr', 'start', 'end']), dfs)

    # Calculate mean and standard deviation
    out['mean'] = out[df.columns[3:]].mean(axis=1)
    out['stdv'] = out[df.columns[3:]].std(axis=1)

    # Write output
    out.to_csv(out_mean, sep='\t', columns=['chr', 'start', 'end', 'mean'], header=False, index=False, na_rep='.')
    out.to_csv(out_stdv, sep='\t', columns=['chr', 'start', 'end', 'stdv'], header=False, index=False, na_rep='.')

calculate_stats(
    list(snakemake.input['files']),
    snakemake.output['mean_gz'],
    snakemake.output['stdv_gz'],
    snakemake.params['sample']
)
