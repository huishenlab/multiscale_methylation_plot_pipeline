import pandas as pd
import sys

# Add a fourth column to an input BED file of the chromosome length
# Restrict to canonical chromosomes at the same time
def add_chrom_length(in_bed, in_ref, outfn):
    # Load data
    bed = pd.read_table(in_bed, header=None, names=['chr', 'start', 'end'])
    fai = pd.read_table(in_ref, header=None, names=['chr', 'length'], usecols=[0, 1])

    # Replace chromosome name with chromosome length
    set_length = fai.set_index('chr').to_dict()['length']

    bed['length'] = bed['chr']
    bed = bed.replace({'length': set_length})

    # Filter out non-canonical chromosomes and contigs
    filter_non_canonical = bed['chr'].str.contains('^chr[0-9XY]*$')
    bed = bed[filter_non_canonical]

    # Sort output
    bed = bed.sort_values(['chr', 'start', 'end'], ascending=(True, True, True))

    bed.to_csv(snakemake.output[0], sep='\t', columns=['chr', 'start', 'end', 'length'], header=False, index=False)

add_chrom_length(snakemake.input['bed'], snakemake.input['ref'], snakemake.output['bed'])
